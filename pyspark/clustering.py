import dask.dataframe as dd
import dask.array as da
from tdigest import TDigest
import pandas as pd
import numpy as np
from dask import delayed

# Step 1: Load Data
df = dd.read_parquet('path/to/parquet_files/*')

# Step 2: Optimize Data Types
df['PCA_1'] = df['PCA_1'].astype('float32')
df['PCA_2'] = df['PCA_2'].astype('float32')
df['PCA_3'] = df['PCA_3'].astype('float32')

# Step 3: Compute Global Percentile Thresholds for PCA_1
pca1_array = df['PCA_1'].to_dask_array(lengths=True)
percentiles_1 = da.percentile(pca1_array, q=[i * 2 for i in range(1, 50)]).compute()
percentiles_1 = [-float('inf')] + list(percentiles_1) + [float('inf')]

# Step 4: Assign bucket_1
def assign_bucket_pca1(pca1, thresholds):
    return np.searchsorted(thresholds, pca1, side='right') - 1

thresholds_1 = percentiles_1  # Already computed
df['bucket_1'] = df['PCA_1'].map_partitions(
    lambda s: s.apply(lambda x: assign_bucket_pca1(x, thresholds_1)),
    meta=('bucket_1', 'int16')
)

df = df.persist()

# Step 5: Compute TDigest for PCA_2 within each bucket_1
def compute_tdigest_pca2(df_group):
    digest = TDigest()
    for value in df_group['PCA_2']:
        digest.update(value)
    return pd.Series({'tdigest_pca2': digest})

tdigest_pca2_df = df.groupby('bucket_1').apply(compute_tdigest_pca2, meta={'tdigest_pca2': 'object'}).compute()

# Step 6: Extract and Store PCA_2 Thresholds
percentiles_2 = {}
for bucket, row in tdigest_pca2_df.iterrows():
    digest = row['tdigest_pca2']
    percentiles = [digest.percentile(q) for q in range(2, 100, 2)]
    percentiles_2[bucket] = [-float('inf')] + percentiles + [float('inf')]

# Step 7: Assign bucket_2
def get_pca2_thresholds(bucket):
    return percentiles_2.get(bucket, [-float('inf'), float('inf')])

df['thresholds_2'] = df['bucket_1'].map(
    lambda x: get_pca2_thresholds(x),
    meta=('thresholds_2', 'object')
)

def assign_bucket_pca2_row(row):
    return np.searchsorted(row['thresholds_2'], row['PCA_2'], side='right') - 1

df['bucket_2'] = df.map_partitions(
    lambda partition: partition.apply(assign_bucket_pca2_row, axis=1),
    meta=('bucket_2', 'int16')
)

df = df.persist()

# Step 8: Compute TDigest for PCA_3 within each (bucket_1, bucket_2)
def compute_tdigest_pca3(df_group):
    digest = TDigest()
    for value in df_group['PCA_3']:
        digest.update(value)
    return pd.Series({'tdigest_pca3': digest})

tdigest_pca3_df = df.groupby(['bucket_1', 'bucket_2']).apply(compute_tdigest_pca3, meta={'tdigest_pca3': 'object'}).compute()

# Step 9: Extract and Store PCA_3 Thresholds
percentiles_3 = {}
for (bucket1, bucket2), row in tdigest_pca3_df.iterrows():
    digest = row['tdigest_pca3']
    percentiles = [digest.percentile(q) for q in range(2, 100, 2)]
    percentiles_3[(bucket1, bucket2)] = [-float('inf')] + percentiles + [float('inf')]

# Step 10: Assign bucket_3
def get_pca3_thresholds(row):
    return percentiles_3.get((row['bucket_1'], row['bucket_2']), [-float('inf'), float('inf')])

df['thresholds_3'] = df.apply(
    lambda row: get_pca3_thresholds(row),
    axis=1,
    meta=('thresholds_3', 'object')
)

def assign_bucket_pca3_row(row):
    return np.searchsorted(row['thresholds_3'], row['PCA_3'], side='right') - 1

df['bucket_3'] = df.map_partitions(
    lambda partition: partition.apply(assign_bucket_pca3_row, axis=1),
    meta=('bucket_3', 'int16')
)

df = df.persist()

# Step 11: Compute cluster_id
df['cluster_id'] = df['bucket_1'] * 2500 + df['bucket_2'] * 50 + df['bucket_3']
df['cluster_id'] = df['cluster_id'].astype('int32')

# Optional: Drop threshold columns to save space
df = df.drop(['thresholds_2', 'thresholds_3'], axis=1)

# Step 12: Write to Parquet in Partitioned Files
df.to_parquet('path/to/output_clusters/', write_index=False, partition_on=['bucket_1', 'bucket_2'])
