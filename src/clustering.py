from pyspark.sql import SparkSession
from pyspark.ml.feature import Bucketizer

# Increase memory allocation
spark = SparkSession.builder \
    .appName("HierarchicalBinning") \
    .config("spark.driver.memory", "64g") \
    .config("spark.executor.memory", "64g") \
    .getOrCreate()


df = spark.read.parquet("/home/afloresep/work/chelombus/data/230M/output/*")
print(df.count())


# Define the percentiles for 25 bins
percentiles = [i / 25.0 for i in range(1, 25)]

# Compute approximate quantiles for PCA_1
pca1_thresholds = df.approxQuantile("PCA_1", percentiles, 0.01)

# Create splits for Bucketizer
pca1_splits = [-float("inf")] + pca1_thresholds + [float("inf")]

# Initialize Bucketizer
bucketizer_pca1 = Bucketizer(
    splits=pca1_splits,
    inputCol="PCA_1",
    outputCol="bin_PCA1_temp"
)

# Transform the DataFrame
df = bucketizer_pca1.transform(df)

# Convert bin_PCA1_temp to IntegerType and adjust to start from 0
from pyspark.sql.functions import col
df = df.withColumn("bin_PCA1", col("bin_PCA1_temp").cast("integer"))

# Drop the temporary column
df = df.drop("bin_PCA1_temp")

from pyspark.sql.window import Window
from pyspark.sql.functions import ntile

# Create a window specification for PCA_2 within each bin_PCA1
window_pca2 = Window.partitionBy("bin_PCA1").orderBy("PCA_2")

# Assign bins for PCA_2, starting from 0
df = df.withColumn("bin_PCA2", ntile(25).over(window_pca2) - 1)

# Create a window specification for PCA_3 within each bin_PCA1 and bin_PCA2
window_pca3 = Window.partitionBy("bin_PCA1", "bin_PCA2").orderBy("PCA_3")

# Assign bins for PCA_3, starting from 0
df = df.withColumn("bin_PCA3", ntile(25).over(window_pca3) - 1)

df.select("smiles", "PCA_1", "PCA_2", "PCA_3", "cluster_id").show()

# Import necessary functions
from pyspark.sql.functions import col


# Save clusters based on first dimension
for i in range(25):
    df_cluster1 = df.filter(col("bin_PCA1") == i)
    df_cluster1_pd = df_cluster1.toPandas()
    df_cluster1_pd.to_csv(f"/home/afloresep/work/chelombus/data/230M/clustered_output/cluster{i}.csv", index=False)

from pyspark.sql import SparkSession
from pyspark.sql.functions import avg, sqrt, pow, col, row_number, broadcast
from pyspark.sql.window import Window

# Step 1: Compute average PCA values for each cluster
cluster_centers = df.groupBy("cluster_id").agg(
    avg("PCA_1").alias("avg_PCA_1"),
    avg("PCA_2").alias("avg_PCA_2"),
    avg("PCA_3").alias("avg_PCA_3")
).cache()

# Step 2: Join the cluster centers back to the original DataFrame
df_with_centers = df.join(broadcast(cluster_centers), on="cluster_id")

# Step 3: Calculate the Euclidean distance to the cluster center
df_with_distance = df_with_centers.withColumn(
    "distance",
    sqrt(
        pow(col("PCA_1") - col("avg_PCA_1"), 2) +
        pow(col("PCA_2") - col("avg_PCA_2"), 2) +
        pow(col("PCA_3") - col("avg_PCA_3"), 2)
    )
)

# Step 4: Find the closest molecule in each cluster
window_spec = Window.partitionBy("cluster_id").orderBy("distance")

df_with_rank = df_with_distance.withColumn(
    "rn",
    row_number().over(window_spec)
)

closest_molecules = df_with_rank.filter(col("rn") == 1).select(
    "smiles",
    "PCA_1",
    "PCA_2",
    "PCA_3",
    "cluster_id"
)

# Step 5: Create the new_cluster_dataframe
new_cluster_dataframe = closest_molecules

cluster_representatives = new_cluster_dataframe.toPandas()

cluster_representatives.to_csv('cluster_representatives.csv', index=False)