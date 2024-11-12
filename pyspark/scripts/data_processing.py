from pyspark.sql import SparkSession
from pyspark.ml.feature import Bucketizer
from pyspark.sql.functions import col, concat_ws
from pyspark.sql.window import Window
from pyspark.sql.functions import ntile

# Initialize Spark session with increased memory
spark = SparkSession.builder \
    .appName("HierarchicalBinning") \
    .config("spark.driver.memory", "64g") \
    .config("spark.executor.memory", "64g") \
    .getOrCreate()

# Read the Parquet files from the shared data directory
input_path = "/data/parquet_files/*.parquet"
output_path = "/shared_volume/clustering_result.parquet"

df = spark.read.parquet(input_path)

# Define percentiles for 25 bins
percentiles = [i / 25.0 for i in range(1, 25)]
pca1_thresholds = df.approxQuantile("PCA_1", percentiles, 0.01)
pca1_splits = [-float("inf")] + pca1_thresholds + [float("inf")]

# Bucketize PCA_1
bucketizer_pca1 = Bucketizer(splits=pca1_splits, inputCol="PCA_1", outputCol="bin_PCA1_temp")
df = bucketizer_pca1.transform(df).withColumn("bin_PCA1", col("bin_PCA1_temp").cast("integer")).drop("bin_PCA1_temp")

# Define windows for PCA_2 and PCA_3 binning
window_pca2 = Window.partitionBy("bin_PCA1").orderBy("PCA_2")
window_pca3 = Window.partitionBy("bin_PCA1", "bin_PCA2").orderBy("PCA_3")

# Assign bins for PCA_2 and PCA_3
df = df.withColumn("bin_PCA2", ntile(25).over(window_pca2) - 1)
df = df.withColumn("bin_PCA3", ntile(25).over(window_pca3) - 1)

# Create cluster identifier
df = df.withColumn("cluster_id", concat_ws("_", "bin_PCA1", "bin_PCA2", "bin_PCA3"))

# Select the relevant columns
df_final = df.select("smiles", "PCA_1", "PCA_2", "PCA_3", "cluster_id")

# Save the DataFrame as a Parquet file in the shared volume
df_final.write.mode("overwrite").parquet(output_path)

# Stop the Spark session
spark.stop()
