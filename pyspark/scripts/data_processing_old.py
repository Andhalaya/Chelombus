# File to do the clustering in PySpark instead of using the clustering pipeline. This will feed the Clickhouse database. 


from pyspark.sql import SparkSession
from pyspark.ml.feature import Bucketizer
from pyspark.sql.functions import col



# Increase memory allocation
spark = SparkSession.builder \
    .appName("HierarchicalBinning") \
    .config("spark.driver.memory", "64g") \
    .config("spark.executor.memory", "64g") \
    .getOrCreate()



df = spark.read.parquet("/home/afloresep/work/chelombus/data/230M/output/*")
df.count()


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

from pyspark.sql.functions import concat_ws

# Combine bin columns to form a cluster identifier
df = df.withColumn("cluster_id", concat_ws("_", "bin_PCA1", "bin_PCA2", "bin_PCA3"))

df.select("smiles", "PCA_1", "PCA_2", "PCA_3", "cluster_id").show()