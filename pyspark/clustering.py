from pyspark.sql import SparkSession
from pyspark.sql.functions import col, concat_ws, ntile, lit
from pyspark.sql.window import Window
import os

# Initialize SparkSession with optimized configurations
spark = SparkSession.builder \
    .appName("GridClusteringJob") \
    .master("local[8]") \
    .config("spark.driver.memory", "50g") \
    .config("spark.executor.memory", "50g") \
    .config("spark.sql.shuffle.partitions", "2000") \
    .config("spark.executor.extraJavaOptions", "-XX:+UseG1GC") \
    .config("spark.local.dir", "/mnt/10tb_hdd/pyspark_tmp") \
    .getOrCreate()

input_dir = "/mnt/10tb_hdd/results/" # Here's where the first 50 buckets are 
output_dir = "/mnt/samsung_2tb/mixed_data/results/clickhouse/"

# Loop over all bucket files
for bucket_num in range(50):
    # Construct the input file path
    input_file = os.path.join(input_dir, f"bucket{bucket_num}.parquet")
    
    # Read the Parquet file
    df = spark.read.parquet(input_file)
    
    # Repartition the DataFrame for better parallelism
    df = df.repartition(2000)
    
    # First Level of Clustering on PCA_2
    window_pca2 = Window.orderBy("PCA_2")
    
    # Assign bins for PCA_2
    df = df.withColumn("bin_PCA2", ntile(50).over(window_pca2) - 1)
    
    # Repartition based on 'bin_PCA2' before the next window function
    df = df.repartition("bin_PCA2")
    
    # Second Level of Clustering on PCA_3
    window_pca3 = Window.partitionBy("bin_PCA2").orderBy("PCA_3")
    
    # Assign bins for PCA_3
    df = df.withColumn("bin_PCA3", ntile(50).over(window_pca3) - 1)
    
    # Convert bin_PCA2 and bin_PCA3 to integer type
    df = df.withColumn("bin_PCA2", col("bin_PCA2").cast("int"))
    df = df.withColumn("bin_PCA3", col("bin_PCA3").cast("int"))
    
    # Create cluster_id with the prefix based on the bucket number
    df = df.withColumn("cluster_id", concat_ws("_", lit(str(bucket_num)), col("bin_PCA2"), col("bin_PCA3")))
    
    # Drop intermediate bin columns
    df = df.drop("bin_PCA2", "bin_PCA3")
    
    # Define the output path for the current bucket
    # Option 1: Save each bucket's results in a separate subdirectory
    bucket_output_dir = os.path.join(output_dir, f"bucket{bucket_num}")
    
    # Write the final DataFrame to Parquet files with Snappy compression
    df.write.mode("overwrite").option('compression', 'snappy').parquet(bucket_output_dir)
    
    print(f"Finished processing bucket{bucket_num}.parquet")

