from pyspark.sql import SparkSession
from pyspark.ml.feature import Bucketizer
from pyspark.sql.functions import col, concat_ws
from pyspark.sql.window import Window
from pyspark.sql.functions import ntile

# optimized_clustering_script.py

from pyspark.sql import SparkSession
from pyspark.sql.functions import col, ntile, concat_ws
from pyspark.sql.window import Window
from pyspark.ml.feature import QuantileDiscretizer
from pyspark.sql.functions import col, ntile, concat_ws
from pyspark.sql.types import IntegerType

from pyspark.sql import SparkSession
from pyspark.ml.feature import QuantileDiscretizer
from pyspark.sql.functions import col, concat_ws, ntile
from pyspark.sql.window import Window

def main():
    # Initialize SparkSession with the new temporary directory
    spark = SparkSession.builder \
        .appName("GridClusteringJob") \
        .master("local[8]") \
        .config("spark.driver.memory", "56g") \
        .config("spark.executor.memory", "56g") \
        .getOrCreate()

    df = spark.read.parquet("/mnt/samsung_2tb/mixed_data/results/test/*.parquet")
    df.count()

    # Step 1: Assign Clusters Based on PCA_1
    discretizer_pca1 = QuantileDiscretizer(
        numBuckets=50, inputCol="PCA_1", outputCol="cluster_pca1", relativeError=0.01
    )
    df = discretizer_pca1.fit(df).transform(df)
    df = df.withColumn("cluster_pca1", col("cluster_pca1").cast(IntegerType()))

    # Step 2: Assign Clusters Based on PCA_2 within PCA_1 Clusters
    window_pca2 = Window.partitionBy("cluster_pca1").orderBy("PCA_2")
    df = df.withColumn("cluster_pca2", ntile(50).over(window_pca2))
    df = df.withColumn("cluster_pca2", col("cluster_pca2").cast(IntegerType()))

    # Step 3: Assign Clusters Based on PCA_3 within PCA_1 and PCA_2 Clusters
    window_pca3 = Window.partitionBy("cluster_pca1", "cluster_pca2").orderBy("PCA_3")
    df = df.withColumn("cluster_pca3", ntile(50).over(window_pca3))
    df = df.withColumn("cluster_pca3", col("cluster_pca3").cast(IntegerType()))

    # Step 4: Create Final Cluster ID and Select Relevant Columns
    df = df.withColumn(
        "cluster_id",
        concat_ws(
            "_",
            col("cluster_pca1"),
            col("cluster_pca2"),
            col("cluster_pca3")
        )
    )

    df_result = df.select("smiles", "PCA_1", "PCA_2", "PCA_3", "cluster_id")

     # Save the Results into 10 Parquet Files
    df_result.coalesce(10).write.mode("overwrite").parquet("path_to_output_directory")

  # Stop Spark Session
    spark.stop()

if __name__ == "__main__":  
    main()