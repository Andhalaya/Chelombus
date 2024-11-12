import sys
import os 
import logging
from pyspark.sql import SparkSession
from pyspark.ml.feature import Bucketizer
from pyspark.sql.functions import col
from pyspark.sql.window import Window
from pyspark.sql.functions import ntile
from pyspark.sql.functions import avg, sqrt, pow, col, row_number, broadcast

#TODO: Complete pipeline for constructing the Apache Dataframe

class ClusterMethod():
    def __init__(self, bin_size: list, output_path: str):
        self.spark = SparkSession.builder \
            .appName("HierarchicalBinning") \
            .config("spark.driver.memory", "64g") \
            .config("spark.executor.memory", "64g") \
            .getOrCreate()
        self.output_path = output_path 

        self.bin_size = bin_size

    def _get_count(self):
        """
        Get number of samples in the Apache Spark dataframe
        """
        pass

    def spark_constructor(self, file_path):
        """
        Method to create a Apache Spark database to perform the clustering of large datastets (i.e. that cannot be allocated entirely in memory)
        param: file_path: the path to the folder where all the parquet files are stored
        """
        database= self.spark.read.parquet(f"{file_path}/*")
        print(database.count())


    def create_clusters(self, dataframe):
        """
        Method to create the clusters from Spark dataframe. The idea is to cluster by sorting the first PCA dimension, then creating n-size bins on that dimension
        And then on each bin, sort by the second PCA dimension, create n-size bins. And so on for every dimension we have.
        At the end we will have n**(number of PCA dimensions) equally sized clusters.
        """

        # Define the percentiles for 25 bins
        percentiles_pca_1= [i / self.bin_size[0] for i in range(1, self.bin_size[0])]

        # Compute approximate quantiles for PCA_1
        pca1_thresholds = dataframe.approxQuantile("PCA_1", percentiles_pca_1, 0.01)

        # Create splits for Bucketizer
        pca1_splits = [-float("inf")] + pca1_thresholds + [float("inf")]

        # Initialize Bucketizer
        bucketizer_pca1 = Bucketizer(
            splits=pca1_splits,
            inputCol="PCA_1",
            outputCol="bin_PCA1_temp"
        )

        # Transform the DataFrame
        dataframe = bucketizer_pca1.transform(dataframe)

        # Convert bin_PCA1_temp to IntegerType and adjust to start from 0
        dataframe = dataframe.withColumn("bin_PCA1", col("bin_PCA1_temp").cast("integer"))

        # Drop the temporary column
        dataframe = dataframe.drop("bin_PCA1_temp")

        # Create a window specification for PCA_2 within each bin_PCA1
        window_pca2 = Window.partitionBy("bin_PCA1").orderBy("PCA_2")

        # Assign bins for PCA_2, starting from 0
        dataframe = dataframe.withColumn("bin_PCA2", ntile(self.bin_size[1]).over(window_pca2) - 1)

        # Create a window specification for PCA_3 within each bin_PCA1 and bin_PCA2
        window_pca3 = Window.partitionBy("bin_PCA1", "bin_PCA2").orderBy("PCA_3")

        # Assign bins for PCA_3, starting from 0
        dataframe = dataframe.withColumn("bin_PCA3", ntile(self.bin_size[2]).over(window_pca3) - 1)

        # dataframe.select("smiles", "PCA_1", "PCA_2", "PCA_3", "cluster_id").show()

        # return new dataframe
        return dataframe 

    def save_first_dimension_cluster(self, dataframe):
        """
        Save into a csv all the points that belongs to each PCA1 cluster. 
        e.g. cluster_10.csv will have all compounds with cluster_id == 19_*_*
        """
        for i in range(self.bin_size[0]):
            df_cluster1 = dataframe.filter(col("bin_PCA1") == i)
            df_cluster1_pd = df_cluster1.toPandas()
            df_cluster1_pd.to_csv(f"{self.output_path}/cluster_{i}.csv", index=False)

    def get_representatives(self, dataframe):
        """
        Method to get compounds that act as representatives of the cluster. 
        This will constitute the dataframe that is used for the first TMAP
        The representative is the molecule that is the closest to the average value for every PCA dimension
        """

        # Step 1: Compute average PCA values for each cluster
        cluster_centers = dataframe.groupBy("cluster_id").agg(
            avg("PCA_1").alias("avg_PCA_1"),
            avg("PCA_2").alias("avg_PCA_2"),
            avg("PCA_3").alias("avg_PCA_3")
        ).cache()

        # Step 2: Join the cluster centers back to the original DataFrame
        df_with_centers = dataframe.join(broadcast(cluster_centers), on="cluster_id")

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

        cluster_representatives.to_csv(f'{self.output_path}/cluster_representatives.csv', index=False)