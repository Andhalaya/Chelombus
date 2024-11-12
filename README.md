docker-compose.yml: Defines and runs multi-container Docker applications.

frontend/: Contains files related to the frontend service.

pyspark/: Contains the PySpark data processing code and Dockerfile.

clickhouse/: Contains configurations for the ClickHouse database (optional if using default settings).

chelombus/: Contains existing GitHub repository code, related to generating the clustering, including its own README.md, LICENSE, and requirements.txt.

shared_volume/: A directory used for sharing generated HTML files between services.

data/: A directory for data files or scripts, These are point folders to my local folders. 

.github/workflows/: Contains CI/CD pipeline configurations for GitHub Actions.