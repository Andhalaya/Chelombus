# Use the official PySpark image or a Python image with Spark installed
FROM bitnami/spark:latest

# Set working directory
WORKDIR /app

# Copy requirements (if any)
# COPY requirements.txt .
# RUN pip install --no-cache-dir -r requirements.txt

# Copy your PySpark scripts
COPY scripts/ /app/scripts/
