# Use Python 3.7.12 slim as the base image
FROM python:3.7.12-slim
 

# Install Miniconda and build tools
RUN apt-get update && apt-get install -y wget bzip2 build-essential gcc && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Add Conda to the PATH
ENV PATH="/opt/conda/bin:$PATH"

# Copy the environment.yml file and create the Conda environment
COPY chelombus.yml /app/chelombus.yml
RUN conda env create -f /app/chelombus.yml -p /opt/conda/envs/chelombus

# Activate the environment using conda run
SHELL ["conda", "run", "-p", "/opt/conda/envs/chelombus", "/bin/bash", "-c"]

# Set the working directory
WORKDIR /app

# Copy the rest of the application files
COPY . /app

# Set the default command
CMD ["conda", "run", "-p", "/opt/conda/envs/chelombus", "bash"]
