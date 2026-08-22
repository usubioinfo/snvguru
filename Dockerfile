FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Copy dependency files first to cache dependencies layer
COPY environment.yml .

# Create the Conda/Mamba environment and clean up cache to save space
RUN mamba env create -f environment.yml && mamba clean -afy

# Set PATH to use the conda environment's binaries by default
ENV PATH /opt/conda/envs/snvguru/bin:$PATH

# Copy package source files
COPY . .

# Install the SNVGuru Python package
RUN pip install .

# Define the entrypoint script
ENTRYPOINT ["snvguru"]
