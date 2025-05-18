FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app

# Install RDKit and other dependencies
RUN conda install -y -c conda-forge rdkit python=3.9 \
    && pip install streamlit py3Dmol pillow

# Let Render know which port to use
ENV PORT=10000
EXPOSE 10000

# Start the Streamlit app
CMD streamlit run app.py --server.port=$PORT --server.address=0.0.0.0
