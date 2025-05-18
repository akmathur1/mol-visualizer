FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy your code
COPY . /app

# Create environment
RUN conda install -y -c conda-forge rdkit python=3.9 \
    && pip install streamlit py3Dmol pillow

# Expose the default Streamlit port
EXPOSE 10000

# Run the Streamlit app
CMD ["streamlit", "run", "app.py", "--server.port=10000", "--server.address=0.0.0.0"]
