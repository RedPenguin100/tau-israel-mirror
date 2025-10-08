FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy environment file and create conda env
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Copy and install requirements
COPY requirements.txt .
RUN /bin/bash -c "source activate aso_design && pip install -r requirements.txt"

# Copy and install aso_gen package
COPY aso_gen ./aso_gen
RUN /bin/bash -c "source activate aso_design && cd aso_gen && pip install . && cd .."

# Create directories for downloads
RUN mkdir -p /app/models /app/data

# Copy BOTH download scripts
COPY download_files.sh .
COPY download_files.py .
RUN chmod +x download_files.sh

# ===== DOWNLOAD FILES USING BOTH SCRIPTS =====
# Shell script runs first
RUN bash download_files.sh

# Then Python script runs
RUN /bin/bash -c "source activate aso_design && python download_files.py"

# Copy the rest of your application code
COPY . .

# Expose port
EXPOSE 8000

# Start command (files already downloaded at build time)
CMD ["/bin/bash", "-c", "source activate aso_design && uvicorn main:app --host 0.0.0.0 --port ${PORT:-8000}"]