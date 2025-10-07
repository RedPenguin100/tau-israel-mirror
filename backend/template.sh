#!/bin/bash
#PBS -l select=2:ncpus=16:mem=16gb
# for running qsub -q tamirQ (or other queue) "$job_file"
hostname

module load miniconda/miniconda3-2023-environmentally
conda activate /tamir2/dinsaadon/miniconda3/envs/aso_design

cd /tamir2/dinsaadon/backend_igem
python aso_generator.py 
