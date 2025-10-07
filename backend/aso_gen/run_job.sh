#!/bin/bash
#PBS -l select=1:ncpus=2:mem=4gb
hostname

module load miniconda/miniconda3-2023-environmentally
conda activate /tamir2/dinsaadon/miniconda3/envs/aso_design

cd /tamir2/dinsaadon/aso_gen/ASOdesign
python ./run_pipe.py 
