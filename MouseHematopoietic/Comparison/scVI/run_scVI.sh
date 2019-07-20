#!/bin/sh

# Conduct analysis by scanorama
cd /scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/scVI
# Activate conda environment
source activate tensorflow_env
# Python version 3.7.3
python scVI_log.py

# Summarize the results of scVI
R-351 --vanilla --slave < summarize_scVI.R 