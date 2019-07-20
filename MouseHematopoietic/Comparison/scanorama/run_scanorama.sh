#!/bin/sh

# Conduct analysis by scanorama
cd /scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/scanorama
# Activate conda env
source activate tensorflow_gpuenv
# Python version 3.7.3
python scanorama_log.py

# Summarize the results of scanorama
R-351 --vanilla --slave < summarize_scanorama.R 