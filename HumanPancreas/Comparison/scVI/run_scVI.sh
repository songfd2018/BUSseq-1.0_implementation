#!/bin/sh

# Conduct analysis by scVI
# Activate conda environment
source activate tensorflow_env
# Python version 3.7.3
python scVI_log.py

# Summarize the results
R --vanilla --slave < summarize_scVI.R
