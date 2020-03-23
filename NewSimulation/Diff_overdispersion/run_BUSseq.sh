#!/bin/sh

cd Diff_overdispersion
# simulate data
R --vanilla --slave < simulate_data_Diff_overdispersion.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K5" and the posterior inference in the folder "Inference_K5" 
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 2 -K 5 -i 4000 -o 1000 -s 9123 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 2 -K 5 -i 4000 -b 2000 -c 8

# Summarize the results of BUSseq
R --vanilla --slave < summarize_BUSseq.R
