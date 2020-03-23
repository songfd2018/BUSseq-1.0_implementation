#!/bin/sh
cd /scratch/data01/BUSseq_cpp/Revision/GitHub_sub/SimulationPlus

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K8" and the posterior inference in the folder "Inference_K6"
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 1 -K 5 -i 8000 -o 2000 -s 9123 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 5 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
R-351 --vanilla --slave < summarize_BUSseq.R
