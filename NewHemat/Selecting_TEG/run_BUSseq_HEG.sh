#!/bin/sh

cd /scratch/data01/BUSseq_cpp/Revision/GitHub_sub/NewHemat/Selecting_TEG

# filter out the common TEG
R --vanilla --slave < Hemat_TEG_filtering.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K5" and the posterior inference in the folder "Inference_K5" 
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./  -r./RawCountData/ -p hemat -v 3 -K 6 -i 8000 -o 2000 -s 372 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./  -r./RawCountData/ -p hemat -v 3 -K 6 -i 8000 -b 4000 -c 8

# compare the clustering results
R-351 --vanilla --slave < Gene_filter_comparison_hemat.R