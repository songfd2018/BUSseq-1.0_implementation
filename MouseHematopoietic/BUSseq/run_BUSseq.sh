#!/bin/sh
cd /scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K6" and the posterior inference in the folder "Inference_K6" 
../../BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p hemat -v 1 -K 6 -i 8000 -o 2000 -s 4367 -c 8
../../BUSseq_omp/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p hemat -v 1 -K 6 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
cd /scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/BUSseq
R-351 --vanilla --slave < summarize_BUSseq.R 