#!/bin/sh
cd /scratch/data01/BUSseq_cpp/Revision/GitHub_sub/LUAD

R-351 --vanilla --slave < data_preprocessing.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX"
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 2 -i 8000 -o 2000 -s 4629 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 2 -i 8000 -b 4000 -c 8

/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 3 -i 8000 -o 2000 -s 2548 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 3 -i 8000 -b 4000 -c 8

/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 4 -i 8000 -o 2000 -s 8880 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 4 -i 8000 -b 4000 -c 8

/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 5 -i 8000 -o 2000 -s 924 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 5 -i 8000 -b 4000 -c 8

/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p LUAD -v 1 -K 6 -i 8000 -o 2000 -s 798 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p LUAD -v 1 -K 6 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
R-351 --vanilla --slave < summarize_BUSseq.R
