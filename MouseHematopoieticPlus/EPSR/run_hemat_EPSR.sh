#!/bin/sh
cd EPSR

# run another BUSseq MCMC chain
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r../RawCountData/ -p hemat -v 1 -K 6 -i 8000 -o 2000 -s 372 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p hemat -v 1 -K 6 -i 8000 -b 4000 -c 8

# Calculate EPSR
R-351 --vanilla --slave < EPSR_analysis.R
