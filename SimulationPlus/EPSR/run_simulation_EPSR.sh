#!/bin/sh
cd /scratch/data01/BUSseq_cpp/Revision/GitHub_sub/SimulationPlus/EPSR/

# run another BUSseq MCMC chain
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq -d./ -r../RawCountData/ -p simulation -v 1 -K 5 -i 4000 -o 1000 -s 2187 -c 8
/scratch/data01/BUSseq_cpp/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 1 -K 5 -i 4000 -b 2000 -c 8

# Calculate EPSR
R-351 --vanilla --slave < EPSR_analysis.R
