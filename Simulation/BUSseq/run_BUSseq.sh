#!/bin/sh

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX", where X = 3 - 10, and the posterior inference in the folder "Inference_KX" 
cd /your/working/directory/BUSseq_implementation-1.0/Simulation/

# K = 5
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 4475 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 3
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 7970 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 4
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 4870 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 6
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 579 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 7
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 3757 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 8
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 5340 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 9
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 6794 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# K = 10
the/source/code/directory/BUSseq_omp/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -o 2000 -s 184 -c 8
the/source/code/directory/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p simulation -v 1 -K 3 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
cd /your/working/directory/BUSseq_implementation-1.0/Simulation/BUSseq
R --vanilla --slave < summarize_BUSseq.R