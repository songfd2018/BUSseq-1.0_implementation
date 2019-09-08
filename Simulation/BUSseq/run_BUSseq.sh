#!/bin/sh

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX" 
cd /your/working/directory/BUSseq_implementation-1.0/Simulation/
# K = 3
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 3 -i 8000 -o 2000 -s 6007 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 3 -i 8000 -b 4000 -c 8
# K = 4
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 4 -i 8000 -o 2000 -s 2914 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 4 -i 8000 -b 4000 -c 8
# K = 5
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 5 -i 8000 -o 2000 -s 7787 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 5 -i 8000 -b 4000 -c 8
# K = 6
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 6 -i 8000 -o 2000 -s 1710 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 6 -i 8000 -b 4000 -c 8
# K = 7
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 7 -i 8000 -o 2000 -s 1374 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 7 -i 8000 -b 4000 -c 8
# K = 8
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 8 -i 8000 -o 2000 -s 1841 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 8 -i 8000 -b 4000 -c 8
# K = 9
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 9 -i 8000 -o 2000 -s 4408 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 9 -i 8000 -b 4000 -c 8
# K = 10
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p simulation -v 4 -K 10 -i 8000 -o 2000 -s 6796 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p pancreas -v 4 -K 10 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
R --vanilla --slave < ./BUSseq/summarize_BUSseq.R