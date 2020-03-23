#!/bin/sh
cd Extreme_zero

# simulate data
R --vanilla --slave < simulate_data_extreme_zero.R

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_KX" and the posterior inference in the folder "Inference_KX" 
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 4 -K 3 -i 4000 -o 1000 -s 8592 -c 8
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 4 -K 3 -i 4000 -b 2000 -c 8

/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 4 -K 4 -i 4000 -o 1000 -s 9956 -c 8
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 4 -K 4 -i 4000 -b 2000 -c 8

/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 4 -K 5 -i 4000 -o 1000 -s 5493 -c 8
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 4 -K 5 -i 4000 -b 2000 -c 8

/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 4 -K 6 -i 4000 -o 1000 -s 9646 -c 8
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 4 -K 6 -i 4000 -b 2000 -c 8

/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 4 -K 7 -i 4000 -o 1000 -s 5657 -c 8
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 4 -K 7 -i 4000 -b 2000 -c 8

/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq -d./ -r./RawCountData/ -p simulation -v 4 -K 8 -i 4000 -o 1000 -s 8883 -c 8
/lustre/project/Stat/s1155082896/BUSseq/BUSseq_omp/BUSseq_inference -d./ -r./RawCountData/ -p simulation -v 4 -K 8 -i 4000 -b 2000 -c 8

# Summarize the results of BUSseq
R --vanilla --slave < summarize_BUSseq.R
