#!/bin/sh

# Run BUSseq to generate the posterior sampling in the folder "MCMC_sampling_K6" and the posterior inference in the folder "Inference_K6"
cd /your/working/directory/BUSseq_implementation-1.0/MouseHematopoietic/
../../BUSseq-1.0/BUSseq -d./BUSseq/ -r./RawCountData/ -p hemat -v 1 -K 6 -i 8000 -o 2000 -s 4367 -c 8
../../BUSseq-1.0/BUSseq_inference -d./BUSseq/ -r./RawCountData/ -p hemat -v 1 -K 6 -i 8000 -b 4000 -c 8

# Summarize the results of BUSseq
cd /your/working/directory/BUSseq_implementation-1.0/MouseHematopoietic/BUSseq
R --vanilla --slave < summarize_BUSseq.R