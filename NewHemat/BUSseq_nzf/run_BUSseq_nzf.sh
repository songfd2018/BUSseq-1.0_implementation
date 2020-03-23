#!/bin/sh

# compile the source code
make

# run BUSseq with no zero inflation 
cd /scratch/data01/BUSseq_cpp/Revision/GitHub_sub/NewHemat/BUSseq_nzf/
BUSseq_nd -d./  -r./RawCountData/ -p hemat -v 2 -K 6 -i 8000 -o 2000 -s 6232 -c 8
BUSseq_inference_nd -d./  -r./RawCountData/ -p hemat -v 2 -K 6 -i 8000 -b 4000 -c 8
