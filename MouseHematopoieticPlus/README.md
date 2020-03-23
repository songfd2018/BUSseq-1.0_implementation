# Generate the Figures of Hematopoietic Study.

This folder contains scripts that generate results for the hematopoietic study in the manuscript.

## Mimic the real data

BUSseq tries to mimics the data generating procedure of scRNA-seq experiments, we first evaluate how BUSseq recapitulates the properties of real scRNAseq data.

First, download the raw count data of two studies from [the OneDrive link](https://mycuhk-my.sharepoint.com/:u:/g/personal/1155082896_link_cuhk_edu_hk/EW-UIYqDLIRNk8DwWP823AUBKmeh_c9__Rs-7YrtOO34zA?e=PJlW1u) by password `BUSseq2019`. Next, please run `run_BUSseq.sh` to apply BUSseq on the hematopoietic dataset and analyze the posterior inference, including corrected-count-data generation.  

### Dropout rate

BUSseq can detect the existence of dropout events automatically. To calcluate dropout rates in the hematopoietic study, please run `DropoutRate\Dropout_rate_hemat.R` to output the batch-speicfic observed zero rates and estimated dropout rates.

### Posterior predictive check 

To compare BUSseq with its simplified form (without considering zero inflation), we conduct posterior prediction check in terms of zero rates to check which model mimics the observed scRNA-seq data better. More specifically, we run 8,000 iterations with the first 4,000 iterations as burn-ins in the MCMC algorithm, so we generated *J=8,000-4,000=4,000* replicated datasets. For each generated replicate dataset, we calculated the zero rates of each batch. Finally, we averaged the zero rates over all *J* iterations to calculate the posterior mean of the zero rate of each batch and compared it with the corresponding observed zero rate. To calculate the posterior mean of the zero rate of each batch, please run `PPC/PPC_hemat_original_zerorate.R` to generate the zero rates of the whole dataset and each batch for all *J* iterations in the txt file `PPC/ZeroRate_hemat_v1.txt`. For further comparison with no-zero-inflation model, please refer to the `README` file in the directory `../NewHemat`.

### Mean-variance trend

In order to explore how well a model recapitulates the properties of real scRNA-seq data, Zappia et al. [[1]](#1) suggested to simulate a scRNA-seq dataset and compare its properties with those of real data. However, it is inapproriate to conduct this comparsion when batch effects exist. Without loss of generality, here we draw the mean-variance trends for the second batch *b=2* of the hematopoietic study as illustrations. Please run 

## Covergence of MCMC alogrithm

### Acceptance rate

### EPSR

## More Comparison

### Unfied clustering method

### Silhouette on PCA

## Others

### Slingshot

### Reproduce results of MNN

## References
<a id="1">[1]</a> 
Zappia, L., Phipson, B. and Oshlack, A., 2017. Splatter: simulation of single-cell RNA sequencing data. Genome biology, 18(1), p.174.