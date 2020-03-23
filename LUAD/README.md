# Generate the Figures for the Lung Adenocarcinoma (LUAD) study

This folder contains scripts that generate figures for the human pancreatic study in the manuscript. To reproduce all the results, please

1. Run `data_preprocessing` to conduct preprocessing on the LUAD dataset. 
   - conduct BUSseq and obtain the posterior sampling and inference (according to the BIC criterion, we select K = 8)
   - write out the list of identified intrinsic genes and the corrected read count matrix into `x_corrected.txt`;
   - implement batch effect correction to obtain the corrected read count matrix;
   - calculate ARI between the estimated cell types by BUSseq and true cell type labels;
   - draw the t-SNE and PCA plots of BUSseq correction as well as the expression levels of marker genes on the t-SNE plot and save them in the folder `BUSseq/Image`;
   - save the ARI and t-SNE coordinates as `BUSseq/BUSseq_results.RData`;

2. Run `Comparison/liger/run_liger.R` to apply [LIGER](https://github.com/MacoskoLab/liger) to the pancreatic dataset.
   - calculate ARI between the estimated cell types by LIGER and true cell type labels;
   - draw the t-SNE and PCA plots of LIGER correction and save them in the folder `Comparison/liger/Image`;
   - save the ARI and t-SNE coordinates as `Comparison/liger/liger_results.RData`;

3. Run `Comparison/MNN/run_MNN.R` to apply [MNN](https://github.com/MarioniLab/MNN2017) to the pancreatic dataset.
   - calculate ARI between the estimated cell types by MNN and true cell type labels;
   - draw the t-SNE and PCA plots of MNN correction and save them in the folder `Comparison/MNN/Image`;
   - save the ARI and t-SNE coordinates as `Comparison/MNN/MNN_results.RData`;

4. Run `Comparison/scanorama/run_scanorama.sh` to apply [Scanorama](https://github.com/brianhie/scanorama) to the pancreatic dataset.
   - calculate ARI between the estimated cell types by Scanorama and true cell type labels;
   - draw the t-SNE and PCA plots of Scanorama correction and save them in the folder `Comparison/scanorama/Image`;
   - save the ARI and t-SNE coordinates as `Comparison/scanorama/scanorama_result.RData`;

5. Run `Comparison/scVI/run_scVI.sh` to apply [scVI](https://github.com/YosefLab/scVI) to the pancreatic dataset.
   - calculate ARI between the estimated cell types by scVI and true cell type labels;
   - draw the t-SNE and PCA plots of scVI correction and save them in the folder `Comparison/scVI/Image`;
   - save the ARI and t-SNE coordinates as `Comparison/scVI/scVI_results.RData`;

6. Run `Comparison/Seurat/run_Seurat.R` to apply [Seurat](https://satijalab.org/seurat/) to the pancreatic dataset.
   - calculate ARI between the estimated cell types by Seurat and true cell type labels;
   - draw the t-SNE and PCA plots of Seurat correction and save them in the folder `Comparison/Seurat/Image`;
   - save the ARI and t-SNE coordinates as `Comparison/Seurat/Seurat_results.RData`;

7. Run `Comparison/ZINBWaVE/run_ZINBWaVE.R` to apply [ZINBWaVE](https://github.com/drisso/zinbwave) to the pancreatic dataset
   - calculate ARI between the estimated cell types by ZINBWaVE and true cell type labels;
   - draw the t-SNE and PCA plots of ZINBWaVE correction and save them in the folder `Comparison/ZINBWaVE/Image`;
   - save the ARI and t-SNE coordinates as `Comparison/ZINBWaVE/ZINBWaVE_results.RData`;

8. Run `collect_evaluation.R` to
   - collect the ARIs of all methods, calculate Silhouette coefficients and save them in the `Results` folder; 
   - generate t-SNE and PCA plots of uncorrected count matrix and save them in the `Image` folder.
   - draw the boxplot of Silhouette coefficients and save them in the `Image` folder.
