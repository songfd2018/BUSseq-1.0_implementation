# Generate the Figures of Mouse Hematopoietic Study

This folder contains scripts that generate results for the mouse hematopoietic study in the manuscript. To reproduce all the results, please

1. Run `BUSseq/run_BUSseq.sh` to apply BUSseq to the hematopoietic dataset. 
   - conduct BUSseq and obatin the posterior sampling and infernece (we select K = 8, according the BIC criterion)
   - write out the list of identified intrinsic genes
   - implement batch effect correction to obtain the corrected read count matrix;
   - calculate ARI between the estimated cell types by BUSseq and true cell type labels;
   - draw the t-SNE and PCA plots of BUSseq correction and save them in the folder `BUSseq/Image`;
   - save the workerspace as `BUSseq/BUSseq_workspace.RData`;

2. Run `Comparison/liger/run_liger.R` to apply [LIGER](https://github.com/MacoskoLab/liger) to the hematopoietic dataset.
   - calculate ARI between the estimated cell types by liger and true cell type labels;
   - draw the t-SNE and PCA plots of liger correction and save them in the folder `Comparison/liger/Image`;
   - save the workerspace as `Comparison/liger/liger_workspace.RData`;

3. Run `Comparison/MNN/run_MNN.R` to apply [MNN](https://github.com/MarioniLab/MNN2017) to the hematopoietic dataset.
   - calculate ARI between the estimated cell types by MNN and true cell type labels;
   - draw the t-SNE and PCA plots of MNN correction and save them in the folder `Comparison/MNN/Image`;
   - save the workerspace as `Comparison/MNN/MNN_workspace.RData`;

4. Run `Comparison/scanorama/run_scanorama.sh` to apply [Scanorama](https://github.com/brianhie/scanorama) to the hematopoietic dataset.
   - calculate ARI between the estimated cell types by scanorama and true cell type labels;
   - draw the t-SNE and PCA plots of scanorama correction and save them in the folder `Comparison/scanorama/Image`;
   - save the workerspace as `Comparison/scanorama/scanorama_workspace.RData`;

5. Run `Comparison/scVI/run_scVI.sh` to apply [scVI](https://github.com/YosefLab/scVI) to the hematopoietic dataset.
   - calculate ARI between the estimated cell types by scVI and true cell type labels;
   - draw the t-SNE and PCA plots of scVI correction and save them in the folder `Comparison/scVI/Image`;
   - save the workerspace as `Comparison/scVI/scVI_workspace.RData`;

6. Run `Comparison/Seurat/run_Seurat.R` to apply [Seurat](https://satijalab.org/seurat/) to the hematopoietic dataset.
   - calculate ARI between the estimated cell types by Seurat and true cell type labels;
   - draw the t-SNE and PCA plots of Seurat correction and save them in the folder `Comparison/Seurat/Image`;
   - save the workerspace as `Comparison/Seurat/Seurat_workspace.RData`;

7. Run `Comparison/ZINBWaVE/run_ZINBWaVE.R` to apply [ZINBWaVE](https://github.com/drisso/zinbwave) to the hematopoietic dataset
   - calculate ARI between the estimated cell types by ZINBWaVE and true cell type labels;
   - draw the t-SNE and PCA plots of ZINBWaVE correction and save them in the folder `Comparison/ZINBWaVE/Image`;
   - save the workerspace as `Comparison/ZINBWaVE/ZINBWaVE_workspace.RData`;

8. Run `collect_evaluation.R` to
   - collect the ARIs of all methods, calculate Silhouette coefficients and save them in the `./Results` folder; 
   - generate t-SNE and PCA plots of uncorrected count matrix and save them in the `./Image/` folder.
   - draw the boxplot of Silhouette coefficients and save it in the `./Image/` folder.