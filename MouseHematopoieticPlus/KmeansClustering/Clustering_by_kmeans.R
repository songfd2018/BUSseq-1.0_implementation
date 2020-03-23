rm(list=ls())
library(mclust) # For ARI
library(cluster) 

####################################################################################
# Calculate the ARI by kmeans clustering with the knowledge of existing 7 clusters #
####################################################################################
rm(list=ls())
# Load the analysis of hematopotietic study
# Working directory
load("../comparison_workspace.RData")
ARI_values
set.seed(123)
ARI_kmeans <- rep(NA, 7)
# BUSseq: output is the corrected count data of the identified intrinsic genes 
kclust_BUSseq <- pam(t(log1p(x_corrected[D.est==1,])), k = 7)
w_BUSseq_kmeans <- kclust_BUSseq$clustering
ARI_kmeans[1] <- adjustedRandIndex(w_BUSseq_kmeans,metadata$CellType)

# LIGER: output is the low-dimensional shared factors
kclust_liger <- pam(liger_hemat@H.norm, k = 7)
w_liger_kmeans <- kclust_liger$clustering
ARI_kmeans[2] <- adjustedRandIndex(w_liger_kmeans,metadata$CellType)

# MNN: output is the corrected count data
kclust_MNN <- pam(t(omat), k = 7)
w_MNN_kmeans <- kclust_MNN$clustering
ARI_kmeans[3] <- adjustedRandIndex(w_MNN_kmeans,metadata$CellType)

# Scanorama : the batch-corrected dataset
kclust_Scanorama <- pam(scanorama_corrected, k = 7)
w_Scanorama_kmeans <- kclust_Scanorama$clustering
ARI_kmeans[4] <- adjustedRandIndex(w_Scanorama_kmeans,metadata$CellType)

# scVI : the batch-corrected dataset
kclust_scVI <- pam(scVI_corrected, k = 7)
w_scVI_kmeans <- kclust_scVI$clustering
ARI_kmeans[5] <- adjustedRandIndex(w_scVI_kmeans,metadata$CellType)

# Seurat : the batch-corrected dataset
kclust_Seurat <- pam(hemat.integrated@reductions$pca@cell.embeddings, k = 7)
w_Seurat_kmeans <- kclust_Seurat$clustering
ARI_kmeans[6] <- adjustedRandIndex(w_Seurat_kmeans,metadata$CellType)

# ZINBWaVE : the batch-corrected dataset
kclust_ZINBWaVE <- pam(zinb_batch@W, k = 7)
w_ZINBWaVE_kmeans <- kclust_ZINBWaVE$clustering
ARI_kmeans[7] <- adjustedRandIndex(w_ZINBWaVE_kmeans,metadata$CellType)

save(ARI_values, ARI_kmeans, file = "ARI_hemat.RData")