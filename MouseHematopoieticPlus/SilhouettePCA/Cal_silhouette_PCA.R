##########################################################################
# Silhouette coefficient of pca by FACS cell type labels of each method #
##########################################################################
library(cluster)
library(reshape2)
setwd("SilhouettePCA")

# load the workspace of BUSseq
load("../BUSseq_workspace.RData")
# Storing the Silhouette coefficients
sil_pca_true <- matrix(NA,N,7)
pc_rank <- 10
colnames(sil_pca_true) <- c("BUSseq", "LIGER", "MNN", "Scanorama", "scVI", "Seurat", "ZINBWaVE")

# 10-dimension PCs
PCs_BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=pc_rank)$x
BUSseq_dist_pca <- dist(PCs_BUSseq)
sil_pca_true[,1] <- silhouette(as.integer(metadata$CellType), dist = BUSseq_dist_pca)[,3]

# load the workspace of liger
load("Workspace/liger_workspace.RData")

PCs_liger <- prcomp(loading_liger, rank=pc_rank)$x
liger_dist_pca <- dist(PCs_liger)
sil_pca_true[,2] <- silhouette(as.integer(metadata$CellType), dist = liger_dist_pca)[,3]

# load the workspace of MNN
load("Workspace/MNN_workspace.RData")

PCs_MNN <- prcomp(t.mnn, rank=pc_rank)$x
MNN_dist_pca <- dist(PCs_MNN)
sil_pca_true[,3] <- silhouette(as.integer(metadata$CellType), dist = MNN_dist_pca)[,3]

# load the workspace of scanorama
load("Workspace/scanorama_workspace.RData")

PCs_scanorama <- prcomp(scanorama_embedding, rank=pc_rank)$x
scanorama_dist_pca <- dist(PCs_scanorama)
sil_pca_true[,4] <- silhouette(as.integer(metadata$CellType), dist = scanorama_dist_pca)[,3]

# load the workspace of scVI
load("Workspace/scVI_workspace.RData")

PCs_scVI <- prcomp(scVI_corrected, rank=pc_rank)$x
scVI_dist_pca <- dist(PCs_scVI)
sil_pca_true[,5] <- silhouette(as.integer(metadata$CellType), dist = scVI_dist_pca)[,3]

# load the workspace of seurat
load("Workspace/Seurat_workspace.RData")

PCs_Seurat <- Seurat_PCA[,1:10]
seurat_dist_pca <- dist(PCs_Seurat)
sil_pca_true[,6] <- silhouette(as.integer(metadata$CellType), dist = seurat_dist_pca)[,3]

# load the workspace of ZINBWaVE
load("Workspace/ZINBWaVE_workspace.RData")

PCs_ZINBWaVE <- prcomp(zinb_batch@W, rank=pc_rank)$x
ZINBWaVE_dist_pca <- dist(PCs_ZINBWaVE)
sil_pca_true[,7] <- silhouette(as.integer(metadata$CellType), dist = ZINBWaVE_dist_pca)[,3]

# Drawing the boxplot plot
colnames(sil_pca_true) <- c("BUSseq", "LIGER", "MNN", "Seurat", "Scanorama", "scVI", "ZINB-WaVE")
sil_com_pca_true<-melt(sil_pca_true)
colnames(sil_com_pca_true) <- c("No_Cell","Method","Silhouette_Coefficient")

jpeg("Image/boxplot_of_hemat_Silhouette_coefs.jpg",width = 900, height = 600,quality = 100)
p <- ggplot(sil_com_pca_true, aes(x=Method, y=Silhouette_Coefficient,fill=Method)) +
  geom_boxplot() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=48),
        axis.text.x =element_text(face = "bold", #color = "#993333", 
                                  size = 36,angle = 90), 
        axis.text.y = element_text(face = "bold", #color = "#993333", 
                                   size = 36),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

save.image("hemat_silhouette_pca.RData")

