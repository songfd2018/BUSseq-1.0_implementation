#Apply BUSseq to the pancreatic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(12345)

########################
# Load Pancreatic Data #
########################
# Working directory
# setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/HumanPancreas/BUSseq")
setwd("/your/working/directory/BUSseq_implementation-1.0/HumanPancreas/BUSseq/")

# Loading pancreatic count data
y_obs <- read.table("../RawCountData/count_data_pancreas_v1.txt")

# Load dimension
dim <- read.table("../RawCountData/dim_pancreas_v1.txt")
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table("../RawCountData/metadata_pancreas_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../RawCountData/gene_list_pancreas_v1.txt",stringsAsFactors = F))

##########################
#load posterior inference#
##########################
K <- 8

# load w_est
w.est <- read.table("Inference_K8/w_est.txt")
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table("Inference_K8/alpha_est.txt")
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table("Inference_K8/beta_est.txt")
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table("Inference_K8/nu_est.txt")
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table("Inference_K8/delta_est.txt")
delta.est <- unlist(delta.est)
plot(delta.est,col=rep(1:B,nb) + 1)


# load gamma_est
gamma.est <- read.table("Inference_K8/gamma_est.txt")
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table("Inference_K8/phi_est.txt")
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table("Inference_K8/pi_est.txt")
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table("Inference_K8/p_est.txt")

# load tau0_est
tau0.est <- read.table("Inference_K8/tau0_est.txt")

# load PPI_est
PPI.est <- read.table("Inference_K8/PPI_est.txt")
D.est <- unlist(read.table("Inference_K8/IG_est.txt"))

# load X_imputed
x_imputed <- read.table("Inference_K8/x_imputed.txt")


####################
# Pathway analysis #
####################
intri_gene <- gene_list[D.est==1]
write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

# Pathway_analysis <- read.table("pathway_analysis_by_DAVID.txt",sep = "\t",stringsAsFactors = FALSE)
# Col_name <- Pathway_analysis[1,]
# colnames(Pathway_analysis) <- Col_name
# Pathway_analysis <- Pathway_analysis[-1,]
# 
# write.csv(Pathway_analysis,file="Pathway_analysis.csv",row.names = FALSE)


############################
# Batch effects correction #
############################
adjusted_values<-function(Truereads, .B, .nb, .N, .G, .K,
                          .alpha,.beta,.nu,.delta,.phi,.w){
  read_corrected <- Truereads
  .w <- .w + 1
  for(b in 1:.B){
    for(i in 1:.nb[b]){
      if(b==1){
        row.index <- i
      }else{
        row.index <- i + sum(.nb[1:(b-1)])
      }
      for(j in 1:.G){
        p_x<-pnbinom(Truereads[j,row.index],size = .phi[j,b], mu = exp(.alpha[j] + .beta[j,.w[row.index]] + .nu[j,b] + .delta[row.index]))
        p_xminus1<-pnbinom(Truereads[j,row.index]-1,size = .phi[j,b], mu = exp(.alpha[j] + .beta[j,.w[row.index]] + .nu[j,b] + .delta[row.index]))
        u <- runif(1,min=p_xminus1,max=p_x)
        u <- min(u,0.999)
        read_corrected[j,row.index] <- qnbinom(u,size =.phi[j,1], mu = exp(.alpha[j] + .beta[j,.w[row.index]]))
      }
      print(paste("Finish the ",row.index,"-th cell",sep=""))
    }
  }
  return(read_corrected)
}

start_time<-Sys.time()
print("Calculate corrected read counts:")
x_corrected<-adjusted_values(x_imputed, B, nb, N, G, K,
                             alpha.est,beta.est,nu.est,delta.est,phi.est,w_BUSseq)
write.table(x_corrected, file = "x_corrected.txt", row.names = FALSE, col.names = FALSE)
end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)

####################
# Pathway analysis #
####################
intri_gene <- gene_list[D.est==1]
write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)


#######
# ARI #
#######
ARI_BUSseq <- adjustedRandIndex(metadata$CellType,w_BUSseq)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
allcolors<-labels2colors(metadata$CellType)#[allsamples])

allcolors[allcolors=="red"]<-"deeppink"
allcolors[allcolors=="yellow"]<-"orange1"#"darkgoldenrod1"

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=20, 
       col=alpha(allcolors[subset],0.6),
       xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

colors4<-brewer.pal(4,"Set3")
batch.cols<-colors4[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=20, 
       col=alpha(batch.cols[subset],0.6),
       xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- "Image/tsne_pancreas_BUSseq_by_celltype.jpeg"
plot_by_celltype(BUSseq_by_celltype, tsne_BUSseq_dist$Y)

BUSseq_by_batch <- "Image/tsne_pancreas_BUSseq_by_batch.jpeg"
plot_by_batch(BUSseq_by_batch, tsne_BUSseq_dist$Y)

##########################################
# Draw the PCA plot on common cell types #
##########################################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)
BUSseq_PCA<- "Image/pca_pancreas_BUSseq.jpeg"
plot_by_celltype(BUSseq_PCA, pca.BUSseq$x, xlab = "PC 1", ylab = "PC 2")

# legend
pdf(file="Image/legend.pdf", width=10, height=8)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=c("GSE81076","GSE85241","GSE86473","E-MTAB-5061"), pch=20, cex=2.5, col=alpha(colors4,0.6), title = "Batch", bty="n")
legend(x=-0, y=4, legend=names(table(metadata$CellType)), pch=20, cex=2.5, col=alpha(c("turquoise","blue","brown","orange1","green","deeppink","black"),0.6), title = "Cell Type", bty="n")
dev.off()

##############################################
# Draw the expression levels of marker genes #
##############################################
# GCG for alpha cell
GCG_index <-  which(gene_list=="GCG")

# INS for beta cells
INS_index <-  which(gene_list=="INS")

# SST for delta cells
SST_index <-  which(gene_list=="SST")

# PPY for gamma cells
PPY_index <-  which(gene_list=="PPY")

# PRSS1 for acinar cells
PRSS1_index <-  which(gene_list=="PRSS1")

# KRT19 for ductal cells
KRT19_index <-  which(gene_list=="KRT19")

marker_genes <- data.frame(Cell_Index = 1:N, Protocol = metadata$Protocol, 
                           PCA_intri_x =  pca.BUSseq$x[,1], PCA_intri_y =  pca.BUSseq$x[,2],
                           Study = metadata$Study, Celltype = metadata$CellType, 
                           GCG = log1p(unlist(x_corrected[GCG_index,])),
                           INS = log1p(unlist(x_corrected[INS_index,])), 
                           SST = log1p(unlist(x_corrected[SST_index,])),
                           PPY = log1p(unlist(x_corrected[PPY_index,])),
                           PRSS1 = log1p(unlist(x_corrected[PRSS1_index,])),
                           KRT19 = log1p(unlist(x_corrected[KRT19_index,]))
                           )

# Higher expression levels of GCG gene in alpha cells
jpeg(filename = "Image/Pancreas_GCG_for_alpha.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=GCG,shape=(Celltype=="Alpha"))) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values = c(5, 16), guide = guide_legend(reverse=FALSE)) +
  guides(shape = guide_legend(order = 1)) + 
  labs(shape = "Is Alpha Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Higher expression levels of INS gene in beta cells
jpeg(filename = "Image/Pancreas_INS_for_beta.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=INS,
                         shape=(Celltype=="Beta"))) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values = c(5, 16), guide = guide_legend(reverse=TRUE)) +
  guides(shape = guide_legend(order = 1)) + 
  labs(shape = "Is Beta Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Higher expression levels of SST gene in delta cells
jpeg(filename = "Image/Pancreas_SST_for_delta.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=SST,
                         shape=(Celltype=="Delta"))) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values = c(5, 16), guide = guide_legend(reverse=TRUE)) +
  guides(shape = guide_legend(order = 1)) + 
  labs(shape = "Is Delta Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Higher expression levels of PPY gene in gamma cells
jpeg(filename = "Image/Pancreas_PPY_for_gamma.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=PPY,
                         shape=(Celltype=="Gamma"))) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values = c(5, 16), guide = guide_legend(reverse=TRUE)) +
  guides(shape = guide_legend(order = 1)) + 
  labs(shape = "Is Gamma Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Higher expression levels of PRSS1 gene in acinar cells
jpeg(filename = "Image/Pancreas_PRSS1_for_acinar.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=PRSS1,
                         shape=(Celltype=="Acinar"))) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values = c(5, 16), guide = guide_legend(reverse=TRUE)) +
  guides(shape = guide_legend(order = 1)) + 
  labs(shape = "Is Acinar Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Higher expression levels of KRT19 gene in ductal cells
jpeg(filename = "Image/Pancreas_KRT19_for_ductal.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=KRT19,
                         shape=(Celltype=="Ductal"))) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values = c(5, 16), guide = guide_legend(reverse=TRUE)) +
  guides(shape = guide_legend(order = 1)) + 
  labs(shape = "Is Ductal Cell") +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

########################################
# Draw the violin plot of marker genes #
########################################
jpeg(filename = paste0("Image/Pancreas_violin_GCG_for_alpha.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, GCG)) + 
  geom_violin(aes(fill = Celltype),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"))#,
#legend.text=element_text(size=18),
#legend.title = element_text(size = 24))
p
dev.off()

jpeg(filename = paste0("Image/Pancreas_violin_INS_for_beta.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, INS)) + 
  geom_violin(aes(fill = Celltype),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"))#,
#legend.text=element_text(size=18),
#legend.title = element_text(size = 24))
p
dev.off()

jpeg(filename = paste0("Image/Pancreas_violin_SST_for_delta.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, SST)) + 
  geom_violin(aes(fill = Celltype),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"))#,
#legend.text=element_text(size=18),
#legend.title = element_text(size = 24))
p
dev.off()

jpeg(filename = paste0("Image/Pancreas_violin_PPY_for_gamma.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, PPY)) + 
  geom_violin(aes(fill = Celltype),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"))#,
#legend.text=element_text(size=18),
#legend.title = element_text(size = 24))
p
dev.off()

jpeg(filename = paste0("Image/Pancreas_violin_PRSS1_for_acinar.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, PRSS1)) + 
  geom_violin(aes(fill = Celltype),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"))#,
#legend.text=element_text(size=18),
#legend.title = element_text(size = 24))
p
dev.off()

jpeg(filename = paste0("Image/Pancreas_violin_KRT19_for_ductal.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, GCG)) + 
  geom_violin(aes(fill = Celltype),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"))#,
#legend.text=element_text(size=18),
#legend.title = element_text(size = 24))
p
dev.off()

# Store the workspace
# save.image("BUSseq_workspace.RData")
save(ARI_BUSseq,tsne_BUSseq_dist,file = "BUSseq_results.RData")