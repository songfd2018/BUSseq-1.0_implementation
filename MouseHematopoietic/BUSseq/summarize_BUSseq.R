#Apply BUSseq to the hematopoietic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

set.seed(123)
###################
# Load hemat Data #
###################
# Working directory
#setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/BUSseq")
setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/BUSseq")

# Loading hemat count data
y_obs <- read.table("../RawCountData/count_data_hemat_v1.txt")

# Load dimension
dim <- read.table("../RawCountData/dim_hemat_v1.txt")
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table("../RawCountData/hemat_metadata.txt")

##########################
#load posterior inference#
##########################
K <- 6

# load w_est
w.est <- read.table("Inference_K6/w_est.txt")
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table("Inference_K6/alpha_est.txt")
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table("Inference_K6/beta_est.txt")
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table("Inference_K6/nu_est.txt")
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table("Inference_K6/delta_est.txt")
delta.est <- unlist(delta.est)
# plot(delta.est,col=rep(1:B,nb) + 1)


# load gamma_est
gamma.est <- read.table("Inference_K6/gamma_est.txt")
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table("Inference_K6/phi_est.txt")
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table("Inference_K6/pi_est.txt")
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table("Inference_K6/p_est.txt")

# load tau0_est
tau0.est <- read.table("Inference_K6/tau0_est.txt")

# load PPI_est
PPI.est <- read.table("Inference_K6/PPI_est.txt")
D.est <- unlist(read.table("Inference_K6/IG_est.txt"))

# load X_imputed
x_imputed <- read.table("Inference_K6/x_imputed.txt")

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
end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)

####################
# Pathway analysis #
####################
load("../RawCountData/hemat_countdata.RData")
gene_set <- rownames(dataA2)
intri_gene <- gene_set[D.est==1]
write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

#######
# ARI #
#######
ARI_BUSseq <- adjustedRandIndex(metadata$CellType,w_BUSseq)

###########################
# Silhouette Coefficients #
###########################
# Silhouette coefficients of the corrected read count matrix of 
# the identified intrinsic genes by BUSseq
BUSseq_dist <- dist(t(log1p(x_corrected[D.est==1,])))
sil_BUSseq <- silhouette(as.integer(w_BUSseq), dist = BUSseq_dist)

sil_BUSseq_true <- silhouette(as.integer(metadata$CellType), dist = BUSseq_dist)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

# coloring for cell types
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", other="grey")
colmatF <- col2rgb(color.legendF) 
colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendA[as.character(metadata$CellType[1:nb[1]])],color.legendF[as.character(metadata$CellType[1:nb[2]+nb[1]])])
first.batch <- rep(c(TRUE, FALSE), nb)

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, allcolors, "black")[subset],
       bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  
  # legend(x=-15, y=-5, legend=names(color.legendF), pch=21, cex=2, col="black", pt.bg=color.legendF, title = "GSE81682", bty="n")
  # legend(x=-11, y=-5, legend=names(color.legendA)[1:3], pch=1, cex=2, col=color.legendA[1:3], pt.bg=color.legendA[1:3], title = "GSE72857", bty="n")
  # 
  dev.off()
}

batchcolor <- c("lightcoral","lavender")
batch_ind <- rep(1:2,nb)
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, batchcolor[batch_ind[subset]], "black"),
       bg=batchcolor[batch_ind[subset]], xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")

  dev.off()
}

celltype_est <- rep(NA,N)
celltype_est[which(w_BUSseq == 0)] <- "other"
celltype_est[which(w_BUSseq == 1)] <- "CMP"
celltype_est[which(w_BUSseq == 2)] <- "GMP"
celltype_est[which(w_BUSseq == 3)] <- "MPP"
celltype_est[which(w_BUSseq == 4)] <- "MEP"
celltype_est[which(w_BUSseq == 5)] <- "LTHSC"

estcolors <- color.legendF[celltype_est]
plot_by_celltype_est<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, estcolors, "black")[subset],
       bg=estcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  
  
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- "Image/tsne_hemat_BUSseq_by_celltype.jpeg"
plot_by_celltype(BUSseq_by_celltype, tsne_BUSseq_dist$Y)

BUSseq_by_batch <- "Image/tsne_hemat_BUSseq_by_batch.jpeg"
plot_by_batch(BUSseq_by_batch, tsne_BUSseq_dist$Y)

BUSseq_by_celltype_est<- "Image/tsne_hemat_BUSseq_by_celltype_est.jpeg"
plot_by_celltype_est(BUSseq_by_celltype_est, tsne_BUSseq_dist$Y)

#####################
# Draw the PCA plot #
#####################
pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)

BUSseq_PCA<- "Image/pca_hemat_BUSseq.jpeg"
plot_by_celltype(BUSseq_PCA, pca.BUSseq$x, xlab = "PC 1", ylab = "PC 2")


BUSseq_PCA_by_batch<- "Image/pca_hemat_BUSseq_by_batch.jpeg"
plot_by_batch(BUSseq_PCA_by_batch, pca.BUSseq$x, xlab = "PC 1", ylab = "PC 2")

BUSseq_PCA_by_celltype_est<- "Image/pca_hemat_BUSseq_by_celltype_est.jpeg"
plot_by_celltype_est(BUSseq_PCA_by_celltype_est, pca.BUSseq$x, xlab = "PC 1", ylab = "PC 2")

##############################################
# Draw the expression levels of marker genes #
##############################################
gene_list <- rownames(dataA2)

# ENSMUSG00000002985 is APOE gene
APOE_index <- which(gene_list=="ENSMUSG00000002985")
# ENSMUSG00000015053 is GATA2 gene
GATA2_index <- which(gene_list=="ENSMUSG00000015053")
# ENSMUSG00000004552 is CTSE gene for MEP
CTSE_index <-  which(gene_list=="ENSMUSG00000004552")

marker_genes <- data.frame(Cell_Index = 1:N, Protocol = metadata$Protocol, 
                           PCA_intri_x =  pca.BUSseq$x[,1], PCA_intri_y =  pca.BUSseq$x[,2],
                           Study = metadata$Study, Celltype = metadata$CellType,
                           Celltype_est = w_BUSseq,
                           APOE = log1p(unlist(x_corrected[APOE_index,])),
                           GATA2 = log1p(unlist(x_corrected[GATA2_index,])), 
                           CTSE = log1p(unlist(x_corrected[CTSE_index,])))

# Expression levels of APOE gene indicating different stages of erythropoiesis
jpeg(filename = "Image/Hmeat_APOE_from_CMP_to_MEP.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=APOE,
                         shape=Celltype)) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Expression levels of GATA2 gene indicating different stages of erythropoiesis
jpeg(filename = "Image/Hmeat_GATA2_from_CMP_to_MEP.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=GATA2,
                         shape=Celltype)) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0)) +
  xlab("PC 1") + ylab("PC 2") +
  theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
    size = 18), #angle = 45),
    axis.text.y = element_text(#face = "bold", color = "blue", 
      size = 18),#, angle = 45))
    axis.title=element_text(size=24,face="bold"),
    legend.text=element_text(size=18),
    legend.title = element_text(size = 24))
dev.off()

# Higher expression levels of CTSE gene in MEPs
jpeg(filename = "Image/Hmeat_CTSE_for_MEP.jpeg",width = 960, height = 720,quality = 100)
ggplot(marker_genes, aes(x=PCA_intri_x, y=PCA_intri_y, colour=CTSE,
                         shape=Celltype)) +
  geom_point(size=3) + theme_classic() +
  scale_colour_continuous(type = "viridis") +
  scale_shape_manual(values=c(18, 19, 24, 15, 11, 103, 7, 23, 0)) +
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
jpeg(filename = paste0("Image/Hmeat_violin_CTSE_for_MEP.jpeg"),width = 960, height = 720,quality = 100)
p <- ggplot(marker_genes, aes(x = Celltype, CTSE)) + 
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

# jpeg(filename = paste0("Image/Hmeat_violin_CTSE_for_MEP_est.jpeg"),width = 960, height = 720,quality = 100)
# p <- ggplot(marker_genes, aes(x = Celltype_est, CTSE)) + 
#   geom_violin(aes(fill = factor(Celltype_est)),show.legend = FALSE,draw_quantiles = c(0.25, 0.50 , 0.75),scale = "width") +
#   theme(axis.text.x = element_text(#face = "bold", color = "#993333", 
#     size = 18), #angle = 45),
#     axis.text.y = element_text(#face = "bold", color = "blue", 
#       size = 18),#, angle = 45))
#     axis.title=element_text(size=24,face="bold"))#,
# #legend.text=element_text(size=18),
# #legend.title = element_text(size = 24))
# p
# dev.off()

# Store the workspace
save.image("BUSseq_workspace.RData")