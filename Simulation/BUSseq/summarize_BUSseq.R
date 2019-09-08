#Apply BUSseq to the simulation study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

# coloring
library(scales)
library(WGCNA)
library(RColorBrewer)

# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv, cex.axis = 2.5)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 3.5, cex =3)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 3.5, cex = 3)
    else mtext(side = 1, KeyValueName, line = 3.5, cex = 3)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y),cex.axis=2.5)
      title(main = "Color Key", cex = 2)
      par(cex = 0.5)
      #mtext(side = 2, "Density", line = 2.3, cex = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy),cex.axis=2.5)
      title(main = "Color Key", cex.main = 3)
      par(cex = 0.5)
      #mtext(side = 2, "Count", line = 2.3, cex = 2)
    }
    else title(main = "Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

set.seed(12345)

#######################
# Load Simulated Data #
#######################
# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/Simulation/BUSseq/")

# Loading pancreatic count data
y_obs <- read.table("../RawCountData/count_data_simulation_v4.txt")

# Load dimension
dim <- read.table("../RawCountData/dim_simulation_v4.txt")
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table("../RawCountData/metadata_simulation_v4.txt")

# Load gene_list
gene_list <- paste0("gene-",1:G)

##########################################################
# Compare the BIC values for different Ks as Figure 3(j) #
##########################################################
proj <- "simulation"
k_sel <-3:10
N_k <- length(k_sel)

BIC.record <- rep(NA,N_k)
names(BIC.record) <- paste0("K=",k_sel)
v <- 4
for(k in k_sel){
  post_dir<-paste("Inference_K",k,sep="")
  if(dir.exists(post_dir)){
    # load BIC
    BIC<- unlist(read.table(paste(post_dir,"/BIC.txt",sep="")))
    BIC.record[k-k_sel[1] + 1] <- BIC[1]
    message(paste("Finish loading the results of ",proj,"_v",v,"_K",k,"!/n",sep=""))
  }else{
    message(paste("The results of ",proj,"_v",v,"_K",k," doesn't exist!/n",sep=""))
  }
}

jpeg("Image/BIC_values.jpeg",width = 540, height = 720)
par(mar=c(5.1,5.1,4.1,2.1))
plot(k_sel,BIC.record,xlab= "K",ylab = "BIC",type="n",cex.axis=2.5,cex.lab=3)
points(k_sel,BIC.record,type="b",pch=19,cex=3)
dev.off()

#####################################
# Load posterior inference of K = 5 #
#####################################
K <- 5

# load w_est
w.est <- read.table("Inference_K5/w_est.txt")
w_BUSseq <- unlist(w.est)

# load alpha_est
# alpha.post <- as.matrix(read.big.matrix("alpha_post.txt",sep=" ",skip=num.burntin,type="double"))
alpha.est <- read.table("Inference_K5/alpha_est.txt")
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table("Inference_K5/beta_est.txt")
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table("Inference_K5/nu_est.txt")
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table("Inference_K5/delta_est.txt")
delta.est <- unlist(delta.est)

# load gamma_est
gamma.est <- read.table("Inference_K5/gamma_est.txt")
gamma.est <- matrix(unlist(gamma.est),B,2)


# load phi_est
phi.est <- read.table("Inference_K5/phi_est.txt")
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table("Inference_K5/pi_est.txt")
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table("Inference_K5/p_est.txt")

# load tau0_est
tau0.est <- read.table("Inference_K5/tau0_est.txt")

# load PPI_est
PPI.est <- read.table("Inference_K5/PPI_est.txt")
D.est <- unlist(read.table("Inference_K5/IG_est.txt"))

# load X_imputed
x_imputed <- read.table("Inference_K5/x_imputed.txt")


################################################
# Prepare intrinsic genes for Pathway analysis #
################################################
intri_gene <- gene_list[D.est==1]
write.table(intri_gene, file = paste0("intrinsic_gene_list.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE)

############################
# Batch effects correction #
############################
adjusted_values = function(Truereads, .B, .nb, .N, .G, .K,
                          .alpha,.beta,.nu,.delta,.phi,.w){
    .w = .w + 1
    .b = unlist(sapply(seq_len(.B), function(b) rep(b, .nb[b])))
    #uniform random numbers are obtained outside of parallel code to handle seed issue.
    global_u = array(runif(.N*.G),c(.G,.N))
    get_corrected_read = function(i){
        b = .b[i]
        # The following are vectors of length .G
        p_x <- pnbinom(Truereads[,i], size = .phi[,b], mu = exp(.alpha + .beta[,.w[i]] + .nu[,b] + .delta[i]))
        p_xminus1 <- pnbinom(Truereads[,i] - 1, size = .phi[,b], mu = exp(.alpha + .beta[,.w[i]] + .nu[,b] + .delta[i]))
        local_u = mapply(function(u) min(u, .999), global_u[,i]*(p_x - p_xminus1) + p_xminus1)
        return(qnbinom(local_u, size = .phi[,1], mu = exp(.alpha + .beta[,.w[i]])))
    }
    return(mapply(get_corrected_read, seq_len(.N)))
    
}

start_time<-Sys.time()
print("Calculate corrected read counts:")
x_corrected<-adjusted_values(x_imputed, B, nb, N, G, K,
                             alpha.est,beta.est,nu.est,delta.est,phi.est,w_BUSseq)
write.table(x_corrected, file = "x_corrected.txt", row.names = FALSE, col.names = FALSE)
end_time<-Sys.time()
running_time<-end_time-start_time
print(running_time)

#######
# ARI #
#######
ARI_BUSseq <- adjustedRandIndex(metadata$celltype,w_BUSseq)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")


# plot colored by cell type
allcolors<-color_by_celltype[metadata$celltype]#[allsamples])

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=alpha(allcolors[subset],0.6)) 
  dev.off()
}

# plot colored by batch
batch.cols<-color_by_batch[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=alpha(batch.cols[subset],0.6)) 
  dev.off()
}

#############################################################
# Draw t-SNE plots by batch and true cell types as Figure 4 #
#############################################################
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- "Image/tsne_simulation_BUSseq_by_celltype.jpeg"
plot_by_celltype(BUSseq_by_celltype, tsne_BUSseq_dist$Y)

BUSseq_by_batch <- "Image/tsne_simulation_BUSseq_by_batch.jpeg"
plot_by_batch(BUSseq_by_batch, tsne_BUSseq_dist$Y)

# #####################
# # Draw the PCA plot #
# #####################
# pca.BUSseq <- prcomp(t(log1p(x_corrected[D.est==1,])), rank=2)
# BUSseq_PCA<- "Image/pca_simulation_BUSseq.jpeg"
# plot_by_celltype(BUSseq_PCA, pca.BUSseq$x, xlab = "PC 1", ylab = "PC 2")

####################
# Draw the Heatmap #
####################
# gene reorder
gene_matched <- cbind(1:G, rep(rep(1:4,each = 50), G/50/4))
gene_matched[501:G,] <- gene_matched[500 + order(gene_matched[501:G,2]),]

# Colors in the heatmap
colorsChoice <- colorRampPalette(c("#F2F2F2","black"))

# The splitting points for log-scale mean expression levels and batch effects
break_effects <- seq(-3, 5, 0.8)

##########################################################################
# The heatmap of the estimated log-scale mean expression levels (Fig 3f) #
log_mu_celltype_est <- alpha.est + beta.est

# Sorting cell-type-specific expression levels 
# in the decreasing order of the proprotions of the first cell type
# because the true cell types are sorted in the same decreasing order
order.est <- order(pi.est[1,],decreasing = T)

jpeg("Image/heatmap_simulation_est_log_mean_expression_levels.jpeg",width = 1080, height = 1440)
heatmap.3(log_mu_celltype_est[gene_matched[,1], order.est],
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_celltype,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(1,0.4,3.6),
          col=colorsChoice(10),breaks = break_effects,symkey = FALSE,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"), trace = "none",
          keysize = 4)
dev.off()

##################################################
# The heatmap of location batch effects (Fig 3g) #
batch_effects_est <- nu.est

jpeg("Image/heatmap_simulation_est_batch_effects.jpeg",width = 1080, height = 1440)
heatmap.3(batch_effects_est[gene_matched[,1],],
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_batch,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(1,0.4,3.6),
          col=colorsChoice(10),breaks = break_effects,symkey = FALSE,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"), cexRow=0.5,trace = "none",
          keysize = 4)
dev.off()

#####################################################################
# The heatmap of log-scale underlying true read count data (Fig 3h) #
log_imputed_count_est <- log1p(x_imputed)

# Batch color bar for count data matrix
color_by_batch_count <- rep(color_by_batch,nb)

# Cell type color bar for count data matrix
color_by_celltype_count <- rep(NA,N)
celltype_est <- unlist(w.est) + 1
for(i in 1:N){
  color_by_celltype_count[i] <- color_by_celltype[which(celltype_est[i]==order.est)]
}


# Upper color bar for batch, lower color bar for cell type
col_annotation<-cbind(color_by_celltype_count,color_by_batch_count)
colnames(col_annotation) <- NULL

# The splitting points for log-scale count data
break_logcount <- seq(0,9,0.9)

jpeg("Image/heatmap_simulation_log_imputed_count.jpeg",width = 1080, height = 1440)
heatmap.3(log_imputed_count_est[gene_matched[,1],],
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(1,0.4,3.6),
          col=colorsChoice(10),breaks = break_logcount,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"), cexRow=0.5,trace = "none",
          keysize = 4)
dev.off()


###############################################################
# The heatmap of log-scale corrected read count data (Fig 3i) #

# Sharing the same double color bar and splitting points 
# as that of the underlying true count data
jpeg("Image/heatmap_simulation_log_corrected_count.jpeg",width = 1080, height = 1440)
heatmap.3(log1p(x_corrected)[gene_matched[,1],],
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(1,0.4,3.6),
          col=colorsChoice(10),breaks = break_logcount,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"), cexRow=0.5,trace = "none",
          keysize = 4)
dev.off()

####################################################################################################################
# The scatter plot of the estimated cell-specific size factors versus the true cell-specific size factors (Fig 3e) #

# Loading the true cell-specific size factors
cell_effect_true_vec <- unlist(read.table("../True_para/delta_syn.txt"))

# Obtaining the estimated cell-specific size factors
cell_effect_est <- delta.est

jpeg("./Image/scatter_simulation_cell_size_factors.jpeg",width = 540, height = 720)
par(mar = c(5.1,6.1,4.1,2.1)) 
plot(cell_effect_true_vec,cell_effect_est,xlab = expression(paste("True ",delta)), ylab = expression(paste("Estimated ",delta)),type="n",ylim = c(-3,2),xlim = c(-3,2),cex.axis = 3, cex.lab = 3)
points(cell_effect_true_vec,cell_effect_est,pch = 1,cex=3)
abline(a=0,b=1,lty=3,cex=3)
dev.off()

# Store the workspace
save.image("BUSseq_workspace.RData")
save(ARI_BUSseq,tsne_BUSseq_dist,file = "BUSseq_results.RData")