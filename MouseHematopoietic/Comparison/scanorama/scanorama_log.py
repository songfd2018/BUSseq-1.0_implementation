#ssh cn024
#cd /scratch/project/project24/busseq/scanorama
# conda activate tensorflow_gpuenv
#ipython

import numpy as np
import scanorama

simulation_data = np.transpose(np.genfromtxt("../../BUSseq_omp/data/count_data_simulation_v1.txt", delimiter = " ",dtype=int))
hemat_data = np.transpose(np.genfromtxt("../../data/count_data_hemat_v1.txt", delimiter = " ",dtype=int))
pancreas_data = np.transpose(np.genfromtxt("../../data/count_data_pancreas_v1.txt", delimiter = " ",dtype=int))

simulation_gene_list = ["gene_"+str(i) for i in range(1,10001)]
simulation_gene_list = [simulation_gene_list for i in range(0,4)]

hemat_gene_list = ["gene_"+str(i) for i in range(1,3471)]
hemat_gene_list = [hemat_gene_list for i in range(0,2)]

pancreas_gene_list = ["gene_"+str(i) for i in range(1,2481)]
pancreas_gene_list = [pancreas_gene_list for i in range(0,4)]

simulation_data = [simulation_data[0:150, 0:10000], simulation_data[150:300, 0:10000],
                    simulation_data[300:450, 0:10000], simulation_data[450:600, 0:10000]]
hemat_data = [hemat_data[0:2729, 0:3470], hemat_data[2729:4649, 0:3470]]
pancreas_data = [pancreas_data[0:1006, 0:2480],pancreas_data[1006:3337, 0:2480],
                pancreas_data[3337:4932, 0:2480],pancreas_data[4932:7095, 0:2480]]

simulation_integrated, simulation_corrected, simulation_genes = scanorama.correct(simulation_data, simulation_gene_list, return_dimred=True)
hemat_integrated, hemat_corrected, hemat_genes = scanorama.correct(hemat_data, hemat_gene_list, return_dimred=True)
pancreas_integrated, pancreas_corrected, pancreas_genes = scanorama.correct(pancreas_data, pancreas_gene_list, return_dimred=True)

#Please note that *_corrected is of type "scipy.sparse.csr.csr_matrix".
for b in range(4):
    np.savetxt("simulation_integrated_batch"+str(b+1)+".txt", simulation_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("simulation_corrected_batch"+str(b+1)+".txt", simulation_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("simulation_genes.txt", simulation_genes, fmt="%s", delimiter="\t")
for b in range(2):
    np.savetxt("hemat_integrated_batch"+str(b+1)+".txt", hemat_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("hemat_corrected_batch"+str(b+1)+".txt", hemat_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("hemat_genes.txt", hemat_genes, fmt="%s", delimiter="\t")
for b in range(4):
    np.savetxt("pancreas_integrated_batch"+str(b+1)+".txt", pancreas_integrated[b], fmt = '%10.9f', delimiter="\t")
    np.savetxt("pancreas_corrected_batch"+str(b+1)+".txt", pancreas_corrected[b].toarray(), fmt = '%10.9f', delimiter="\t")
np.savetxt("pancreas_genes.txt", pancreas_genes, fmt="%s", delimiter="\t")

#R part
library(data.table)

simulation_corrected = NULL
for (b in 1:4)
    simulation_corrected = rbind(simulation_corrected, fread(paste0("simulation_corrected_batch",b,".txt"),data.table=F))
simulation_genes = read.table("simulation_genes.txt",header=F)
colnames(simulation_corrected) = simulation_genes[,1]

hemat_corrected = NULL
for (b in 1:2)
    hemat_corrected = rbind(hemat_corrected, fread(paste0("hemat_corrected_batch",b,".txt"),data.table=F))
hemat_genes = read.table("hemat_genes.txt",header=F)
colnames(hemat_corrected) = hemat_genes[,1]

pancreas_corrected = NULL
for (b in 1:4)
    pancreas_corrected = rbind(pancreas_corrected, fread(paste0("pancreas_corrected_batch",b,".txt"),data.table=F))
pancreas_genes = read.table("pancreas_genes.txt",header=F)
colnames(pancreas_corrected) = pancreas_genes[,1]

simulation_kmeans = kmeans(simulation_corrected, 5)
hemat_kmeans = kmeans(hemat_corrected, 6)
pancreas_kmeans = kmeans(pancreas_corrected, 8)

write.table(simulation_kmeans$cluster,file = "simulation_clusters.txt",sep="\t",row.names = F,col.names=F)
write.table(hemat_kmeans$cluster,file = "hemat_clusters.txt",sep="\t",row.names = F,col.names=F)
write.table(pancreas_kmeans$cluster,file = "pancreas_clusters.txt",sep="\t",row.names = F,col.names=F)