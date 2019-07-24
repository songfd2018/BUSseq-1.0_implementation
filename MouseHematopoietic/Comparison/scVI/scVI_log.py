#Please note that the pwd has to be ./ while data is in ./data/

# Python code starts here
import os
import numpy as np
import torch
import matplotlib.pyplot as plt

from scvi.models import SCANVI, VAE
from scvi.inference import UnsupervisedTrainer, JointSemiSupervisedTrainer, SemiSupervisedTrainer
from scvi.dataset.csv import CsvDataset
from scvi.dataset.dataset import GeneExpressionDataset

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

import numpy.random as random
import pandas as pd
import scanpy as sc
import louvain

#from umap import UMAP #This is only for plots

############
##pancreas##
############
os.chdir("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/HumanPancreas/Comparison/scVI")
if not os.path.exists("./data"):
    os.mkdir("./data")
pancreas_full=pd.read_csv("count_data_pancreas_v1.txt",sep=" ",header=None)

pancreas_full.index = ["gene_"+str(i) for i in range(1,2481)]
pancreas_full.columns = ["sample_"+str(i) for i in range(1,7096)]

pancreas_full.iloc[:,0:1006].to_csv("./data/count_data_pancreas_v1_batch1.csv",sep=",")
pancreas_full.iloc[:,1006:3337].to_csv("./data/count_data_pancreas_v1_batch2.csv",sep=",")
pancreas_full.iloc[:,3337:4932].to_csv("./data/count_data_pancreas_v1_batch3.csv",sep=",")
pancreas_full.iloc[:,4932:7095].to_csv("./data/count_data_pancreas_v1_batch4.csv",sep=",")

# write count data into desired format

pancreas_batch_1=CsvDataset("count_data_pancreas_v1_batch1.csv", new_n_genes = 2480)
pancreas_batch_2=CsvDataset("count_data_pancreas_v1_batch2.csv", new_n_genes = 2480)
pancreas_batch_3=CsvDataset("count_data_pancreas_v1_batch3.csv", new_n_genes = 2480)
pancreas_batch_4=CsvDataset("count_data_pancreas_v1_batch4.csv", new_n_genes = 2480)

pancreas_data = GeneExpressionDataset.concat_datasets(pancreas_batch_1,pancreas_batch_2,
                                                        pancreas_batch_3,pancreas_batch_4)

pancreas_vae = VAE(pancreas_data.nb_genes, n_batch=pancreas_data.n_batches, n_labels=pancreas_data.n_labels,
                n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

pancreas_trainer = UnsupervisedTrainer(pancreas_vae, pancreas_data, train_size=0.9)

pancreas_trainer.train(n_epochs=100)

pancreas_full = pancreas_trainer.create_posterior(pancreas_trainer.model, pancreas_data, indices=np.arange(len(pancreas_data)))
pancreas_latent, pancreas_batch_indices, pancreas_labels = pancreas_full.sequential().get_latent()
pancreas_batch_indices = pancreas_batch_indices.ravel()

np.savetxt("scVI_pancreas_v1_latent.txt", pancreas_latent, fmt="%10.9f",delimiter="\t")

pancreas_adata_latent = sc.AnnData(pancreas_latent)
sc.pp.neighbors(pancreas_adata_latent, use_rep='X', n_neighbors=30, metric='minkowski')
sc.tl.louvain(pancreas_adata_latent, partition_type=louvain.ModularityVertexPartition, use_weights=False)
pancreas_clusters = pancreas_adata_latent.obs.louvain.values.to_dense().astype(int)

np.savetxt("scVI_pancreas_v1_clusters.txt", pancreas_clusters, fmt="%d",delimiter="\t")

pancreas_de_res, pancreas_de_clust = pancreas_full.one_vs_all_degenes(cell_labels=pancreas_clusters, n_samples=10000, 
                                           M_permutation=10000, output_file=False,
                                           save_dir="./data", filename='pancreas_harmonized_clusterDE',
                                           min_cells=1)

pancreas_de_genes = pancreas_de_res[0].loc[np.abs(pancreas_de_res[0].bayes1)>3].index.values
for b in range(1,2):
    pancreas_de_genes = np.append(pancreas_de_genes, pancreas_de_res[b].loc[np.abs(pancreas_de_res[b].bayes1)>3].index.values)

np.savetxt("scVI_pancreas_v1_de_genes.txt",pancreas_de_genes, fmt="%s", delimiter="\t")