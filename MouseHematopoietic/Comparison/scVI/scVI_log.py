#This has to be run on Central Cluster (or a system with Centos 7/ RHEL 7)
# ssh chpc-cn001
# conda activate tensorflow_env
#cd scVI
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

batch_1=CsvDataset("count_data_simulation_v1_batch1.csv", new_n_genes = 10000)
batch_2=CsvDataset("count_data_simulation_v1_batch2.csv", new_n_genes = 10000)
batch_3=CsvDataset("count_data_simulation_v1_batch3.csv", new_n_genes = 10000)
batch_4=CsvDataset("count_data_simulation_v1_batch4.csv", new_n_genes = 10000)

data = GeneExpressionDataset.concat_datasets(batch_1,batch_2,batch_3,batch_4)

vae = VAE(data.nb_genes, n_batch=data.n_batches, n_labels=data.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

trainer = UnsupervisedTrainer(vae, data, train_size=0.9)

trainer.train(n_epochs=100)

# Save result
#torch.save(trainer.model.state_dict(),'data/simulation_vae_harmonization.pkl')
#torch.save(hemat_trainer.model.state_dict(),'data/hemat_vae_harmonization.pkl')

# Load result
#trainer.model.load_state_dict(torch.load('data/simulation_vae_harmonization.pkl'))
#trainer.model.eval()

full = trainer.create_posterior(trainer.model, data, indices=np.arange(len(data)))
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("scVI_simulation_latent.txt", latent, fmt="%10.9f",delimiter="\t")

#Get clusters
adata_latent = sc.AnnData(latent)
sc.pp.neighbors(adata_latent, use_rep='X', n_neighbors=30, metric='minkowski')
sc.tl.louvain(adata_latent, partition_type=louvain.ModularityVertexPartition, use_weights=False)
clusters = adata_latent.obs.louvain.values.to_dense().astype(int)

np.savetxt("scVI_simulation_clusters.txt", clusters, fmt="%d",delimiter="\t")

#Get differentiallyexpressed genes
de_res, de_clust = full.one_vs_all_degenes(cell_labels=clusters, n_samples=10000, 
                                           M_permutation=10000, output_file=False,
                                           save_dir="./data", filename='Harmonized_ClusterDE',
                                           min_cells=1)

de_genes = de_res[0].loc[np.abs(de_res[0].bayes1)>3].index.values
for b in range(1,2):
    de_genes = np.append(de_genes, de_res[b].loc[np.abs(de_res[b].bayes1)>3].index.values)

np.savetxt("scVI_simulation_de_genes.txt",de_genes, fmt="%s", delimiter="\t")

#de_res_stim, de_clust_stim = full.within_cluster_degenes(cell_labels=clusters,
#                                                         states=data.batch_indices.ravel() == 1,
#                                                         output_file=False, batch1=[0], batch2=[1],
#                                                         save_dir="./data", filename='Harmonized_Simulation_DE',
#                                                         min_cells=1)

#hemat_de_res_stim, hemat_de_clust_stim = hemat_full.within_cluster_degenes(cell_labels=hemat_clusters,
#                                                        states=hemat_data.batch_indices.ravel() == 1,
#                                                        output_file=False, batch1=[0], batch2=[1],
#                                                        save_dir="./data", filename='Harmonized_hemat_DE',
#                                                        min_cells=1)

#########
##hemat##
#########
hemat_batch_1=CsvDataset("count_data_hemat_v1_batch1.csv", new_n_genes = 3470)
hemat_batch_2=CsvDataset("count_data_hemat_v1_batch2.csv", new_n_genes = 3470)

hemat_data = GeneExpressionDataset.concat_datasets(hemat_batch_1,hemat_batch_2)

hemat_vae = VAE(hemat_data.nb_genes, n_batch=hemat_data.n_batches, n_labels=hemat_data.n_labels,
                n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

hemat_trainer = UnsupervisedTrainer(hemat_vae, hemat_data, train_size=0.9)

hemat_trainer.train(n_epochs=100)

hemat_full = hemat_trainer.create_posterior(hemat_trainer.model, hemat_data, indices=np.arange(len(hemat_data)))
hemat_latent, hemat_batch_indices, hemat_labels = hemat_full.sequential().get_latent()
hemat_batch_indices = hemat_batch_indices.ravel()

np.savetxt("scVI_hemat_v1_latent_0716.txt", hemat_latent, fmt="%10.9f",delimiter="\t")

hemat_adata_latent = sc.AnnData(hemat_latent)
sc.pp.neighbors(hemat_adata_latent, use_rep='X', n_neighbors=30, metric='minkowski')
sc.tl.louvain(hemat_adata_latent, partition_type=louvain.ModularityVertexPartition, use_weights=False)
hemat_clusters = hemat_adata_latent.obs.louvain.values.to_dense().astype(int)

np.savetxt("scVI_hemat_v1_clusters_0716.txt", hemat_clusters, fmt="%d",delimiter="\t")

hemat_de_res, hemat_de_clust = hemat_full.one_vs_all_degenes(cell_labels=hemat_clusters, n_samples=10000, 
                                           M_permutation=10000, output_file=False,
                                           save_dir="./data", filename='hemat_harmonized_clusterDE',
                                           min_cells=1)

hemat_de_genes = hemat_de_res[0].loc[np.abs(hemat_de_res[0].bayes1)>3].index.values
for b in range(1,2):
    hemat_de_genes = np.append(hemat_de_genes, hemat_de_res[b].loc[np.abs(hemat_de_res[b].bayes1)>3].index.values)

np.savetxt("scVI_hemat_v1_de_genes_0716.txt",hemat_de_genes, fmt="%s", delimiter="\t")

############
##pancreas##
############

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