import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

outcome_path = '../scRNA_data/preprocessed_data/'
adata_pbmc_pca = ad.read_h5ad(outcome_path + 'pbmc_pca.h5ad')

# Let's compute the neighborhood graph of cells using the PCA representation of the count matrix
sc.pp.neighbors(adata_pbmc_pca, n_neighbors=10, n_pcs=10)
adata_pbmc_pca
# AnnData object with n_obs × n_vars = 2638 × 1838
#     obs: 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'
#     var: 'gene_ids', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std'
#     uns: 'hvg', 'log1p', 'pca', 'neighbors'
#     obsm: 'X_pca'
#     varm: 'PCs'
#     obsp: 'distances', 'connectivities'

# Let's embed the neighborhood graph using 2-D UMAP
# pip3 install leidenalg
sc.tl.leiden(adata_pbmc_pca)
sc.tl.paga(adata_pbmc_pca, groups='leiden') # Key for categorical in adata.obs. Default: The first present key of 'leiden' or 'louvain'
sc.pl.paga(adata_pbmc_pca, plot=False) 
sc.tl.umap(adata_pbmc_pca, init_pos='paga')
