import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

outcome_path = '../scRNA_data/preprocessed_data/'
adata_pbmc_filter = ad.read_h5ad(outcome_path + 'pbmc_filter.h5ad')

# Let's build the principal component analysis
sc.tl.pca(adata_pbmc_filter, svd_solver='arpack') # arpack: for the ARPACK wrapper in SciPy
sc.pl.pca(adata_pbmc_filter, color='CST3') # Colour the data points by gene expression of CST3
sc.pl.highest_expr_genes(adata_pbmc_filter, n_top=20)
sc.pl.pca(adata_pbmc_filter, color='CD37')
sc.pl.pca(adata_pbmc_filter, color='CTSS')

# Let's view the contribution of each PC to the total variance in the data. 
# We reply on this information to decide how many PCs is going to be used to compute the neighborhood relations of cells
sc.pl.pca_variance_ratio(adata_pbmc_filter, log=True)

adata_pbmc_filter.write(outcome_path + 'pbmc_pca.h5ad')
