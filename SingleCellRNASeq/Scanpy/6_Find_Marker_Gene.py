import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

outcome_path = '../scRNA_data/preprocessed_data/'
adata_pbmc_leiden = ad.read_h5ad(outcome_path + 'pbmc_leiden.h5ad')
adata_pbmc_leiden.uns['log1p']["base"] = None
adata_pbmc_leiden.uns['log1p']
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
# https://github.com/scverse/scanpy/issues/2239

# Let's compute the ranking for the highly differential genes per cluster vs the rest clusters 
# By default, the .raw attribute of AnnData is used if it has been initialized
# Approach 1: use t test to examine the differential genes
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_pbmc_leiden, n_genes=25, sharey=False)

# Approach 2: use wilcoxon test to examine the differential genes
sc.settings.verbosity = 2 
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', method='wilcoxon') # Recommend wilcoxon over t test
sc.pl.rank_genes_groups(adata_pbmc_leiden, n_genes=25, sharey=False)
# Consider MAST, limma, DESeq2, and, diffxpy

adata_pbmc_leiden.write(outcome_path + 'pbmc_wilcoxon.h5ad')

# Approach 3: use logistic regression to examine the differential genes
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata_pbmc_leiden, n_genes=25, sharey=False)

# Find the overlapping marker genes detected by all approaches
