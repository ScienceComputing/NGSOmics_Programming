# https://peps.python.org/pep-0008/#imports
# https://stackoverflow.com/questions/53014306/error-15-initializing-libiomp5-dylib-but-found-libiomp5-dylib-already-initial
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import numpy as np
import pandas as pd

sc.settings.set_figure_params(dpi=300, facecolor='white')
sc.settings.verbosity = 3
sc.logging.print_header()
outcome_path = '../preprocessed_data/pbmc.h5ad'

adata_pbmc = sc.read_10x_mtx(
    '../scRNA_data/filtered_gene_bc_matrices/hg19/',  # Specify the directory that contains the .mtx feature-barcode matrix
    var_names='gene_symbols', # Alternatively, we can use gene_ids
    cache=True) # Write an h5ad cache file to speedup reading next time 
adata_pbmc
# AnnData object with n_obs × n_vars = 2700 × 32738
#    var: 'gene_ids'

adata_pbmc.var_names_make_unique()
adata_pbmc
# AnnData object with n_obs × n_vars = 2700 × 32738
#    var: 'gene_ids'
