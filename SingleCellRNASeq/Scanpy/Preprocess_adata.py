import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

outcome_path = '../scRNA_data/preprocessed_data/pbmc_v0.h5ad'
adata_pbmc = ad.read_h5ad(outcome_path)

# Display the genes that contribute the largest proportion of counts within each individual cell, spanning all the cells.
sc.pl.highest_expr_genes(adata_pbmc, n_top=30)
