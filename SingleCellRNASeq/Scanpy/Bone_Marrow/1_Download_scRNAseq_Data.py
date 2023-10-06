import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams

sc.settings.verbosity = 3  # verbosity: hints (3)
sc.logging.print_versions()
file_path = '../result/'
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(4, 4), facecolor='white') 

adata = sc.datasets.paul15()
adata
# AnnData object with n_obs × n_vars = 2730 × 3451
#     obs: 'paul15_clusters'
#     uns: 'iroot'

adata.X = adata.X.astype('float64') # Data type conversion from float32 to float64
adata.write(file_path + 'bone_marrow_raw.h5ad')
