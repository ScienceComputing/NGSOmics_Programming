import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams

sc.settings.verbosity = 3  # verbosity: hints (3)
sc.logging.print_versions()
file_path = '../result/bone_marrow.h5ad'
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(4, 4), facecolor='white') 
