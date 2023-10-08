# pip install anndata==0.9.2
# pip install git+https://github.com/theislab/scgen.git
# pip install pertpy
# pip install scanpy
# pip install numba
# pip install llvmlite

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import pertpy as pt
import scgen # pip install git+https://github.com/theislab/scgen.git

adata = pt.dt.kang_2018()
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata

adata.obs.label
# index
# AAACATACATTTCC-1    ctrl
# AAACATACCAGAAA-1    ctrl
# AAACATACCATGCA-1    ctrl
# AAACATACCTCGCT-1    ctrl
# AAACATACCTGGTA-1    ctrl
#                     ...
# TTTGCATGCCTGAA-2    stim
# TTTGCATGCCTGTC-2    stim
# TTTGCATGCTAAGC-2    stim
# TTTGCATGGGACGA-2    stim
# TTTGCATGTCTTAC-2    stim
# Name: label, Length: 24673, dtype: category
# Categories (2, object): ['ctrl', 'stim']
adata.obs.rename({"label": "condition"}, axis=1, inplace=True)
adata.obs["condition"].replace({"ctrl": "control", "stim": "stimulated"}, inplace=True)
adata.obs.condition
# index
# AAACATACATTTCC-1       control
# AAACATACCAGAAA-1       control
# AAACATACCATGCA-1       control
# AAACATACCTCGCT-1       control
# AAACATACCTGGTA-1       control
#                        ...
# TTTGCATGCCTGAA-2    stimulated
# TTTGCATGCCTGTC-2    stimulated
# TTTGCATGCTAAGC-2    stimulated
# TTTGCATGGGACGA-2    stimulated
# TTTGCATGTCTTAC-2    stimulated
# Name: condition, Length: 24673, dtype: category
# Categories (2, object): ['control', 'stimulated']

# Count the cell types
adata.obs.cell_type.value_counts()

# Exclude CD4T cells from the training data (adata_train) to mimic a real-world scenario where a specific population isn't captured in an experiment.
adata_train = adata[~((adata.obs["cell_type"] == "CD4 T cells") & (adata.obs["condition"] == "stimulated"))].copy()
# ~ is the logical NOT operator, which negates the combined condition. It returns True for rows that do not satisfy the combined condition.

cd4_train_stim = adata[((adata.obs["cell_type"] == "CD4 T cells") & (adata.obs["condition"] == "stimulated"))].copy()
cd4_train_stim.to_df() # View the cell x gene count matrix in CD4+T cells with perturbations
