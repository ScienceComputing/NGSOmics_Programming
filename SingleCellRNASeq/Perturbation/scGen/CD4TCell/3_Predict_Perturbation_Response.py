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
import scgen # https://scgen.readthedocs.io/en/stable/api/reference/scgen.SCGEN.html#scgen.SCGEN
file_path = 'gdrive/MyDrive/Perturbation/'
adata = pt.dt.kang_2018()
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata # (16798, 15706)
cd4_train_stim = adata[((adata.obs["cell_type"] == "CD4 T cells") & (adata.obs["condition"] == "stimulated"))]

# Load the trained model
model = scgen.SCGEN.load(file_path + 'single_cell_data/perturbation_model.pt')

# Make the prediction
pred, delta = model.predict(
    ctrl_key="control", stim_key="stimulated", celltype_to_predict="CD4 T cells"
)
pred # !Predicted gene expression per gene in each cell

# Annotate the predicted CD4+T cells with the label "predicted stimulated" to distinguish them from ground truth cells
pred.obs["condition"] = "predicted stimulated"
ctrl_adata = adata[((adata.obs["cell_type"] == "CD4 T cells") & (adata.obs["condition"] == "control"))]
