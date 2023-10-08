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
import scgen
file_path = 'gdrive/MyDrive/Perturbation/'

# Create the model 
# Recap: 
# adata_train = adata[~((adata.obs["cell_type"] == "CD4 T cells") & (adata.obs["condition"] == "stimulated"))].copy()
adata_train.shape # (18995, 15706)
scgen.SCGEN.setup_anndata(adata_train, batch_key="condition", labels_key="cell_type")
model = scgen.SCGEN(adata_train, n_hidden=800, n_latent=100, n_layers=2)
"""
n_hidden: number of nodes per hidden layer
n_latent: dimensionality of the latent space
n_layers: number of hidden layers used for encoder and decoder NNs
"""

# Train the model
model.train(max_epochs=100, batch_size=32, early_stopping=True, early_stopping_patience=25)
"""
max_epochs: the maximum number of iterations the model is allowed to update its parameters; the higher values training epochs will take more computation time but might help better results
batch_size: number of samples (individual cells) the model sees to update its parameters;t the lower numbers usually help better results
early_stopping: which enables the model to stop the training if its results are not improved after `early_stopping_patience` training epochs; the early stopping mechanism prevents potential overfitting of the training data, which can lead to poor generalization to unseen cell populations
"""

type(model) # scgen._scgen.SCGEN

# Save and reload the model
model.save(file_path + 'single_cell_data/perturbation_model.pt', overwrite=True, save_anndata=True)
model = scgen.SCGEN.load(file_path + 'single_cell_data/perturbation_model.pt')

# Show the learned latent space
model.get_latent_representation()
model.get_latent_representation().shape # (18995, 100)

# Visualize the representation of data learned by the model using UMAP
adata_train.obsm["scgen"] = model.get_latent_representation() # get_latent_representation() returns a 100-dimensional vector for each cell
adata_train.obsm["scgen"]
# array([[ -7.4319363 ,   6.0347595 ,   4.360623  , ...,  -4.1190267 ,
#           0.4719472 ,   5.999382  ],
#        [ -2.1332135 ,   0.66588956,  -1.8832861 , ...,   3.388834  ,
#          -0.79437214,  -1.1460924 ],
#        [  6.881926  ,  -0.81465465,  -2.5517473 , ...,   3.8420978 ,
#           3.7350404 ,  -7.246197  ],
#        ...,
#        [ -1.1007156 ,   1.9041696 ,   0.51787966, ...,   6.500449  ,
#          -1.1392326 ,   3.870847  ],
#        [ -4.348127  ,  -3.089072  ,   1.0366997 , ...,   4.119344  ,
#          -1.5113178 ,  -0.7055234 ],
#        [-10.115226  ,   4.760014  ,   4.996005  , ...,  -9.193241  ,
#           3.9584777 ,   4.2281528 ]], dtype=float32)

# Retrive the latent space of gene expression learned by the encoder for each cell
adata_train.obsm["scgen"] = model.get_latent_representation() # get_latent_representation() returns a 100-dimensional vector for each cell
adata_train.obsm["scgen"].shape

# Visualize the learned representation of the gene expression space by the model
sc.pp.neighbors(adata_train, use_rep="scgen") # use_rep="scgen" : here we use the latent space leaned by the encoder as the compressed representation of the gene expression space 
sc.tl.umap(adata_train)
sc.pl.umap(adata_train, color=["condition", "cell_type"], wspace=0.4, frameon=False)
