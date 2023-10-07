import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import pertpy as pt
import scgen

# Create and save the model 
# Recap: adata_train = adata[~((adata.obs["cell_type"] == "CD4 T cells") & (adata.obs["condition"] == "stimulated"))].copy()
model = scgen.SCGEN(adata_train, n_hidden=800, n_latent=100, n_layers=2)
"""
n_hidden: number of nodes per hidden layer
n_latent: dimensionality of the latent space
n_layers: number of hidden layers used for encoder and decoder NNs
"""
model.save("../model/perturbation.pt", overwrite=True)

# Train the model
model.train(max_epochs=100, batch_size=32, early_stopping=True, early_stopping_patience=25)
"""
max_epochs: the maximum number of iterations the model is allowed to update its parameters; the higher values training epochs will take more computation time but might help better results
batch_size: number of samples (individual cells) the model sees to update its parameters;t the lower numbers usually help better results
early_stopping: which enables the model to stop the training if its results are not improved after `early_stopping_patience` training epochs; the early stopping mechanism prevents potential overfitting of the training data, which can lead to poor generalization to unseen cell populations
"""

# Visualize the representation of data learned by the model using UMAP
adata_train.obsm["scgen"] = model.get_latent_representation() # get_latent_representation() returns a 100-dimensional vector for each cell
