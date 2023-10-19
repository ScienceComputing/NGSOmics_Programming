# Citation
# Melsted, P., Booeshaghi, A.S. et al. Modular and efficient pre-processing of single-cell RNA-seq. bioRxiv (2019). doi:10.1101/673285 
# Wolf, F. A., Angere, P. and Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology (2018). doi:10.1186/s13059-017-1382-0
from scipy.io import mmread
import anndata
from sklearn.decomposition import TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os
import seaborn as sns
import itertools

from scipy.sparse import csr_matrix
matplotlib.rcParams.update({'font.size': 22})

import sys
sra_id = sys.argv[1]
sra_list = sra_id.split(",")

# Define a density plot function
from scipy.interpolate import interpn

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    sc = ax.scatter( x, y, c=z, **kwargs )
    return sc

import warnings
warnings.filterwarnings("ignore")
color = itertools.cycle(sns.color_palette("rocket", 5))
mtx_dic = {}

print("\n~~~Visulizing the low-dimensional representation of cells...~~~")
for sra_id, color_i in zip(sra_list, color):

    # Import cells x genes matrix
    mtx = mmread("../result_b/{}/counts_unfiltered/cells_x_genes.mtx".format(sra_id))
    
    # Create sparse matrix representation of the count matrix
    mtx_dic[sra_id] = csr_matrix(mtx)

    # Set the layout
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))    

    print ("Left plot shows 2D PCA for: " + sra_id)
    # Represent the cells in 2D using PCA
    # Perform SVD
    tsvd = TruncatedSVD(n_components=2)
    tsvd.fit(mtx)
    X = tsvd.transform(mtx)
    # Plot the cells in the 2D PCA projection
    ax1.scatter(X[:,0], X[:,1], alpha=0.5, c=color_i)

    print ("Middle plot shows the density PCA plot for: " + sra_id)
    # Plot the density display for 2D PCA plot
    x = X[:,0]
    y = X[:,1]
    sc = density_scatter(x, y, ax=ax2, cmap=sns.color_palette("rocket_r", as_cmap=True))
    fig.colorbar(sc, ax=ax2)

    print ("Right plot shows 3rd and 4th subspaces for: " + sra_id)
    # Explore more PCA subspsaces
    n_components = 4 # Number of dimensions for the PCA reduction 
    dimension_A = 3 # 3rd subspace
    dimension_B = 4 # 4th subspace
    # Perform SVD
    tsvd = TruncatedSVD(n_components)
    tsvd.fit(mtx)
    X = tsvd.transform(mtx)
    ax3.scatter(X[:,dimension_A-1], X[:,dimension_B-1], alpha=0.5, c=color_i)
    
    ax1.axis("off")
    ax2.axis("off")
    ax3.axis("off")

    plt.savefig("../result_b/figure/PCA_QC_{}.png".format(sra_id))
    print("\n")


print("\n~~~Visulizing the genes detected as a function of UMI counts...~~~")
for sra_id, color_i in zip(sra_list, color):
    fig, ax = plt.subplots(figsize=(4, 4))
    mtx = mtx_dic[sra_id]
    
    ax.scatter(np.asarray(mtx.sum(axis=1))[:,0], np.asarray(np.sum(mtx>0, axis=1))[:,0], color=color_i, alpha=0.01)
    ax.set_title(sra_id)
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim((0.5, 4500))
    ax.set_ylim((0.5,2000))

    plt.savefig("../result_b/figure/gene_UMIcount_{}.png".format(sra_id))
    print("Complete the plotting for: " + sra_id + "\n")

print("\n~~~Visulizing the knee plot...~~~")
cutoff = 200
for sra_id, color_i in zip(sra_list, color):
    mtx = mtx_dic[sra_id]
    knee = np.sort((np.array(mtx.sum(axis=1))).flatten())[::-1]
    cell_set = np.arange(len(knee))
    num_cells = cell_set[knee > cutoff][::-1][0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))    
    
    ax1.loglog(knee, range(len(knee)), linewidth=5, color=color_i)
    ax1.axvline(x=cutoff, linewidth=3, color="k")
    ax1.axhline(y=num_cells, linewidth=3, color="k")
    ax1.set_title(sra_id + " \n(transposed)")
    ax1.set_xlabel("UMI Counts")
    ax1.set_ylabel("Set of Barcodes")

    ax2.loglog(range(len(knee)), knee, linewidth=5, color=color_i) # Reverse the position of knee and range(len(knee))
    ax2.axvline(x=cutoff, linewidth=3, color="k")
    ax2.axhline(y=num_cells, linewidth=3, color="k")
    ax2.set_title(sra_id + " \n(normal)")
    ax2.set_xlabel("Set of Barcodes")
    ax2.set_ylabel("UMI Counts")
    
    plt.grid(True, which="both")
    plt.savefig("../result_b/figure/kneeplot_{}.png".format(sra_id))
    print("Complete the plotting for: " + sra_id)
    print(f"{num_cells:,.0f} cells passed the {cutoff} UMI threshold for: " + sra_id + "\n")

mtx_filtered_dic = {}
for sra_id in sra_list:
    mtx = mtx_dic[sra_id]
    row_mask = np.asarray(mtx.sum(axis=1)>30).reshape(-1)
    col_mask = np.asarray(mtx.sum(axis=0)>0).reshape(-1)
    mtx_filtered_dic[sra_id] = mtx[row_mask,:][:,col_mask]
    print("\n~~~Filtering the cells x genes count by a threshold for: " + sra_id + "...~~~")
    print("Before filtering ... ")
    print(mtx_dic)
    print("\n")
    print("After filtering ... ")
    print(mtx_filtered_dic)