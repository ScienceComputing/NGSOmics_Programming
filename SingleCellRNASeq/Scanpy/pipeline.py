import argparse
import scanpy as sc
import numpy as np

print('*****************************************')
print('  Single-Cell RNA Sequencing Analysis Tool')
print('*****************************************\n')

# Creating an argument parser
parser = argparse.ArgumentParser(description="scRNA-seq Analysis Tool")

# Adding arguments
parser.add_argument("-i", "--input", required=True,
                    help="Path to the input file (e.g., a .h5ad file containing scRNA-seq data)")

parser.add_argument("-n", "--normalization", type=str, default='CPM',
                    choices=['Yes', 'None'],
                    help="Use the normalization or not")

parser.add_argument("-c", "--cluster", action='store_true',
                    help="Perform clustering analysis")

parser.add_argument("-d", "--dim_reduction", type=str, default='PCA',
                    choices=['PCA', 'tSNE', 'UMAP'],
                    help="Dimensionality reduction technique (PCA, tSNE, UMAP)")

parser.add_argument("-g", "--gene_list", nargs='*',
                    help="List of genes to include in the analysis (optional)")

parser.add_argument("-o", "--output", default='output.h5ad',
                    help="Output file name (e.g., a .h5ad file)")

parser.add_argument("-t", "--threads", type=int, default=1,
                    help="Number of threads to use for parallel processing (default: 1)")

# Parse arguments
args = parser.parse_args()

# Use args to access the arguments
input_file = args.input
normalization_method = args.normalization
do_clustering = args.cluster
dim_reduction_technique = args.dim_reduction
gene_list = args.gene_list
output_file = args.output
num_threads = args.threads

# Process the scRNA-seq data using Scanpy workflow
def process_scrna_data(input_file, normalization, clustering, dim_reduction, gene_list, output_file, num_threads):
    # Load the data
    print("Loading data")
    adata = sc.read(input_file)

    # Basic Preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Normalization
    if normalization != 'None':
        print("Normalizing data")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Selecting genes if provided
    if gene_list:
        adata = adata[:, gene_list]

    # Dimensionality Reduction
    print("Running PCA")
    sc.tl.pca(adata, svd_solver='arpack')

    if dim_reduction == 'tSNE':
        print("Running tSNE")
        sc.tl.tsne(adata, n_jobs=num_threads)
    elif dim_reduction == 'UMAP':
        print("Running UMAP")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)

    # Clustering
    if clustering:
        print("Performing clustering...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.louvain(adata)

    # Saving the result
    print("Saving results")
    adata.write(output_file)

    # Optionally, we can output plots
    if dim_reduction == 'tSNE':
        sc.pl.tsne(adata, save='tsne_plot.png')
    elif dim_reduction == 'UMAP':
        sc.pl.umap(adata, save='umap_plot.png')

# Call the processing function with parsed arguments
process_scrna_data(input_file=args.input, 
                   normalization=args.normalization, 
                   clustering=args.cluster, 
                   dim_reduction=args.dim_reduction, 
                   gene_list=args.gene_list, 
                   output_file=args.output,
                   num_threads=args.threads)

# In the shell
# python3 pipeline.py -i scRNA.h5ad -n Yes -d UMAP -o output.h5ad -t 30
