import argparse
import scanpy as sc
import numpy as np
import subprocess

print('*****************************************')
print('  Single-Cell RNA Sequencing Analysis Tool')
print('*****************************************\n')

# Set up the argument parser
parser = argparse.ArgumentParser(description="scRNA-seq Analysis Tool")

parser.add_argument("fastq", nargs=2, help="Path to the paired-end FASTQ files")

parser.add_argument("-s", "--sample_id", required=True, help="Sample ID for Cell Ranger")

parser.add_argument("-r", "--ref_genome", required=True, help="Path to the reference genome for Cell Ranger")

parser.add_argument("-n", "--normalization", type=str, default='Yes',
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
fastq_files=args.fastq
sample_id=args.sample_id
ref_genome=args.ref_genome
input_file=args.input
normalization_method=args.normalization
do_clustering=args.cluster
dim_reduction_technique=args.dim_reduction
gene_list=args.gene_list
output_file=args.output
num_threads=args.threads

# Process the scRNA-seq data using Scanpy workflow
def run_fastqc(fastq_files):
    print("Running FastQC...")
    subprocess.run(["fastqc", *fastq_files])

def run_trim_galore(fastq_files):
    print("Running Trim Galore...")
    trimmed_fastq_files = []
    for fastq in fastq_files:
        output_fastq = fastq.replace('.fastq', '_trimmed.fastq')
        subprocess.run(["trim_galore", "--paired", fastq, output_fastq])
        trimmed_fastq_files.append(output_fastq)
    return trimmed_fastq_files

def run_cellranger(trimmed_fastq_files, sample_id, ref_genome):
    print("Running Cell Ranger...")
    subprocess.run([
        "cellranger", "count", 
        "--id", sample_id, 
        "--transcriptome", ref_genome, 
        "--fastqs", ",".join(trimmed_fastq_files), 
        "--sample", sample_id
    ])

def process_count_data(input_file, normalization, clustering, dim_reduction, gene_list, output_file, num_threads):
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


# Call the functions with parsed arguments
run_fastqc(fastq_files=args.fastq)
trimmed_fastq = run_trim_galore(fastq_files=args.fastq)
run_cellranger(trimmed_fastq, 
               sample_id=args.sample_id, 
               ref_genome=args.ref_genome)
input_file_path=args.sample_id + "/outs"
process_count_data(input_file=input_file_path, 
                   normalization=args.normalization, 
                   clustering=args.cluster, 
                   dim_reduction=args.dim_reduction, 
                   gene_list=args.gene_list, 
                   output_file=args.output,
                   num_threads=args.threads)

# In the shell
# python pipeline.py PBMC_R1.fastq PBMC_R2.fastq -s PBMC_count -r ../data/refdata-gex-GRCh38-2020-A -n Yes -d UMAP -o output.h5ad -t 30
