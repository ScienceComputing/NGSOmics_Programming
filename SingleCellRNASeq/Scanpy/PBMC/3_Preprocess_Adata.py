import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

outcome_path = '../scRNA_data/preprocessed_data/'
adata_pbmc = ad.read_h5ad(outcome_path + 'pbmc_v0.h5ad')

# Display the genes that contribute the largest proportion of total counts within each individual cell, spanning all the cells.
sc.pl.highest_expr_genes(adata_pbmc, n_top=30)

# Let's perform the basic quality control
sc.pp.filter_cells(adata_pbmc, min_genes=200) # We set 200 as the minimum number of genes expressed required for a barcode to pass filtering
sc.pp.filter_genes(adata_pbmc, min_cells=3) # We set 3 as the minimum number of barcodes expressed required for a gene to pass filtering

# Compute the distribution of mitochondrial genes
adata_pbmc.var['mt'] = adata_pbmc.var_names.str.startswith('MT-') 
# adata_pbmc = sc.read_10x_mtx(..., var_names='gene_symbols', ...)
# for mouse datasets, the prefix of mitochondrial genes is typically lower case: mt-
sc.pp.calculate_qc_metrics(adata_pbmc, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) 
# inplace: whether to place calculated metrics in adata’s .obs and .var
adata_pbmc.var
#                       gene_ids  n_cells     mt  n_cells_by_counts  mean_counts  pct_dropout_by_counts  total_counts
# AL627309.1     ENSG00000237683        9  False                  9     0.003333              99.666667           9.0
# AP006222.2     ENSG00000228463        3  False                  3     0.001111              99.888889           3.0
# RP11-206L10.2  ENSG00000228327        5  False                  5     0.001852              99.814815           5.0
# RP11-206L10.9  ENSG00000237491        3  False                  3     0.001111              99.888889           3.0
# LINC00115      ENSG00000225880       18  False                 18     0.006667              99.333333          18.0
# ...                        ...      ...    ...                ...          ...                    ...           ...
# AC145212.1     ENSG00000215750       16  False                 16     0.006667              99.407407          18.0
# AL592183.1     ENSG00000220023      323  False                323     0.134815              88.037037         364.0
# AL354822.1     ENSG00000215615        8  False                  8     0.002963              99.703704           8.0
# PNRC2-1        ENSG00000215700      110  False                110     0.042963              95.925926         116.0
# SRSF10-1       ENSG00000215699       69  False                 69     0.025926              97.444444          70.0

adata_pbmc.obs
#                   n_genes  n_genes_by_counts  total_counts  total_counts_mt  pct_counts_mt
# AAACATACAACCAC-1      781                779        2419.0             73.0       3.017776
# AAACATTGAGCTAC-1     1352               1352        4903.0            186.0       3.793596
# AAACATTGATCAGC-1     1131               1129        3147.0             28.0       0.889736
# AAACCGTGCTTCCG-1      960                960        2639.0             46.0       1.743085
# AAACCGTGTATGCG-1      522                521         980.0             12.0       1.224490
# ...                   ...                ...           ...              ...            ...
# TTTCGAACTCTCAT-1     1155               1153        3459.0             73.0       2.110436
# TTTCTACTGAGGCA-1     1227               1224        3443.0             32.0       0.929422
# TTTCTACTTCCTCG-1      622                622        1684.0             37.0       2.197150
# TTTGCATGAGAGGC-1      454                452        1022.0             21.0       2.054795
# TTTGCATGCCTCAC-1      724                723        1984.0             16.0       0.806452

# Use the violin plot/scatter plot to visualize the following quality control measures:
# (1) the number of genes expressed with at least 1 count per barcode; 
# (2) the total counts per barcode; 
# (3) the percentage of counts in mitochondrial genes among total counts per barcode
sc.pl.violin(adata_pbmc, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.5, multi_panel=True)


sc.pl.scatter(adata_pbmc, x='total_counts', y='pct_counts_mt') # (2) vs (3)
sc.pl.scatter(adata_pbmc, x='total_counts', y='n_genes_by_counts') # (2) vs (1)

# Filter out cells that express too many mitochondrial genes or too high total counts/numbers of genes expressed with at least 1 count
adata_pbmc = adata_pbmc[adata_pbmc.obs.pct_counts_mt < 5, :]
adata_pbmc = adata_pbmc[adata_pbmc.obs.n_genes_by_counts < 2500, :]

# Normalize the total count to 10000 reads per barcode
sc.pp.normalize_total(adata_pbmc, target_sum=1e4)
adata_pbmc.obs["total_counts"].describe()
# adata_pbmc.obs.describe()
adata_pbmc[:,:].to_df()
#                   AL627309.1  AP006222.2  RP11-206L10.2  RP11-206L10.9  ...  AL592183.1  AL354822.1  PNRC2-1  SRSF10-1
# AAACATACAACCAC-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# AAACATTGAGCTAC-1         0.0         0.0            0.0            0.0  ...    2.039568    0.000000      0.0       0.0
# AAACATTGATCAGC-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# AAACCGTGCTTCCG-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# AAACCGTGTATGCG-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# ...                      ...         ...            ...            ...  ...         ...         ...      ...       ...
# TTTCGAACTCTCAT-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# TTTCTACTGAGGCA-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# TTTCTACTTCCTCG-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# TTTGCATGAGAGGC-1         0.0         0.0            0.0            0.0  ...    0.000000    9.784736      0.0       0.0
# TTTGCATGCCTCAC-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0

# Do the natural logarithm (log(x+1)) on the count in the barcode x gene count matrix
sc.pp.log1p(adata_pbmc)
adata_pbmc[:,:].to_df()
#                   AL627309.1  AP006222.2  RP11-206L10.2  RP11-206L10.9  ...  AL592183.1  AL354822.1  PNRC2-1  SRSF10-1
# AAACATACAACCAC-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# AAACATTGAGCTAC-1         0.0         0.0            0.0            0.0  ...    1.111715    0.000000      0.0       0.0
# AAACATTGATCAGC-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# AAACCGTGCTTCCG-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# AAACCGTGTATGCG-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# ...                      ...         ...            ...            ...  ...         ...         ...      ...       ...
# TTTCGAACTCTCAT-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# TTTCTACTGAGGCA-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# TTTCTACTTCCTCG-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0
# TTTGCATGAGAGGC-1         0.0         0.0            0.0            0.0  ...    0.000000    2.378132      0.0       0.0
# TTTGCATGCCTCAC-1         0.0         0.0            0.0            0.0  ...    0.000000    0.000000      0.0       0.0

# Identify the highly variable genes
adata_pbmc[:,:].to_df().describe()
#         AL627309.1   AP006222.2  ...      PNRC2-1     SRSF10-1
# count  2638.000000  2638.000000  ...  2638.000000  2638.000000
# mean      0.005436     0.001847  ...     0.066251     0.040642
# std       0.093517     0.055012  ...     0.324988     0.255053
# min       0.000000     0.000000  ...     0.000000     0.000000
# 25%       0.000000     0.000000  ...     0.000000     0.000000
# 50%       0.000000     0.000000  ...     0.000000     0.000000
# 75%       0.000000     0.000000  ...     0.000000     0.000000
# max       1.882199     1.828519  ...     2.867758     2.498723
sc.pp.highly_variable_genes(adata_pbmc, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_pbmc)

# Attach the .raw attribute of adata_pbmc to the normalized and logarithmized raw count for differential testing and visualizations of gene expression
adata_pbmc.raw = adata_pbmc
adata_pbmc.raw
# <anndata._core.raw.Raw object at 0x14df5c400>
# adata_pbmc.raw.to_adata() # Return to the AnnData
adata_pbmc.shape
# (2638, 13714)
adata_pbmc.write(outcome_path + 'pbmc_raw.h5ad')

# Select highly variable genes
adata_pbmc_filter = adata_pbmc[:, adata_pbmc.var.highly_variable]
adata_pbmc_filter.shape
# (2638, 1838) # So we remove (13714 - 1838) non highly variable genes

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
sc.pp.regress_out(adata_pbmc_filter, ['total_counts', 'pct_counts_mt'])

# Normalize each gene's variance to a unit scale, and truncate any values that surpass 10 standard deviations.
sc.pp.scale(adata_pbmc_filter, max_value=10)
adata_pbmc_filter.write(outcome_path + 'pbmc_filter.h5ad')
