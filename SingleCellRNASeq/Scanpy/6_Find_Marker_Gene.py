import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

outcome_path = '../scRNA_data/preprocessed_data/'
adata_pbmc_leiden = ad.read_h5ad(outcome_path + 'pbmc_leiden.h5ad')
adata_pbmc_leiden.uns['log1p']["base"] = None
adata_pbmc_leiden.uns['log1p']
# if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
# KeyError: 'base'
# https://github.com/scverse/scanpy/issues/2239

# Let's compute the ranking for the highly differential genes between each cluster and the rest of clusters 
# By default, the .raw attribute of AnnData is used if it has been initialized
# Approach 1: use t test to examine the differential genes
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', method='t-test') # groups=clusters
sc.pl.rank_genes_groups(adata_pbmc_leiden, n_genes=25, sharey=False)

# Approach 2: use wilcoxon rank-sum test to examine the differential genes between each cluster and the rest of clusters 
sc.settings.verbosity = 2 
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', method='wilcoxon') # Recommend wilcoxon rank-sum test over t test
sc.pl.rank_genes_groups(adata_pbmc_leiden, n_genes=25, sharey=False)
# Consider MAST, limma, DESeq2, and, diffxpy

adata_pbmc_leiden.write(outcome_path + 'pbmc_wilcoxon.h5ad')

# Approach 3: use logistic regression to examine the differential genes between each cluster and the rest of clusters 
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata_pbmc_leiden, n_genes=25, sharey=False)

# Find the overlapping marker genes detected by all approaches
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 
                'LYZ', 'CD14', 'LGALS3', 'S100A8', 'GNLY', 
                'NKG7', 'KLRB1', 'FCGR3A', 'MS4A7', 'FCER1A', 
                'CST3', 'PPBP']

# Show the names of top 10 ranked genes per cluster 0, 1, …, 12 in a dataframe
adata_pbmc_wilcoxon = sc.read(outcome_path + 'pbmc_wilcoxon.h5ad')
pd.DataFrame(adata_pbmc_wilcoxon.uns['rank_genes_groups']['names']).head(10)
#         0         1       2      3       4       5       6       7      8      9      10        11     12
# 0     LYZ      CD74   RPS12  RPS12    LDHB    IL32    LST1    NKG7   CCL5   CCL5     LTB  HLA-DPA1    PF4
# 1  S100A9     CD79A    RPS6  RPS3A     LTB    CD3D  FCER1G    GNLY   NKG7   GZMK    LDHB  HLA-DPB1   SDPR
# 2  S100A8   HLA-DRA   RPS14   RPS6   RPS25     B2M    AIF1    GZMB   CST7   NKG7    RPSA   HLA-DRA  GNG11
# 3  TYROBP     CD79B    RPS3  RPL32   RPS27     LTB   COTL1    CTSW    B2M   IL32    IL7R  HLA-DRB1   PPBP
# 4    FCN1  HLA-DPB1   RPS25  RPL13   RPL30   HLA-A    FTH1    PRF1  HLA-C   GZMA   RPS12  HLA-DQA1   NRGN
# 5     FTL  HLA-DQA1   RPL32  RPS14  RPS27A    LDHB  IFITM3    GZMA   GZMA   CTSW  EEF1A1      CD74   GPX1
# 6    CST3     MS4A1   RPS27   RPL3    JUNB    CD3E    SAT1    CST7   GZMH    B2M    CD3D      CST3  SPARC
# 7  S100A6      CD37   RPL31   RPL9   RPS29     CD2     FTL   HLA-C   CTSW   LYAR    RPL3    FCER1A   TPM4
# 8  LGALS2  HLA-DQB1   RPL13   RPS3  RPS15A   HINT1    PSAP     B2M  HLA-A   IL7R    TPT1  HLA-DRB5  RGS18
# 9    FTH1  HLA-DRB1  RPL13A  RPLP2    TPT1  MYL12A  FCGR3A  FGFBP2   CD3D  HLA-C    AQP3   HLA-DMA  CALM3

# Show the names and p values of top 10 ranked genes per cluster 0, 1, …, 12 in a dataframe
marker_gene_result = adata_pbmc_wilcoxon.uns['rank_genes_groups']
clusters = marker_gene_result['names'].dtype.names
# marker_gene_result['names']
pd.DataFrame(
    {cluster_i + '_' + key[:1]: marker_gene_result[key][cluster_i]
    for cluster_i in clusters for key in ['names', 'pvals']}).head(10)
marker_gene_result['names'] 
# array([('LYZ', 'CD74', 'RPS12', 'RPS12', 'LDHB', 'IL32', 'LST1', 'NKG7', 'CCL5', 'CCL5', 'LTB', 'HLA-DPA1', 'PF4'),
#        ('S100A9', 'CD79A', 'RPS6', 'RPS3A', 'LTB', 'CD3D', 'FCER1G', 'GNLY', 'NKG7', 'GZMK', 'LDHB', 'HLA-DPB1', 'SDPR'),
#        ('S100A8', 'HLA-DRA', 'RPS14', 'RPS6', 'RPS25', 'B2M', 'AIF1', 'GZMB', 'CST7', 'NKG7', 'RPSA', 'HLA-DRA', 'GNG11'),
#        ...,
#        ('RPS27', 'S100A6', 'CYBA', 'S100A4', 'HLA-DRA', 'HLA-DRB1', 'RPL13A', 'RPL28', 'RPLP1', 'CTSS', 'CD74', 'RPL21', 'RPL11'),
#        ('RPS27A', 'TMSB4X', 'HLA-DRA', 'CD74', 'CD74', 'HLA-DRA', 'RPL13', 'RPL18A', 'FOS', 'HLA-DRA', 'HLA-DRB1', 'RPS27', 'MALAT1'),
#        ('MALAT1', 'S100A4', 'CD74', 'CYBA', 'CYBA', 'TYROBP', 'RPL3', 'RPL32', 'LTB', 'TMSB10', 'FTL', 'MALAT1', 'RPL10')],
#       dtype=[('0', 'O'), ('1', 'O'), ('2', 'O'), ('3', 'O'), ('4', 'O'), ('5', 'O'), ('6', 'O'), ('7', 'O'), ('8', 'O'), ('9', 'O'), ('10', 'O'), ('11', 'O'), ('12', 'O')])
marker_gene_result['pvals']
# array([(1.25611944e-236, 4.25455920e-184, 2.71837261e-53, 4.55104827e-78, 3.49421355e-40, 1.44724928e-37, 6.63436256e-102, 2.08391191e-94, 3.04081566e-77, 2.73479191e-58, 5.16909620e-18, 1.77755156e-20, 4.72288553e-10),
#        (7.42159504e-236, 1.10313255e-171, 5.01050757e-53, 1.73561469e-72, 1.68019515e-32, 8.08588915e-30, 3.23239121e-098, 2.10541489e-90, 1.23109815e-70, 5.30654399e-55, 8.54755604e-16, 2.09727623e-20, 4.73389941e-10),
#        (6.99633918e-230, 1.32339233e-167, 9.84567145e-48, 1.29939068e-71, 2.47067938e-28, 2.67230101e-27, 7.80036583e-098, 3.09593235e-86, 4.56539804e-62, 8.18241509e-41, 2.21544352e-14, 4.24385539e-19, 4.73389941e-10),
#        ...,
#        (6.96403538e-141, 4.80493563e-083, 4.46765523e-34, 1.01324255e-43, 2.55486113e-27, 1.05275323e-13, 1.94026931e-064, 9.36473849e-52, 4.21431478e-17, 5.63287405e-12, 2.21434405e-08, 1.78459332e-08, 1.07832197e-09),
#        (9.21984913e-156, 1.81088101e-095, 2.92769162e-41, 2.95414983e-47, 2.70687326e-28, 7.29194735e-16, 5.43221156e-066, 8.85632098e-56, 1.99666906e-20, 1.68154041e-15, 5.72293847e-09, 1.99398075e-09, 8.91801689e-10),
#        (3.38925469e-189, 5.32504861e-110, 6.37439763e-47, 2.07536660e-60, 7.45014635e-36, 6.21686238e-16, 1.66082644e-067, 8.47132409e-59, 3.89698090e-26, 7.95971981e-31, 3.29421043e-10, 8.10599435e-13, 5.68796812e-10)],
#       dtype=[('0', '<f8'), ('1', '<f8'), ('2', '<f8'), ('3', '<f8'), ('4', '<f8'), ('5', '<f8'), ('6', '<f8'), ('7', '<f8'), ('8', '<f8'), ('9', '<f8'), ('10', '<f8'), ('11', '<f8'), ('12', '<f8')])

# Use wilcoxon rank-sum test to examine the differential genes between a cluster (e.g., cluster 0) and another single cluster (e.g., cluster 1) 
sc.tl.rank_genes_groups(adata_pbmc_leiden, 'leiden', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata_pbmc_leiden, groups=['0'], n_genes=25)

# Use the violin plot to view the gene expression levels of top 10 ranked differentially expressed genes between the cluster 0 and the cluster 1
sc.pl.rank_genes_groups_violin(adata_pbmc_leiden, groups='0', n_genes=10)

# Use the violin plot to view the gene expression levels of top 10 ranked differentially expressed genes between the cluster 0 and the rest of clusters
adata_pbmc_wilcoxon = sc.read(outcome_path + 'pbmc_wilcoxon.h5ad')
sc.pl.rank_genes_groups_violin(adata_pbmc_wilcoxon, groups='0', n_genes=10)

# Use the violin plot to view the gene expression level of a particular gene among clusters
sc.pl.violin(adata_pbmc_wilcoxon, ['LYZ', 'S100A9', 'S100A8'], groupby='leiden')

# Mark the cell types
