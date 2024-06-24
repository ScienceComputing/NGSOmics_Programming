# Reference: https://github.com/ScienceComputing/NGSOmics_Programming/tree/main/SingleCellRNASeq/Anndata/Combine_Data.py

import anndata as ad
data_s1 = ad.read_h5ad('output_s1/counts_unfiltered/adata.h5ad')
data_s2 = ad.read_h5ad('output_s2/counts_unfiltered/adata.h5ad')
data_s1
data_s2
data_s1.X
data_s2.X
data_s1.obs.head() # Take a look at first few rows
data_s2.obs.head()
data_s1.var.head() # Take a look at first few columns
data_s2.var.head()

data_s12 = ad.concat({"data_s1": data_s1, "data_s2": data_s2}, join='outer', index_unique='-', label="data_origin")
# The shared names of obs can be made unique by appending the relevant key using the index_unique argument

data_s12.obs.head()
#                          data_origin
# barcode
# AAACCTGAGAAACCAT-data_s1     data_s1
# AAACCTGAGAAACCGC-data_s1     data_s1
# AAACCTGAGAAACCTA-data_s1     data_s1
# AAACCTGAGAAAGTGG-data_s1     data_s1
# AAACCTGAGAACAATC-data_s1     data_s1
