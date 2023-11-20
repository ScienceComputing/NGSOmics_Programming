# Convert a specific column (metadata) within the obs attribute of an AnnData object (adata) into a categorical data type
adata.obs[column_name] = adata.obs[column_name].astype('category')
