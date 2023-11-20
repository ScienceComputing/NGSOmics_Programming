# Before conducting differential gene expression (DGE) analysis, it's essential to aggregate the expression values specific to each cell type within an individual. 
# This aggregation can be achieved by using either a sum, mean, or random effect for each individual. 
# This process, known as pseudobulk generation, is critical for addressing correlations that exist within a sample.

import warnings
warnings.filterwarnings("ignore") # Suppress warnings
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import random
import sc_toolbox
import logging
sc.settings.verbosity = 0 # Configure scanpy to not print any output messages

def count_entries(df, *args):
    """Return a dictionary with counts of occurrences as value for each key."""
    cols_count = {}
    for col_name in args:
        col = df[col_name]
        for entry in col:
            if entry in cols_count.keys():
                cols_count[entry] += 1
            else:
                cols_count[entry] = 1
    return cols_count

"""
Prepare the dataset
"""
# This is a preprocessed data that filter cells with less than 200 genes and genes found in less than 3 cells
adata = ad.read_h5ad('adata_v2.h5ad')
adata
adata[0:100,0:100].to_df()
np.max(adata.X)
np.min(adata.X)
adata.obs
adata.layers['counts'] = adata.X.copy()
adata
# Treat sm_name as a categorical variable - cat, 
# get the distinct categories - categories, 
# and calculate the number of categories - len
len(adata.obs.sm_name.cat.categories)
adata.obs['condition'] = ['perturbated' if sm_i == 'Belinostat' else 'control' for sm_i in adata.obs['sm_name']]

# Convert the condition into categories
adata.obs.condition = adata.obs.condition.astype('category')
adata.obs.condition
adata.obs.condition.cat.categories
count_entries(adata.obs, 'condition')
count_entries(adata.obs, 'plate_name')
df = adata.obs
crosstab = pd.crosstab(df['condition'], df['plate_name'], margins=False)
print(crosstab)
con_1 = (df['condition'] == 'perturbated') & (df['plate_name'].isin(['plate_2', 'plate_3', 'plate_5']))
con_2 = (df['condition'] == 'control') & (df['plate_name'].isin(['plate_0', 'plate_1', 'plate_4']))
bdata = adata[con_1 | con_2]

df_2 = bdata.obs
pd.crosstab(df_2['condition'], df_2['plate_name'], margins=False)
pd.crosstab(df_2['cell_type'], df_2['plate_name'], margins=False)

cdata = bdata[bdata.obs.cell_type.isin(['B cells', 'NK cells', 'T cells CD4+'])]
df_3 = cdata.obs
pd.crosstab(df_3['cell_type'], df_3['plate_name'], margins=False)
pd.crosstab(df_3['condition'], df_3['plate_name'], margins=False)

"""
Create the pseudobulk from the single-cell RNA dataset
"""
# For each plate, one pseudobulk sample per cell type is created. This is done by aggregating the cell of each type using the summation of gene expression within that cell type.
cdata.obs['sample'] = [
    f'{p}_{c}' for p, c in zip(cdata.obs.plate_name, cdata.obs.condition)
]

cdata.obs['sample']

# Replace spaces with _; remove +
cdata.obs['cell_type'] = [ct.replace(" ", "_") for ct in cdata.obs['cell_type']]
cdata.obs['cell_type'] = [ct.replace("+", "") for ct in cdata.obs['cell_type']]
count_entries(cdata.obs, 'cell_type')

# Force the following variables to be categorical
cat_var = ['plate_name', 'condition', 'sample', 'cell_type']
cdata.obs[cat_var] = cdata.obs[cat_var].astype('category')

# Original function
NUM_OF_CELL_PER_DONOR = 30 # Filter out donors that have fewer than 30 cells for the specified population

def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key='sample',
    condition_key='condition',
    cell_identity_key='cell_type',
    obs_to_keep=[], # Which additional metadata to keep, e.g. gender, age, etc.
    replicates_per_patient=1,
):
    """
    Aggregate and filter single-cell data from an AnnData object based on specified cell identity, 
    with additional donor-based filtering and data aggregation.

    This function filters cells based on the specified cell identity, then aggregates data for each donor,
    and further divides the data into specified numbers of replicates per donor. It then aggregates 
    this data (summing gene expression) and retains specified metadata.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object containing the single-cell dataset.
    cell_identity : str
        The cell identity to filter on. Only cells with this identity are considered.
    donor_key : str, optional
        The key in adata.obs used to identify different donors. Defaults to 'sample'.
    condition_key : str, optional
        The key in adata.obs used to identify different conditions. Defaults to 'condition'.
    cell_identity_key : str, optional
        The key in adata.obs used to identify cell types. Defaults to 'cell_type'.
    obs_to_keep : list of str, optional
        A list of additional metadata fields to keep from adata.obs. Defaults to an empty list.
    replicates_per_patient : int, optional
        The number of replicates to create for each donor. Defaults to 1.

    Returns:
    --------
    AnnData
        An aggregated AnnData object with summed gene expression and retained metadata, 
        filtered for specified cell identity and donors meeting the cell count criteria.

    Notes:
    ------
    - The function filters out donors with fewer than 30 cells for the specified cell identity.
    - It assumes that the entire filtered dataset fits into memory.
    - Random shuffling is used for creating replicates, which introduces a random element to the output.

    Example:
    --------
    aggregated_data = aggregate_and_filter(adata, 'T-cell', obs_to_keep=['age', 'gender'])
    """
    
    # Filter the adata object to include only cells of a specified type 
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # Calculate the number of cells per donor 
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    # Filter out donors with fewer cells than NUM_OF_CELL_PER_DONOR
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print('Drop the following samples:')
        print(donors_to_drop)

    # Initialize a DataFrame df to store aggregated data
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])
    # Categorize donor_key
    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype('category')
    # Iterate over each donor (which is represented as a categorical level)
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f'\tProcess donor {i+1} out of {len(donors)}...', end='\r')

        # Check if the donor is not in the filtered list
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            # Split the data for each donor into replicates
            indices = list(adata_donor.obs_names) # Extract the unique identifiers of each cell in an AnnData object (adata_donor) and store them in a Python list
            random.shuffle(indices) # Randomly shuffle (reorder) the elements in the indices list in-place 
            indices = np.array_split(np.array(indices), replicates_per_patient) # Split an array (or list) of indices into multiple subarrays or sublists, where the number of subarrays is determined by the replicates_per_patient variable

            # For each replicate, it aggregates data (sums gene expression and keeps additional metadata)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx] # rep_idx is a sublist of cell identifiers for a specific replicate; if replicate per donor is 1, then rep_idx is a whole list of cell identifiers for that donor
                agg_dict = {gene: 'sum' for gene in adata_replicate.var_names} # The 'sum' specifies the aggregation method to compute the sum of gene expression values when aggregating data for each gene within the adata_replicate 
                for obs in obs_to_keep:
                    agg_dict[obs] = 'first' # For the metadata variables associated with each replicate, it should take the first value encountered
                # Build a df with all genes, donor and metadata
                df_donor = pd.DataFrame(adata_replicate.X.A) # adata_replicate.X typically contains the matrix of gene expression values, and .A is used to access the underlying NumPy array representation of the data
                df_donor.index = adata_replicate.obs_names # The row (index) labels of the df_donor DataFrame are set to the observation names from adata_replicate
                df_donor.columns = adata_replicate.var_names # Set the column names of the df_donor DataFrame to the gene names (variable names) from adata_replicate
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep]) # Add additional metadata columns to the df_donor DataFrame 
                # Aggregate data within a specific donor 
                df_donor = df_donor.groupby(donor_key).agg(agg_dict) # Apply the specified aggregation methods (summation) to each replicate belonging to the same donor
                df_donor[donor_key] = donor # A new column named donor_key is added to the df_donor DataFrame, and its value is set to the donor identifier
                df.loc[f"donor_{donor}_{i}"] = df_donor.loc[donor] # It assigns the aggregated data for the current donor (identified by donor) to a new row in df. The row label is constructed using an f-string, which includes the donor identifier (donor) and a unique identifier (i) for the current replicate
    print('\n')
    # Convert the aggregated DataFrame back into an AnnData object for further analysis
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names) # It assigns the columns from the DataFrame df that are not included in adata_cell_pop.var_names as the observation metadata. In other words, it includes any metadata columns from df that are not related to gene expression
    )
    return adata_cell_pop

obs_to_keep = ['condition', 'cell_type', 'sample']
cdata.X = cdata.layers["counts"].copy()
cdata.obs['cell_type']

cell_categories = cdata.obs["cell_type"].cat.categories
cdata_pb = None
for i, cell_type in enumerate(cell_categories):
    print(
        f'Processing {cell_type} ({i+1} out of {len(cell_categories)})...'
    )
    cdata_cell_type = aggregate_and_filter(cdata, cell_type, obs_to_keep=obs_to_keep)
    if cdata_pb is None:
        cdata_pb = cdata_cell_type
    else:
        cdata_pb = cdata_pb.concatenate(cdata_cell_type)

# View the pseudobulk
cdata_pb.to_df()
