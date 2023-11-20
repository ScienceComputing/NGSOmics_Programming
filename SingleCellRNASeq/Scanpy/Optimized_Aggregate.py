# Optimized the function aggregate_and_filter
# Avoid Copying Data: The copy() operation is removed to save memory.
# Efficient Donor Filtering: The donor filtering is streamlined with vectorized operations.
# DataFrame Initialization: The DataFrame is pre-initialized with the necessary columns.
# Loop Over Donors: The loop iterates directly over donors_to_keep.
# Aggregation: The aggregation is done using numpy operations, which are generally faster.
# Data Appending: Data is appended row-wise to the DataFrame, which is more efficient than building a DataFrame in a loop.

NUM_OF_CELL_PER_DONOR = 30  # Filter out plates that have fewer than 30 cells for the specified population

def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key='sample',
    condition_key='condition',
    cell_identity_key='cell_type',
    obs_to_keep=[],  # Which additional metadata to keep, e.g. gender, age, etc.
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
    
    # Filter by cell identity without copying the data
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity]

    # Efficient donor filtering
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_keep = size_by_donor[size_by_donor > NUM_OF_CELL_PER_DONOR].index

    # Filter the adata to include only the donors to keep
    adata_cell_pop = adata_cell_pop[adata_cell_pop.obs[donor_key].isin(donors_to_keep)]

    # Initialize DataFrame
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep, donor_key])

    # Processing each donor
    for donor in donors_to_keep:
        print(f'\tProcessing donor {donor}...', end='\r')
        adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]

        # Split data for each donor into replicates
        indices = list(adata_donor.obs_names)
        random.shuffle(indices)
        indices = np.array_split(np.array(indices), replicates_per_patient)

        for rep_idx in indices:
            adata_replicate = adata_donor[rep_idx]
            # Aggregate data
            agg_data = np.sum(adata_replicate.X, axis=0)  # Sum gene expression
            additional_data = adata_replicate.obs[obs_to_keep].iloc[0].to_list()  # Get additional data

            # Append to DataFrame
            df_row = np.append(agg_data, additional_data)
            df_row = np.append(df_row, donor)
            df.loc[len(df)] = df_row

    print('\nDone.')

    # Create AnnData object from the df
    adata_aggregated = sc.AnnData(
        X=df[adata_cell_pop.var_names].values, 
        obs=pd.DataFrame(df[obs_to_keep + [donor_key]].values, columns=obs_to_keep + [donor_key])
    )
    return adata_aggregated
