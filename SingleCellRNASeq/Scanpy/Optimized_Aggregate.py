# TD: Optimize the function aggregate_and_filter
# Avoid Copying Data: The copy() operation is removed to save memory.
# Efficient Donor Filtering: The donor filtering is streamlined with vectorized operations.
# DataFrame Initialization: The DataFrame is pre-initialized with the necessary columns.
# Loop Over Donors: The loop iterates directly over donors_to_keep.
# Aggregation: The aggregation is done using numpy operations, which are generally faster.
# Data Appending: Data is appended row-wise to the DataFrame, which is more efficient than building a DataFrame in a loop.
