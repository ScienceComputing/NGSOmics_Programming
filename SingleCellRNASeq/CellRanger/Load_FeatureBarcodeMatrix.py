import scipy.io, csv, gzip, os, pandas as pd
print(os.getcwd()) # https://github.com/microsoft/vscode-python/issues/12173
mat_path = "./data/CTRL/"
mat = scipy.io.mmread(os.path.join(mat_path, "matrix.mtx.gz"))
print(mat)

# View the transcript ID, for example, ENSMUSG00000095523
feature_path = os.path.join(mat_path, "features.tsv.gz")
feature_id = [row[0] for row in csv.reader(gzip.open(feature_path, mode="rt"), delimiter="\t")] # The feature ID and name are stored in the first and second column of the unzipped features.tsv.gz file
# mode="rt" open the file for reading as a text file
print(feature_id)
# CMD approach: gzip -cd ./data/CTRL/features.tsv.gz | awk '{print $1}' | head -n 6 # Take a look at the first 6 rows of the transcript ID 

# View the gene symbol, for example, Pdzd8
gene_symbol = [row[1] for row in csv.reader(gzip.open(feature_path, mode="rt"), delimiter="\t")]
print(gene_symbol)
# CMD approach: gzip -cd ./data/CTRL/features.tsv.gz | awk '{print $2}' | head -n 6 # Take a look at the first 6 rows of the gene symbols

# View the feature types, for example, Gene Expression/Antibody Capture/CRISPR Guide Capture/Multiplexing Capture/CUSTOM
feature_type = [row[2] for row in csv.reader(gzip.open(feature_path, mode="rt"), delimiter="\t")]
print(feature_type)
# CMD approach: gzip -cd ./data/CTRL/features.tsv.gz | awk '{print $3}' | head -n 6 # Take a look at the first 6 rows of the feature types

# View the barcode, for example, TTTGTCATCGTTTAGG-1
# All detected barcodes vs cell-associated barcodes
barcode_path = os.path.join(mat_path, "barcodes.tsv.gz")
barcode = [row[0] for row in csv.reader(gzip.open(barcode_path, mode="rt"), delimiter="\t")]
print(barcode)

# Covert the sparse feature-barcode matrix to the pandas DataFrame and label rows and columns
# Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column)
count_mat = pd.DataFrame.sparse.from_spmatrix(mat)
count_mat.columns = barcode # Assign the new column names
count_mat.insert(loc=0, column="feature_id", value=feature_id) # new "feature_id" column will be inserted as the first column in the DataFrame
count_mat.insert(loc=0, column="gene", value=gene_symbol)
count_mat.insert(loc=0, column="feature_type", value=feature_type)

# This shows information on each of the columns, such as the data type and number of missing values
count_mat.info()
# <class 'pandas.core.frame.DataFrame'>
# RangeIndex: 27998 entries, 0 to 27997
# Columns: 737282 entries, gene to TTTGTCATCTTTCCTC-1
# dtypes: Sparse[int64, 0](737280), object(2)
# memory usage: 244.0+ MB
    
# View matrix
print(count_mat.head())
# Save the pandas DataFrame as a CSV
count_mat.to_csv("mex_matrix.csv", index=False)

