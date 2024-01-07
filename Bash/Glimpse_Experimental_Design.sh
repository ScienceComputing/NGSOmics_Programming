# Print the meta information excluding the cell type of treatment1
grep -v treatment1 ~/scRNAseq/cell_metadata.txt 

# How many cell types in this meta dataset?
cut -f 2 ~/scRNAseq/cell_metadata.txt | sort -ufb | wc -l
