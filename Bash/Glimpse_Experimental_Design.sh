# Print the meta information excluding the cell type of treatment1
grep -v treatment1 ~/scRNAseq/cell_metadata.txt 

file=~/scRNAseq/cell_metadata.txt
grep -v treatment1 $file

# How many cell types are present in this meta dataset?
cut -f 2 ~/scRNAseq/cell_metadata.txt | sort -ufb | wc -l
