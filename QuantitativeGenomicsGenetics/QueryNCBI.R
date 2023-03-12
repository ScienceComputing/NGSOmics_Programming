# BiocManager::install("KEGGREST")

library(KEGGREST)
listDatabases()
# [1] "pathway"  "brite"    "module"   "ko"       "genome"  
# [6] "vg"       "ag"       "compound" "glycan"   "reaction"
# [11] "rclass"   "enzyme"   "disease"  "drug"     "dgroup"  
# [16] "environ"  "genes"    "ligand"   "kegg"   

# Find immune-related genes/drugs
# Ref: https://www.kegg.jp/kegg/docs/keggapi.html
keggFind("genes", "immune")
keggFind("drug", "immune")

# Get the identifiers of the genes/drugs
keggFind("genes", "immune") |> names()
immune.gene.id <- keggFind("genes", "immune") |> names()
keggFind("drug", "immune") |> names()

# Find all ligands that interact with ATP
keggFind("ligand", "ATP")

# Find immune-related pathways
path.id <- keggFind("pathway", "immune") |> names()

# Check the genes involving in each immune-related pathway in **humans (hsa)**
for (id.i in path.id){
  id.new <- gsub("path:map", "hsa", id.i)
  res <- keggGet(id.new)
  cat("Pathway Name", res[[1]]$NAME, "\n") 
  print(res[[1]]$GENE[1:10])
}

# Find more details on immune-related genes
immune.gene <- keggGet(immune.gene.id)

# Let us see more information stored in the `immune.gene` object
names(immune.gene[[1]])
# [1] "ENTRY"       "NAME"        "PRODUCT"     "FORMULA"    
# [5] "EXACT_MASS"  "MOL_WEIGHT"  "SEQUENCE"    "SOURCE"     
# [9] "CLASS"       "REMARK"      "EFFICACY"    "COMMENT"    
# [13] "TARGET"      "METABOLISM"  "INTERACTION" "STR_MAP"    
# [17] "OTHER_MAP"   "BRITE"       "DBLINKS"     "ATOM"       
# [21] "BOND"  

# View the sequence of the first immune-related gene
# keggGet(dbentries = immune.gene.id[1], option = "sequence")
immune.gene[[1]]$SEQUENCE

# Generate a more complex search term
paste("genomics", "proteomics", "transcriptomics", sep = " AND ")

# More info: https://www.bioconductor.org/packages/devel/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html
