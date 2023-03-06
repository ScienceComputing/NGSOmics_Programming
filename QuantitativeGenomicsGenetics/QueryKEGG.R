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
keggFind("drug", "immune") |> names()

# Find all ligands that interact with ATP
keggFind("ligand", "ATP")
