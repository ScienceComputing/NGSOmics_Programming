library(biomaRt)

# Define a function to convert Ensembl IDs to NIH Gene IDs
ensembl_to_nih <- function(ensembl_ids) {
  # Connect to the Ensembl BioMart database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Retrieve the NIH Gene IDs corresponding to the input Ensembl IDs
  results <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Convert the results to a named vector for easy lookup
  nih_ids <- as.character(results$entrezgene_id)
  names(nih_ids) <- results$ensembl_gene_id
  
  # Return the NIH Gene IDs corresponding to the input Ensembl IDs
  return(nih_ids[ensembl_ids])
}

# Example usage
ensembl_ids <- c("ENSG00000198821", "ENSG00000244734", "ENSG00000141510")
nih_ids <- ensembl_to_nih(ensembl_ids)
nih_ids
