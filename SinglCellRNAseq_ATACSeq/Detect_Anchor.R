##------Identify anchors between scRNA-seq and scATAC-seq data------
# Quantify gene activity
# Estimate the transcriptional activity of each gene by quantifying ATAC-seq counts in the 2 kb-upstream region and gene body
gene_activities <- GeneActivity(pbmc_atac, features = VariableFeatures(pbmc_rna))

# Allocate gene activities to a new assay
pbmc_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)

# Normalize gene activities
DefaultAssay(pbmc_atac) <- "ACTIVITY"
pbmc_atac <- NormalizeData(pbmc_atac)
pbmc_atac <- ScaleData(pbmc_atac, features = rownames(pbmc_atac))

# Detect anchors
transfer_anchors <- FindTransferAnchors(reference = pbmc_rna, 
                                        query = pbmc_atac, 
                                        features = VariableFeatures(object = pbmc_rna),
                                        reference.assay = "RNA", 
                                        query.assay = "ACTIVITY", 
                                        reduction = "cca")
