library(Seurat)

# Load the dataset
sc_data <- Read10X(data.dir = "~/scRNA_analysis/data/filtered_gene_bc_matrices/hg19/")")

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq", min.cells = 3, min.features = 200)

# QC and filter out low-quality cells
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)

# Scale the data
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Perform linear dimensional reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Cluster the cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run non-linear dimensional reduction (t-SNE or UMAP)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Save results
saveRDS(seurat_obj, file = "seurat_obj.rds")
