##------Attach the libraries and data------
# install.packages('Seurat')
# devtools::install_github('satijalab/seurat-data')
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("biovizBase")
library(SeuratData)
options(timeout = 1000) # positive integer. The timeout for some Internet operations, in seconds. Default 60 (seconds) but can be set from environment variable R_DEFAULT_INTERNET_TIMEOUT.
InstallData("pbmcMultiome.SeuratData")
library(Seurat)
packageVersion("Seurat")
# [1] '5.0.3'
library(Signac) # Analysis of Single-Cell Chromatin Data
library(biovizBase)
library(EnsDb.Hsapiens.v86)
# This package loads an SQL connection to a database containing annotations from Ensembl
# This data package was made from resources at Ensembl on Thu May 18 16:32:27 2017 and based on the 86
library(ggplot2)
library(cowplot)
library(Matrix)

pbmc_rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc_atac <- LoadData("pbmcMultiome", "pbmc.atac")

##------Convert pbmc_rna[["RNA"]] into another class type "Assay5"------
pbmc_rna[["RNA"]] <- as(pbmc_rna[["RNA"]], Class = "Assay5")
class(pbmc_rna[["RNA"]])

##------Repeat QC steps------
pbmc_rna_qc <- subset(pbmc_rna, seurat_annotations != "filtered") # seurat_annotations is a column in the meta.data of pbmc_rna dataset
# pbmc_rna@meta.data[["seurat_annotations"]]
pbmc_atac_qc <- subset(pbmc_atac, seurat_annotations != "filtered")

##------scRNAseq analysis------
pbmc_rna_norm <- NormalizeData(pbmc_rna_qc, normalization.method = "LogNormalize")
# LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p
pbmc_rna_vf <- FindVariableFeatures(selection.method = "vst", pbmc_rna_norm)
# vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
pbmc_rna_scale <- ScaleData(pbmc_rna_vf)
# ScaleData: scales and centers features in the dataset. 
pbmc_rna_pca <- RunPCA(pbmc_rna_scale, npcs = 50)
pbmc_rna_umap <- RunUMAP(pbmc_rna_pca, dims = 1:30)

##------scATACseq analysis------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations@seqnames
# factor-Rle of length 3021151 with 25 runs
# Lengths:  92459  65287 276555 135207 190802 141516 182677 182178 116305 ... 133785  64670 107112   6155  51812 113543  33885     26
# Values :     X      20     1      6      3      7      12     11     4  ...     5      22     10     Y      18     15     21     MT
# Levels(25): X 20 1 6 3 7 12 11 4 17 2 16 8 19 9 13 14 5 22 10 Y 18 15 21 MT

seqlevelsStyle(annotations) <- "UCSC" # Rename the seqlevels of an object according to a given style
annotations@seqnames
# factor-Rle of length 3021151 with 25 runs
# Lengths:  92459  65287 276555 135207 190802 141516 182677 182178 116305 ... 133785  64670 107112   6155  51812 113543  33885     26
# Values :  chrX   chr20  chr1   chr6   chr3   chr7   chr12  chr11  chr4  ...  chr5   chr22  chr10  chrY   chr18  chr15  chr21  chrM 
# Levels(25): chrX chr20 chr1 chr6 chr3 chr7 chr12 chr11 chr4 chr17 chr2 ... chr9 chr13 chr14 chr5 chr22 chr10 chrY chr18 chr15 chr21 chrM

genome(annotations) <- "hg38" # Set the genome as hg38 for a ChromatinAssay
Annotation(pbmc_atac) <- annotations # Set the annotation for a ChromatinAssay

# Exclude the first dimension typically linked to sequencing depth
pbmc_atac <- RunTFIDF(pbmc_atac) # Run term frequency inverse document frequency (TF-IDF) normalization on a matrix
pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = "q0") 
# Find most frequently observed features; 
# Cutoff for feature to be included in the VariableFeatures for the object. This can be a percentile specified as 'q' followed by the minimum percentile, for example 'q5' to set the top 95% most common features as the VariableFeatures for the object. 
# "q0" include all features.
pbmc_atac <- RunSVD(pbmc_atac)
pbmc_atac <- RunUMAP(pbmc_atac, 
                     reduction = "lsi", 
                     dims = 2:30, 
                     reduction.name = "umap.atac", 
                     reduction.key = "atacUMAP_")

# Visualize the results from scRNAseq and scATACseq
# Predict annotations for the scATAC-seq cells
pbmc_rna_viz <- DimPlot(pbmc_rna_umap, 
                        group.by = "seurat_annotations", # seurat_annotations is a column in the meta.data of pbmc_rna dataset
                        label = TRUE) + 
  NoLegend() + 
  ggtitle("RNA")
# Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class

pbmc_atac_viz <- DimPlot(pbmc_atac, 
                         group.by = "orig.ident", 
                         label = FALSE) + 
  NoLegend() + 
  ggtitle("ATAC")

pbmc_rna_viz + pbmc_atac_viz

final_plot <- (pbmc_rna_viz + pbmc_atac_viz) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
ggsave(filename = "EDA_scRNAseq_scATACseq.jpg", 
       height = 8, 
       width = 11, 
       plot = final_plot,
       dpi = 300)
