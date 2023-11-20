# Quantitative Biology

This repository houses coding practice, assignment/competition solutions based on the materials from a variety of computational biology/bioinformatics courses, workshops, technical manuals, academic articles, and others. 

## Featured contents

* [Single cell RNA-seq data analysis](#Analyze-single-cell-RNA-seq-data)
* [Bulk RNA-seq data analysis](#Analyze-bulk-RNA-seq-data)
* [Multi-omics idea](HighLevelIdea_MultiOmics.md)


## Technical aspect
### Analyze single cell RNA-seq data
- A mini scRNA-seq [pipeline](SingleCellRNASeq/Scanpy/pipeline.py)
- If given raw `bcl` files, we [convert them to fastq files](FastQC/bcl_to_fastq.sh)
- As inputs are `fastq` files, we can ...
  - Run FastQC to [evaluate sequence quality and content](FastQC/Run_FastQC.sh)
  - Use Trim Galore to [trim reads](FastQC/Trim_Read.sh) if we spot low-quality base calls/adaptor contamination
  - Re-run FastQC to [re-evaluate sequence quality and content](FastQC/Run_FastQC.sh)
  - If single-cell RNA-seq data is generated from the plate-based protocol, we can ...
    - Use STAR to perform alignment and FeatureCounts to generate the count matrix
  - Else if single-cell RNA-seq data is generated from the droplet-based protocol, we can ...
    - Use kb-python package to perform [psuedo sequence alignment and generate the count matrix](SingleCellRNASeq/kb-python)
    - Use Cell Ranger pipelines to perform [sequence alignment and generate the count matrix](SingleCellRNASeq/CellRanger/cellranger_count.sh)
- After having the `feature-barcode matrices` at hand, we can ...
  - Use Scanpy workflow to perform [quality assurance, cell clustering, marker gene detection](SingleCellRNASeq/Scanpy/PBMC), and [trajectory inference](SingleCellRNASeq/Scanpy/Bone_Marrow)
  - Use Seurat workflow to perform [quality assurance, cell clustering, and marker gene detection](SingleCellRNASeq/Seurat/SkinCell.Rmd)
  - Use Bioconductor packages to [perform single cell RNA-Seq data analysis](SingleCellRNASeq/Bioconductor/BioconductorSkinCell.Rmd)
  - Use scGen to model the [perturbation responses](SingleCellRNASeq/Perturbation/scGen)  

<hr>

### Analyze bulk RNA-seq data

  - [Recommend] Use STAR to [align the reads](BulkRNASeq/STAR_Align.sh)
  - Use Rsubread to [align the reads](BulkRNASeq/AlignmentCountingTCell.Rmd)
    - **Why align?** To pinpoint the specific location on the human genome from which our reads originated
  - Use Qualimap to perform [quality assurance](BulkRNASeq/Qualimap_QC.sh) on the aligned reads
  - Use GenomicAlignments for aligned reads to [obtain the gene-level or exon-level quantification](BulkRNASeq/AlignmentCountingTCell.Rmd)
  - Use featureCounts for aligned reads to [count the fragments](BulkRNASeq/featureCounts.sh)
  - [Recommend] Use Salmon for unaligned reads to [obtain the transcript-level quantification](BulkRNASeq/Salmon_quant.sh)
    - **Why unalign?** To speed up the counting process of reads
    - Next step: Use tximport to aggregate transcript-level quantification to the gene level
  - [Perform differential gene expression analysis](BulkRNASeq/DEAnalysisTCell.Rmd)
  - [Perform principal component analysis, heatmap, and clustering](BulkRNASeq/PCAHeatmapClusteringTissue.Rmd)
  - [Perform gene set enrichment analysis](BulkRNASeq/GeneSetTCell.Rmd)

<hr>


### Analyze ATAC-seq data
  
  - [Align the fastq files](ATACSeq/AlignFASTQ.Rmd)
  - [Perform post-alignment processing](ATACSeq/PostAlignment.Rmd)
  - [Perform ATAC-seq quality assurance using ATACseqQC](ATACSeq/ATACseqQC.Rmd)
  - [Evaluate the transcriptional start site signal](ATACSeq/EvaluateTSS.Rmd)
  - [Call peaks with quality control](ATACSeq/CallPeak.Rmd)
  - [Perform differential and enrichment analysis with peaks](ATACSeq/DifferentialAnalysis.Rmd)

<hr>

### Analyze integrated single cell RNA-seq and ATAC-seq data

*Under Active Construction*

<hr>

### Analyze proteomics data

- A quick start from [loading an online spectrum, performing peak quality control, annotating peaks, to visualizing the annotated peaks](Proteomics/spectrum_utils/0_Quick_Start.py)

<hr>

## Conceptual aspect

- [High level multi-omics idea](HighLevelIdea_MultiOmics.md)
- [Case-control design](CaseControl_Design.md)
- [Graph neural networks in biomedical studies](Reference/Graph_Neural_Network.md)

