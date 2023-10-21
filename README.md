# Quantitative Biology

This repository houses coding practice, assignment/competition solutions based on the materials from a variety of computational biology/bioinformatics courses, workshops, technical manuals, academic articles, and others. 

## Technical aspect
**Analyze single cell RNA-Seq data**
- If we are given raw `bcl` files, we need to [convert them to `fastq` files](FastQC/bcl_to_fastq.sh).
- As our inputs are `fastq` files, we can ...
  - Run FastQC to [evaluate sequence quality and content](FastQC)
  - Use kb-python package to perform [psuedo sequence alignment and generate the count matrix](SingleCellRNASeq/kb-python)
  - Use Cell Ranger pipelines to perform [sequence alignment and generate the count matrix](SingleCellRNASeq/CellRanger/cellrangercount.sh)
- After having the `feature-barcode matrices` at hand, we can ...
  - Use Scanpy workflow to perform [preprocessing, cell clustering, marker gene detection](SingleCellRNASeq/Scanpy/PBMC), and [trajectory inference](SingleCellRNASeq/Scanpy/Bone_Marrow)
  - Use Seurat workflow to perform [quality assurance, clustering, and marker gene detection](SingleCellRNASeq/SeuratSkinCell.Rmd)
  - Use Bioconductor packages to [orchestrate single cell RNA-Seq data analysis](SingleCellRNASeq/BioconductorSkinCell.Rmd)
  - Use scGen to model the [perturbation responses](SingleCellRNASeq/Perturbation/scGen)  

<hr>

**Analyze bulk RNA-Seq data**

  - [Align and count the reads](BulkRNASeq/AlignmentCountingTCell.Rmd)
  - [Perform differential gene expression analysis](BulkRNASeq/DEAnalysisTCell.Rmd)
  - [Perform principal component analysis, heatmap, and clustering](BulkRNASeq/PCAHeatmapClusteringTissue.Rmd)
  - [Perform gene set enrichment analysis](BulkRNASeq/GeneSetTCell.Rmd)

<hr>


**Analyze ATAC-Seq data**
  
  - [Align the FASTQ files](ATACSeq/AlignFASTQ.Rmd)
  - [Perform post-alignment processing](ATACSeq/PostAlignment.Rmd)
  - [Perform ATAC-Seq quality assurance using ATACseqQC](ATACSeq/ATACseqQC.Rmd)
  - [Evaluate the transcriptional start site signal](ATACSeq/EvaluateTSS.Rmd)
  - [Call peaks with quality control](ATACSeq/CallPeak.Rmd)

<hr>

**Analyze proteomics data**

*Under Active Construction*
- A quick start from [loading an online spectrum, performing peak quality control, annotating peaks, to visulizing the annotated peaks](Proteomics/spectrum_utils/0_Quick_Start.py)

<hr>

## Conceptual aspect

[**High-level scientific idea**](HighLevelIdea_MultiOmics.md)



