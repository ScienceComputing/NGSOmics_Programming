# Quantitative Biology (Under Active Construction)

This repository houses various coding practice, learning notes, assignment/competition solutions based on the materials from various computational biology/bioinformatics courses, workshops, technical manuals, academic articles, and others. 

**Analyze single cell RNA-Seq data**

- As our inputs are `fastq` files, we can ...
  - Use kb-python package to perform [psuedo sequence alignment and generate the count matrix](SingleCellRNASeq/kb-python)
  - Use Cell Ranger pipelines to perform sequence alignment and generate the count matrix
- After having the `feature-barcode matrices` at hand, we can ... 
  - Use Seurat workflow to perform [quality assurance, clustering, and marker gene detection](SingleCellRNASeq/SeuratSkinCell.Rmd)
  - Use Bioconductor packages to [orchestrate single cell RNA-Seq data analysis](SingleCellRNASeq/BioconductorSkinCell.Rmd)
  - Use Scanpy workflow to perform [preprocessing, cell clustering, marker gene detection](SingleCellRNASeq/Scanpy/PBMC), and [trajectory inference](SingleCellRNASeq/Scanpy/Bone_Marrow)
  - Use scGen to model the [perturbation responses](SingleCellRNASeq/Perturbation/scGen)  
  - Spatial transcriptomics

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

**[Multimodal data integration](Note_MultimodalDataIntegration.md)**

Multimodal integration of radiology, pathology and genomics for prediction of response to PD-(L)1 blockade in patients with non-small cell lung cancer. Nat Cancer 3, 1151–1164 (2022). https://doi.org/10.1038/s43018-022-00416-8

<br>

**Expression quantitative trait locus analysis using single cell RNA-Seq data**

Bryois, J., Calini, D., Macnair, W. et al. Cell-type-specific cis-eQTLs in eight human brain cell types identify novel risk genes for psychiatric and neurological disorders. Nat Neurosci 25, 1104–1112 (2022). https://doi.org/10.1038/s41593-022-01128-z

