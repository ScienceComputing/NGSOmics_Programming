# Quantitative Biology

This repository houses coding practice, assignment/competition solutions based on the materials from a variety of computational biology/bioinformatics courses, workshops, technical manuals, academic articles, and others. 

## Featured contents

* [Single cell RNA-seq data analysis](#Analyze-single-cell-RNA-seq-data)
* [Bulk RNA-seq data analysis](#Analyze-bulk-RNA-seq-data)
* [COVID-19 RNA-seq data resources](https://github.com/ScienceComputing/COVID-19-RNA-Seq-datasets)
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
    - Use kb-python package to perform [pseudo sequence alignment and generate the count matrix](SingleCellRNASeq/kb-python)
    - Use Cell Ranger pipelines to perform [sequence alignment and generate the count matrix](SingleCellRNASeq/CellRanger/cellranger_count.sh)
- After having the `feature-barcode matrices` at hand, we can ...
  - Use Scanpy workflow to perform [quality assurance, cell clustering, marker gene detection](SingleCellRNASeq/Scanpy/PBMC), and [trajectory inference](SingleCellRNASeq/Scanpy/Bone_Marrow)
  - Use Seurat workflow to perform [quality assurance, cell clustering, and marker gene detection](SingleCellRNASeq/Seurat/SkinCell.Rmd)
  - Use Bioconductor packages to [perform single cell RNA-Seq data analysis](SingleCellRNASeq/Bioconductor/BioconductorSkinCell.Rmd)
  - Generate [pseudobulk](SingleCellRNASeq/Scanpy/Pseudobulk.py), which aggregates the gene expression levels specific to each cell type within an individual
  - Perform pseudobulk-based [differentially gene expression analysis](SingleCellRNASeq/Scanpy/scRNAseq_DE_Part1.ipynb)
  - Use GSEApy to evaluate if a predefined set of genes shows statistically significant and consistent variations between two biological conditions
  - Use scGen to model [perturbation responses](SingleCellRNASeq/Perturbation/scGen)  

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

  - Run [ENCODE ATAC-seq pipeline](https://github.com/ScienceComputing/atac-seq-pipeline/blob/master/README.md) to perform alignment, quality assurance, peaking calling, and signal track generation
  - If we're interested in inspecting every step in each analytical phase, or even leveraging advanced/unique features of other tools that the current pipeline ignores, ...
  - For alignment and post-alignment phases, we can ...
    - Use Rsubread or Rbowtie2 to [align the fastq files relative to hg19/hg38/hs1](ATACSeq/AlignFASTQ.Rmd)
    - Use GenomicAlignments and GenomicRanges to perform post-alignment processing including [reading properly paired reads, estimating MapQ scores/insert sizes, reconstructing the full-length fragment, and others](ATACSeq/PostAlignment.Rmd)
    - Use ATACseqQC to perform [comprehensive ATAC-seq quality assurance](ATACSeq/ATACseqQC.Rmd)
  - For TSS analysis phase, we can ...
    - Use soGGi to [assess the transcriptional start site signal](ATACSeq/EvaluateTSS.Rmd) in the nucleosome-free open region
  - For peaking calling phase, we can ...
    - Use MACS2 and ChIPQC to [call peaks in the nucleosome-free open region, and perform quality assurance](ATACSeq/CallPeak.Rmd)
    - Or use Genrich to call peaks in the nucleosome-free open region
    - Or use MACS3/MACSr (R wrapper of MACS3) to [call peaks in the nucleosome-free open region](ATACSeq/CallPeak.Rmd)
    - Use ChIPseeker to [annotate peak regions with genomic features](ATACSeq/CallPeak.Rmd)
  - For other types of downstream analyses, we can ...
    - Use rGREAT to [functionally interpret the peak regions based on the GO database](ATACSeq/FunctionalAnalysis.Rmd) 
    - Use GenomicRanges and GenomicAlignments to [select and count non-redundant peaks](ATACSeq/DifferentialAnalysis.Rmd)
    - Use DESeq2/DESeq2-based DiffBind and ChIPseeker to [perform differential analysis with gene annotations](ATACSeq/DifferentialAnalysis.Rmd)
    - Use clusterProfiler to [perform enrichment analysis of differential peak regions](ATACSeq/DifferentialAnalysis.Rmd)
  - [Search and visualize motifs](ATACSeq/Search_Visualize_Motif.Rmd)
  - [Map peaks to motifs](ATACSeq/IdentifyMotif.Rmd)
  - [Analyze differences in motifs across conditions](ATACSeq/Detect_Difference_Motif.Rmd)

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

