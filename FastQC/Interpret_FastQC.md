## Table of Contents

* [Basic Statistics](#Basic-Statistics)
* [Per Base Sequence Quality](#Per-Base-Sequence-Quality)
* [Per Tile Sequence Quality](#Per-Tile-Sequence-Quality)
* [Per Sequence Quality Scores](#Per-Sequence-Quality-Scores)


## Basic Statistics
- High-quality single-cell data typically exhibit minimal instances of low-quality sequences and frequently have consistent sequence lengths. 
- The GC content should closely align with the overall GC composition of the genome or transcriptome specific to the species being sequenced. 

## Per Base Sequence Quality
- The horizontal axis denotes the read positions, while the vertical axis displays the quality scores. 
- High-quality single-cell data: all the yellow boxes which represent the inter-quantile range of position quality scores, should align with the green region. Likewise, the whiskers which represent the 10th and 90th percentiles of the distribution, should also fall within the green area.
- The quality scores tend to decrease as we move along the body of the reads, causing certain base calls at the last positions to enter the orange region (indicating reasonably good quality calls). This drop in quality is often due to the common reduction in the signal-to-noise ratio in most sequencing-by-synthesis methods. Nevertheless, it's important that the **boxes remain outside the red area** (representing calls of poor quality).
- If we observe the poor quality calls, we may consider performing the quality trimming of reads.

## Per Tile Sequence Quality
The per tile sequence quality plot indicates the variation in read quality across each sequenced tile in the flowcell. Hotter colors represent greater deviations in quality. Ideally, the plot should display uniform blue, indicating consistent quality across all tiles. However, if some tiles show poor quality (indicated by Hotter  colors), it may suggest temporary issues during sequencing, like bubbles or smudges in the flowcell.

## Per Sequence Quality Scores
The per sequence quality score plot displays each read's average quality score distribution in the file. The x-axis represents average quality scores, while the y-axis indicates the occurrence frequency of each score. Ideally, this plot should exhibit a single peak towards the end, signifying good data quality. The emergence of additional peaks may suggest quality issues in a subset of the reads.
