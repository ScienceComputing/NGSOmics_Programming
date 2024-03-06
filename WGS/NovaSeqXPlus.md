- XLEAP-SBS Chemistry with novel polymerase
  - faster incorporation
  - higher accuracy: large majority of bases have quality scores >= 40
  - more stable to thermal stress, stable in liquid
- Dual flow cell processing up to 16T
- 25B flow cell
  - handle 26 billion reads
  - run time: 48 hours
  - 64 genomes per flow cell (8 per lane; 8 lanes per flow cell)
  - output per FC: 3T
- 10B flow cell (for single cell projects)
  - handle 10 billion reads
  - cycles: 100
  - run time: 18 hours
  - cells per flow cell: 400000 cells (25K reads per cell)
- Accurate (*F1 score*) secondary analysis with [DRAGEN graph](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/In/Informatics/DRAGEN/GraphMapper_fDG.htm) and machine learning 

```
total_reads_per_flow_cell = 26000 # in million
expected_reads_per_sample = 40 # in million
number_of_samples_per_flow_cell = total_reads_per_flow_cell / expected_reads_per_sample = 26000 / 40 = 650
```
