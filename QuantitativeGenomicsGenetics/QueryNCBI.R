# install.packages("rentrez")

library(rentrez)

# Find the available databases, and return a summary of each of those databases
e.dbs <- entrez_dbs() 
for (ele in e.dbs){
  print(entrez_db_summary(ele))
  cat("\n") 
}

# DbName: pubmed
# MenuName: PubMed
# Description: PubMed bibliographic record
# DbBuild: Build-2023.03.09.09.47
# Count: 35398194
# LastUpdate: 2023/03/09 09:47 
# 
# DbName: protein
# MenuName: Protein
# Description: Protein sequence record
# DbBuild: Build230304-1526m.1
# Count: 1131717045
# LastUpdate: 2023/03/06 17:14 
# 
# DbName: nuccore
# MenuName: Nucleotide
# Description: Core Nucleotide db
# DbBuild: Build230305-1015m.1
# Count: 587074925
# LastUpdate: 2023/03/07 20:21 
# 
# DbName: ipg
# MenuName: Identical Protein Groups
# Description: Identical Protein Groups DB
# DbBuild: Build230228-0203.1
# Count: 555827569
# LastUpdate: 2023/03/04 22:24 
# 
# DbName: nuccore
# MenuName: Nucleotide
# Description: Core Nucleotide db
# DbBuild: Build230305-1015m.1
# Count: 587074925
# LastUpdate: 2023/03/07 20:21 
# 
# DbName: structure
# MenuName: Structure
# Description: Three-dimensional molecular model
# DbBuild: Build230309-0100.1
# Count: 202329
# LastUpdate: 2023/03/09 02:30 
# 
# DbName: genome
# MenuName: Genome
# Description: Genomic sequences, contigs, and maps
# DbBuild: Build230307-1905.1
# Count: 99186
# LastUpdate: 2023/03/07 20:03 
# 
# DbName: annotinfo
# MenuName: AnnotInfo
# Description: Annotinfo Database
# DbBuild: Build230309-0015.1
# Count: 1715
# LastUpdate: 2023/03/09 00:44 
# 
# DbName: assembly
# MenuName: Assembly
# Description: Genome Assembly Database
# DbBuild: Build230309-0840.1
# Count: 1660548
# LastUpdate: 2023/03/09 11:35 
# 
# DbName: bioproject
# MenuName: BioProject
# Description: BioProject Database
# DbBuild: Build230309-0620.1
# Count: 668858
# LastUpdate: 2023/03/09 07:28 
# 
# DbName: biosample
# MenuName: BioSample
# Description: BioSample Database
# DbBuild: Build230309-0301m.1
# Count: 31289726
# LastUpdate: 2023/03/09 07:35 
# 
# DbName: blastdbinfo
# MenuName: BlastdbInfo
# Description: BlastdbInfo Database
# DbBuild: Build230308-1422.1
# Count: 22793812
# LastUpdate: 2023/03/08 15:59 
# 
# DbName: books
# MenuName: Books
# Description: Books Database
# DbBuild: Build230309-0325.1
# Count: 1115591
# LastUpdate: 2023/03/09 04:30 
# 
# DbName: cdd
# MenuName: Conserved Domains
# Description: Conserved Domain Database
# DbBuild: Build220919-0918.1
# Count: 64234
# LastUpdate: 2022/09/20 11:22 
# 
# DbName: clinvar
# MenuName: ClinVar
# Description: ClinVar Database
# DbBuild: Build230306-0240.1
# Count: 2208280
# LastUpdate: 2023/03/06 09:47 
# 
# DbName: gap
# MenuName: dbGaP
# Description: dbGaP Data
# DbBuild: Build230222-0855m.1
# Count: 363716
# LastUpdate: 2023/02/22 09:24 
# 
# DbName: gapplus
# MenuName: GaPPlus
# Description: Internal Genotypes and Phenotypes database
# DbBuild: Build170929-0435.1
# Count: 136796
# LastUpdate: 2017/09/29 04:56 
# 
# DbName: grasp
# MenuName: grasp
# Description: grasp Data
# DbBuild: Build150126-1400.1
# Count: 7862970
# LastUpdate: 2015/01/26 16:10 
# 
# DbName: dbvar
# MenuName: dbVar
# Description: dbVar records
# DbBuild: Build230308-2325.1
# Count: 7744904
# LastUpdate: 2023/03/09 05:03 
# 
# DbName: gene
# MenuName: Gene
# Description: Gene database
# DbBuild: Build230307-2355m.1
# Count: 65252800
# LastUpdate: 2023/03/08 16:04 
# 
# DbName: gds
# MenuName: GEO DataSets
# Description: GEO DataSets
# DbBuild: Build230308-1822.1
# Count: 5794797
# LastUpdate: 2023/03/08 21:37 
# 
# DbName: geoprofiles
# MenuName: GEO Profiles
# Description: Genes Expression Omnibus
# DbBuild: Build160819-1300.288
# Count: 128414055
# LastUpdate: 2023/03/07 11:48 
# 
# DbName: homologene
# MenuName: HomoloGene
# Description: HomoloGene Database
# DbBuild: Build140512-1105.4
# Count: 141268 
# 
# DbName: medgen
# MenuName: MedGen
# Description: Medgen Database
# DbBuild: Build230309-0747.1
# Count: 213613
# LastUpdate: 2023/03/09 08:20 
# 
# DbName: mesh
# MenuName: MeSH
# Description: MeSH Database
# DbBuild: Build230309-0315.1
# Count: 353125
# LastUpdate: 2023/03/09 03:42 
# 
# DbName: nlmcatalog
# MenuName: NLM Catalog
# Description: NLM Catalog Database
# DbBuild: Build230309-0615.1
# Count: 1630028
# LastUpdate: 2023/03/09 06:56 
# 
# DbName: omim
# MenuName: OMIM
# Description: OMIM records
# DbBuild: Build230309-0305.1
# Count: 28178
# LastUpdate: 2023/03/09 03:51 
# 
# DbName: orgtrack
# MenuName: Orgtrack
# Description: Orgtrack Database
# DbBuild: Build230309-0800.1
# Count: 7945
# LastUpdate: 2023/03/09 08:44 
# 
# DbName: pmc
# MenuName: PMC
# Description: PubMed Central
# DbBuild: Build230308-1045m.1
# Count: 8900398
# LastUpdate: 2023/03/08 19:13 
# 
# DbName: popset
# MenuName: PopSet
# Description: PopSet sequence record
# DbBuild: Build230309-0427m.1
# Count: 395368
# LastUpdate: 2023/03/09 06:25 
# 
# DbName: proteinclusters
# MenuName: Protein Clusters
# Description: Protein Cluster record
# DbBuild: Build171204-1005.1
# Count: 1137329
# LastUpdate: 2017/12/04 13:20 
# 
# DbName: pcassay
# MenuName: PubChem BioAssay
# Description: PubChem BioAssay Database
# DbBuild: Build230303-0546.1
# Count: 1506716
# LastUpdate: 2023/03/03 08:03 
# 
# DbName: protfam
# MenuName: Protein Family Models
# Description: protfam DB
# DbBuild: Build230306-2330.1
# Count: 163907
# LastUpdate: 2023/03/07 11:43 
# 
# DbName: pccompound
# MenuName: PubChem Compound
# Description: PubChem Compound Database
# DbBuild: Build230308-0245m.1
# Count: 113992523
# LastUpdate: 2023/03/08 21:50 
# 
# DbName: pcsubstance
# MenuName: PubChem Substance
# Description: PubChem Substance Database
# DbBuild: Build230308-0245m.1
# Count: 301525869
# LastUpdate: 2023/03/08 14:46 
# 
# DbName: seqannot
# MenuName: SeqAnnot
# Description: SeqAnnot Database
# DbBuild: Build230309-0804.1
# Count: 387014
# LastUpdate: 2023/03/09 09:35 
# 
# DbName: snp
# MenuName: SNP
# Description: Single Nucleotide Polymorphisms
# DbBuild: Build221118-1625.1
# Count: 1121739543
# LastUpdate: 2022/11/22 11:07 
# 
# DbName: sra
# MenuName: SRA
# Description: SRA Database
# DbBuild: Build230309-0703m.1
# Count: 26884609
# LastUpdate: 2023/03/09 11:13 
# 
# DbName: taxonomy
# MenuName: Taxonomy
# Description: Taxonomy db
# DbBuild: Build230309-0810.1
# Count: 2626678
# LastUpdate: 2023/03/09 09:19 
# 
# DbName: biocollections
# MenuName: Biocollections
# Description: Biocollections db
# DbBuild: Build230309-0315.1
# Count: 8497
# LastUpdate: 2023/03/09 03:52 
# 
# DbName: gtr
# MenuName: GTR
# Description: GTR Database
# DbBuild: Build230309-0800.1
# Count: 76600
# LastUpdate: 2023/03/09 08:46 

# Find which fields we can search by in a database
f.genome <- entrez_db_searchable("genome") 
f.gene <- entrez_db_searchable("gene")
f.genome
f.gene
gene.search <- entrez_search(db = "gene", 
                             term = "(OPN1LW[GENE]) AND (Homo sapiens[ORGN])")

library(XML)
xmlrec <- entrez_fetch(db = "gene", 
                       id = "OPN1LW", 
                       rettype = "xml", 
                       parsed = T) 
xmldf <- xmlToDataFrame(xmlrec)
names(xmldf)
xmldf


