library(clusterProfiler)
library(data.table)
library(ggplot2)
data(geneList, package = "DOSE")
geneList |> head()
# 4312     8318    10874    55143    55388      991 
# 4.572613 4.514594 4.418218 4.144075 3.876258 3.677857 

##------Extract Entrez gene ID of a subset of gene interest------
gene.id <- names(geneList[abs(geneList) > 2])
head(gene.id)
# [1] "4312"  "8318"  "10874" "55143" "55388" "991" 

##------Construct the GO term------
## groupGO() function is designed for gene classification based on GO distribution at a specific level. 
# One of "MF", "BP", and "CC" subontologies
# MF: Molecular Function
# BP: Biological Process
# CC: Cellular Component
ggo <- groupGO(gene = gene.id,
               OrgDb = "org.Hs.eg.db",
               ont = "MF",
               level = 3,
               readable = T) # T: gene IDs map to gene symbols
head(ggo)
#                    ID                               Description Count GeneRatio                                             geneID
# GO:0000146 GO:0000146              microfilament motor activity     1     1/207                                              MYH11
# GO:0003777 GO:0003777                microtubule motor activity     8     8/207 KIF23/CENPE/KIF18A/KIF11/KIFC1/KIF18B/KIF20A/KIF4A
# GO:0061791 GO:0061791                     GTPase motor activity     0     0/207                                                   
# GO:0140605 GO:0140605 proton motive force-driven motor activity     0     0/207                                                   
# GO:0004133 GO:0004133      glycogen debranching enzyme activity     0     0/207                                                   
# GO:0009975 GO:0009975                          cyclase activity     0     0/207        

##------Translate biological ID------
## Ref: https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html#id-convert
## Example
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
##  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
## [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
## [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
## [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [26] "UNIPROT"

x <- "ENSG00000211445"
ids <- bitr(x, fromType = "ENSEMBL", 
            toType = c("ENTREZID", "SYMBOL"), 
            OrgDb = "org.Hs.eg.db")
# 'select()' returned 1:1 mapping between keys and columns
# ids
# ENSEMBL ENTREZID SYMBOL
# 1 ENSG00000211445     2878   GPX3

##------GO over-representation analysis - enrichGO()------
ego <- enrichGO(gene          = gene.id,
                universe      = names(geneList),
                OrgDb         = "org.Hs.eg.db",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
#                   ID                 Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue
# GO:0008017 GO:0008017         microtubule binding    18/190 200/11755 3.623923e-09 1.674252e-06 1.583082e-06
# GO:0015631 GO:0015631             tubulin binding    18/190 278/11755 5.680781e-07 9.551120e-05 9.031020e-05
# GO:0003777 GO:0003777  microtubule motor activity     8/190  46/11755 6.202026e-07 9.551120e-05 9.031020e-05
# GO:0042379 GO:0042379  chemokine receptor binding     8/190  51/11755 1.412537e-06 1.631480e-04 1.542639e-04
# GO:0008009 GO:0008009          chemokine activity     7/190  38/11755 2.131359e-06 1.969376e-04 1.862135e-04
# GO:0003774 GO:0003774 cytoskeletal motor activity     9/190  83/11755 7.307109e-06 5.626474e-04 5.320088e-04
#                                                                                                                      geneID Count
# GO:0008017 S100A9/KIF23/CENPE/S100A8/DLGAP5/SKA1/NUSAP1/TPX2/KIF18A/BIRC5/KIF11/PRC1/KIFC1/KIF18B/KIF20A/KIF4A/CCDC170/MAPT    18
# GO:0015631 S100A9/KIF23/CENPE/S100A8/DLGAP5/SKA1/NUSAP1/TPX2/KIF18A/BIRC5/KIF11/PRC1/KIFC1/KIF18B/KIF20A/KIF4A/CCDC170/MAPT    18
# GO:0003777                                                               KIF23/CENPE/KIF18A/KIF11/KIFC1/KIF18B/KIF20A/KIF4A     8
# GO:0042379                                                              CXCL10/CXCL13/CXCL11/CXCL9/CCL18/CCL8/CXCL14/CX3CR1     8
# GO:0008009                                                                     CXCL10/CXCL13/CXCL11/CXCL9/CCL18/CCL8/CXCL14     7
# GO:0003774                                                         KIF23/CENPE/KIF18A/KIF11/KIFC1/KIF18B/KIF20A/KIF4A/MYH11     9

ego.dt <- as.data.table(ego)

## Visualize the results
clusterProfiler::dotplot(ego, showCategory = 9) +
  theme(axis.text.y = element_text(size = 7))
## X-axis: gene ratio: no. of DEGs that fall into the GO term set / total no. of DEGs

## Visualize enriched GO terms as a directed acyclic graph
goplot(ego)

## Use other gene IDs
gene.id.ENSEMBL <- bitr(gene.id, fromType = "ENTREZID",
                        toType = "ENSEMBL",
                        OrgDb = "org.Hs.eg.db")

ego2 <- enrichGO(gene          = gene.id.ENSEMBL$ENSEMBL,
                 OrgDb         = "org.Hs.eg.db",
                 keyType       = "ENSEMBL",
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego2)        

#                   ID                                                  Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue
# GO:0008017 GO:0008017                                          microtubule binding    24/222 285/20758 6.145324e-15 2.379662e-12 2.033582e-12
# GO:0032395 GO:0032395                               MHC class II receptor activity    14/222  67/20758 1.006199e-14 2.379662e-12 2.033582e-12
# GO:0015234 GO:0015234                  thiamine transmembrane transporter activity     8/222  11/20758 2.421074e-14 3.817227e-12 3.262079e-12
# GO:0015220 GO:0015220                   choline transmembrane transporter activity     8/222  14/20758 4.286313e-13 5.068565e-11 4.331432e-11
# GO:0015651 GO:0015651 quaternary ammonium group transmembrane transporter activity     8/222  18/20758 6.019994e-12 4.260522e-10 3.640904e-10
# GO:1901474 GO:1901474                     azole transmembrane transporter activity     8/222  18/20758 6.019994e-12 4.260522e-10 3.640904e-10
#                                                                                                                                                                                                                                                                                                                                                                                                     geneID
# GO:0008017 ENSG00000163220/ENSG00000137807/ENSG00000138778/ENSG00000143546/ENSG00000126787/ENSG00000154839/ENSG00000262634/ENSG00000137804/ENSG00000088325/ENSG00000121621/ENSG00000089685/ENSG00000138160/ENSG00000198901/ENSG00000237649/ENSG00000233450/ENSG00000056678/ENSG00000204197/ENSG00000186185/ENSG00000112984/ENSG00000090889/ENSG00000120262/ENSG00000186868/ENSG00000276155/ENSG00000277956
# GO:0032395                                                                                                                                                                 ENSG00000196735/ENSG00000225890/ENSG00000233192/ENSG00000231526/ENSG00000206301/ENSG00000236418/ENSG00000257473/ENSG00000231823/ENSG00000223793/ENSG00000232062/ENSG00000228284/ENSG00000225103/ENSG00000206305/ENSG00000237541
# GO:0015234                                                                                                                                                                                                                                                                 ENSG00000204385/ENSG00000235336/ENSG00000203463/ENSG00000232180/ENSG00000231479/ENSG00000229077/ENSG00000206378/ENSG00000228263
# GO:0015220                                                                                                                                                                                                                                                                 ENSG00000204385/ENSG00000235336/ENSG00000203463/ENSG00000232180/ENSG00000231479/ENSG00000229077/ENSG00000206378/ENSG00000228263
# GO:0015651                                                                                                                                                                                                                                                                 ENSG00000204385/ENSG00000235336/ENSG00000203463/ENSG00000232180/ENSG00000231479/ENSG00000229077/ENSG00000206378/ENSG00000228263
# GO:1901474                                                                                                                                                                                                                                                                 ENSG00000204385/ENSG00000235336/ENSG00000203463/ENSG00000232180/ENSG00000231479/ENSG00000229077/ENSG00000206378/ENSG00000228263
#            Count
# GO:0008017    24
# GO:0032395    14
# GO:0015234     8
# GO:0015220     8
# GO:0015651     8
# GO:1901474     8

##------GO gene set enrichment analysis - gseGO()------
# GSEA analysis requires a ranked gene list, which contains three features:
# numeric vector: fold change or other type of numerical variable
# named vector: every number has a name, the corresponding gene ID
# sorted vector: number should be sorted in decreasing order
ego3 <- gseGO(geneList     = geneList, # Notice here we use the full list
              OrgDb        = "org.Hs.eg.db",
              ont          = "MF",
              minGSSize    = 100, # Minimal size of each geneSet for analyzing
              maxGSSize    = 500, # Maximal size of genes annotated for testing
              pvalueCutoff = 0.05,
              verbose      = FALSE)

##------KEGG enrichment analysis - enrichKEGG()------
# search kegg organism, listed in https://www.genome.jp/kegg/catalog/org_list.html
search_kegg_organism("hsa", by = "kegg_code")
# kegg_code           scientific_name               common_name
# 1          hsa              Homo sapiens                     human
# 7964      hsai      Halorubrum salinarum      Halorubrum salinarum
# 7980      hsal Haloterrigena salifodinae Haloterrigena salifodinae
data(geneList, package = "DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = "hsa",
                 pvalueCutoff = 0.05)
head(kk)

## KEGG pathway gene set enrichment analysis
kk2 <- gseKEGG(geneList     = geneList,
               organism     = "hsa",
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

##------KEGG module over-representation analysis------
# KEGG Module is a collection of manually defined function units. 
# In some situation, KEGG Modules have a more straightforward interpretation.
mkk <- enrichMKEGG(gene = gene,
                   organism = "hsa",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
mkk@result |> View()    

##------KEGG module gene set enrichment analysis------
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = "hsa",
                 pvalueCutoff = 1)
mkk2@result |> View()   

##------Visualize enriched KEGG pathways------
browseKEGG(mkk, "M00938") # pathway ID
# Explore selected KEGG pathway. 
# Differentially expressed genes that are enriched in the selected pathway will be highlighted.


# The following example illustrates how to visualize the “hsa04110” pathway, which was enriched in our previous analysis.
library(pathview)
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene = max(abs(geneList)), cpd = 1))
hsa04110
