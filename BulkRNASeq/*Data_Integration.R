library(dplyr)
tables <- c(sprintf("../data/Table%d.csv", 1:3), sprintf("../data/Table%d.rds", 4:5))
load_table <- function(table) {
  if (grepl('\\.csv$', table)) { read.csv(table) } else if (grepl('\\.rds$', table)) { readRDS(table) }
}
list2env(x = setNames(object = lapply(tables, load_table), nm = paste0('table_', 1:5)), envir = .GlobalEnv)
names(table_3)
any(table_3$PF_READS > table_3$TOTAL_READS)
hist(table_3$PF_ALIGNED_BASES/table_3$PF_READS_ALIGNED) 
any(table_3$PF_READS_ALIGNED != table_3$READS_ALIGNED_IN_PAIRS)
class(table_3$STRAND_BALANCE); summary(table_3$STRAND_BALANCE)

# View the dimension of each dataframe
sapply(X = paste0('table_', 1:5), FUN = function(table) dim(get(table)), simplify = F) 

# Explore the relationships between sorts of ids from each dataframe
full_join(table_1, table_2, by = c('specimenID' = 'Sampleid')) |> View()
setdiff(table_2$projid, table_4$id)
setdiff(table_2$projid, table_5$projid)
table_3$sample_id |> unique() |> length()
table_4$id |> unique() |> length()
setdiff(table_4$id |> unique(), table_5$projid)

# Create a datatable including specified columns
library(data.table)
data_table <- table_1 |> dplyr::filter(name != '280_129.bam') |> 
  left_join(y = table_2, by = c('specimenID' = 'Sampleid')) |> 
  left_join(y = table_3 |> dplyr::filter(CATEGORY == 'FIRST_OF_PAIR'), 
            by = c('specimenID' = 'sample_id')) |> 
  left_join(y = table_4, by = c('projid' = 'id')) |> 
  select(c(specimenID, individualID, projid, 
           RINcontinuous, PCT_PF_READS, PCT_PF_READS_ALIGNED, PCT_READS_ALIGNED_IN_PAIRS, PCT_READS_ALIGNED_IN_PAIRS, PF_ALIGNED_BASES, PF_READS, PF_READS_ALIGNED, TOTAL_READS,
           visit, age_at_visit, cogdx, sex, educ)) |>
  rename(RIN = RINcontinuous, `age at visit` = age_at_visit, diagnosis = cogdx, education = educ) 

data_table <- as.data.table(data_table) 
class(data_table) 
head(data_table)

# Recode the diagnosis variable as a new column
data_table[,diagnosis_2:=gsub('^T1.*', 1, 
                              gsub('^T2.*', 2,
                                   gsub('^T3.*', 3,
                                        gsub('^Other.*', 4, diagnosis)))) |> as.numeric()]
table(data_table$diagnosis_2, useNA = 'ifany')

# Create a table only including participants’ initial baseline visit information
data_table_base <- data_table[visit==0]
head(data_table_base)
dim(data_table) # longitudinal table
dim(data_table_base) # cross-sectional table
round(table(data_table_base$sex, useNA = 'ifany')/nrow(data_table_base) * 100, 0)
mean(data_table_base$education) |> round(2)
sd(data_table_base$education) |> round(2)
mean(data_table_base$`age at visit`) |> round(2)
data_table_last <- data_table |> arrange(visit) |> group_by(projid) |> slice_tail()
mean(data_table_last$`age at visit`) |> round(2)

# View the distributions of first and last 6 gene expression values
library(gridExtra)
library(ggplot2)
library(tidyr)
library(factoextra)
make_hist <- function(gene_data, gene_id, fill_palette) {
  ggplot(data = get(gene_data), 
         mapping = aes(x = get(gene_id))) + 
    geom_histogram(bins = 30, fill = fill_palette) + 
    theme_bw() +
    labs(y = 'Frequency', x = 'Gene expression', title = paste0('Distribution of ', gene_id)) + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
}
plot_array_1 <- mapply(FUN = make_hist, gene_data = 'table_5', gene_id = names(table_5)[2:7], 
                       fill_palette = c('#E18727FF', '#EFC000FF', '#868686FF', '#CD534CFF', '#0072B5FF', '#20854EFF'), 
                       SIMPLIFY = F)
grid.arrange(do.call('arrangeGrob', c(plot_array_1, ncol = 2))) 
plot_array_2 <- mapply(FUN = make_hist, gene_data = 'table_5', gene_id = names(table_5)[96:101], 
                       fill_palette = c('#E18727FF', '#EFC000FF', '#868686FF', '#CD534CFF', '#0072B5FF', '#20854EFF'), 
                       SIMPLIFY = F)
grid.arrange(do.call('arrangeGrob', c(plot_array_2, ncol = 2)))
mean_sd <- table_5[,-1] |>
  summarize(across(where(is.numeric), list(mean = ~ mean(., na.rm = T), sd = ~ sd(., na.rm = T)), .names = "{col}_{fn}")) |>
  pivot_longer(cols = everything(), names_to = c("Gene id", ".value"), names_sep = "_")
hist(mean_sd$mean) 
hist(mean_sd$sd) 

# Check whether sampleids are within their own column
any(rownames(table_5) != table_5$projid) 
x <- colSums(is.na(table_5)); x[x!=0] 
table_5_imputed <- table_5 |> 
  select(!projid) |> 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = T), .))) # Alter: use KNN to impute
# https://www.rdocumentation.org/packages/impute/versions/1.46.0/topics/impute.knn
# impute.knn uses $k$-nearest neighbors in the space of genes to impute missing expression values

# Build PCA and draw a PCA plot
pca_data <- prcomp(x = table_5_imputed, center = T, scale. = T) 
pca_data$x |> View() 
(pca_plot <- fviz_pca_ind(pca_data, geom.ind = c('point', 'text'), repel = T, title = 'PCA on Table_5') + 
    labs(caption = 'Each point represents each sample with the projid label.'))
ggsave(filename = 'PCA_Plot.jpg', height = 9, width = 6,  plot = pca_plot, dpi = 300)

# Translate ENSG ids to gene symbols
library(clusterProfiler)
library(org.Hs.eg.db)
(ensg_to_symbol <- bitr(colnames(table_5)[-1], fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db'))
