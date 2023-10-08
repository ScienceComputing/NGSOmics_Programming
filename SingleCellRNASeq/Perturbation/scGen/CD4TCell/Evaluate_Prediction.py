# Put predicted CD4+T cells among control populations, control CD4+T cells, and true CD4+T cells with perturbations into one object
eval_adata = ctrl_adata.concatenate(cd4_train_stim, pred)
eval_adata.shape # (16798, 15706); 16789 cells; 15706 genes
eval_adata.obs.condition.value_counts() # The number of predicted CD4+T cells should equal the number of control CD4+T cells

# To evaluate the predictions qualitatively, we embed the latent space of gene expression in all 3 types of CD4+T cells via PCA
sc.tl.pca(eval_adata)
sc.pl.pca(eval_adata, color="condition", frameon=False)

# To evaluate the predictions quantitatively,
# (1) we measure the correlation between mean gene expression across all genes of the predicted CD4+T control cells vs real CD4+T cells, both types of cells under perturbations
# (2) we measure the correlation between mean gene expression across top 100 DEGs of predicted CD4+T control cells and real CD4+T cells, both types of cells under perturbations

# Retrieve DEGs between CD4+T control cells and real CD4+T perturbed cells
cd4t_adata = adata[adata.obs["cell_type"] == "CD4 T cells"]
sc.tl.rank_genes_groups(cd4t_adata, groupby="condition", method="wilcoxon")
de_gene_names = cd4t_adata.uns["rank_genes_groups"]["names"]["stimulated"]
# https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
de_gene_names
# array(['ISG15', 'IFI6', 'ISG20', ..., 'EEF1A1', 'FTH1', 'RGCC'],
#       dtype=object)

r2_value = model.reg_mean_plot(
    eval_adata,
    axis_keys={"x": "predicted stimulated", "y": "stimulated"},
    gene_list=de_gene_names[:10],
    top_100_genes=de_gene_names,
    labels={"x": "predicted", "y": "ground truth"},
    show=True,
    legend=False,
)

# View the distribution of the top 3 DEGs (ISG15, IFI6, ISG20) between CD4+T control cells and real CD4+T perturbed cells  
sc.pl.violin(eval_adata, keys=["ISG15", "IFI6", "ISG20"], groupby="condition", multi_panel=True)
