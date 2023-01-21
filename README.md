# scNMF
NMF analysis in tumor single-cell data (R scripts)
Data Input: Seurat Object



```
all_metas = list()

for ( i in unique(cell_subgroups)){
  mini = subset(ss2, cells = colnames(ss2)[cell_subgroups == i ])
  min_cells = round(ncol(mini)*0.02) ## 2% of cells
  k_val = as.numeric(comps[names(comps) == i ])
  i = stringr::str_replace(i, "/", "_") ## avoid / bc it messes filepaths
  print(i)
  
datExpr = nmf_preprocess(seurat_obj = mini,
               genes,
               counts_slot = "scale.data",
               n_cells = min_cells)

#estim_rank = nmf_rank(datExpr,
#         data_dir = data_dir,
#         data_prefix = i,
#         k_range = 2:12,
#         n_run = 15,
#         cell_labels = mini$Sample)

print(paste0("K-value = ", k_val))
res = nmf_run(datExpr,
          data_dir = data_dir,
          data_prefix = i,
          k = k_val,
          n_run = 150)
 
meta_list = nmf_features(nmf_res = res,
            datExpr = datExpr,
            data_dir = data_dir,
            data_prefix = i,
            feature_method = NULL,
            MT_filt = TRUE,
            topn = 50)

metagenes = meta_list[[1]][[1]]
names(metagenes) = paste0(i,"-", toupper(letters[1:length(metagenes)]))

## collect all metagenes
all_metas = c(all_metas, metagenes)

dir.create( paste0(data_dir, "GS/"), showWarnings = FALSE)
meta_functions = annotate_metagenes(metagenes, bp_db="GO_Biological_Process_2021")
if (length(meta_functions$Subgroup) > 0 ){
  meta_functions$Subgroup = i
  plot_metagenes(meta_functions, plot_dir = paste0(data_dir, "GS/"), suffix = "_GO")
}

meta_functions = annotate_metagenes(metagenes, bp_db="PanglaoDB_Augmented_2021")
if (length(meta_functions$Subgroup) > 0 ){
meta_functions$Subgroup = i
plot_metagenes(meta_functions, plot_dir = paste0(data_dir, "GS/"), suffix = "_PanglaoDB")
  }
}

saveRDS(all_metas, paste0(data_dir, "All_metagenes.rds"))


```
