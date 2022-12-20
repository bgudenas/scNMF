# scNMF
NMF analysis in tumor single-cell data
Data Input: Seurat Object


## Set up data
```{r}
GO = readRDS("/home/bgudenas/Annots/GSEA/MisgDB_c2_c5_c6_c8_h1.rds")
genes = unique(unlist(GO))
seurat_obj = readRDS("./MB_ss2_SO.rds")

## filter genes to only those with known function
genes = genes[genes %in% rownames(seurat_obj)]

cell_samples = seurat_obj$Sample
cell_subgroup = seurat_obj$Subgroup

data_dir = "../test_raw_counts/"
dir.create(data_dir)

```


## Run NMF per group (Tumor subgroup or Batch)
```{r}
for ( i in unique(cell_subgroups)){
  print(i)
  mini = subset(seurat_obj, cells = colnames(seurat_obj)[cell_subgroups == i ])
  
datExpr = nmf_preprocess(seurat_obj = mini,
               genes,
               counts_slot = "counts",
               n_cells = 50)

estim_rank = nmf_rank(datExpr,
         data_dir = data_dir,
         data_prefix = i,
         k_range = 2:8,
         n_run = 20,
         cell_labels = mini$Sample)

res = nmf_run(datExpr,
         data_dir = data_dir,
         data_prefix = i,
         k = 3,
         n_run = 30)

meta_list = nmf_features(nmf_res = res,
            datExpr,
            data_dir = data_dir,
            data_prefix = i,
            feature_method = NULL,
            topn = 30)
}



```