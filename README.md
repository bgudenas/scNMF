# scNMF
NMF analysis in tumor single-cell data (R scripts)
Data Input: Seurat Object


### Load data
```
library(Seurat)
source("~/src/scNMF/singlecell_nmf.R")
source("~/src/scNMF/annotate_nmf_metagenes.R")

ss2 = readRDS("../Data/snPB_ss2_TPM_so.rds") ## seurat object

map = readRDS("~/Annots/Annotables/hg38.rds")
map = map[map$gene_biotype == "protein_coding", ]
#genes = genes[genes %in% map$external_gene_name]
genes = map$external_gene_name
genes = genes[genes %in% rownames(ss2)]

cell_samples = ss2$Sample
cell_subgroups =  ss2$Subgroup ## defines subset of cells NMF is run on

data_dir = "../NMF/Raw_Scaled_TPM_r6/"
dir.create(data_dir)
```

## intialize rank (k) (can finetune after rank estimation)
```
comps = c("PB-MYC/FOXR2" = 6,
          "PPTID" = 6,
          "PB-miRNA1" = 6,
          "PB-RB" = 6,
          "PC" = 6,
          "PB-miRNA2" = 6,
          "PAT" = 6,
          "PTPR" = 6)
```

## Loop through cell_subgroups and run NMF on each subset
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

## creates scatterplot for rank estimation
estim_rank = nmf_rank(datExpr,
         data_dir = data_dir,
         data_prefix = i,
         k_range = 2:12,
         n_run = 15,
         cell_labels = mini$Sample)

print(paste0("K-value = ", k_val))
res = nmf_run(datExpr,
          data_dir = data_dir,
          data_prefix = i,
          k = k_val,
          n_run = 150)

## creates NMF metagenes
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


### functionally annotate metagenes

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
