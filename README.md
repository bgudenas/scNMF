# scNMF
### NMF analysis in tumor single-cell data

#### Workflow
NMF analysis is run for all cells within a given sample to identify intra-tumoral metagenes. Next, the metagenes of a single sample are scored across cells from other samples that belong to the same tumor / subgroup. Unsupervised clustering through dynamic tree cutting is performed to merge sample-specific metagenes into conserved subgroup-level metagenes. These metagenes are then functionally annotated using the enrichR API to define functional enrichment in MsigDB biological processes.

#### Dependecies
 ```R
 #R/4.0.2
library(NMF)
library(GSVA)
library(stringr)
library(enrichR)
library(ggplot2)
```

#### parameter defaults
```R
genes = rownames(seurat_obj) ## only analyze the provided gene list (could be genes in MsigDB or passing an expression filter)
nrun = 2 # how many iterations to run NMF per sample (set to 2 for fast testing but ~200 for accurate results)
k = 3 # how many metagenes to identify -- through metagene merging can be more or less than this at subgroup level
topn = 30 # how many genes should be selected to characterize the metagenes, i.e. 30 or 50...

## Currently uses normalized counts as input -- seurat_obj@assays$RNA@data
```

#### Running scNMF
```R
source("singlecell_nmf.R")
source("annotate_nmf_metagenes.R")

meta_list = NMF_samples(seurat_obj = seurat_obj, 
                        genes = genes, 
                        cell_samples = cell_samples, 
                        nrun = 200, k = 4)
                        
metagenes = NMF_consensus(seurat_obj = seurat_obj,
                          meta_list = meta_list,
                          cell_samples = cell_samples,
                          cell_subgroups = cell_subgroups,
			                       topn=30)

## Functional annotation by EnrichR
meta_functions = annotate_metagenes(metagenes)
## Plotting results in ./Figures
plot_metagenes(meta_functions)


```
