# scNMF
### NMF analysis in tumor single-cell data

#### Workflow
NMF analysis is run for all cells within a given sample to identify intra-tumoral metagenes. Next, the metagenes of a single sample are scored across cells from other samples that belong to the same tumor / subgroup. Unsupervised clustering through dynamic tree cutting is performed to merge sample-specific metagenes into conserved subgroup-level metagenes. These metagenes are then functionally annotated using the enrichR API to define functional enrichment in MsigDB biological processes.

#### Dependecies
 ```R
library(NMF)
library(GSVA)
library(stringr)
library(enrichR)
library(ggplot2)
```
