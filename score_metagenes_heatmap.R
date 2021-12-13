
seurat_obj = readRDS("../../snPB_integrated4_renorm.rds")
cell_subgroups = seurat_obj$Subgroup
cell_samples = seurat_obj$Sample

score_genesets = function(seurat_obj, 
                          metagenes, 
                          cell_subgroups, 
                          cell_samples, 
                          plotdir="Figures/Heatmaps", 
                          method="plage"){
  library(pheatmap)
  for ( i in unique(cell_subgroups)){
  print(i)
  cell_filt = cell_subgroups == i
  df = as.matrix(seurat_obj@assays$RNA@data[ ,cell_filt] )
  gene_filt = rowSums(df >= 1) >= 100
  df = df[gene_filt, ]
  subgroup_metagenes = metagenes[grepl(i, names(metagenes))]
  
  ssGSEA = gsva(df, gset.idx.list = subgroup_metagenes, method=method, min.sz=8, max.sz=300, ssgsea.norm = FALSE)
  ssGSEA = t(scale(t(ssGSEA)))
  ssGSEA[ssGSEA > 3 ] = 3
  ssGSEA[ssGSEA < -3 ] = -3
  cols = colorRampPalette(colors = c("blue4","blue","white","red","red4") )(100)
  if (nrow(ssGSEA) >= 2){
  p1 = pheatmap(ssGSEA,
                show_rownames = TRUE,
                show_colnames = FALSE,
                color = cols,
                annotation_col = data.frame(row.names = colnames(df), "Sample" = cell_samples[cell_filt] ),
                main = i )
  } else {
    store_colnames = colnames(ssGSEA)
    ssGSEA = t(matrix(ssGSEA[ ,order(ssGSEA[1, ], decreasing = TRUE)]))
    colnames(ssGSEA) = store_colnames
    rownames(ssGSEA) = i
    p1 = pheatmap(ssGSEA,
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE,
                  color = cols,
                  annotation_col = data.frame(row.names = colnames(df), "Sample" = cell_samples[cell_filt] ),
                  main = i )
  }
  
  i = stringr::str_replace(i, "\\/", "_")
  ptitle = file.path(".", plotdir, paste0("metagenes_", i, ".pdf"))
  pdf(ptitle, width = 12, height = 8)
  print(p1)
  dev.off()

  }
}