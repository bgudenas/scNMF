
nmf_preprocess = function(seurat_obj,
                          genes,
                          counts_slot = "scale.data",
                          n_cells = 50) {
  # @seurat_obj processed seurat object
  # @genes eligible gene list for NMF (ie, proteing coding or genes with functional annotations)
  stopifnot(sum(genes %in% rownames(seurat_obj)) > 1000) ## at least 1k genes remain after filter
    
    #seurat_obj = ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE, vars.to.regress = "nCount_RNA")
    # print("regress nCount")
    
    datExpr = as.matrix(seurat_obj@assays$RNA@counts )
    datExpr = datExpr[genes, ]
    #var_filt = matrixStats::rowVars(datExpr) > quantile(matrixStats::rowVars(datExpr), 0.10) ## remove bottom 10% of low variable genes
    var_filt = rowSums(datExpr > 0) >= n_cells ## expressed in at least n_cells
    genes = genes[genes %in% (rownames(datExpr)[var_filt]) ]
    
    if (counts_slot == "scale.data"){
      message("Scaled counts ---------")
      if ( sum(seurat_obj@assays$RNA@counts[ ,1]) > 800000 ){
        print("TPM detected ----")
        seurat_obj = ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE )
      } else {
        "Counts detected ----"
      seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      seurat_obj = ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE, vars.to.regress = "nCount_RNA")
      }
      datExpr = as.matrix(seurat_obj@assays$RNA@scale.data )
      #datExpr = datExpr + abs(min(datExpr)) ## add pseudocount to avoid negative values for scaled data
      datExpr[datExpr < 0 ] = 0
      print("no less than 0")
      print(quantile(datExpr))
      datExpr = datExpr[rownames(datExpr) %in% genes, ]
      print(dim(datExpr))
    } else if  (counts_slot == "data"){
      message("Normalized counts ---------")
      seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      datExpr = as.matrix(seurat_obj@assays$RNA@data )
      datExpr = datExpr[rownames(datExpr) %in% genes, ]
      print(quantile(datExpr))
      print(dim(datExpr))
    } else if  (counts_slot == "counts"){
      message("Raw counts ---------")
      datExpr = as.matrix(seurat_obj@assays$RNA@counts )
      datExpr = datExpr[rownames(datExpr) %in% genes, ]
      print(quantile(datExpr))
      print(dim(datExpr))
    }

stopifnot(sum(is.na(datExpr)) == 0) ## no NAs present
return(datExpr)
  
}






nmf_rank = function(datExpr,
                    data_dir = NULL,
                    data_prefix = NULL,
                    k_range = 2:8,
                    n_run = 30,
                    cell_labels = NULL) {
  ## how to select optimal (k), where cophenetic starts decreasing or RSS inflection
  library(NMF)
  stopifnot(!is.null(data_dir)) # specifcy outputs
  stopifnot(!is.null(data_prefix)) # specifcy output prefix
  stopifnot(!is.null(cell_labels)) # cell labels for consensus heatmap
  
  out_ranks = file.path(data_dir, paste0(data_prefix, "_NMF_ranks.rds"))
  if (!file.exists(out_ranks)){
    message("Computing NMF rnaks -----")
    estim.r <- nmf(datExpr, k_range, nrun=n_run, .opt='vp8', seed=123456)
    saveRDS(estim.r, out_ranks)
  } else {
    message("Reading in saved results -----")
    estim.r = readRDS(out_ranks)
  }
  
  if(requireNamespace("Biobase", quietly=TRUE)){
    
  pdf(file.path(data_dir, paste0(data_prefix, "_NMF_ranks_scatterplot.pdf")), width = 14, height = 14, pointsize = 10)
    print(plot(estim.r, cex = 0.5 ))
  dev.off()
  
  pdf(file.path(data_dir, paste0(data_prefix, "_NMF_ranks_heatmap.pdf")), width = 8, height = 10)
    consensusmap(estim.r, annCol=cell_labels, labCol=NA, labRow=NA)
    dev.off()
  }
  
  return(estim.r)
}




nmf_run = function(datExpr,
                   data_dir = NULL,
                   data_prefix = NULL,
                   k = 3,
                   n_run = 100) {
  library(NMF)
  
  out_file = file.path(data_dir, paste0(data_prefix, "_NMF_result.rds"))
  if (!file.exists(out_file)){
  message("STARTING NMF ---------")
  res <- nmf(datExpr,
             rank = k,
             nrun = n_run,
             .opt='vp8',
             seed=123456 )
  message("NMF FINISHED ---------")
  saveRDS(res, out_file)
  } else {
    message("LOADING NMF ---------")
   res = readRDS(out_file)
  }
  return(res)
}



nmf_features = function(nmf_res,
                        datExpr,
                        data_dir = NULL,
                        data_prefix = NULL,
                        feature_method = NULL,
                        topn = 30,
                        MT_filt = FALSE){
  library(NMF)
  meta_list = list()
  metagenes = extractFeatures(nmf_res, topn)
  message(paste0("Features (Topn) =", topn))
  if (!is.null(feature_method)){
    print(paste0("FEATURE EXTRACTION =", feature_method))
    metagenes = extractFeatures(nmf_res, method = feature_method)
  } else {
    metagenes = extractFeatures(nmf_res, 200) ## arbitrarily large num
  }
  for (j in 1:length(metagenes)){
    metagene_list =  rownames(datExpr)[metagenes[[j]]]
    if (MT_filt == TRUE){
      metagene_list = metagene_list[!grepl("^MT-", metagene_list)]
      metagene_list = metagene_list[!grepl("^MRP", metagene_list)]
      metagene_list = metagene_list[!grepl("^RPS", metagene_list)]
      metagene_list = metagene_list[!grepl("^RPL", metagene_list)]
      metagene_list = metagene_list[1:topn]
    }
    metagene_list = metagene_list[1:topn]
    metagenes[[j]] =metagene_list
  }
  meta_list[[i]]$gene_set = metagenes
  
  gweights = basis(nmf_res)
  meta_list[[i]]$weights = gweights
  
  saveRDS(meta_list, file.path(data_dir, paste0(data_prefix, "_NMF_metamatrix.rds")))
return(meta_list)
} 
  