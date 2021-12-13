## Author: Brian Gudenas
## Goal: perform NMF analysis on tumor single-cell data to find metagenes characterizing intra-subgroup variability
## Logic: Tumor single-cell data has large inter-sample variability therefore I identify NMF metagenes within each sample
# to avoid inter-sample bias, then i do an unsupervised merging of all sample NMF components within a tumor subgroup

NMF_samples = function(seurat_obj, genes, cell_samples, nrun=2, k=4) {
  library(NMF)
  # @seurat_obj processed seurat object
  # @genes eligible gene list for NMF (ie, proteing coding or genes with functional annotations)
  # @cell_samples  list of samples to subset for NMF within each sample
  # nrun = number of iterations to run nmf (default = 2 for testing) should be over 100 for real stuff
  stopifnot(length(cell_samples) == ncol(seurat_obj)) ## sample vec == number of cells
  stopifnot(sum(genes %in% rownames(seurat_obj)) > 1000) ## at least 1k genes remain after filter
  
# Run NMF to find metagenes per sample ------------------------------------
# meta_list is a list = length(unique(samples)) where the "gene_set" is stored and the gene "weights" are stored
  meta_list = list()
  for ( i in unique(cell_samples)){
    datExpr = as.matrix(seurat_obj@assays$RNA@data[ ,cell_samples == i ]  )
    datExpr = datExpr[rowSums(datExpr >= 1) >= 7, ] ## arbitrary filter to remove lowly expressed genes ( at least 1 norm count in 5 cells)
    datExpr = datExpr[rownames(datExpr) %in% genes, ]
    print(dim(datExpr))
    #  ptm <- proc.time()
    res <- nmf(datExpr, rank = k, nrun = nrun, seed=123456 )
    #  proc.time() - ptm
    metagenes = extractFeatures(res)
    for (j in 1:length(metagenes)){
        metagenes[[j]] = rownames(datExpr)[metagenes[[j]]]
        
        }
    meta_list[[i]]$gene_set = metagenes
    
    gweights = basis(res)
    #colnames(gweights) = paste("_", 1:k, sep="")
    meta_list[[i]]$weights = gweights
    print(paste0(Sys.time(), "  finished ---- ", i ))
  }
stopifnot(length(unique(cell_samples)) == length(meta_list))
return(meta_list)
}


# Define consensus metagenes ----------------------------------------------

NMF_consensus = function(seurat_obj, meta_list, cell_samples, cell_subgroups, topn=30){
  library(GSVA)
  # @seurat_obj processed seurat object
  # @meta_list list of sample NMF programs - output from NMF_samples
  # @cell_subgroups list of cell subgroups  to group cells for consensus analysis
  stopifnot(length(cell_subgroups) == ncol(seurat_obj))
  stopifnot(length(cell_samples) == ncol(seurat_obj))
  
  consensus_programs = list()
  for (i in unique(cell_subgroups)){
    print(i)
    GO = c()
    ## get two things, list of all NMF genes
    subgroup_samples = unique(cell_samples[ cell_subgroups == i ])
    subgroup_genesets = lapply(meta_list[subgroup_samples], "[[", 1 )
    for (samp_index in 1:length(subgroup_genesets)){
      sample_name = names(subgroup_genesets)[samp_index]
      samp_nmf_comps = subgroup_genesets[[sample_name]]
      names(samp_nmf_comps) = paste0(sample_name, "_", 1:length(samp_nmf_comps ))
      GO = c(GO, samp_nmf_comps)
    }
    if  ( all(is.na(unlist(subgroup_genesets))) ) { ##if no metagenes were found -- skip
     print(paste0( i," -- No genesets found" ))
      next
    }
# Score cells within each subgroup for metagenes --------------------------
    df = as.matrix(seurat_obj@assays$RNA@data[ ,cell_subgroups == i ]  )
    ssGSEA = gsva(df, gset.idx.list = GO, method="plage", min.sz=8, max.sz=300, ssgsea.norm = FALSE)
    ## How to merge?!? could do dynamic tree cut or just cor based
    hords = hclust(as.dist(1-cor(t(ssGSEA), method = "spearman")), method = "ward.D2")
    clusts = dynamicTreeCut::cutreeDynamicTree(hords, minModuleSize = 2)
    consensus_assignments = cbind(hords$labels, clusts) ## matrix 1/ 2 cols (samp_nmf_comps, cluster assignments)
    colnames(consensus_assignments) = c("samp_nmf_comps", "cluster")
    consensus_assignments = as.data.frame(consensus_assignments)
    final_metagenes = merge_clusters(meta_list, GO, consensus_assignments, topn=topn)
    consensus_programs[[i]] = final_metagenes
  }
  final = c()
  for (subgroup  in names(consensus_programs)){
      for (gs in 1:length(consensus_programs[[subgroup]])){
        new_name = paste0(subgroup, "_", gs)
        final[new_name] = consensus_programs[[subgroup]][gs]
      }
  }
  ## test consensus gene sets dont contain duplicate genes
  for (gs in final){
    stopifnot( sum(duplicated(gs)) == 0  ) 
    
  }
  return(final)
}


rename_columns = function(gene_weight_list){
  ## rename columns of gene weight matrices to sample_componentNumber
  samples = names(gene_weight_list)
  for ( i in 1:length(gene_weight_list)){
    colnames(gene_weight_list[[i]]) = paste0(samples[i], "_", 1:ncol(gene_weight_list[[i]]))
  }
  return(gene_weight_list) 
}


# Merge sample NMF components into consensus subgroup metagenes -----------

merge_clusters = function(meta_list, GO, consensus_assignments, topn=50){
  ## merge clusters
  final_metagenes = c()
  for (i in unique(consensus_assignments$cluster)){
    samples = unique(stringr::str_remove(consensus_assignments$samp_nmf_comps, "_[0-9]+"))
    stopifnot(samples > 0)
    
    components = consensus_assignments$samp_nmf_comps[consensus_assignments$cluster == i ]
    merged_metagene = unique(unlist(GO[names(GO) %in% components]))
    gene_weight_list = lapply(meta_list[names(meta_list) %in% samples], "[[", 2)
    gene_weight_list = rename_columns(gene_weight_list)
    
    ave_weights = c()
    for (j in 1:length(merged_metagene)){
      ## Loop through each gene in merged_metagene, extract weights from relevant comps and finally take median as metric
        gene = merged_metagene[j]
        g_weight_vec =c()
        for (df_index in 1:length(gene_weight_list)){
          sel_gene = 0
          if (gene %in% rownames(gene_weight_list[[df_index]])){
            sel_gene = gene_weight_list[[df_index]][gene, ]
            sel_gene = sel_gene[names(sel_gene) %in% components]
          }
          g_weight_vec = c(g_weight_vec, sel_gene)
          
        }
      ave_weights = c(ave_weights, median(g_weight_vec))
      
    }
    names(ave_weights) = merged_metagene
    ave_weights = sort(ave_weights, decreasing = TRUE)
    stopifnot( length(ave_weights) ==  length(merged_metagene))
    final_consensus_cluster = names(ave_weights)[1:topn]
    final_consensus_cluster = final_consensus_cluster[!is.na(final_consensus_cluster)]
    if (i != "0"){ ## 0 is the discard module from dynamic tree cut
      final_metagenes[[i]] = final_consensus_cluster
    }
  }
return(final_metagenes)
}
#final = merge_clusters(meta_list, GO, consensus_assignments, topn=30)
# TESTING -----------------------------------------------------------------
#map = readRDS("../../Annots/Annotables/hg38.rds")
# GO = fgsea::gmtPathways("../../Annots/GSEA/c5.bp.v7.0.symbols.gmt")
# genes = unique(unlist(GO))
# seurat_obj = readRDS("./Data/ss2_SO.rds")
# cell_samples = seurat_obj$Sample
# cell_subgroups = seurat_obj$Subgroup ## tumor groups
# 
# 
# meta_list = NMF_samples(seurat_obj = seurat_obj, 
#                         genes = genes, 
#                         cell_samples = cell_samples, 
#                         nrun = 2, k = 4)
# meta_list = readRDS("./")
#metagenes = NMF_consensus(seurat_obj = seurat_obj,
#                          meta_list = meta_list,
#                          cell_samples = cell_samples,
#                          cell_subgroups = cell_subgroups)
#
