## Author: Brian Gudenas
## Goal: Function to annotate NMF meta-programs
## Logic: Perform hyper-geometric test to functionally annotate NMF metaprograms
#metagenes = readRDS("./test_data/metagenes_testdata.rds")

annotate_metagenes = function(metagenes, bp_db=  "GO_Biological_Process_2021", pval = 0.05){
library(enrichR)
library(stringr)
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
websiteLive <- TRUE
if (is.null(dbs)) websiteLive <- FALSE
stopifnot(websiteLive == TRUE)
stopifnot(is.list(metagenes) == TRUE)

meta_functions = c()
for (gs in 1:length(metagenes)){
  gene_list = unlist(metagenes[gs])
  enriched <- enrichr(gene_list, bp_db)
  Sys.sleep(1)
  
  enriched = enriched[[bp_db]]
  overlap_num = as.numeric(unlist(lapply(str_split(enriched$Overlap, "\\/"), "[[", 1)))
  enriched = enriched[enriched$Adjusted.P.value <= pval & overlap_num > 2, ] ## overlap should be at least 3 genes
  if (nrow(enriched) > 0 ){
      enriched$metagene = names(metagenes[gs])
      if (nrow(enriched) > 0 ) meta_functions = rbind(meta_functions, enriched)
      }
}
meta_functions$Subgroup =  as.character(unlist(lapply(str_split(meta_functions$metagene, "_"), "[[", 1)))
return(meta_functions)
}

#meta_functions = annotate_metagenes(metagenes)

plot_metagenes = function(meta_functions, plot_dir="./Figures/", suffix = "_GO"){
  dir.create(plot_dir, showWarnings = FALSE)
 library(ggplot2)
  library(dplyr)
  th <- theme(text = element_text(size=12, family = "Helvetica" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
  df  = meta_functions
  df$Term = stringr::str_trim(stringr::str_remove(df$Term, "\\(GO:[0-9]+\\)"))
  df$Subgroup = factor(df$Subgroup, levels = unique(df$Subgroup))
  df = df[order(df$Odds.Ratio, decreasing = TRUE), ]
  df = df[order(df$Subgroup), ]
  #df$Term = factor(df$Term, levels = unique(df$Term))
  
  plist = list()
  for ( i in unique(df$Subgroup )) {
    print(paste0("Plotting-----",i))
    sdf = df[df$Subgroup == i, ]
   # sdf$Term = droplevels(sdf$Term)
    sdf = sdf[order(sdf$metagene), ]
    sdf = sdf %>% 
      group_by(metagene) %>% 
      top_n(n = 12, wt = Odds.Ratio) %>% ungroup()
    sdf$Term = factor(sdf$Term, levels = unique(sdf$Term))
    g1 = ggplot(sdf, aes(x = Odds.Ratio, y= Term, fill = metagene )) +
      geom_bar(stat = "identity") +
      theme_bw() +
      th + 
      facet_wrap( ~Subgroup ) +
      ggtitle(i)
    ggsave(g1, device="pdf", filename = paste0(plot_dir, i,suffix, ".pdf") , width=10, height=8)
    plist[[i]] = g1
  }
  ## Plot All (Subgroup level)
  # df$Term = factor(df$Term, levels = unique(df$Term))
  # g1 = ggplot(df, aes(x = Odds.Ratio, y= Term, fill = Subgroup )) +
  #   geom_bar(stat = "identity") +
  #   theme_bw() +
  #   th
  # ggsave(g1, device="pdf", filename = paste0(plot_dir, "All_Subgroups",".pdf") , width=16, height=12)
  # 
  # return(g1)
}


#plot_metagenes(meta_functions)



plot_overlaps = function(megagenes, plot_dir){
  ## make heatmap of NMF metagene overlaps
  library(pheatmap)
  ndims = length(names(metagenes))
  adjmat = data.frame(matrix(nrow = ndims, ncol = ndims, data = 0))
  colnames(adjmat) = rownames(adjmat) = names(metagenes)
  
  
  for ( i in 1:nrow(adjmat)){
    gene_set_a_name = names(adjmat)[i]
    gene_set_a_genes = metagenes[[ gene_set_a_name]]
    for (j in 1:ncol(adjmat)){
      gene_set_b_name = names(adjmat)[j]
      gene_set_b_genes = metagenes[[ gene_set_b_name]]
      
      overlaps = sum(gene_set_a_genes %in% gene_set_b_genes)
      jaccards = overlaps/length(unique(c(gene_set_a_genes, gene_set_b_genes)))
      adjmat[i, j] = jaccards
    }
  }
  cols = colorRampPalette(colors = c("white","red","red4","black") )(100)
  p1 = pheatmap(adjmat,
                color = cols)
  pdf( paste0(plot_dir, "Metagene_overlap_heatmap.pdf"), width=8, height = 8)
  print(p1)
  dev.off()
  
  return(p1)
}
