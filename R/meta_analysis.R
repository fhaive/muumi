#' Batch-adjust and visualize multi-study expression data
#'
#' This function takes in an expression matrix, sample labels, and batch labels,
#' performs batch effect adjustment using the pamr package, and visualizes the
#' principal component analysis (PCA) before and after batch correction.
#' 
#' @param expr_mat A dataframe with genes on the rows and samples in the columns.
#' @param samples_label A factor of samples labels in the same order as the samples in expr_mat columns
#' @param batch_labels A factor of labels indicating the batches (the studies where the samples are coming from) 
#' 
#' @return A batch-adjusted expression matrix of the same dimension of expr_mat
#' 
#' @importFrom pamr pamr.batchadjust
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats prcomp
#' @examples
#' \dontrun{
#' # Example Usage:
#' adjusted_data <- multi_studies_adjust(expr_mat=expr_mat, samples_label=samples_label, batch_labels=batch_labels)
#' }
#' @export
multi_studies_adjust <- function(expr_mat, samples_label, batch_labels){
  mylist <- list(x=as.matrix(expr_mat), y=as.factor(samples_label), batchlabels=as.factor(batch_labels))
  adjusted_mat <- pamr::pamr.batchadjust(data = mylist)
  
  table_transpose<-as.data.frame(t(expr_mat))
  df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
  p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=batch_labels))
  p <- p+geom_point(size=3)+guides(color = guide_legend(title = "Batch labels"))+
    ggtitle("Before batch correction")+
    theme_bw()+
    theme(legend.key.size = unit(5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16, face="bold"), #change legend title font size
          legend.text = element_text(size=15),
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(color = "grey20", size=14),
          axis.text.y = element_text(color = "grey20", size = 14),
          axis.title=element_text(size=14,face="bold"))
  
  
  table_transpose2<-as.data.frame(t(adjusted_mat$x))
  df_pca2 <- prcomp(table_transpose2, center = TRUE, scale. = TRUE)
  p2<-ggplot(as.data.frame(df_pca2$x), aes(x=PC1, y=PC2, color=batch_labels))
  p2 <- p2+geom_point(size=3)+ guides(color = guide_legend(title = "Batch labels"))+
    ggtitle("After batch correction")+
    theme_bw()+
    theme(legend.key.size = unit(5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text( size=16, face="bold"), #change legend title font size
          legend.text = element_text(size=15),
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(color = "grey20", size=14),
          axis.text.y = element_text(color = "grey20", size = 14),
          axis.title=element_text(size=14,face="bold"))
  
  require(gridExtra)
  grid.arrange(p, p2, ncol=2)
  
  return(adjusted_mat)
}


#' Calculate and rank effect sizes from a meta-dataframe
#'
#' This function calculates effect sizes based on p-values in a meta-dataframe,
#' and ranks the genes by their effect sizes.
#'
#' @param meta_dataframe A dataframe with genes on the rows (as rownames) and samples in the columns. The columns should contain (adjusted) p-values deriving from gene-based statistical testing.
#'
#' @return A gene list ranked on the base of the effect size.
#'
#' @importFrom esc effect_sizes
#'
#' @examples
#' 
#' \dontrun{
#' # Example Usage:
#' ranked_effect_sizes <- calc_effect_size_rank(meta_dataframe)
#' }
#'
#' @export
calc_effect_size_rank<- function(meta_dataframe) {
  #index_pval<-which(grepl('pval',colnames(meta_dataframe)))
  #meta_data_pval<-meta_dataframe[,index_pval]
  #colnames(meta_data_pval)<-colnames(meta_dataframe)[index_pval]
  #heatmap(as.matrix(meta_data_pval))
  tmp_data<-meta_dataframe
  pval<-as.vector(rowMeans(tmp_data))
  tmp <- data.frame(
    pval = pval,
    n =rep(ncol(meta_dataframe),length(pval)),
    studyname = rownames(meta_dataframe)
  )
  effect_sizes_values<-esc::effect_sizes(tmp, p = pval, totaln = n, study = studyname, fun = "chisq")
  #check is study label is kept, in case add colnames(meta_data_log)
  effect_size_rank<-effect_sizes_values[,c('study','es')]
  rownames(effect_size_rank)<-effect_size_rank$study
  effect_size_rank$study<-NULL
  #ranked final list 
  effect_size_rank<-effect_size_rank[order(-effect_size_rank$es), , drop = FALSE]
  return(effect_size_rank)
}

#' Calculate p-value-based rank using Fisher-based method from a meta-dataframe for meta-analysis
#'
#' This function calculates p-value-based rank using a Fisher-based method using meta-dataframe as an input,
#' and ranks the genes based on this p-value rank.
#'
#' @param meta_dataframe A dataframe with genes on the rows (as rownames) and samples in the columns. The columns should contain (adjusted) p-values deriving from gene-based statistical testing.
#'
#' @return Ranked genes based on p-value based rank.
#'
#' @importFrom metap sumlog
#'
#' @examples
#' \dontrun{
#' # Example Usage:
#' ranked_pvalues <- calc_pvalue_based_rank(meta_dataframe)
#' }
#'
#' @export
calc_pvalue_based_rank<-function(meta_dataframe){
  #retrieve data
  #index_pval<-which(grepl('pval',colnames(meta_dataframe)))
  #meta_data_pval<-as.data.frame(meta_dataframe[,index_pval])
  #colnames(meta_data_pval)<-colnames(meta_dataframe)[index_pval]
  #for each gene combine the p-values by the sum of logs method
  fisher_based_res<-list()
  for(i in 1:length(rownames(meta_dataframe))){
    metap<-metap::sumlog(meta_dataframe[i,])
    fisher_based_res[[i]]<-metap$p
  } 
  fisher_based_pvalues<-as.data.frame(fisher_based_res)
  colnames(fisher_based_pvalues)<-rownames(meta_dataframe)
  fisher_based_pvalues<-t(fisher_based_pvalues)
  fisher_based_pvalues<-as.data.frame(fisher_based_pvalues)
  colnames(fisher_based_pvalues)<-"pValue"
  #ranked final list 
  fisher_based_pvalues<- fisher_based_pvalues[order(fisher_based_pvalues$pValue), , drop = FALSE]
  return(fisher_based_pvalues)
}

 
#' Calculate and the genes using Rank Product method from a meta-dataframe
#'
#' This function calculates gene ranks using the Rank Product method from a meta-dataframe.
#'
#' @param meta_dataframe dataframe with genes on the rows (as rownames) and samples in the columns. The columns should contain (adjusted) p-values deriving from gene-based statistical testing.
#' @param class A numeric vector indicating the class of each sample.  In the two class unpaired case, the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). For one class data, the label for each sample should be 1.
#' @param origin A character vector containing containing the origin labels of the samples.
#'
#' @return A ranked gene list using the Rank Product method.
#'
#' @importFrom RankProd RP.advance
#'
#' @examples
#' \dontrun{
#' # Example Usage:
#' ranked_pvalues <- calc_rank_base_rank(meta_dataframe, class, origin)
#' }
#'
#' @export
calc_rank_base_rank <- function(meta_dataframe, class, origin){
  index_pval_adj <- which(grepl('_adj_pval',colnames(meta_dataframe)))
  meta_data_adpval <- as.data.frame(meta_dataframe[,index_pval_adj])
  colnames(meta_data_adpval) <- colnames(meta_dataframe)[index_pval_adj]
  #origin contains the  labels for different studies
  #origin <- gsub(pattern = "_adj_pval",colnames(meta_data_adpval),replacement = "")
  #origin <- gsub("\\_.*","",origin)
  o <- origin
  RP_advance_out <- RankProd::RP.advance(data = meta_data_adpval, cl = class, origin = o, calculateProduct =T)
  #ranked  final list
  ranks_based_pvalues<-as.data.frame(RP_advance_out$pval)
  rownames(ranks_based_pvalues)<-rownames(meta_data_adpval)
  rank_1<-rownames(ranks_based_pvalues[order(abs(ranks_based_pvalues$`class1 < class2`)), , drop = FALSE])
  rank_2<-rownames(ranks_based_pvalues[order(abs(ranks_based_pvalues$`class1 > class2`)), , drop = FALSE])
  borda_list<-list()
  borda_list[[1]]<-rank_1
  borda_list[[2]]<-rank_2
  outputBorda<-Borda(borda_list)
  output_borda<-outputBorda$TopK$mean
  return(output_borda)
}


#' Computes a gene rank based on an ensembl of metanalysis methods, including effect size, p-value and rank product.
#' 
#' @importFrom RankProd RP.advance
#' @importFrom metap sumlog
#' @importFrom esc effect_sizes
#' 
#' @param meta_dataframe A dataframe with genes on the rows (as rownames) and samples in the columns. The columns should contain adjusted p-values deriving from gene-based statistical testing.
#' @param method Statistical method(s) to be included in the ensembl metanalysis.
#' @param class a vector containing the class labels of the samples. In the two class unpaired case, the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). For one class data, the label for each sample should be 1.
#' @param origin a vector containing the origin labels of the samples.
#' @param metric One statistical metric between "median" and "mean".
#' @return A gene list ranked on the base of the methods chosen for the metanalysis.
#' @examples 
#' \dontrun {
#' run_ensembl_metanalysis(meta_dataframe)
#' }
#' @export
run_ensembl_metanalysis <- function(meta_dataframe, method = c("effect_size", "pvalue", "rank_product"), class, origin, metric = "median") {
  if (length(method) < 2) {
    stop("Error: choose at least two methods among effect size, pvalue and rank product!")
  }
  
  data <- list()
  
  if ("effect_size" %in% method) {
    es <- calc_effect_size_rank(meta_dataframe)
    data[['Effect_size']] <- rownames(es)
  }
  if ("pvalue" %in% method) {
    pval <- calc_pvalue_based_rank(meta_dataframe)
    data[['Fisher_test']] <- rownames(pval)
  }
  if ("rank_product" %in% method) {
    rankprod <- calc_rank_base_rank(meta_dataframe, class, origin)
    data[['Rank_Prod']] <- rownames(rankprod)
  }
  
  outputBorda <- TopKLists::Borda(data)
  
  result <- switch(metric,
                   "median" = as.data.frame(outputBorda$TopK$median),
                   "mean" = as.data.frame(outputBorda$TopK$mean),
                   stop("Unknown metric specified.")
  )
  
  return(result)
}


#' Computes a GSEA of a gene rank against gene sets (provided through a GMT file).
#' 
#' @importFrom fgsea gmtPathways
#' @importFrom fgsea fgsea
#' @importFrom fgsea plotGseaTable
#' 
#' @param gene_list A ranked gene list.
#' @param gmt_file a vector containing the class labels of the samples. In the two class unpaired case, the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). For one class data, the label for each sample should be 1.
#' @examples 
#' \dontrun {
#' compute_gsea(gene_list, gmt_file = "c2.cp.reactome.v7.2.symbols.gmt")
#' }
#' @export
compute_gsea <- function(gene_list, gmt_file, no_permutations=10000){
  
  if (is.numeric(gene_list)==FALSE) {
    genes<-gene_list[,1]
    ranked_gene_list<-rev(seq(1, length(genes)))
    names(ranked_gene_list)<-genes
  }
  
  if (is.numeric(gene_list)==TRUE) {
    ranked_gene_list<-gene_list
    
  }
  
  mypath <- fgsea::gmtPathways(gmt_file)
  fgRes <- fgsea::fgsea(pathways = mypath, 
                        stats = ranked_gene_list,
                        minSize = 10,
                        #maxSize = 600,
                        nperm = no_permutations)
  
  topPathways <- head(fgRes$pathway[order(fgRes$pval)], 15)
  print(fgsea::plotGseaTable(pathways = mypath[topPathways], stats = ranked_gene_list, fgseaRes = fgRes, gseaParam = 0.5))
  return(fgRes)
}


#' Computes a threshold of the gene rank based on the enrichment score of the GSEA.
#' 
#' @importFrom fgsea calcGseaStat
#' 
#' @param gene_list result from the meta_analysis genelist or ranked genelist
#' @param fgsea_res fgsea object deriving from fgsea results or the compute_gsea function.
#' @param background gmt_file or Whole gene set object deriving from the fgsea::gmtPathways(gmt_file) function.
#' @examples 
#' \dontrun {
#' compute_gsea_thresh(gene_list, fgsea_res, background)
#' }
#' @export
compute_gsea_thresh <- function(geneList, fgsea_res, background){
  
  if (is.numeric(geneList)==FALSE) {
    genes<-geneList[,1]
    ranked_gene_list<-rev(seq(1, length(genes)))
    names(ranked_gene_list)<-genes
  }
  
  if (is.numeric(geneList)==TRUE) {
    ranked_gene_list<-geneList
    
  }
  
  
  if (class(background)=="list"){
      background<-background
  }
  
  if (class(background)=="character"){
    background<- fgsea::gmtPathways(background)
  }
  
  gseaParam=1
  stats <- ranked_gene_list
  fgsea_res <- fgsea_res[order(fgsea_res$pval, decreasing = FALSE),]
  max_vec <- c()
  sigpath <- which(fgsea_res$pval<0.05)
  for(i in 1:length(sigpath)){
    
    pathway <- background[[fgsea_res$pathway[i]]]
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway, 
                                   returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    #diff <- (max(tops) - min(bottoms))/8
    
    max_vec <- c(max_vec, which(names(ranked_gene_list) %in% names(gseaRes$tops)[which(gseaRes$tops==max(gseaRes$tops))]))
  }
  
  median_max_vec<-median(max_vec)
  
  return(median_max_vec)
}
