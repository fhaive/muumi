#' Perform pathway enrichment analysis from modules using Reactome.
#'
#' This function performs pathway enrichment analysis using Reactome for gene modules.
#' It takes an input of modules (either as communities or integer vector), maps gene symbols
#' to ENTREZ IDs using `clusterProfiler::bitr`, and enriches pathways using `ReactomePA::enrichPathway`.
#' Results can be exported either as individual CSV files for each module (if `layout = "single"`)
#' or combined into a single CSV file (`Pathway_results_overall.csv`) if `layout = "overall"`.
#'
#'
#' @importFrom igraph membership
#' @importFrom clusterProfiler bitr
#' @importFrom ReactomePA enrichPathway
#' @param modules An object containing gene modules. Should be either igraph object of class 'communities' (i.e. from get_modules function)
#'or an integer vector (i.e. from get_most_probable_community function) representing modules.
#' @param geneID The type of gene ID used in input. Default is "SYMBOL".
#' @param pval_cutoff The p-value cutoff for pathway enrichment. Default is 0.05.
#' @param outPath Path where the output files will be exported.
#' @param layout The layout type for exporting results. Either "overall" or "single". Default is "overall".
#' @return A data frame containing enriched pathways and associated statistics.
#' @export
#' @examples
#' # Example usage:
#' get_reactome_from_modules(most_probable_communities, geneID="SYMBOL", pval_cutoff=0.05,
#'                           outPath="/path/to/output", layout="overall")
get_reactome_from_modules <- function(modules, geneID = "SYMBOL", pval_cutoff = 0.05, outPath, layout = "overall") {
  # Create directory if it doesn't exist
  if (!file.exists(outPath)) {
    dir.create(outPath, recursive = TRUE)
  }
  cat("The files will be exported in ", outPath, "\n")
  
  # Helper function to process each module
  process_module <- function(mod, members, layout, outPath) {
    x <- names(members[members == mod])
    if (length(x) > 10) {
      eg <- clusterProfiler::bitr(x, fromType = geneID, toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      print(head(eg))
      sigpath <- ReactomePA::enrichPathway(gene = eg$ENTREZID, pvalueCutoff = pval_cutoff, readable = TRUE)
      sigpath <- as.data.frame(sigpath)
      print(head(sigpath))
      if (nrow(sigpath) > 0) {
        sigpath$Module <- mod
        sigpath <- sigpath[order(sigpath$p.adjust), ]  # Sort by p-value within the module
        if (layout == "overall") {
          return(sigpath)
        } else if (layout == "single") {
          write.csv(sigpath, file = file.path(outPath, paste0("Significant_enriched_pathways_module_", mod, ".csv")), 
                    row.names = FALSE, quote = TRUE)
        }
      }
    }
    return(NULL)
  }
  
  # Main processing
  sigpath.overall <- data.frame()
  
  if (class(modules) == "communities" || is.integer(modules)) {
    members <- if (class(modules) == "communities") igraph::membership(modules) else modules
    unique_modules <- sort(unique(members))  # Ensure modules are sorted numerically
    
    for (mod in unique_modules) {
      result <- process_module(mod, members, layout, outPath)
      if (!is.null(result) && layout == "overall") {
        sigpath.overall <- rbind(sigpath.overall, result)
      }
    }
    
    if (layout == "overall" && nrow(sigpath.overall) > 0) {
      write.csv(sigpath.overall, file = file.path(outPath, "Pathway_results_overall.csv"), 
                row.names = FALSE, quote = TRUE)
    }
  } else {
    stop("Invalid input: modules should be of class 'communities' or an integer vector.")
  }
  
  return("Analysis completed!!!")
}

#' Generate bubble plot of enriched pathways from gene modules.
#'
#' This function generates a bubble plot of enriched pathways using `clusterProfiler::compareCluster`
#' and `clusterProfiler::dotplot` for gene modules provided as input. It maps gene symbols to ENTREZ IDs
#' using `clusterProfiler::bitr`, and performs pathway enrichment analysis using `clusterProfiler::compareCluster`.
#'
#' @param modules An object containing gene modules. Should be either igraph object of class 'communities' (i.e. from get_modules function)
#'or an integer vector (i.e. from get_most_probable_community function) representing modules.
#' @param geneID The type of gene ID used in input. Default is "SYMBOL".
#' @return A plot showing enriched pathways.
#' @export
#' @examples
#' # Example usage:
#' get_bubbleplot_from_pathways(most_probable_communities, geneID="SYMBOL")
#' @import clusterProfiler
#' @importFrom clusterProfiler compareCluster dotplot
get_bubbleplot_from_pathways <- function(modules, geneID="SYMBOL") {
  
  lst <- list()  # Initialize an empty list to store results
  
  if (class(modules) == "communities") {
    members <- igraph::membership(modules)
    for(mod in 1:length(modules)){
      x <- names(members[members == mod])  # Get gene names for the current module
      if(length(x) > 10) {
        convgenes <- clusterProfiler::bitr(x, fromType=geneID, toType="ENTREZID", OrgDb=org.Hs.eg.db)
        print(head(convgenes))
        lst[[mod]] <- convgenes$ENTREZID  # Store ENTREZIDs in the list
      }
    }
    
    names(lst) <- seq_along(lst)
    lst[sapply(lst, is.null)] <- NULL
    
    res <- clusterProfiler::compareCluster(lst, fun="enrichPathway")
    print(clusterProfiler::dotplot(res))
    return(res)
  }
  
  if (class(modules) == "integer") {
    members <- modules
    for(mod in 1:length(unique(modules))){
      x <- names(members[members == mod])  # Get gene names for the current module
      if(length(x) > 10) {
        convgenes <- clusterProfiler::bitr(x, fromType=geneID, toType="ENTREZID", OrgDb=org.Hs.eg.db)
        print(head(convgenes))
        lst[[mod]] <- convgenes$ENTREZID  # Store ENTREZIDs in the list
      }
    }
    
    names(lst) <- seq_along(lst)
    lst[sapply(lst, is.null)] <- NULL
    
    res <- clusterProfiler::compareCluster(lst, fun="enrichPathway")
    print(clusterProfiler::dotplot(res))
    return(res)
  }
  
  # If neither "communities" nor "integer" class, return an error
  stop("Invalid input. Expected 'communities' or 'integer' class for modules.")
}

