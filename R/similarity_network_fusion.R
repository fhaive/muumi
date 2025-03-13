#' Integrate Gene Expression Matrices Using Similarity Network Fusion (SNF)
#'
#' This function integrates multiple gene expression matrices using the Similarity Network Fusion (SNF) algorithm. 
#' The input matrices should have genes on the rows and samples on the columns.
#'
#' @param list_of_gene_expression_matrices A list of gene expression matrices. The matrices should have genes on the rows and samples on the columns. At least two matrices should be provided.
#' @param K An integer specifying the number of neighbors in the K-nearest neighbors part of the algorithm. Default is 20. For more details, see the SNFtool documentation.
#' @param sigma A numeric value specifying the variance for the local model. Default is 0.5. For more details, see the SNFtool documentation.
#' @param t An integer specifying the number of iterations for the diffusion process. Default is 20. For more details, see the SNFtool documentation.
#'
#' @return A fused affinity matrix computed by the SNF algorithm. If the input validation fails (less than two matrices or no common genes), the function returns -1.
#'
#' @export
#' 
snf_based_integration = function(list_of_gene_expression_matrices, K = 20, sigma = 0.5, t = 20){
  
  if(length(list_of_gene_expression_matrices)<2){
    print("Error: At least two matrices should be provided!")
    return(-1)
  }
  
  # identify common genes between all the co-expression matrices
  common_genes = rownames(list_of_gene_expression_matrices[[1]])
  for (i in 2:length(list_of_gene_expression_matrices)) {
    common_genes = intersect(common_genes, rownames(list_of_gene_expression_matrices[[i]]))
  }
  
  if (length(common_genes) == 0) {
    print("Error: no common genes in the matrices")
    return(-1)
  }
  
  dist_list = list()
  affinity_matrix_list = list()
  for (i in 1:length(list_of_gene_expression_matrices)) {
    list_of_gene_expression_matrices[[i]] = list_of_gene_expression_matrices[[i]][common_genes,]
    dist_list[[i]] = as.matrix(dist( list_of_gene_expression_matrices[[i]]))
    affinity_matrix_list[[i]] = affinityMatrix(dist_list[[i]] , K = K, sigma = sigma)
    
  }
  
  names(dist_list) = names(list_of_gene_expression_matrices)
  names(affinity_matrix_list) = names(list_of_gene_expression_matrices)
  W = SNF(affinity_matrix_list, K, t)
  
  return(W)
  
}


#' Prune Similarity Network Fusion (SNF) Network
#'
#' This function prunes an SNF network by setting edges above a specified percentile threshold to 1
#' and returning a weighted adjacency matrix where values correspond to the original edge weights.
#'
#' @param snf_network A square, symmetric numeric matrix representing the SNF network.
#' @param percentile_thr A numeric value between 0 and 1 specifying the percentile threshold. Edges with weights above this percentile are retained in the pruned network.
#' @return A pruned weighted adjacency matrix of the same dimensions as `snf_network`, with values corresponding to retained edge weights.
#' @examples
#' # Generate an example SNF network
#' snf_network <- matrix(runif(100, min = 0, max = 1), nrow = 10)
#' snf_network <- (snf_network + t(snf_network)) / 2  # Ensure symmetry
#' result <- prune_snf_network(snf_network, percentile_thr = 0.75)
#' @export
prune_snf_network <- function(snf_network, percentile_thr = 0.75) {
  
  if (!is.matrix(snf_network)) {
    stop("Input `snf_network` must be a matrix.")
  }
  if (nrow(snf_network) != ncol(snf_network)) {
    stop("Input `snf_network` must be a square matrix.")
  }
  if (!is.numeric(snf_network)) {
    stop("Input `snf_network` must contain numeric values.")
  }
  if (!isSymmetric(snf_network)) {
    stop("Input `snf_network` must be symmetric.")
  }
  if (!is.numeric(percentile_thr) || length(percentile_thr) != 1 || percentile_thr < 0 || percentile_thr > 1) {
    stop("`percentile_thr` must be a single numeric value between 0 and 1")
  }
  
  new_mat <- matrix(data = 0, nrow = nrow(snf_network), ncol = ncol(snf_network), dimnames = list(rownames(snf_network), colnames(snf_network)))
  
  for (i in seq_len(nrow(snf_network))) {
    thr <- quantile(snf_network[i, ], probs = percentile_thr, na.rm = TRUE)
    new_mat[i, snf_network[i, ] > thr] <- 1
    new_mat[, i][snf_network[, i] > thr] <- 1
  }
  if (isSymmetric.matrix(new_mat)) {
    message("Returning symmetric matrix...")
  } else {
    warning("Resulting matrix is not symmetric.")
  }
  
  new_mat_weighted <- new_mat * snf_network
  
  return(new_mat_weighted)
}


#' Generate a Ranked Edge List from an SNF Result Matrix
#'
#' This function generates a ranked edge list from an SNF (Similarity Network Fusion) result matrix.
#' The edge list is derived from the upper triangle of the matrix, listing gene pairs and their similarity values.
#' Only edges with positive similarity values are retained and ranked in descending order by similarity.
#'
#' @param snf_result A square, symmetric numeric matrix representing the SNF similarity network with named rows and columns.
#' @return A named numeric vector where each element represents the similarity value between a gene pair, with names formatted as "gene1_gene2".
#' The vector is ranked in descending order by similarity.
#' @examples
#' # Create an example SNF result matrix
#' snf_result <- matrix(runif(100, min = 0, max = 1), nrow = 10)
#' snf_result <- (snf_result + t(snf_result)) / 2  # Ensure symmetry
#' rownames(snf_result) <- colnames(snf_result) <- paste0("Gene", 1:10)
#' ranked_edges <- get_ranked_edge_list(snf_result)
#' @export
get_ranked_edge_list <- function(snf_result) {
  
  if (!is.matrix(snf_result)) {
    stop("Input `snf_result` must be a matrix.")
  }
  if (nrow(snf_result) != ncol(snf_result)) {
    stop("Input `snf_result` must be a square matrix.")
  }
  if (!is.numeric(snf_result)) {
    stop("Input `snf_result` must contain numeric values.")
  }
  if (!isSymmetric(snf_result)) {
    stop("Input `snf_result` must be symmetric.")
  }
  if (is.null(rownames(snf_result)) || is.null(colnames(snf_result))) {
    stop("Input `snf_result` must have row and column names.")
  }
  
  num_genes <- nrow(snf_result)
  gene_names <- rownames(snf_result)
  
  edges <- which(upper.tri(snf_result, diag = FALSE), arr.ind = TRUE)
  edge_list <- data.frame(
    gene1 = gene_names[edges[, 1]],
    gene2 = gene_names[edges[, 2]],
    similarity = snf_result[edges]
  )
  
  edge_list <- edge_list[edge_list$similarity > 0, ]
  edge_list <- edge_list[order(-edge_list$similarity), ]
  edge_list_to_esea <- edge_list$similarity
  names(edge_list_to_esea) <- paste(edge_list$gene1, edge_list$gene2, sep = "_")
  
  return(edge_list_to_esea)
}
