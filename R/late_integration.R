#' Aggregate Community Labels from Multiple Views
#'
#' This function aggregates community labels from multiple views (or datasets) into a unified set of community labels using non-negative matrix factorization (NMF).
#'
#' @param community_labels_list A list of community labels for different views.
#' @param n_aggregate_communities Number of aggregated communities to generate.
#' @param max_iter Maximum number of iterations for NMF optimization. Default is 1000.
#' @param learning_rate Learning rate for the NMF optimizer. Default is 0.001.
#' @return A list containing:
#' \describe{
#'   \item{aggregated_membership}{Transposed matrix of the H factor from NMF, representing the aggregated community memberships.}
#'   \item{view_contributions}{Matrix of the contributions of each view to each aggregated community.}
#' }
#' @examples
#' \dontrun{
#' # Example usage:
#' x <- aggregate_communities(list('view1'=cl1, 'view2'=cl2), 3, max_iter = 50000)
#' print(x)
#' }
#' @export
aggregate_communities <- function(community_labels_list, n_aggregate_communities, max_iter=1000, learning_rate=0.001){
  membership <- get_stacked_membership_matrix(community_labels_list)
  result <- nmfbin::nmfbin(
    membership,
    n_aggregate_communities,
    init="nndsvd",
    optimizer="gradient",
    learning_rate=learning_rate,
    loss_fun="logloss",
    max_iter=max_iter
  )
  return(list(
    "aggregated_membership"=t(result$H),
    "view_contributions"=compute_view_contributions(result$W)
  ))
}

#' Get the most probable community from aggregated communities.
#'
#' This function determines the most probable community for each gene or node from an input
#' of aggregated community memberships. It selects the community with the highest membership
#' probability for each gene or node and orders the results.
#'
#' @param aggregated_communities the output of function aggregate_communities containing the aggregated membership probabilities.
#' @return A named integer vector where each value represents the most probable community for a gene or node,
#'         and the names are the gene or node identifiers.
#' @export
#' @examples
#' # Example usage:
#' most_probable_communities <- get_most_probable_community(aggregated_communities)
get_most_probable_community<-function(aggregated_communities){
  most_probable_communities<- apply(aggregated_communities$aggregated_membership, 1, which.max)
  most_probable_communities_ordered<-most_probable_communities[order(most_probable_communities)]
}


#' Create Stacked Membership Matrix from Community Labels List
#'
#' This function constructs a stacked membership matrix from the community labels of different views.
#'
#' @param community_labels_list A list of community labels for different views.
#' @return A matrix where rows represent genes and columns represent communities across all views.
#' @export
get_stacked_membership_matrix <- function(community_labels_list){
  unique_genes = get_unique_genes(community_labels_list)
  n_genes = length(unique_genes)
  
  stacked_membership = c()
  for(n in names(community_labels_list)){
    cl = community_labels_list[[n]]
    unique_communities <- sort(unique(unlist(cl)))
    communities_names <- sapply(
      unique_communities,
      function(x) { paste0(n, "_comm", as.character(x))}
    )
    n_communities <- length(unique_communities)
    memb = matrix(0, nrow=n_genes, ncol=n_communities)
    dimnames(memb) <- list(unique_genes, communities_names)
    for(g in names(cl)){
      memb[g, cl[[g]]] = 1
    }
    stacked_membership = rbind(stacked_membership, t(memb))
  }
  return(stacked_membership)
}

#' Extract Unique Genes from Community Labels List
#'
#' This function extracts a list of unique genes from the input list of community labels.
#'
#' @param community_labels_list A list of community labels for different views.
#' @return A vector of unique gene names.
#' @export
get_unique_genes <- function(community_labels_list){
  all_genes = c()
  for(m in community_labels_list){
    all_genes = c(all_genes, names(m))
  }
  return(sort(unique(all_genes)))
}

#' Compute Contributions of Each View to Aggregated Communities
#'
#' This function calculates the contributions of each view (original datasets) to each of the aggregated communities.
#'
#' @param W Matrix representing the contributions of the original community memberships.
#' @return A matrix representing the contributions of each view to each aggregated community.
#' @export
compute_view_contributions <- function(W){
  grouped <- stringr::str_split(rownames(W), "_", simplify=TRUE)[,1]
  unique_groups = sort(unique(grouped))
  norms <- colSums(W)
  contributions <- matrix(0, length(unique_groups), ncol(W))
  rownames(contributions) <- unique_groups
  for(g in unique_groups){
    sub_m <- W[which(grouped == g), ]
    sub_norm <- colSums(sub_m)
    contributions[g, ] <- sub_norm / norms
  }
  return(contributions)
}
