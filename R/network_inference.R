#' Infer correlation matrix from gene expression table by using mutual information.
#'
#' calculate_correlation_matrix uses the MINET package to create correlation matrix by mutual information method. User can specify
#' the inference algorithms, correlation calculation methods and discretization methods to create mutiple combinations of parameters 
#' for multiple runs of minet(). Multiple inferences are unified by taking median to create a consensus matrix.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom minet minet
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach getDoParWorkers getDoParName %dopar%
#' @importFrom plyr aaply laply
#' @importFrom stats quantile median
#' @importFrom utils capture.output
#'
#' @param gx_table Gene expression table as a data frame.
#' @param iMethods Vector of valid inference algorithms for MINET package.
#' @param iEst Vector of valid correlation methods for MINET package.
#' @param iDisc Vector of valid discretization methods for MINET package.
#' @param ncores Number of cores for running instances of MINET in parallel default:2.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A binary symmetrix matrix representing the median of mutual information correlation computed across various MINET combinations
#' @examples
#' \dontrun{
#' calculate_correlation_matrix(gx_table=gene_expression.df,
#' iMethods=c("clr","aracne","mrnet","mrnetb"),
#' iEst=c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
#' iDisc=c("none","equalfreq","equalwidth","globalequalwidth"),
#' ncores=12,
#' debug_output=TRUE
#' )
#' }
#' @keywords internal
#' @export
calculate_correlation_matrix <- function(gx_table, iMethods, iEst, iDisc, ncores=2, debug_output=FALSE, updateProgress=NULL){
  parList <- list()
  out.GX.list<- list()
  out.tmp.list<- list()
  out.list<- list()
  tGX <- t(gx_table)
  stdGX <- scale(tGX)
  
  if(is.null(iMethods)){
    stop("Please Select atLeast One Inference Algorithm!")
  }
  
  if(is.null(iEst)){
    stop("Please Select At Least One Correlation!")
  }
  
  if(is.null(iDisc)){
    stop("Please Select At Least One Discretization Method!")
  }
  
  cntr <- 0
  for(i in 1:length(iMethods)){
    #print(methods2[i])
    for(j in 1:length(iEst)){
      #print(paste("-",est.opt2[j]))
      for(k in 1:length(iDisc)){
        cntr <- cntr + 1
        parList[[cntr]] <- list()
        parList[[cntr]][["mt"]] <- iMethods[i]
        parList[[cntr]][["est"]] <- iEst[j]
        parList[[cntr]][["disc"]] <- iDisc[k]
      }
    }
  }
  #print(paste0("List of MINET Combinations: ", parList))
  
  if(ncores > parallel::detectCores()){
    ncores <- parallel::detectCores()
  }
  
  cl <- parallel::makeCluster(ncores, outfile="")
  doParallel::registerDoParallel(cl)
  print(paste("Total Cores: ", parallel::detectCores(), sep=""))
  print(paste("Total Workers: ", foreach::getDoParWorkers(), sep=""))
  print(paste("DoPar Name: ",  foreach::getDoParName(), sep=""))
  
  print(paste0("Before For Each, ", "Number of Combinations:", length(parList)))
  utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-log.txt", append=TRUE)
  utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-to-run.txt", append=TRUE)
  utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-completed.txt", append=TRUE)
  #out.tmp.list <- foreach(i=1:length(parList)) %do% {
  out.tmp.list <- foreach::foreach(i=1:length(parList)) %dopar% {
    #capture.output(print("Start For Each"), file="minet-log.txt", append=T)
    mt <- parList[[i]]$mt
    est <- parList[[i]]$est
    disc <- parList[[i]]$disc
    #print(parList[[i]]$mt)
    
    if((est == "mi.empirical" | est == "mi.mm" | est == "mi.shrink" | est == "mi.sg") & disc == "none") {
      miMat <- -1
      miMatName <- "np"
    }
    else{
      if(debug_output==TRUE)
        utils::capture.output(print(paste("----",mt,est,disc, sep="__")), file="minet-log.txt", append=TRUE)
      
      utils::capture.output(print(paste0("Iteration-", i, ": ", mt, "-", est, "-", disc)), file="minet-to-run.txt", append=TRUE)
      #print("Before Updating text")
      #if (is.function(updateProgress)){
      #	text <- paste("MINET: ", mt, "-", est, "-", disc, "-", sep="")
      #	print("Updating text")
      #	updateProgress(detail = text)
      #}
      
      ptm <- proc.time()
      miMat <- minet::minet(stdGX, method=mt, estimator=est, disc=disc)
      utils::capture.output(print(paste0("Iteration-", i, ", ", mt, "-", est, "-", disc, ": ", "MINET Execution Time - ", round(proc.time() - ptm)[3], " sec")), file="minet-completed.txt", append=TRUE)
      #capture.output(print(proc.time() - ptm), file="minet-log.txt", append=T)
      miMatName <- paste(mt,est,disc,sep="__")
    }
    out.list <- list("mat"=miMat, "name"=miMatName)
  }
  utils::capture.output(print("For Each Finished, Stopping Cluster..."), file="minet-log.txt", append=TRUE)
  parallel::stopCluster(cl)
  
  for(i in 1:length(out.tmp.list)){
    out.GX.list[[i]] <- out.tmp.list[[i]]$mat
    names(out.GX.list)[[i]] <- out.tmp.list[[i]]$name
  }
  
  #if(debug_output==TRUE)
  #print(paste("Got minet list of lists --- ", str(out.GX.list)))
  
  #removing ones having no values
  #llGX <- out.GX.list[-(which(names(out.GX.list)=="np"))]
  idx <- which(names(out.GX.list)=="np")
  if(length(idx)>0){
    llGX <- out.GX.list[-idx]
  }else{
    llGX <- out.GX.list
  }
  
  if(length(llGX)>1){
    print("Making Median Matrix...")
    arr1<-abind::abind(llGX,along=3)
    matGX <- apply(arr1,c(1,2),median)
  }else{
    print("Only one matrix, not computing median!")
    matGX <- llGX[[1]]
  }
  
  print("Return Matrix")
  return(matGX)
}



#' Create edge ranked inference matrix by combining information from different inference algorithms with the help of internal function calculate_correlation_matrix().
#'
#' get_ranked_consensus_matrix uses the internal function calculate_correlation_matrix() to get a single consensus matrix per inference algorithm. User can specify
#' the inference algorithms, correlation calculation methods and discretization methods, a combination of parameters will be created per inference algorithm to run calculate_correlation_matrix(), 
#' this will generate a consensus matrix per inference algorithm. The consesus matrices from different inference algorithms are used to create a single binary matrix by rank based selection of edges.
#'
#' @importFrom TopKLists Borda
#'
#' @param gx_table Gene expression table as a data frame.
#' @param iMethods Vector of valid inference algorithms for MINET package.
#' @param iEst Vector of valid correlation methods for MINET package.
#' @param iDisc Vector of valid discretization methods for MINET package.
#' @param ncores Number of cores for running instances of MINET in parallel default:2.
#' @param edge_selection_strategy How to select top ranked edges default:default. By default selectis top edges until all the nodes have at least one edge.
#' @param topN Top n percentage edges to select if edge_selection_strategy is 'top' default:10.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A symmetrix matrix with median edge ranks representing the edge rank based consensus from different inference algorithms.
#' @examples
#' \dontrun{
#' get_ranked_consensus_matrix(gx_table=gene_expression.df,
#' iMethods=c("clr","aracne","mrnet","mrnetb"),
#' iEst=c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
#' iDisc=c("none","equalfreq","equalwidth","globalequalwidth"),
#' ncores=12,
#' edge_selection_strategy="default",
#' topN=10,
#' debug_output=TRUE
#' )
#' }
#' @keywords internal
#' @export
#get_ranked_consensus_binary_matrix <- function(gx_table, iMethods, iEst, iDisc, ncores=2, debug_output=FALSE, updateProgress=NULL){
get_ranked_consensus_matrix <- function(gx_table=NULL, iMethods=NULL, iEst=NULL, iDisc=NULL, ncores=2, matList=NULL, mat_weights="rank", ensemble_strategy="minet", debug_output=FALSE, updateProgress=NULL){
  mat_ll <- list()
  ranked_edges_ll <- list()
  
  if(length(grep("minet" ,ensemble_strategy))>0){
    mthdCount <- 1
    totalMthds <- length(iMethods)+1
    for(mthd in iMethods){
      if (is.function(updateProgress)) {
        text <- paste0("'", mthd, "' Median")
        value <- mthdCount / totalMthds
        updateProgress(detail = text, value = value)
      }
      mthdCount <- mthdCount + 1
      
      print(paste0("Calculate correlation matrix for method : ", mthd))
      mat_ll[[mthd]] <- calculate_correlation_matrix(gx_table=gx_table, iMethods=mthd, iEst=iEst, iDisc=iDisc, ncores=ncores)
      
      print(paste0("Get ranked edges for method : ", mthd))
      mat_ll[[mthd]][lower.tri(mat_ll[[mthd]], diag=TRUE)] <- NA
      
      edge_df <- as.data.frame(as.table(mat_ll[[mthd]]))
      #print("Minet edge_df before:")
      #print(str(edge_df))
      edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
      edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
      #print("Minet edge_df after:")
      #print(str(edge_df))
      
      ranked_edges_ll[[mthd]] <- edge_df[order(edge_df$weight, decreasing=TRUE), "edge"]
    }
  }
  
  if(length(grep("user", ensemble_strategy))>0){
    matList <- as.list(matList)
    for(i in c(1:length(matList))){
      itmName <- paste0("user_", i)
      matList[[i]][lower.tri(matList[[i]], diag=TRUE)] <- NA
      
      edge_df <- as.data.frame(as.table(matList[[i]]))
      #print("User edge_df before:")
      #print(str(edge_df))
      edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
      edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
      hasNulls <- FALSE
      if(mat_weights=="rank"){
        nullIdx <- which(edge_df$weight==0)
        if(length(nullIdx)>0){
          edge_df_null <- edge_df[nullIdx,]
          edge_df <- edge_df[-nullIdx,]
          hasNulls <- TRUE
        }
      }
      #print("User edge_df after:")
      #print(str(edge_df))
      
      ordOption <- ifelse(mat_weights=="rank", FALSE, TRUE)
      #print(paste0("ordOption : ", ordOption))
      edge_df <- edge_df[order(edge_df$weight,decreasing=ordOption),]
      if(hasNulls){
        edge_df <- as.data.frame(rbind(edge_df, edge_df_null), stringsAsFactors=FALSE)
      }
      #ranked_edges_ll[[itmName]] <- edge_df[order(edge_df$weight, decreasing=ordOption), "edge"]
      ranked_edges_ll[[itmName]] <- edge_df$edge
    }
  }
  #print("ranked_edges_ll:")
  #print(length(ranked_edges_ll))
  #print(names(ranked_edges_ll))
  #print(lapply(ranked_edges_ll, length))
  
  if (is.function(updateProgress)) {
    updateProgress(detail = "Consensus Binary", value = 1)
  }
  
  print("Perform Borda on list of list of ranked edges.")
  borda_res <- TopKLists::Borda(ranked_edges_ll)
  #Borda.plot(borda_res)
  
  print("Get a consensus binary matrix by selecting the most significant ranked edges from median rank of Borda result.")
  if(length(mat_ll)>0){
    rank_mat <- mat_ll[[1]]
  }else{
    rank_mat <- matList[[1]]
  }
  rank_mat[,] <- 0
  #print("rank_mat")
  #print(dim(rank_mat))
  
  median_list <- borda_res$TopK$median
  #Kendall.plot(ranked_edges_ll, median_list)
  
  #if(edge_selection_strategy=="default"){
  #	#input_genes <- dim(gx_table)[1]
  #	input_genes <- nrow(rank_mat)
  #	genes <- NULL
  #	total_genes <- 0
  #	cutoffIdx <- NULL
  #	for(i in c(1:length(median_list))){
  #		if(total_genes<input_genes){
  #			local_genes <- strsplit(median_list[i], ";")[[1]]
  #			rank_mat[local_genes[1],local_genes[2]] <- i
  #			rank_mat[local_genes[2],local_genes[1]] <- i
  #			genes[local_genes[1]] <- 1
  #			genes[local_genes[2]] <- 1
  #			total_genes <- length(genes)
  #		}else{
  #			cutoffIdx <- i-1
  #			break
  #		}
  #	}
  #}else if(edge_selection_strategy=="top"){
  #	cutOff <- round((as.numeric(topN)*length(median_list))/100)
  #	for(i in c(1:length(median_list))){
  #		if(i<=cutOff){
  #			local_genes <- strsplit(median_list[i], ";")[[1]]
  #			rank_mat[local_genes[1],local_genes[2]] <- i
  #			rank_mat[local_genes[2],local_genes[1]] <- i
  #		}else{
  #			break
  #		}
  #	}
  #}
  
  #print("median_list:")
  #print(str(median_list))
  #print(length(median_list))
  #print(head(median_list))
  for(i in c(1:length(median_list))){
    local_genes <- strsplit(median_list[i], ";")[[1]]
    rank_mat[local_genes[1],local_genes[2]] <- i
    rank_mat[local_genes[2],local_genes[1]] <- i
  }
  
  print("Rank matrix computed, returning!")
  return(rank_mat)
}


#' Get edge rank list and binary inference matrix from edge rank matrix computed by get_ranked_consensus_matrix().
#'
#' parse_edge_rank_matrix parses the edge rank matrix created by using the internal function get_ranked_consensus_matrix_matrix() to get a ranked edge list and a binary matrix.
#'
#' @importFrom TopKLists Borda
#'
#' @param edge_rank_matrix A symmetrix matrix with edge ranks as weight.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A list containing a vector of consensus edge ranks and a binary symmetrix matrix representing the edge rank matrix.
#' @examples
#' \dontrun{
#' parse_edge_rank_matrix <- function(edge_rank_matrix, debug_output=FALSE, updateProgress=NULL)
#' }
#' @keywords internal
#' @export
parse_edge_rank_matrix <- function(edge_rank_matrix, edge_selection_strategy="default", mat_weights="rank", topN=10, debug_output=FALSE, updateProgress=NULL){
  bin_mat <- edge_rank_matrix
  #idx <- which(bin_mat>0)
  #if(length(idx)>0){
  #        print("Creating binary matrix...")
  #        bin_mat[which(bin_mat>0)] <- 1
  #}
  bin_mat[,] <- 0
  
  print("Getting edge list ordered by rank...")
  rank_matrix <- edge_rank_matrix
  rank_matrix[lower.tri(rank_matrix, diag=TRUE)] <- NA
  edge_df <- as.data.frame(as.table(rank_matrix))
  edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
  edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
  edge_df <- edge_df[which(edge_df$weight>0),]
  ordOption <- ifelse(mat_weights=="rank", FALSE, TRUE)
  print(paste0("ordOption : ", ordOption))
  print(class(edge_df$weight))
  edge_rank <- edge_df$edge[order(edge_df$weight, decreasing=ordOption)]
  print(paste0("edge rank before:", length(edge_rank)))
  print(head(edge_rank))
  
  if(edge_selection_strategy=="default"){
    #input_genes <- nrow(rank_matrix)
    input_genes <- length(unique(unlist(strsplit(edge_rank, ";"))))
    print(paste0("input genes : ", input_genes))
    genes <- NULL
    total_genes <- 0
    cutoffIdx <- NULL
    for(i in c(1:length(edge_rank))){
      if(total_genes<input_genes){
        local_genes <- strsplit(edge_rank[i], ";")[[1]]
        bin_mat[local_genes[1],local_genes[2]] <- 1
        bin_mat[local_genes[2],local_genes[1]] <- 1
        genes[local_genes[1]] <- 1
        genes[local_genes[2]] <- 1
        total_genes <- length(genes)
        cutoffIdx <- i
      }else{
        print(paste0("Cutoff before break:", cutoffIdx))
        break
      }
    }
  }else if(edge_selection_strategy=="top"){
    cutoffIdx <- round((as.numeric(topN)*length(edge_rank))/100)
    for(i in c(1:length(edge_rank))){
      if(i<=cutoffIdx){
        local_genes <- strsplit(edge_rank[i], ";")[[1]]
        bin_mat[local_genes[1],local_genes[2]] <- 1
        bin_mat[local_genes[2],local_genes[1]] <- 1
      }else{
        break
      }
    }
  }
  print(paste0("Cutoff:", cutoffIdx))
  edge_rank <- edge_rank[c(1:cutoffIdx)]
  print(paste0("edge rank after:", length(edge_rank)))
  
  res_ll <- list(bin_mat=bin_mat, edge_rank=edge_rank)
  
  print("Returning!")
  return(res_ll)
}


#' Create iGraph object from a symmetrical adjacency matrix, annotate it with centrality attributes and return the annotated iGraph.
#'
#' @importFrom igraph graph.adjacency vertex_attr betweenness closeness degree eigen_centrality transitivity edge_attr list.vertex.attributes
#'
#' @param adj_mat Binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @return Annotated igraph object representing the binary consensus matrix computed by get_ranked_consensus_binary_matrix() function.
#' @examples
#' \dontrun{
#' get_iGraph(adj_mat=binary.inference.matrix)
#' }
#' @keywords internal
#' @export
get_iGraph <- function(adj_mat){
  if(is.null(adj_mat)){
    stop("Input Adjacency Matrix is NULL!")
  }
  
  if(!isSymmetric(adj_mat)){
    print("Matrix is not symmetric. Getting symmetric matrix!")
    #adj_mat <- get_symmetric_matrix(adj_mat)
    adj_mat <- pmax(adj_mat, t(adj_mat))
  }
  
  iG <- igraph::graph.adjacency(adj_mat, mode="undirected", weighted=NULL)
  
  igraph::vertex_attr(iG, name="betweenness") <- as.vector(igraph::betweenness(iG))
  igraph::vertex_attr(iG, name="closeness") <- as.vector(igraph::closeness(iG, normalized=TRUE))
  igraph::vertex_attr(iG, name="degree") <- as.vector(igraph::degree(iG))
  igraph::vertex_attr(iG, name="eigenvector") <- as.vector(igraph::eigen_centrality(iG)$vector)
  igraph::vertex_attr(iG, name="cc") <- igraph::transitivity(iG, type="local", isolates="zero")
  
  igraph::vertex_attr(iG, name="color") <- "lightgray"
  igraph::vertex_attr(iG, name="highlightcolor") <- "darkgray"
  igraph::edge_attr(iG, name="color") <- "lightgray"
  igraph::edge_attr(iG, name="highlightcolor") <- "darkgray"
  
  print(igraph::list.vertex.attributes(iG))
  iG
}


#' Annotate iGraph object with centrality attributes and return the annotated iGraph.
#'
#' @importFrom igraph vertex_attr betweenness closeness degree eigen_centrality transitivity edge_attr list.vertex.attributes list.edge.attributes
#'
#' @param iG igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @return Annotated igraph object representing the binary consensus matrix computed by get_ranked_consensus_binary_matrix() function.
#' @examples
#' \dontrun{
#' annotate_iGraph(iG=inferred.iGraph)
#' }
#' @keywords internal
#' @export
annotate_iGraph <- function(iG){
  igraph::vertex_attr(iG, name="betweenness") <- as.vector(igraph::betweenness(iG))
  igraph::vertex_attr(iG, name="closeness") <- as.vector(igraph::closeness(iG, normalized=TRUE))
  igraph::vertex_attr(iG, name="degree") <- as.vector(igraph::degree(iG))
  igraph::vertex_attr(iG, name="eigenvector") <- as.vector(igraph::eigen_centrality(iG)$vector)
  igraph::vertex_attr(iG, name="cc") <- igraph::transitivity(iG, type="local", isolates="zero")
  
  vertex_attr_list <- igraph::list.vertex.attributes(iG)
  edge_attr_list <- igraph::list.edge.attributes(iG)
  
  if(!("color" %in% vertex_attr_list)){
    igraph::vertex_attr(iIG, name="color") <- "lightgray"
  }
  if(!("highlightcolor" %in% vertex_attr_list)){
    igraph::vertex_attr(iIG, name="highlightcolor") <- "darkgray"
  }
  
  if(!("color" %in% edge_attr_list)){
    igraph::edge_attr(iIG, name="color") <- "lightgray"
  }
  if(!("highlightcolor" %in% edge_attr_list)){
    igraph::edge_attr(iIG, name="highlightcolor") <- "darkgray"
  }
  
  #print(igraph::list.vertex.attributes(iG))
  return(iG)
}


#' Get top ranked candidates computed by Borda on the provided list of attributes and the cutoff.
#'
#' get_ranked_gene_list uses the annotation associated with each vertex to rank them and generated a list of vertices names ordered by rank. 
#' By default the annotations calculated by INfORM are "betweenness", "cc", "degree", "eccentricity", "closeness" & "eigenvector". The vertices are 
#' ranked by each annotation separately nad then the ranks are unified by means of Borda(). User must provide the annotations to use in ranking scheme, 
#' if the user has custom annotations such as "score" for differential gene expression then it can also be used for ranking.
#'
#' @importFrom igraph get.vertex.attribute
#' @importFrom TopKLists Borda
#' @importFrom utils head
#'
#' @param iGraph igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @param rank_list_attr Vector of network attributes/scores to use for generating a combined gene rank
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @return vector of genes ordered by rank, based on the selected attributes associated to each gene.
#' @examples
#' \dontrun{
#' get_ranked_gene_list(iGraph=inferred.igraph,
#' rank_list_attr=c("betweenness",
#' "cc",
#' "degree",
#' "eccentricity",
#' "closeness",
#' "eigenvector",
#' "score"),
#' debug_output=FALSE
#' )
#' }
#' @keywords internal
#' @export
get_ranked_gene_list <- function(iGraph, rank_list_attr, debug_output=FALSE){
  if(is.null(rank_list_attr)){
    stop("List of attributes for ranking is NULL!")
  }
  
  attrOrdMat <- list()
  for(a in 1:length(rank_list_attr)){
    val_ord=TRUE
    
    if(debug_output==TRUE)
      print(paste("Using attribute: ", rank_list_attr[a], sep=""))
    
    attrValueVector <- igraph::get.vertex.attribute(iGraph, rank_list_attr[a])
    
    if(rank_list_attr[a] == "score")
      attrValueVector <- abs(attrValueVector)
    
    attrOrdList <- cbind(igraph::V(iGraph)$name, attrValueVector)[order(attrValueVector, decreasing=val_ord)]
    
    if(debug_output==TRUE)
      print(utils::head(attrOrdList))
    
    attrOrdMat[[a]] <- attrOrdList
  }
  
  attrBorda <- TopKLists::Borda(attrOrdMat)
  
  if(debug_output==TRUE)
    print(utils::head(attrBorda$TopK, n=10))
  
  return(attrBorda$TopK$median)
}


#' Get modules from the main igraph.
#'
#' Extract modules from the main graph by using a specified community detection algorithm from igraph.
#'
#' @importFrom igraph cluster_walktrap cluster_spinglass cluster_louvain cluster_fast_greedy
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param method community detection method from igraph.
#' @return communities object cotaining communities identified from igraph object by using a specific method
#' @examples
#' \dontrun{
#' get_modules(iGraph=main.graph, method="walktrap")
#' }
#' @keywords internal
#' @export
get_modules <- function(iGraph, method="walktrap")
{
  switch(method,
         "walktrap" = {
           optimalStep <- c(2:10)[which.max(sapply(c(2:10), function(x){igraph::modularity(igraph::cluster_walktrap(iGraph, step=x))}))]
           igraph::cluster_walktrap(iGraph, step=optimalStep)
         },
         "spinglass" = igraph::cluster_spinglass(iGraph),
         "louvain" = igraph::cluster_louvain(iGraph),
         "greedy" = igraph::cluster_fast_greedy(iGraph)
  ) 
}
