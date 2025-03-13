#' Run ESEA Analysis
#'
#' This function performs ESEA (Edge Set Enrichment Analysis) on a set of edges, using their correlation scores and a provided pathway edge database. It computes enrichment scores, nominal p-values, FDR-adjusted p-values, and provides both summary results and detailed per-pathway outputs.
#'
#' @param EdgeCorScore A numeric vector of correlation scores or rank indexes (one score per edge) in decreasing order. 
#' @param pathwayEdge.db A list (or similar structure) containing pathway edge definitions. Each element should be a tab-delimited string where the first two fields are the pathway name and the pathway type (e.g., "kegg"), followed by edge identifiers.
#' @param weighted.score.type Numeric. The weighting scheme for calculating enrichment scores. Common values are \code{0} (unweighted), \code{1} (moderately weighted; default), or \code{2} (strongly weighted).
#' @param pathway Character. Specifies which pathway type to filter from \code{pathwayEdge.db}. Possible options include \code{"kegg"}, \code{"reactome"}, \code{"biocarta"}, \code{"nci"}, \code{"spike"}, \code{"huamncyc"}, or \code{"panther"}. Using \code{"null"} means no filtering by pathway type.
#' @param gs.size.threshold.min, gs.size.threshold.max Numeric. The minimum and maximum sizes of an edge set (gene set) to be considered. Defaults are \code{1} and \code{1000}, respectively.
#' @param reshuffling.type Character. The method to use for permutations:
#'   \itemize{
#'     \item \code{"edge.labels"} (default) - permutes the edge labels
#'     \item \code{"gene.labels"} - permutes gene-level scores
#'   }
#' @param nperm Integer. The number of permutations to run for p-value and FDR estimation. Default is \code{100}.
#' @param p.val.threshold Numeric. Nominal p-value threshold for reporting. Default is \code{-1} (disabled).
#' @param FDR.threshold Numeric. FDR q-value threshold for reporting. Default is \code{0.05}.
#' @param topgs Integer. How many top gene sets to report from each end of the enrichment distribution. Default is \code{1}.
#'
#' @details
#' The function calculates enrichment scores (ES) for each edge set in \code{pathwayEdge.db}, using the provided \code{EdgeCorScore} vector. It then carries out permutations to determine nominal p-values and FDR q-values. The output separates enrichment results into high-ranking (positive ES) and low-ranking (negative ES) sets.
#'
#' @return A list with two main components:
#' \item{SUMMARY.RESULTS}{A list of two data frames:
#'   \itemize{
#'     \item \code{SUMMARY.RESULTS.High_portion_of_the_rank}
#'     \item \code{SUMMARY.RESULTS.Low_portion_of_the_rank}
#'   }
#'   These summarize the enrichment results, including ES, normalized ES, p-values, FDR, and other metrics.
#' }
#' \item{Pathway results}{A list of data frames (one per significantly enriched pathway), each containing:
#'   \itemize{
#'     \item Edge IDs
#'     \item Edge correlation scores
#'     \item RES (running enrichment score)
#'     \item CORE_ENRICHMENT flags indicating whether each edge is in the "leading edge" subset
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Generate random scores for 100 edges
#' set.seed(123)
#' EdgeCorScore <- rnorm(100)
#'
#' # Create a mock pathway edge database
#' pathwayEdge.db <- list(
#'   "pathway1\tkegg\tedge1\tedge2",
#'   "pathway2\tkegg\tedge3\tedge4\tedge5"
#' )
#'
#' # Run the ESEA analysis
#' results <- compute_esea(
#'   EdgeCorScore = EdgeCorScore,
#'   pathwayEdge.db = pathwayEdge.db,
#'   weighted.score.type = 1,
#'   pathway = "kegg",
#'   gs.size.threshold.min = 1,
#'   gs.size.threshold.max = 1000,
#'   reshuffling.type = "edge.labels",
#'   nperm = 100,
#'   p.val.threshold = 0.05,
#'   FDR.threshold = 0.05,
#'   topgs = 1
#' )
#' # Inspect the results
#' results$SUMMARY.RESULTS$SUMMARY.RESULTS.High_portion_of_the_rank
#' results$SUMMARY.RESULTS$SUMMARY.RESULTS.Low_portion_of_the_rank
#' }
#'
#' @export
compute_esea <- function(
    EdgeCorScore,
    pathwayEdge.db,
    weighted.score.type = 1, 
    pathway = "kegg", # "kegg" or "reactome", "biocarta", "nci", "spike", "huamncyc", "panther"
    gs.size.threshold.min =1, 
    gs.size.threshold.max = 1000,
    reshuffling.type = "edge.labels", # or "gene.labels"
    nperm = 100, 
    p.val.threshold=-1,
    FDR.threshold = 0.05, 
    topgs =1
) 
{   
  print("Running ESEA Analysis...")
  #   input edge set database	
  edge.labels <- names(EdgeCorScore) #pvalue?
  temp<-pathwayEdge.db
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }
  
  max.size.G <- max(temp.size.G)      
  gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    edge.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    if(pathway == "null" ) {
      edge.set.name <- gs.line[1] 
      edge.set.desc <- gs.line[2] 
      edge.set.tags <- vector(length = edge.set.size, mode = "character")
      for (j in 1:edge.set.size) {
        edge.set.tags[j] <- gs.line[j + 2]
      } 
      existing.set <- is.element(edge.set.tags, edge.labels)
      set.size <- length(existing.set[existing.set == T])
      if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
      temp.size.G[gs.count] <- set.size
      gs[gs.count,] <- c(edge.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
      temp.names[gs.count] <- edge.set.name
      temp.desc[gs.count] <- edge.set.desc
      gs.count <- gs.count + 1
    } else  if (gs.line[2]==pathway)	{
      edge.set.name <- gs.line[1] 
      edge.set.desc <- gs.line[2] 
      edge.set.tags <- vector(length = edge.set.size, mode = "character")
      for (j in 1:edge.set.size) {
        edge.set.tags[j] <- gs.line[j + 2]
      } 
      existing.set <- is.element(edge.set.tags, edge.labels)
      set.size <- length(existing.set[existing.set == T])
      if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
      temp.size.G[gs.count] <- set.size
      gs[gs.count,] <- c(edge.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
      temp.names[gs.count] <- edge.set.name
      temp.desc[gs.count] <- edge.set.desc
      gs.count <- gs.count + 1
    }
  }
  
  Ng <- gs.count - 1
  print('intersection')
  print(Ng)
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric") 
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng] 
  size.G <- temp.size.G[1:Ng]
  
  rm(temp)
  gc()
  N<-length(EdgeCorScore)
  
  Obs.indicator <- matrix(nrow= Ng, ncol=N)
  Obs.RES <- matrix(nrow= Ng, ncol=N)
  
  Obs.ES <- vector(length = Ng, mode = "numeric")
  Obs.arg.ES <- vector(length = Ng, mode = "numeric")
  Obs.ES.norm <- vector(length = Ng, mode = "numeric")
  
  # ESEA methodology	
  ########################### start of auxiliary functions ###########################
  ESEA.EnrichmentScore <- function(edge.list, edge.set, weighted.score.type = 1, correl.vector = NULL) {  
    
    tag.indicator <- sign(match(edge.list, edge.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
    no.tag.indicator <- 1 - tag.indicator 
    N <- length(edge.list) 
    Nh <- length(edge.set) 
    Nm <-  N - Nh 
    if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
    norm.tag    <- 1.0/sum.correl.tag
    norm.no.tag <- 1.0/Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > - min.ES) {
      #      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      #      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
  }
  ###########################################################	
  ESEA.EnrichmentScore2 <- function(edge.list, edge.set, weighted.score.type = 1, correl.vector = NULL) {  
    
    N <- length(edge.list) 
    Nh <- length(edge.set) 
    Nm <-  N - Nh 
    
    loc.vector <- vector(length=N, mode="numeric")
    peak.res.vector <- vector(length=Nh, mode="numeric")
    valley.res.vector <- vector(length=Nh, mode="numeric")
    tag.correl.vector <- vector(length=Nh, mode="numeric")
    tag.diff.vector <- vector(length=Nh, mode="numeric")
    tag.loc.vector <- vector(length=Nh, mode="numeric")
    
    loc.vector[edge.list] <- seq(1, N)
    tag.loc.vector <- loc.vector[edge.set]
    
    tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
    
    if (weighted.score.type == 0) {
      tag.correl.vector <- rep(1, Nh)
    } else if (weighted.score.type == 1) {
      tag.correl.vector <- correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
    } else if (weighted.score.type == 2) {
      tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
    } else {
      tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
    }
    
    norm.tag <- 1.0/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1.0/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
    
    return(list(ES = ES))
  }
  ########################end of auxiliary functions###############################
  
  obs.s2n <- vector(length=N, mode="numeric")
  signal.strength <- vector(length=Ng, mode="numeric")
  tag.frac <- vector(length=Ng, mode="numeric")
  edge.frac <- vector(length=Ng, mode="numeric")
  
  obs.s2n<-as.vector(EdgeCorScore)
  obs.index <- order(obs.s2n, decreasing=T)            
  obs.s2n   <- sort(obs.s2n, decreasing=T) 
  
  obs.edge.labels <- edge.labels[obs.index]  	
  obs.edge.list2 <- obs.index
  for (i in 1:Ng) {
    
    edge.set <- gs[i,gs[i,] != "null"]
    edge.set2 <- vector(length=length(edge.set), mode = "numeric")
    edge.set2 <- match(edge.set, edge.labels)
    ESEA.results <- ESEA.EnrichmentScore(edge.list=obs.edge.list2, edge.set=edge.set2, weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
    Obs.ES[i] <- ESEA.results$ES
    Obs.arg.ES[i] <- ESEA.results$arg.ES
    Obs.RES[i,] <- ESEA.results$RES
    Obs.indicator[i,] <- ESEA.results$indicator
    if (Obs.ES[i] >= 0) {  # compute signal strength
      tag.frac[i] <- sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.G[i]
      edge.frac[i] <- Obs.arg.ES[i]/N
    } else {
      tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
      edge.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
    }
    signal.strength[i] <- tag.frac[i] * (1 - edge.frac[i]) * (N / (N - size.G[i]))
  }
  
  phi <- matrix(nrow = Ng, ncol = nperm)
  
  obs.phi <- matrix(nrow = Ng, ncol = nperm)
  
  for (i in 1:nperm) {
    obs.phi[ ,i] <- Obs.ES
  }
  
  if (reshuffling.type == "gene.labels") { # reshuffling phenotype labels
    for (i in 1:Ng) {
      edge.set <- gs[i,gs[i,] != "null"]
      edge.set2 <- vector(length=length(edge.set), mode = "numeric")
      edge.set2 <- match(edge.set, edge.labels)
      for (r in 1:nperm) {
        s2n<-rnorm(length(EdgeCorScore),mean=mean(EdgeCorScore),sd=sd(EdgeCorScore))
        
        index <- order(s2n, decreasing=T)            
        s2n   <- sort(s2n, decreasing=T) 
        edge.list2 <- index
        ESEA.results <- ESEA.EnrichmentScore2(edge.list=edge.list2, edge.set=edge.set2, weighted.score.type=weighted.score.type, correl.vector=s2n )   
        phi[i, r] <- ESEA.results$ES
        
      }
      gc()
    }		
  } else if (reshuffling.type == "edge.labels") { # reshuffling gene labels
    for (i in 1:Ng) {
      edge.set <- gs[i,gs[i,] != "null"]
      edge.set2 <- vector(length=length(edge.set), mode = "numeric")
      edge.set2 <- match(edge.set, edge.labels)
      for (r in 1:nperm) {
        edge.list2<-sample(1:N)
        ESEA.results <- ESEA.EnrichmentScore2(edge.list=edge.list2, edge.set=edge.set2, weighted.score.type=weighted.score.type, correl.vector=obs.s2n )   
        phi[i, r] <- ESEA.results$ES
        
      }
      gc()
    }
  }
  # Find nominal p-values       
  p.vals <- matrix(0, nrow = Ng, ncol = 2)
  for (i in 1:Ng) {
    if (Obs.ES[i] >= 0) {
      p.vals[i, 1] <-  sum(phi[i,] >= Obs.ES[i])/length(phi[i,])
      p.vals[i, 1] <-  signif(p.vals[i, 1], digits=5)
    } else {
      p.vals[i, 1] <-  sum(phi[i,] <= Obs.ES[i])/length(phi[i,])
      p.vals[i, 1] <-  signif(p.vals[i, 1], digits=5)
    }
  }
  # Compute FDRs 
  p.vals[, 2]<-p.adjust(p.vals[, 1],method="fdr")
  # Rescaling normalization for each edge set null
  for (i in 1:Ng) {
    pos.phi <- Obs.ES[i]
    neg.phi <- Obs.ES[i]
    for (j in 1:nperm) {
      if(is.na(phi[i, j])) next()
      if (phi[i, j] >= 0) {
        pos.phi <- c(pos.phi, phi[i, j]) 
      } else {
        neg.phi <- c(neg.phi, phi[i, j]) 
      }
    }
    
    pos.m <- mean(pos.phi)
    neg.m <- mean(abs(neg.phi))
    
    if (Obs.ES[i] >= 0) {
      Obs.ES.norm[i] <- Obs.ES[i]/pos.m
    } else {
      Obs.ES.norm[i] <- Obs.ES[i]/neg.m
    }
  }
  
  Obs.ES.index <- order(Obs.ES.norm, decreasing=T)
  
  # Produce results report
  p.val <- vector(length=Ng, mode="numeric")
  FDR.val <- vector(length=Ng, mode="numeric")
  Obs.ES <- signif(Obs.ES, digits=5)
  Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
  p.val <- signif(p.vals[,1], digits=4)
  FDR.val <- signif(p.vals[,2], digits=4)
  signal.strength <- signif(signal.strength, digits=3)
  tag.frac <- signif(tag.frac, digits=3)
  edge.frac <- signif(edge.frac, digits=3)
  
  report <- data.frame(cbind(gs.names, gs.desc, size.G, Obs.ES, Obs.ES.norm, p.val, FDR.val, tag.frac, edge.frac, signal.strength))
  names(report) <- c("GS", "SOURCE","SIZE" , "ES", "NES", "NOM p-val", "FDR q-val",  "Tag \\%", "Edge \\%", "Signal")
  report2 <- report
  report.index2 <- order(Obs.ES.norm, decreasing=T)
  for (i in 1:Ng) {
    report2[i,] <- report[report.index2[i],]
  }   
  report3 <- report
  report.index3 <- order(Obs.ES.norm, decreasing=F)
  for (i in 1:Ng) {
    report3[i,] <- report[report.index3[i],]
  }   
  phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
  phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
  report.phen1 <- report2[1:phen1.rows,]
  report.phen2 <- report3[1:phen2.rows,]
  result1<-list("SUMMARY.RESULTS.High_portion_of_the_rank"=report.phen1,"SUMMARY.RESULTS.Low_portion_of_the_rank"=report.phen2)
  # Produce report for each pathway passing the nominal p-value,FDR q-value  or the top topgs in each side
  n<-1
  result2<-list()
  
  for (i in 1:Ng) {
    if(is.na(p.vals[i, 1]) || is.na(p.vals[i, 2])) next
    if ((p.vals[i, 1] <= p.val.threshold) ||
        (p.vals[i, 2] <= FDR.threshold) ||
        (is.element(i, c(Obs.ES.index[1:topgs], Obs.ES.index[(Ng - topgs + 1): Ng])))) {
      kk <- 1
      edge.number <- vector(length = size.G[i], mode = "character")
      edge.names <- vector(length = size.G[i], mode = "character")
      
      edge.list.loc <- vector(length = size.G[i], mode = "numeric")
      core.enrichment <- vector(length = size.G[i], mode = "character")
      edge.s2n <- vector(length = size.G[i], mode = "numeric")
      edge.RES <- vector(length = size.G[i], mode = "numeric")
      rank.list <- seq(1, N)
      
      if (Obs.ES[i] >= 0) {
        set.k <- seq(1, N, 1)
        phen.tag <- "High_part_of_the_rank"
        loc <- match(i, Obs.ES.index)
      } else {
        set.k <- seq(N, 1, -1)
        phen.tag <- "Low_part_of_the_rank"
        loc <- Ng - match(i, Obs.ES.index) + 1
      }
      
      for (k in set.k) {
        if (Obs.indicator[i, k] == 1) {
          edge.number[kk] <- kk
          edge.names[kk] <- obs.edge.labels[k]
          
          edge.list.loc[kk] <- k
          edge.s2n[kk] <- signif(obs.s2n[k], digits=3)
          edge.RES[kk] <- signif(Obs.RES[i, k], digits = 3)
          if (Obs.ES[i] >= 0) {
            core.enrichment[kk] <- ifelse(edge.list.loc[kk] <= Obs.arg.ES[i], "YES", "NO")
          } else {
            core.enrichment[kk] <- ifelse(edge.list.loc[kk] > Obs.arg.ES[i], "YES", "NO")
          }
          kk <- kk + 1
        }
      }
      
      edge.report <- data.frame(cbind(edge.number, edge.names,  edge.list.loc, edge.s2n, edge.RES, core.enrichment))
      names(edge.report) <- c("#", "EdgeID",  "List Loc", "EdgeCorScore", "RES", "CORE_ENRICHMENT")
      filename <- paste( gs.names[i], ".", phen.tag,  sep="", collapse="")
      result2[[n]]<-edge.report
      names(result2)[n]<-filename
      n<-n+1
    }
  }
  
  result<-list("SUMMARY.RESULTS"=result1,"Pathway results"=result2)
  
  return(result)
}
