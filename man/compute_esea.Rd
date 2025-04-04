% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esea.R
\name{compute_esea}
\alias{compute_esea}
\title{Run ESEA Analysis}
\usage{
compute_esea(
  EdgeCorScore,
  pathwayEdge.db,
  weighted.score.type = 1,
  pathway = "kegg",
  gs.size.threshold.min = 1,
  gs.size.threshold.max = 1000,
  reshuffling.type = "edge.labels",
  nperm = 100,
  p.val.threshold = -1,
  FDR.threshold = 0.05,
  topgs = 1
)
}
\arguments{
\item{EdgeCorScore}{A numeric vector of correlation scores or rank indexes (one score per edge) in decreasing order.}

\item{pathwayEdge.db}{A list (or similar structure) containing pathway edge definitions. Each element should be a tab-delimited string where the first two fields are the pathway name and the pathway type (e.g., "kegg"), followed by edge identifiers.}

\item{weighted.score.type}{Numeric. The weighting scheme for calculating enrichment scores. Common values are \code{0} (unweighted), \code{1} (moderately weighted; default), or \code{2} (strongly weighted).}

\item{pathway}{Character. Specifies which pathway type to filter from \code{pathwayEdge.db}. Possible options include \code{"kegg"}, \code{"reactome"}, \code{"biocarta"}, \code{"nci"}, \code{"spike"}, \code{"huamncyc"}, or \code{"panther"}. Using \code{"null"} means no filtering by pathway type.}

\item{gs.size.threshold.min, }{gs.size.threshold.max Numeric. The minimum and maximum sizes of an edge set (gene set) to be considered. Defaults are \code{1} and \code{1000}, respectively.}

\item{reshuffling.type}{Character. The method to use for permutations:
\itemize{
\item \code{"edge.labels"} (default) - permutes the edge labels
\item \code{"gene.labels"} - permutes gene-level scores
}}

\item{nperm}{Integer. The number of permutations to run for p-value and FDR estimation. Default is \code{100}.}

\item{p.val.threshold}{Numeric. Nominal p-value threshold for reporting. Default is \code{-1} (disabled).}

\item{FDR.threshold}{Numeric. FDR q-value threshold for reporting. Default is \code{0.05}.}

\item{topgs}{Integer. How many top gene sets to report from each end of the enrichment distribution. Default is \code{1}.}
}
\value{
A list with two main components:
\item{SUMMARY.RESULTS}{A list of two data frames:
\itemize{
\item \code{SUMMARY.RESULTS.High_portion_of_the_rank}
\item \code{SUMMARY.RESULTS.Low_portion_of_the_rank}
}
These summarize the enrichment results, including ES, normalized ES, p-values, FDR, and other metrics.
}
\item{Pathway results}{A list of data frames (one per significantly enriched pathway), each containing:
\itemize{
\item Edge IDs
\item Edge correlation scores
\item RES (running enrichment score)
\item CORE_ENRICHMENT flags indicating whether each edge is in the "leading edge" subset
}
}
}
\description{
This function performs ESEA (Edge Set Enrichment Analysis) on a set of edges, using their correlation scores and a provided pathway edge database. It computes enrichment scores, nominal p-values, FDR-adjusted p-values, and provides both summary results and detailed per-pathway outputs.
}
\details{
The function calculates enrichment scores (ES) for each edge set in \code{pathwayEdge.db}, using the provided \code{EdgeCorScore} vector. It then carries out permutations to determine nominal p-values and FDR q-values. The output separates enrichment results into high-ranking (positive ES) and low-ranking (negative ES) sets.
}
\examples{
\dontrun{
# Example usage:
# Generate random scores for 100 edges
set.seed(123)
EdgeCorScore <- rnorm(100)

# Create a mock pathway edge database
pathwayEdge.db <- list(
  "pathway1\tkegg\tedge1\tedge2",
  "pathway2\tkegg\tedge3\tedge4\tedge5"
)

# Run the ESEA analysis
results <- compute_esea(
  EdgeCorScore = EdgeCorScore,
  pathwayEdge.db = pathwayEdge.db,
  weighted.score.type = 1,
  pathway = "kegg",
  gs.size.threshold.min = 1,
  gs.size.threshold.max = 1000,
  reshuffling.type = "edge.labels",
  nperm = 100,
  p.val.threshold = 0.05,
  FDR.threshold = 0.05,
  topgs = 1
)
# Inspect the results
results$SUMMARY.RESULTS$SUMMARY.RESULTS.High_portion_of_the_rank
results$SUMMARY.RESULTS$SUMMARY.RESULTS.Low_portion_of_the_rank
}

}
