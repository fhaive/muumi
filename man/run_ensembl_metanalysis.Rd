% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_analysis.R
\name{run_ensembl_metanalysis}
\alias{run_ensembl_metanalysis}
\title{Computes a gene rank based on an ensembl of metanalysis methods, including effect size, p-value and rank product.}
\usage{
run_ensembl_metanalysis(
  meta_dataframe,
  method = c("effect_size", "pvalue", "rank_product"),
  class,
  origin,
  metric = "median"
)
}
\arguments{
\item{meta_dataframe}{A dataframe with genes on the rows (as rownames) and samples in the columns. The columns should contain adjusted p-values deriving from gene-based statistical testing.}

\item{method}{Statistical method(s) to be included in the ensembl metanalysis.}

\item{class}{a vector containing the class labels of the samples. In the two class unpaired case, the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). For one class data, the label for each sample should be 1.}

\item{origin}{a vector containing the origin labels of the samples.}

\item{metric}{One statistical metric between "median" and "mean".}
}
\value{
A gene list ranked on the base of the methods chosen for the metanalysis.
}
\description{
Computes a gene rank based on an ensembl of metanalysis methods, including effect size, p-value and rank product.
}
\examples{
\dontrun {
run_ensembl_metanalysis(meta_dataframe)
}
}
