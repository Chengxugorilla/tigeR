% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_biomk.R
\name{compare_biomk}
\alias{compare_biomk}
\title{Compare the AUC of signature between datasets}
\usage{
compare_biomk(
  SE = NULL,
  Signature = NULL,
  method = "Average_mean",
  PT_drop = TRUE,
  show.val = TRUE,
  val.size = 2
)
}
\arguments{
\item{SE}{a SummarizedExperiment object for which you want to calculate the Signature Score.}

\item{Signature}{a gene vector represents user-defined signature for Immunotherapy response. If NULL, the function will only calculate 23 built-in signatures in tigeR.}

\item{method}{the method for calculating gene set scores which has several options: "Average_mean", "Weighted_mean", or "GSVA". The method can be set to NULL if the length of the parameter geneSet is 1. This means that if you are working with only one gene, the specific calculation method may not be applicable or necessary.}

\item{PT_drop}{if TRUE, only untreated patient will be used for model training.}

\item{show.val}{if TRUE, the value will be show edin the heatplot.}

\item{val.size}{an integer represents the font size of AUC value.}
}
\value{
\describe{generate a heatmap of signature AUC of datasets to compare the performance of
  biomarkers of interest with built-in biomarkers.}
}
\description{
Generate a heatmap of signature auc of datasets.
}
\examples{
compare_biomk(
   SE = list(MEL_GSE78220, MEL_GSE93157),
   Signature = ipt[1:72],
   method = "Weighted_mean",
   PT_drop = FALSE)
}
