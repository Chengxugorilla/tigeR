% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_name_standardization.R
\name{cell_name_standardization}
\alias{cell_name_standardization}
\title{Standardize the cell names in TME deconvolution}
\usage{
cell_name_standardization(SE)
}
\arguments{
\item{SE}{a SummarizedExperiment object containing cell fraction data generated by deconv_TME().}
}
\value{
\describe{return a SummarizedExperiment object which rownames have been standardized.}
}
\description{
Standarize the cell names in TME deconvolution.
}
