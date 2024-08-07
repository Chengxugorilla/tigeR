% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_Model.R
\name{build_SURV_Model}
\alias{build_SURV_Model}
\title{Build Lasso-cox prediction model for Survival}
\usage{
build_SURV_Model(SE, Signature, rmBE = FALSE, PT_drop, lambda = 1, ...)
}
\arguments{
\item{SE}{an SummarizedExperiment(SE) object or a list consists of SE objects. The colData of SE objects must contain response information.}

\item{Signature}{an gene set you interested in}

\item{rmBE}{whether remove batch effect between different data set using internal Combat method}

\item{PT_drop}{If TRUE, only Untreated patient will be use for model training.}

\item{lambda}{the lambda value for Lasso.}

\item{...}{the arguments}
}
\description{
Generate a Lasso-cox model.
}
