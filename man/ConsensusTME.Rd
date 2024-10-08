% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconv_TME.R
\name{ConsensusTME}
\alias{ConsensusTME}
\title{ConsensusTME deconvolution}
\usage{
ConsensusTME(SE, cancer = "SKCM", ...)
}
\arguments{
\item{SE}{a SummarizedExperiment object contains the bulk RNA-seq dataset that you want to use for deconvolution and obtaining its cell fraction.}

\item{cancer}{string passed to indicate which TCGA cancer type samples are most similar to. N.B samples of different cancer types should be run seperately. Available cancer types: "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP","LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM".}

\item{...}{the arguments}
}
\description{
use ConsensusTME to predict TME
}
