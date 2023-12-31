% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBinByGroupMatrix.R
\name{getBinByGroupMatrix}
\alias{getBinByGroupMatrix}
\title{Get the BinByGroup matrix}
\usage{
getBinByGroupMatrix(
  binByCell,
  cell_anno,
  group_col = "group",
  n_sample = 1000,
  seed = 1000
)
}
\arguments{
\item{binByCell}{BinByCellMatrix}

\item{cell_anno}{Cell meta with cell group information}

\item{group_col}{The column name of the cell group in cell_anno}

\item{n_sample}{Number of cells sampled in each group}

\item{seed}{A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
reproduce results downstream.}
}
\value{
binByGroup_norm
}
\description{
Get the BinByGroup matrix
}
