% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_event_times.R
\name{sample_event_times}
\alias{sample_event_times}
\title{Sample event times.}
\usage{
sample_event_times(x = NULL, chronun_matrix, times, bigmatrix = NULL)
}
\arguments{
\item{x}{An integer specifying the index of a vector to enable easy
vectorization and/or parallelization, or NULL (default) if not using a
vectorized or parallel function call (e.g., apply() or parApply()). If
you're not using a parallel apply function (e.g., parApply), ignore this
parameter.}

\item{chronun_matrix}{A matrix containing discrete estimates describing
chronological uncertainty. Each column should contain density estimates
for a single event and the rows should each refer to discrete times.}

\item{times}{A vector of possible event times.}

\item{bigmatrix}{A character vector containing a path pointing to a
'bigmemory' matrix descriptor file, or NULL (default).}
}
\value{
A vector containing a sample of probable event times, if not calling
this function with a parallel function (e.g., parApply). If using something
like \code{apply()} or 'parApply()', this function will return a matrix. If the
\code{bigmemory} argument is not NULL, this function will return nothing and
instead write the output to the relevant big matrix file.
}
\description{
Sample event times.
}
