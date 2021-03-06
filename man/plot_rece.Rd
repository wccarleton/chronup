% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_rece.R
\name{plot_rece}
\alias{plot_rece}
\title{Plot radiocarbon-dated event count ensemble.}
\usage{
plot_rece(
  c14_dates,
  nsamples,
  verbose = TRUE,
  BP = TRUE,
  resolution = -1,
  use_ggplot2 = FALSE,
  axis_x_res = 100,
  axis_y_res = 1
)
}
\arguments{
\item{c14_dates}{A matrix with two columns. The first column should contain
mean uncalibrated radiocarbon dates and the second should contain the
associated errors.}

\item{nsamples}{Integer number of sampled sequences to plot.}

\item{verbose}{Logical. The plot takes longer to create for more
dates/samples. So, this logical parameter indicates whether the user wants
to see messages indicating progress.}

\item{BP}{Logical. Use Before Present (BP) timescale?}

\item{resolution}{Desired temporal resolution in years.}

\item{use_ggplot2}{Logical. Should the package ggplot2 be used for plotting?}

\item{axis_x_res}{The resolution for the plotted x axis (time). Not used if
use_ggplot2 == T.}

\item{axis_y_res}{The resolution for the plotted y axis (count). Not used if
use_ggplot2 == T.}
}
\value{
Invisibly returns a list with three elements. The first is an
radiocarbon-dated event count ensemble (rece), which is a matrix containing
probable count sequences where each column contains one probable sequence.
The second element is a vector of times at which counts have been sampled.
These time-stamps correspond to the rows of the rece. The last element is
the return value of the chronup::plot_count_ensemble function, which is NULL
if use_ggplot2 == FALSE or a ggplot2 object if use_ggplot2 == T.
}
\description{
Plot radiocarbon-dated event count ensemble.
}
