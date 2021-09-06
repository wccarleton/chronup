#' Plot radiocarbon-dated event count ensemble.
#'
#' @param c14_dates A matrix with two columns. The first column should contain
#'  mean uncalibrated radiocarbon dates and the second should contain the
#'  associated errors.
#' @param nsamples Integer number of sampled sequences to plot.
#' @param BP Logical. Use Before Present (BP) timescale?
#' @param resolution Desired temporal resolution in years.
#' @param use_ggplot2 Logical. Should the package ggplot2 be used for plotting?
#' @param verbose Logical. The plot takes longer to create for more
#'  dates/samples. So, this logical parameter indicates whether the user wants
#'  to see messages indicating progress.
#' @param axis_x_res The resolution for the plotted x axis (time). Not used if
#'  use_ggplot2 == T.
#' @param axis_y_res The resolution for the plotted y axis (count). Not used if
#'  use_ggplot2 == T.
#' @return Invisibly returns a list with three elements. The first is an
#'  radiocarbon-dated event count ensemble (rece), which is a matrix containing
#'  probable count sequences where each column contains one probable sequence.
#'  The second element is a vector of times at which counts have been sampled.
#'  These time-stamps correspond to the rows of the rece. The last element is
#'  the return value of the chronup::plot_count_ensemble function, which is NULL
#'  if use_ggplot2 == FALSE or a ggplot2 object if use_ggplot2 == T.

plot_rece <- function(c14_dates,
                    nsamples,
                    verbose = TRUE,
                    BP = TRUE,
                    resolution = -1,
                    use_ggplot2 = FALSE,
                    axis_x_res = 100,
                    axis_y_res = 1){
    calibrated_dates <- chronup::build_c14_matrix(c14_dates,
                                                    BP = BP,
                                                    resolution = resolution)
    times <- calibrated_dates$time_range[1]:calibrated_dates$time_range[2]
    if(verbose){
        message("Sampling count sequences...")
    }
    event_count_samples <- pbsapply(1:nsamples,
        chronup::sample_event_counts,
        chronun_matrix = calibrated_dates$chronun_matrix,
        times = times,
        BP = BP)
    if(verbose){
        message("Tabulating count frequencies...")
    }
    event_count_freqs <- t(apply(event_count_samples,
                                1,
                                chronup::tabulate_freqs,
                                n_events=dim(c14_dates)[1]))
    max_count <- chronup::find_max_count(event_count_freqs)
    p <- chronup::plot_count_ensemble(event_count_freqs[,2:max_count],
                        times = times,
                        use_ggplot2 = use_ggplot2,
                        axis_x_res = axis_x_res,
                        axis_y_res = axis_y_res)
    return(invisible(list(rece = event_count_freqs[,2:max_count],
                            times = times,
                            p = p)))
}
