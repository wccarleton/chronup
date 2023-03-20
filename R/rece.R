#' Build and optionally plot a radiocarbon-dated event count ensemble (RECE).
#'
#' @param c14_dates A matrix with two columns. The first column should contain
#'  mean uncalibrated radiocarbon dates and the second should contain the
#'  associated errors.
#' @param nsamples Integer number of sampled sequences to plot.
#' @param breaks Vector of temporal breakpoints to use for sampling. This is 
#'  the analytical/plot resolution, not the calibration resolution. It also
#'  defines the temporal boundaries of the RECE.
#' @param BP Logical. Use Before Present (BP) timescale?
#' @param cal_resolution Desired temporal resolution for calibration in years.
#' @param verbose Logical. The plot takes longer to create for more
#'  dates/samples. So, this logical parameter indicates whether the user wants
#'  to see messages indicating progress.
#' @param plot_it Logical. Should a RECE be plotted (default TRUE).
#' @param time_axis_index Integer value determining which elements of 
#'  breaks vector will be included in the x-axis ticks and labels for the 
#'  optional plot. The axis ticks will refer to time bin midpoints and be drawn 
#'  and labelled at seq(1, length(breaks) - 1, time_axis_label_index).
#' @param do_parallel Should the count sequence sampling be parallelized? The
#'  'parallel' package will need to be installed.
#' @param ncores Integer number of cores to use with parallelization. Default 
#'  is NULL and one fewer than the total available cores will be used.
#' @return Invisibly returns a list with two elements. The first is an
#'  radiocarbon-dated event count ensemble, which is a matrix containing
#'  probable count sequences where each column contains one probable sequence.
#'  The second element is matrix containing frequencies of counts in each time
#'  bin defined by breaks (used in the RECE heat map plot).
#' @export

rece <- function(c14_dates,
            nsamples,
            breaks,
            verbose = TRUE,
            BP = TRUE,
            cal_resolution = -1,
            plot_it = TRUE,
            time_axis_index = 10,
            do_parallel = FALSE,
            ncores = NULL){

    # Calibrate the dates and build a matrix to hold the densities
    calibrated_dates <- chronup::build_c14_matrix(c14_dates,
                                                    BP = BP,
                                                    resolution = cal_resolution)

    # extract the range of of calibrate dates (years bp) and create a sequence
    # of cal dates
    times <- seq(calibrated_dates$time_range[1],
                 calibrated_dates$time_range[2],
                 cal_resolution)

    # sample count sequences form the dates
    if(verbose){
        message("Sampling count sequences...")
    }

    if(do_parallel){
        if(is.null(ncores)){
            ncores <- parallel::detectCores() - 1
        }
        cl <- parallel::makeCluster(ncores)
        parallel::clusterEvalQ(cl,{
                            library(chronup)
                            })
        event_count_samples <- pblapply(cl = cl,
            X = 1:nsamples,
            FUN = chronup::sample_event_counts,
            chronun_matrix = calibrated_dates$chronun_matrix,
            times = times,
            breaks = breaks,
            BP = BP)
        parallel::stopCluster(cl)
    }else{
        event_count_samples <- pblapply(1:nsamples,
            chronup::sample_event_counts,
            chronun_matrix = calibrated_dates$chronun_matrix,
            times = times,
            breaks = breaks,
            BP = BP)
    }
    event_count_samples <- as.matrix(do.call(cbind, 
                                    event_count_samples))

    # build a RECE

    max_count <- max(event_count_samples)
    RECE <- t(apply(event_count_samples + 1,
                1,
                tabulate,
                nbins = max_count + 1))
    
    # turn 0's into NAs for plotting
    RECE[which(RECE == 0)] <- NA

    # dimensions (leaving out 0-counts)
    time_bins <- dim(RECE[, 2:max_count])[1]
    counts <- dim(RECE[, 2:max_count])[2]

    t_delta <- mean(diff(breaks))
    t_labels <- c(breaks - (t_delta / 2))[-1]

    # plotting
    if(plot_it){
        image(x = 1:time_bins,
            y = 1:counts,
            z = RECE[, 2:max_count],
            useRaster = T,
            col = grDevices::hcl.colors(n = 10,
                                        palette = "viridis",
                                        alpha = 0.9),
            axes = FALSE,
            xlab = "Time",
            ylab = "Count")
        axis_x_at <- seq(1, time_bins, time_axis_index)
        axis(1, at = axis_x_at, labels = t_labels[axis_x_at])
        axis(2, at = 1:counts)
    }

    # add row/col names to the RECE for convenience
    rownames(RECE) <- t_labels
    colnames(RECE) <- 0:(ncol(RECE) - 1)

    return(invisible(list(count_ensemble = event_count_samples,
                        RECE = RECE))
}