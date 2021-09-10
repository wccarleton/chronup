#' Simulate event counts.
#'
#' @param process A vector describing the background event-count process (i.e.,
#'  conditional mean of count-based random variable). This will be passed to
#'  R's 'sample()' function as the 'prob' parameter.
#' @param times A vector of event times to sample from.
#' @param nevents Integer number of events with uncertain event-times.
#' @param nsamples Integer number of sample sequences to draw.
#' @param binning_resolution Resolution of the desired time intervals (time
#'  bins). It should be negative if using the BP timescale (i.e, if BP = T).
#' @param BP Logical (default T). Assume a Before Present timescale?
#' @param c14 Logical (default T). Are the events dated with radiocarbon? If
#'  so, the R package 'IntCal' will be used to simulate c14 dates from the
#'  'true' samples of the 'times' vector and then calibrate those dates.
#' @param chronun_matrix A matrix containing discrete estimates describing
#'  chronological uncertainty. Each column should contain density estimates
#'  for a single event and the rows should each refer to discrete times. If c14
#'  = F, then this matrix must be included.
#' @param new_times A vector of times corresponding to the rows in the chronun
#'  matrix.
#' @param bigmatrix A character vector containing a path pointing to a
#'  'bigmemory' matrix descriptor file, or NULL (default).
#' @param parallel Logical (default = T). Use parallel (multiple processors) to
#'  speed up computation?
#' @return A list with 1) a matrix containing probable event count sequences
#'  or a file path pointing to a file-backed bigmatrix (an event count
#'  ensemble); 2) a dataframe containing the true (chronological-error-free)
#'  event-count sample; 3) a vector of time-stamps corresponding to the rows in
#'  the event count ensemble matrix; and 4) a vector of simulated uncalibrated
#'  radiocarbon dates (NULL is c14 == FALSE).
#' @import pbapply graphics
#' @export

simulate_event_counts <- function(process,
                                times,
                                nevents,
                                nsamples,
                                binning_resolution = -1,
                                BP = T,
                                c14 = T,
                                chronun_matrix = NULL,
                                new_times = NULL,
                                bigmatrix = T,
                                parallel = T){

    #check user options for required packages
    #big memory
    if(bigmatrix){
        bigmemory_installed <- requireNamespace("bigmemory", quietly = TRUE)
        if(!bigmemory_installed){
            stop("package 'bigmemory' is required if bigmatrix = T")
        }
    }

    #parallel
    if(parallel){
        parallel_installed <- requireNamespace("parallel", quietly = TRUE)
        if(!parallel_installed){
            stop("package 'parallel' is required if parallel = T")
        }
    }

    #IntCal, for c14 calibration
    if(c14){
        IntCal_installed <- requireNamespace("IntCal", quietly = TRUE)
        if(!IntCal_installed){
            stop("package 'IntCal' is required if c14 = T")
        }
    }

    #check timescale and issue warning
    if(BP & (binning_resolution > 0) ){
        warning("Binning resolution should not be positive if BP = T because time indeces should count down toward the present. You may get unexpected results.")
    }

    #check for chronological error matrix
    if(!c14 & is.null(chronun_matrix)){
        stop("chronun_matrix is required if c14 = FALSE")
    }

    #check for new_times
    if(!is.null(chronun_matrix) & is.null(new_times)){
        stop("If chronun_matrix is supplied, new_times must also be supplied.")
    }

    #set some global parameters
    wd <- paste(getwd(),"/",sep="")

    resolution <- mean(diff(times))

    start <- times[1]

    end <- times[length(times)]

    event_times <- sample(x = times,
                        size = nevents,
                        replace = T,
                        prob = process)

    if(BP){
        breaks <- seq(start, end + binning_resolution, binning_resolution)
    }else{
        breaks <- seq(start, end + binning_resolution, binning_resolution)
    }

    nbins <- length(breaks) - 1

    true_event_counts <- data.frame(Timestamps = breaks[1:nbins],
                                Count = chronup::count_events(x = event_times,
                                                            breaks = breaks,
                                                            BP = BP),
                                Time = c(1:nbins) * abs(binning_resolution))

    if(c14){
        message("Simulating and calibrating c14 dates.")
        simc14 <- t(sapply(event_times, IntCal::calBP.14C))
        c14post <- pblapply(1:nevents,
                            function(x, dates){
                                IntCal::caldist(age = dates[x, 1],
                                                error = dates[x, 2],
                                                BCAD = !BP)},
                                dates = simc14)
        sample_time_range <- range(
                                unlist(
                                    lapply(c14post,
                                        function(x)range(x[, 1])
                                    )
                                )
                            )
        if(BP){
            chronun_matrix <- approx_c14(c14post,
                                        sample_time_range[2],
                                        sample_time_range[1],
                                        resolution)
            new_times <- seq(sample_time_range[2],
                            sample_time_range[1],
                            resolution)

            new_span <- length(new_times)

            new_breaks <- seq(sample_time_range[2],
                            sample_time_range[1] - abs(binning_resolution),
                            binning_resolution)
        }else{
            chronun_matrix <- approx_c14(c14post,
                                        sample_time_range[1],
                                        sample_time_range[2],
                                        resolution)
            new_times <- seq(sample_time_range[1],
                            sample_time_range[2],
                            resolution)

            new_span <- length(new_times)

            new_breaks <- seq(sample_time_range[1],
                            sample_time_range[2] + abs(binning_resolution),
                            binning_resolution)
        }
    }else{
        simc14 <- NULL
        new_span <- length(new_times)
        new_time_range <- range(new_times)
        if(BP){
            new_breaks <- seq(sample_time_range[2],
                            sample_time_range[1] - abs(binning_resolution),
                            binning_resolution)
        }else{
            new_breaks <- seq(sample_time_range[1],
                            sample_time_range[2] + abs(binning_resolution),
                            binning_resolution)
        }
    }

    if(bigmatrix){
        count_ensemble <- bigmemory::filebacked.big.matrix(
                            nrow = new_span,
                            ncol = nsamples,
                            backingpath = wd,
                            backingfile = "count_ensemble_mat",
                            descriptorfile = "count_ensemble_desc")
    }

    message("Simulating event count sequences.")

    if(parallel & bigmatrix){
        ncores <- parallel::detectCores()
        cl <- parallel::makeCluster(ncores - 1)
        parallel::clusterEvalQ(cl,{
                            wd <- getwd()
                            library(chronup)
                            })
        pbapply::pbsapply(
                    cl = cl,
                    X = 1:nsamples,
                    FUN = chronup::sample_event_counts,
                    chronun_matrix = chronun_matrix,
                    times = new_times,
                    breaks = new_breaks,
                    bigmatrix = paste(wd,"count_ensemble_desc",sep=""))
        parallel::stopCluster(cl)
        return(list(
                count_ensemble =  paste(wd,"count_ensemble_desc",sep=""),
                new_times = new_times,
                counts = true_event_counts,
                simc14 = simc14))
    }else if(parallel & !bigmatrix){
        ncores <- parallel::detectCores()
        cl <- parallel::makeCluster(ncores - 1)
        parallel::clusterEvalQ(cl,{
                            wd <- getwd()
                            library(chronup)
                            })
        count_ensemble <- pbapply::pbsapply(
                            cl = cl,
                            X = 1:nsamples,
                            FUN = chronup::sample_event_counts,
                            chronun_matrix = chronun_matrix,
                            times = new_times,
                            breaks = new_breaks,
                            bigmatrix = NULL)
        parallel::stopCluster(cl)
        return(list(count_ensemble = count_ensemble,
                    new_times = new_times,
                    counts = true_event_counts,
                    simc14 = simc14))
    }else if(!parallel & bigmatrix){
        pbapply::pbsapply(
                    X = 1:nsamples,
                    FUN = chronup::sample_event_counts,
                    chronun_matrix = chronun_matrix,
                    times = new_times,
                    breaks = new_breaks,
                    bigmatrix = paste(wd,"count_ensemble_desc",sep=""))
        return(list(count_ensemble =  paste(wd,"count_ensemble_desc",sep=""),
                    new_times = new_times,
                    counts = true_event_counts,
                    simc14 = simc14))
    }else if(!parallel & !bigmatrix){
        count_ensemble <- pbapply::pbsapply(
                            X = 1:nsamples,
                            FUN = chronup::sample_event_counts,
                            chronun_matrix = chronun_matrix,
                            times = new_times,
                            breaks = new_breaks,
                            bigmatrix = NULL)
        return(list(count_ensemble = count_ensemble,
                    new_times = new_times,
                    counts = true_event_counts,
                    simc14 = simc14))
    }
}
