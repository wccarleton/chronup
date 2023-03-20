#' Count events.
#'
#' Counts the number of events from 'x' that fall into time-bins defined by
#' breaks.
#'
#' @param x A numeric vector of event times to be counted.
#' @param breaks A numeric vector of temporal bin edges The bin edges defined 
#'  by this argument will serve as the right-most boundary condition which will 
#'  be closed, i.e., the interval will be left-open and right-closed: ( ].
#' @param BP A logical argument specifying whether the events are on the BP
#'  (before present) timeline.
#' @return A vector containing event counts.
#' @export

count_events <- function(x,
                        breaks,
                        BP){
    n <- length(breaks) - 1
    datum <- breaks[1]
    if(BP){
        datum <- -datum
        x <- -x
    }
    bins <- floor((x - datum) / abs(mean(diff(breaks)))) + 1
    bins_trim <- bins[which(bins >= 1 & bins <= n)]
    return(tabulate(bins_trim))
}
