#' Count events.
#'
#' This is a wrapper for graphics::hist.
#'
#' @param x A numeric vector of event times to be counted.
#' @param breaks A numeric vector of temporal bin edges passed at the `breaks`
#'  argument to graphics::hist. These edges need not define intervals at the
#'  same resolution as the `x` argument (e.g., event times can refer to years
#'  while the bins can refer to decades). But, keep in mind that the bins have
#'  to include all of the possible intervals into which events occurring at
#'  different times can fall. Also, the bin edges defined by this argument will
#'  serve as the right-most boundary condition which will be closed, i.e., the
#'  interval will be left-open and right-closed: ( ].
#' @param BP A logical argument specifying whether the events are on the BP
#'  (before present) timeline. If TRUE (default), the counts vector returned by
#'  graphics::hist is reversed to preserve a rightward/downward direction of
#'  time.
#' @return A vector containing event counts.
#' @export

countEvents <- function(x, breaks, BP = TRUE){
    counts <- hist(
                x,
                breaks = breaks,
                include.lowest = F,
                plot = F)$counts
    if(BP){
        return(rev(counts))
    }else{
        return(counts)
    }
}
