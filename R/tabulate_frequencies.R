#' Tabulate frequency distribution of event counts for a given interval.
#'
#' @param event_count_sample An integer vector containing samples of event
#'  count sequences for a given interval. Each element should refer to the
#'  number of times events fell into the relevant interval.
#' @param n_events Number of events.
#' @return A vector containing frequencies of event counts for each interval.
#'  Each element will refer to the number of times a given count was
#'  sampled for a given time.
#' @export

tabulate_freqs <- function(
                        event_count_sample,
                        n_events){
    breaks <- -1:n_events
    event_count_freqs <- graphics::hist(event_count_sample,
                                            breaks = breaks,
                                            plot = F)$counts
    return(event_count_freqs)
}
