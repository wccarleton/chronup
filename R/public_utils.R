#' Find max count.
#'
#' @param event_count_freqs A matrix where each column refers to the count of
#'  events and the rows refer to times (temporal bins used for counting). The
#'  cells each indicate the number of times a given count--time pair has been
#'  sampled.
#' @return Integer indicating the maximum count sampled from any temporal bin in
#'  the event_count_freqs matrix.
#' @export

find_max_count <- function(event_count_freqs){
    non_zeros <- apply(event_count_freqs,
                        2,
                        function(x)any(x>0))
    return(max(which(non_zeros)))
}
