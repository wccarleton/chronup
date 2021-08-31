#' Sample event counts.
#'
#' @param x An integer specifying the index of a vector to enable easy
#'  vectorization and/or parallelization, or NULL (default) if not using a
#'  vectorized or parallel function call (e.g., apply() or parApply()). If
#'  you're not using a parallel apply function (e.g., parApply), ignore this
#'  parameter.
#' @param ce_matrix A matrix containing discrete estimates describing
#'  chronological uncertainty. Each column should contain density estimates
#'  for a single event and the rows should each refer to discrete times.
#' @param times A vector of possible event times.
#' @param breaks A vector containing (time) bin edges used for counting events.
#'  These edges need not define intervals at the same resolution as the `times`
#'  argument (e.g., times can refer to years while the bins can refer to
#'  decades). But, keep in mind that the bins should include all of the
#'  possible intervals into which events can fall. Also, the bin edges defined
#'  by this argument will serve as the right-most boundary condition which will
#'  be closed, i.e., the interval will be left-open and right-closed: ( ].
#' @param BP Logical (default T). Assume a Before Present timescale?
#' @param bigmatrix A character vector containing a path pointing to a
#'  'bigmemory' matrix descriptor file, or NULL (default).
#' @return A vector containing a sample of probable event counts, if not calling
#'  this function with a parallel function (e.g., parApply). If using something
#'  like `apply()` or 'parApply()', this function will return a matrix. If the
#'  `bigmemory` argument is not NULL, this function will return nothing and
#'  instead write the output to the relevant big matrix file.
#' @export

sample_event_counts <- function(
                        x = NULL,
                        ce_matrix,
                        times,
                        breaks = NULL,
                        BP = T,
                        bigmatrix = NULL){
    # Check user input
    if(dim(ce_matrix)[1] != length(times)){
        stop("times length must be equal to the number of rows in ce_matrix.")
    }
    resolution <- mean(diff(times))
    if(BP & resolution > 0){
        stop("When BP is TRUE, times should be a decreasing vector.")
    }
    if(is.null(breaks)){
        if(BP){
            start <- times[1]
            end <- times[length(times)] - 1
            breaks <- seq(from = start, to = end, by = resolution)
        }else{
            start <- times[1]
            end <- times[length(times)] + 1
            breaks <- seq(from = start, to = end, by = resolution)
        }
    }
    times_sample <- apply(ce_matrix,
                        2,
                        function(j)sample(times, size=1, prob=j))
    count_sample <- chronup::count_events(times_sample, breaks, BP)
    if(!is.null(bigmatrix)){
        if(requireNamespace("bigmemory", quietly = TRUE)){
            m <- bigmemory::attach.big.matrix(bigmatrix)
            m[,x] <- count_sample
            return()
        }else{
            stop("The package 'bigmemory' is required when the 'bigmatrix' argument is not NULL. Install the package to make use of this option")
        }
    }else{
        return(count_sample)
    }
}
