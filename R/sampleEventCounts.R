#' Sample event counts.
#'
#' @param x An integer specifying the index of a vector to enable easy
#'  vectorization and/or parallelization, or NULL (default) if not using a
#'  vectorized or parallel function call (e.g., apply() or parApply()). If
#'  you're not using a parallel apply function (e.g., parApply), ignore this
#'  parameter.
#' @param ceMatrix A matrix containing discrete estimates describing
#'  chronological uncertainty. Each column should contain density estimates
#'  for a single event and the rows should each refer to discrete times.
#' @param times A vector of possible event times.
#' @param breaks A vector containing (time) bin edges used for counting events.
#'  These edges need not define intervals at the same resolution as the `times`
#'  argument (e.g., times can refer to years while the bins can refer to
#'  decades). But, keep in mind that the bins have to include all of the
#'  possible intervals into which events occurring at different times can fall.
#'  Also, the bin edges defined by this argument will serve as the right-most
#'  boundary condition which will be closed, i.e., the interval will be left-
#'  open and right-closed: ( ].
#' @param BP Logical (default T). Assume a Before Present timescale?
#' @param bigmatrix A character vector containing a path pointing to a
#'  'bigmemory' matrix descriptor file, or NULL (default).
#' @return A vector containing a sample of probable event counts, if not calling
#'  this function with a parallel function (e.g., parApply). If using something
#'  like `apply()` or 'parApply()', this function will return a matrix. If the
#'  `bigmemory` argument is not NULL, this function will return nothing and
#'  instead write the output to the relevant big matrix file.
#' @export

sampleEventCounts <- function(
                        x = NULL,
                        ceMatrix,
                        times,
                        breaks,
                        BP = T,
                        bigmatrix = NULL){
    times_sample <- apply(ceMatrix, 2, function(j)sample(times, size=1, prob=j))
    count_sample <- countEvents(times_sample, breaks, BP)
    if(!is.null(bigmatrix)){
        if(requireNamespace("bigmemory", quietly = TRUE)){
            m <- attach.big.matrix(bigmatrix)
            m[,x] <- count_sample
            return()
        }else{
            stop("The package 'bigmemory' is required when the 'bigmatrix' argument is not NULL. Install and load the package to make use of this option")
        }
    }else{
        return(count_sample)
    }
}
