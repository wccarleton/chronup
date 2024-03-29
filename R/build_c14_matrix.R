#' Build C14 Matrix.
#'
#' @param c14_dates A matrix with two columns. The first column should contain
#'  uncalibrated radiocarbon date means and the second should contain the
#'  associated errors.
#' @param BP Logical. Use Before Present (BP) timescale?
#' @param resolution Desired temporal resolution in years.
#' @return A list with two elements. The first element is a chronological
#'  uncertainty (chronun) matrix containing a number of columns equal to the
#'  number of dates (rows) in the c14_matrix. Each column contains a single
#'  calibrated radiocarbon date density. The rows contain samples of the
#'  corresponding calibrated date density. The second element in the returned
#'  list contains the time range spanned by all calibrated dates. This range
#'  and the resolution argument can be used to recover time-stamps for each
#'  column in the chronun matrix.
#' @export

build_c14_matrix <- function(c14_dates,
                            BP = TRUE,
                            resolution = -1){
    nevents <- dim(c14_dates)[1]
    IntCal_installed <- requireNamespace("IntCal", quietly = TRUE)
    if(!IntCal_installed){
        stop("package IntCal is required for calibration")
    }
    message("Calibrating c14 dates...")
    c14post <- pblapply(1:nevents,
                        function(x, dates){
                            IntCal::caldist(age = dates[x, 1],
                                            error = dates[x, 2],
                                            BCAD = !BP)},
                            dates = c14_dates)
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
        time_range <- rev(sample_time_range)
    }else{
        chronun_matrix <- approx_c14(c14post,
                                sample_time_range[1],
                                sample_time_range[2],
                                resolution)
        time_range <- sample_time_range
    }
    return(list(chronun_matrix = chronun_matrix, time_range = time_range))
}
