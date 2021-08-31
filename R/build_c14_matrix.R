build_c14_matrix <- function(
                            c14_dates,
                            BP = TRUE,
                            resolution = -1){
    nevents <- dim(c14_dates)[1]
    #clam, for c14 calibration
    clam_installed <- requireNamespace("clam", quietly = TRUE)
    if(!clam_installed){
        stop("package 'clam' is required for calibration")
    }
    message("Calibrating c14 dates...")
    c14post <- pblapply(
            1:nevents,
            function(x){
                sink(paste(tempdir(),"clam_output.txt",sep=""))
                caldate <- clam::calibrate(
                                    c14_dates[x,1],
                                    c14_dates[x,2],
                                    graph = F,
                                    BCAD = !BP)$calib
                sink()
                return(caldate)
            })
    sample_time_range <- range(
                            unlist(
                                lapply(
                                    c14post,
                                    function(x)range(x[,1])
                                )
                            )
                        )
    if(BP){
        ce_matrix <- approx_c14(
                c14post,
                sample_time_range[2],
                sample_time_range[1],
                resolution)
        time_range <- rev(sample_time_range)
    }else{
        ce_matrix <- approx_c14(
                c14post,
                sample_time_range[1],
                sample_time_range[2],
                resolution)
        time_range <- sample_time_range
    }
    return(list(ce_matrix = ce_matrix, time_range = time_range))
}
