resample <- function(y,
                    time_stamps,
                    breaks,
                    FUN){
    resampled <- c()
    for(j in 2:length(breaks)){
        selection <- time_stamps >= breaks[j - 1] & time_stamps < breaks[j]
        resampled <- c(resampled, FUN(y[selection]))
    }
    return(resampled)
}
