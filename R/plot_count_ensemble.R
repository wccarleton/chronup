#' Plot count ensemble.
#'
#' @param count_ensemble A matrix with comlumns that each refer to a probable
#'  count sequence and rows refer to the count of events at a given time.
#' @param times A vector of times stamps that correspond to the times of the
#'  count_ensemble rows.
#' @param use_ggplot2 Logical. Should the package ggplot2 be used for plotting?
#' @param axis_x_res The resolution for the plotted x axis (time). Not used if
#'  use_ggplot2 == T.
#' @param axis_y_res The resolution for the plotted y axis (count). Not used if
#'  use_ggplot2 == T.
#' @return This plot function returns NULL if use_ggplot2 == FALSE or a ggplot2
#'  object if use_ggplot2 == T.
#' @importFrom rlang .data
#' @export

plot_count_ensemble <- function(count_ensemble,
                                times,
                                use_ggplot2 = FALSE,
                                axis_x_res = 100,
                                axis_y_res = 1){
    # Check for ggplot2
    count_ensemble_na <- count_ensemble
    count_ensemble_na[which(count_ensemble_na == 0)] <- NA
    if(use_ggplot2){
        ggplot2_installed <- requireNamespace("ggplot2", quietly = TRUE)
        if(ggplot2_installed){
            ncols <- dim(count_ensemble_na)[2]
            count_ensemble_df <- as.data.frame(cbind(times,
                                                    count_ensemble_na))
            colnames <- c("x",
                        as.character(1:ncols))
            names(count_ensemble_df) <- colnames
            count_ensemble_df_long <- tidyr::pivot_longer(count_ensemble_df,
                                                    cols = 2:ncols,
                                                    names_to = "y",
                                                    values_to = "frequency")
            p <- ggplot2::ggplot(data = count_ensemble_df_long,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y)) +
                ggplot2::geom_raster(mapping = ggplot2::aes(fill = frequency)) +
                ggplot2::scale_fill_viridis_c(option = "B",
                                    na.value = grDevices::rgb(0, 0, 0, 0),
                                    begin = 0.15,
                                    alpha = 0.9,
                                    trans = "log") +
                ggplot2::labs(x = "Time", y = "Count")
            print(p)
            return(p)
        }else{
            stop("ggplot2 not installed.")
        }
    }else{
        image(x = 1:dim(count_ensemble_na)[1],
            y = 1:dim(count_ensemble_na)[2],
            z = count_ensemble_na,
            useRaster = T,
            col = grDevices::hcl.colors(n = 10, palette = "viridis", alpha = 0.9),
            axes = FALSE,
            xlab = "Time",
            ylab = "Count")
        axis_y_at <- seq(1, dim(count_ensemble_na)[2], axis_y_res)
        axis_x_at <- seq(1, dim(count_ensemble_na)[1], axis_x_res)
        axis(1, at = axis_x_at, labels = times[axis_x_at])
        axis(2, at = axis_y_at)
    }
    return()
}
