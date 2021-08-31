plot_rece <- function(c14dates,
                    nsamples,
                    verbose = TRUE,
                    BP = TRUE,
                    resolution = -1,
                    use_ggplot2 = FALSE,
                    axis_x_res = 100,
                    axis_y_res = 1){
    calibrated_dates <- vroomfondel::build_c14_matrix(c14dates,
                                                    BP = BP,
                                                    resolution = resolution)
    times <- calibrated_dates$time_range[1]:calibrated_dates$time_range[2]
    if(verbose){
        message("Sampling count sequences...")
    }
    event_count_samples <- pbsapply(1:nsamples,
        vroomfondel::sample_event_counts,
        ce_matrix = calibrated_dates$ce_matrix,
        times = times,
        BP = BP)
    if(verbose){
        message("Tabulating count frequencies...")
    }
    event_count_freqs <- t(apply(event_count_samples,
                                1,
                                vroomfondel::tabulate_frequencies,
                                n_events=dim(c14dates)[1]))
    max_count <- vroomfondel::find_max_count(event_count_freqs)
    p <- vroomfondel::plot_count_ensemble(event_count_freqs[,2:max_count],
                        times = times,
                        use_ggplot2 = use_ggplot2,
                        axis_x_res = axis_x_res,
                        axis_y_res = axis_y_res)
    return(invisible(list(rece = event_count_freqs[,2:max_count],
                            times = times,
                            p = p)))
}
