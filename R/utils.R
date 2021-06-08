write_to_big_matrix <- function(j, x, bigmatrix){
    m <- bigmemory::attach.big.matrix(bigmatrix)
    nx <- dim(x)[2]
    column_indeces <- 1:nx + (nx * (j - 1))
    m[, column_indeces] <- x
    return()
}

rep_big_x <- function(x, n, bigmatrix, parallel = T){
    if (parallel){
        ncores <- parallel::detectCores()
        cl <- parallel::makeCluster(ncores - 1)
        parallel::clusterEvalQ(cl,{
                            wd <- getwd()
                            devtools::load_all()
                            })
        pbapply::pbsapply(
                    cl = cl,
                    X = 1:n,
                    FUN = write_to_big_matrix,
                    # arg name conflicts between pbsapply and clusterApply
                    # forced me to use positional passing for 'x' here
                    x,
                    bigmatrix = bigmatrix)
        parallel::stopCluster(cl)
        rm(cl)
    }else{
        pbapply::pbsapply(
                    cl = cl,
                    X = 1:n,
                    FUN = write_to_big_matrix,
                    x = x,
                    bigmatrix = bigmatrix)
    }
    return()
}

approx_c14 <- function(x, t1, t2, r){
    n <- length(x)
    funs <- lapply(x, approxfun)
    y_list <- lapply(1:n, function(j)funs[[j]](seq(t1, t2, r)))
    y_mat <- do.call(cbind, y_list)
    y_mat[which(is.na(y_mat))] <- 0
    return(y_mat)
}
