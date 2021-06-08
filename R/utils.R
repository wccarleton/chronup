write_to_big_matrix <- function(j, x, bigmatrix){
    m <- bigmemory::attach.big.matrix(bigmatrix)
    nx <- dim(x)[2]
    column_indeces <- 1:nx + (nx * (j - 1))
    m[, column_indeces] <- x
    return()
}

rep_big_x <- function(x, n, bigmatrix, parallel = T){
    if (!is.matrix(x)){
        warning("Coercing x to a matrix.")
        x <- as.matrix(x)
    }
    if (parallel){
        ncores <- parallel::detectCores()
        cl <- parallel::makeCluster(ncores - 1)
        parallel::clusterEvalQ(cl,{
                            wd <- getwd()
                            # devtools::load_all()
                            library(vroomfondel)
                            })
        pbapply::pbsapply(cl = cl,
                        X = 1:n,
                        FUN = write_to_big_matrix,
                        # arg name conflict between pbsapply and clusterApply
                        # forces me to use positional argument for 'x' in this
                        # one case.
                        x,
                        bigmatrix = bigmatrix)
        parallel::stopCluster(cl)
        rm(cl)
    }else{
        pbapply::pbsapply(X = 1:n,
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
