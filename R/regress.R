#' Bayesian regression with uncertainty propagation.
#'
#' @description This function runs a Metropolis-Hastings MCMC simulation to
#'  estimate the parameters of a specified regression model with
#'  chronological uncertainty propagation.
#'
#' @param Y Matrix containing integer event-counts. Each column should contain
#'  one probable event-count sequence. The number of columns needs to be equal
#'  to or greater than the number of desired MCMC iterations.
#' @param X A matrix containing the covariates for the model (independent
#'  variables). The matrix will be comprise blocks, each of which contains N
#'  columns where N is the number of covariates including an intercept column
#'  if desired. Thus, the total number of columns in X must be a multiple of
#'  the number of covariates (with intercept) and the number of desired MCMC
#'  iterations.
#' @param model A string specifying the type of model to use. Currently only
#'  Negative-Binomial ('nb') and Poisson ('pois') models are supported.
#' @param startvals A numeric vector of starting values for the MCMC.
#' @param niter (optional) An integer scalar indicating the number of MCMC
#'  iterations if different from the number of columns in Y.
#' @param adapt A numeric vector of starting values for the MCMC
#' @param adapt_amount A scalar containing the ratio by which the proposal
#'  variance is changed during the MCMC adapt step when adapt=T---e.g.,
#'  adapt_amount=0.1 (default) grows/shrinks the variance by 10%.
#' @param adapt_interval A numeric scalar indicating the the n-th mcmc
#'  iteration during which the simulation will attempt to adapt the proposal
#'  variances.
#' @param adapt_window A numeric vector with two elements that define the lower
#'  and upper bounds of the target acceptance value for the MCMC. This is only
#'  used when adapt=T.
#' @param scales A numeric vector of scales for the proposal distributions. The
#'  order is important. Each element refers to the scale of a proposal function
#'  for the model parameters in the following order: (b_1, b_2, ..., b_nX, p_1,
#'  p_2, ..., p_n), where `b` refers to a regression coefficient and `n` refers
#'  to the number of observations in the count sequence.
#' @param priors A numeric vector of parameter values for the model priors---the
#'  order is important and corresponds to order of appearance of the prior
#'  density functions in the source code for the `prior` function.
#' @return If adapt=T, returns a list that includes the mcmc samples,
#'  acceptance rates for all model parameters, and the proposal function scales
#'  for all parameters.
#' @import pbapply progress truncnorm stats utils
#' @export

regress <- function(Y,
                    X,
                    model = "nb",
                    startvals = NULL,
                    niter = NULL,
                    adapt = T,
                    adapt_amount = 0.1,
                    adapt_interval = 10,
                    adapt_window = c(0.21, 0.25),
                    scales = NA,
                    priors = NULL){
    # Check argument validity and set up additional parameters
    is_matrix <- is.matrix(Y) | is.matrix(X)
    bigmemory_installed <- requireNamespace("bigmemory", quietly = TRUE)
    if(bigmemory_installed){
        is_big_matrix <- bigmemory::is.big.matrix(Y) |
                        bigmemory::is.big.matrix(X)
    }else{
        is_big_matrix <- FALSE
    }
    # If neither Y nor X are matrices or big matrix pointers, fail
    if (!is_matrix & !is_big_matrix){
        stop("Y and X must be matrices or big matrix pointers.")
    }
    if (!(model %in% c("nb", "pois"))){
        stop("Invalid model. Must be one of ('nb', 'pois').")
    }
    nX <- dim(X)[2] / dim(Y)[2]
    if ((nX %% 1) != 0){
        alert <- paste("Y and X must have the same number of columns, ",
                        "or X must have an integer multiple of Y's ",
                        "number of columns",
                        sep = "")
        stop(alert)
    }
    if (is.null(niter)){
        niter <- dim(Y)[2]
    }
    N <- dim(Y)[1]
    # the type of model determines the number of parameters in the model, which
    # is important for establishing the matrix intended to contain the mcmc
    # samples and for several steps in the mcmc for loop.
    if (model == "nb"){
        nparams <- nX + N
    }else if (model == "pois"){
        nparams <- nX
    }
    # If the user hasn't provuded any starting values to initialize the mcmc,
    # we can provide some (probably dumb) defaults.
    if (is.null(startvals)){
        alert <- paste("No starting values provided ('startvals == NULL').",
                        "Using defaults.",
                        sep = " ")
        warning(alert)
        if (model == "nb"){
            startvals <- c(rep(0, nX), rep(0.5, N))
        }else if (model == "pois"){
            startvals <- rep(0, nX)
        }
    }
    chain <- array(dim = c(niter + 1, nparams))
    chain[1,] <- startvals
    if (adapt){
        n_adapts <- floor(niter / adapt_interval)
        acceptance <- array(dim = c(n_adapts, nparams))
        scales_matrix <- array(dim = c(n_adapts, nparams))
    }
    # monitor progress
    pb <- progress_bar$new(total = niter)
    for (j in 1:niter){
        pb$tick()
        # adapt step
        if(adapt & (j %% adapt_interval == 0)){
            interval_index <- c(j - adapt_interval + 1):j
            diffs <- apply(chain[interval_index, ], 2, diff)
            acceptance_rate <- 1 - colMeans(diffs == 0)
            acceptance[j / adapt_interval, ] <- acceptance_rate
            scales <- unlist(
                        mapply(
                            adaptScale,
                            acceptance_rate,
                            scales,
                            MoreArgs = list(
                                        adapt_amount = adapt_amount,
                                        adapt_window = adapt_window)
                            )
                        )
            scales_matrix[j / adapt_interval, ] <- scales
        }
        # accept step for main regression params
        proposal_regression <- c(propose_reg(chain[j, 1:nX], scales[1:nX]),
                                chain[j, -c(1:nX)])
        # X is a matrix representing uncertainty in the covariates, so we have
        # to walk over that matrix during the MCMC as well.
        x_cols <- 1:nX + (nX * (j - 1))
        pd_proposal <- posterior(Y[, j],
                                X[, x_cols],
                                proposal_regression,
                                nX,
                                model,
                                priors)
        pd_previous <- posterior(Y[, j],
                                X[, x_cols],
                                chain[j, ],
                                nX,
                                model,
                                priors)
        accept <- exp(pd_proposal - pd_previous)
        if (runif(1) < accept){
            chain[j + 1, ] <- proposal_regression
        }else{
            chain[j + 1, ] <- chain[j, ]
        }
        # if the Negative-Binomial model is used, then we have to estimate a
        # 'p' parameter for every 'y' observation (of which there will be N).
        # Experience has shown that this should be handled in a Gibbs-step as
        # follows.
        if (model == "nb"){
            for (l in 1:N){
                proposal_p <- chain[j + 1, ]
                proposal_p[nX + l] <- propose_p(chain[j, nX + l],
                                                scales[nX + l])
                pd_proposal <- posterior(Y[, j],
                                        X[, x_cols],
                                        proposal_p,
                                        nX,
                                        model,
                                        priors)
                pd_previous <- posterior(Y[, j],
                                        X[, x_cols],
                                        chain[j + 1, ],
                                        nX,
                                        model,
                                        priors)
                accept <- exp(pd_proposal - pd_previous)
                if(runif(1) < accept){
                    chain[j + 1, ] <- proposal_p
                }else{
                    chain[j + 1, ] <- chain[j + 1, ]
                }
            }
        }
    }
    alarm()
    if (adapt){
        return(list(
                samples = chain,
                acceptance = acceptance,
                scales = scales_matrix))
    }else{
        return(chain)
    }
}

likelihood <- function(Y, X, params, nX, model){
    if (model == "nb"){
        pred <- exp(X %*% params[1:nX])
        loglike <- dnbinom(x = Y,
                        size = pred,
                        prob = params[-c(1:nX)],
                        log = T)
        sll <- sum(loglike)
    }else if (model == "pois"){
        pred <- exp(X %*% params[1:nX])
        loglike <- dpois(x = Y,
                        lambda = pred,
                        log = T)
        sll <- sum(loglike)
    }
    return(sll)
}

prior <- function(params, nX, model, priors = NULL){
    if (model == "nb"){
        if(is.null(priors)){
            priors <- c(0, 1000, 1e-7, 1 - 1e-7)
        }
        B_priors <- dnorm(
                        x = params[1:nX],
                        mean = priors[1],
                        sd = priors[2],
                        log = T)
        p_priors <- dunif(
                        x = params[-c(1:nX)],
                        min = priors[3],
                        max = priors[4],
                        log = T)
        sll <- sum(c(B_priors, p_priors))
    }else if (model == "pois"){
        if(is.null(priors)){
            priors <- c(0, 1000)
        }
        B_priors <- dnorm(
                        x = params[1:nX],
                        mean = priors[1],
                        sd = priors[2],
                        log = T)
        sll <- sum(B_priors)
    }
    return(sll)
}

posterior <- function(Y, X, params, nX, model, priors){
    post <- likelihood(Y, X, params, nX, model) +
            prior(params, nX, model, priors)
    return(post)
}

propose_reg <- function(Bs, v){
    nX <- length(Bs)
    Bs <- rnorm(nX, mean = Bs, sd = v)
    return(Bs)
}

propose_p <- function(p, v){
    p <- rtruncnorm(1, a = 1e-7, b = 1 - 1e-7, mean = p, sd = v)
    return(p)
}

adaptScale <- function(acceptance_rate, s, adapt_amount, adapt_window){
    if(acceptance_rate < adapt_window[1]){
        return(s * (1 - adapt_amount))
    }else if(acceptance_rate > adapt_window[2]){
        return(s * (1 + adapt_amount))
    }else {
        return(s)
    }
}
