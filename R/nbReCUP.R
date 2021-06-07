#' MCMC for nbReCUP model.
#'
#' @description This function runs a Metropolis-Hastings MCMC simulation to
#'  estimate the parameters of a negative-binomial regression model with
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

nbReCUP <- function(
                    Y,
                    X,
                    startvals,
                    niter = NULL,
                    adapt = T,
                    adapt_amount = 0.1,
                    adapt_interval = 10,
                    adapt_window = c(0.21, 0.25),
                    scales = NA,
                    priors = NULL){

    #some checking
    nX <- dim(X)[2] / dim(Y)[2]
    print(nX)
    if((nX %% 1) != 0){
        stop("Y and X must have the same number of columns, or X must have an integer multiple of Y's number of columns.")
    }
    if(is.null(niter)){
        niter <- dim(Y)[2]
    }
    pb <- progress_bar$new(total = niter)
    nparams <- length(startvals)
    N <- dim(Y)[1]
    chain <- array(dim = c(niter + 1, nparams))
    chain[1,] <- startvals
    if(adapt){
        n_adapts <- floor(niter / adapt_interval)
        acceptance <- array(dim = c(n_adapts, nparams))
        scales_matrix <- array(dim = c(n_adapts, nparams))
    }
    for(j in 1:niter){
        pb$tick()
        # adapt step
        if(adapt & (j %% adapt_interval == 0)){
            interval_index <- c(j - adapt_interval + 1):j
            diffs <- apply(chain[interval_index, ],2,diff)
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
        proposal_regression <- c(
                        propose_reg(chain[j, 1:nX], scales[1:nX]),
                        chain[j, -c(1:nX)]
                        )
        x_cols <- 1:nX + (nX * (j - 1))
        pd_proposal <- posterior(Y[, j], X[, x_cols], proposal_regression, nX)
        pd_previous <- posterior(Y[, j], X[, x_cols], chain[j, ], nX)
        accept <- exp(pd_proposal - pd_previous)
        if(runif(1) < accept){
            chain[j + 1, ] <- proposal_regression
        }else{
            chain[j + 1, ] <- chain[j, ]
        }
        # accept step for each p separately
        for(l in 1:N){
            proposal_p <- chain[j + 1, ]
            proposal_p[nX + l] <- propose_p(chain[j, nX + l], scales[nX + l])
            pd_proposal <- posterior(Y[, j], X[, x_cols], proposal_p, nX)
            pd_previous <- posterior(Y[, j], X[, x_cols], chain[j + 1, ], nX)
            accept <- exp(pd_proposal - pd_previous)
            if(runif(1) < accept){
                chain[j + 1, ] <- proposal_p
            }else{
                chain[j + 1, ] <- chain[j + 1, ]
            }
        }
    }
    alarm()
    if(adapt){
        return(list(
                samples = chain,
                acceptance = acceptance,
                scales = scales_matrix))
    }else{
        return(chain)
    }
}

likelihood <- function(Y, X, params, nX){
    pred <- exp(X %*% params[1:nX])
    loglike <- dnbinom(x = Y, size = pred, prob = params[-c(1:nX)], log = T)
    sll <- sum(loglike)
    return(sll)
}

prior <- function(params, nX, priors = NULL){
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
    return(sll)
}

posterior <- function(Y, X, params, nX){
    post <- likelihood(Y, X, params, nX) + prior(params, nX)
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
