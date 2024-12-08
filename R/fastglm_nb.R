#### fastglm_nb ####
#' fast generalized linear model fitting for negative binomial model
#'
#' @param x input model matrix. Must be a matrix object
#' @param y numeric response vector of length nobs.
#' @param weights an optional vector of 'prior weights' to be used in the fitting process. Should be a numeric vector.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting.
#' This should be a numeric vector of length equal to the number of cases
#' @param start starting values for the parameters in the linear predictor.
#' @param etastart starting values for the linear predictor.
#' @param mustart values for the vector of means.
#' @param tol threshold tolerance for convergence. Should be a positive real number
#' @param maxit maximum number of IRLS iterations. Should be an integer
#' @param method an integer scalar with value 0 for the column-pivoted QR decomposition, 1 for the unpivoted QR decomposition,
#' @param link specification of which link function to use. Default is log
#' 2 for the LLT Cholesky, or 3 for the LDLT Cholesky
#' @return A list with the elements
#' \item{coefficients}{a vector of coefficients}
#' \item{se}{a vector of the standard errors of the coefficient estimates}
#' \item{rank}{a scalar denoting the computed rank of the model matrix}
#' \item{df.residual}{a scalar denoting the degrees of freedom in the model}
#' \item{residuals}{the vector of residuals}
#' \item{s}{a numeric scalar - the root mean square for residuals}
#' \item{fitted.values}{the vector of fitted values}
#' @export
#' @importFrom MASS theta.ml negative.binomial
#' @examples
#'
#' g <- factor(rep(letters[1:4], each = 10))
#' n <- length(g)
#' x <- model.matrix(~g)
#' b<- rnorm(4)
#'
#' mu <- exp(x%*%b)
#' theta <- 3
#'
#' y <- rnbinom(n = n, size = theta, mu = mu)
#'
#' model <- fastglm_nb(x, y)
#'
fastglm_nb <- function(x,
                       y,
                       weights = rep(1, NROW(y)),
                       offset = rep(0, NROW(y)),
                       start = NULL,
                       etastart = NULL,
                       mustart = NULL,
                       tol = 1e-08,
                       maxit = 25L,
                       method = 0L:3L, 
                       link = "log") {
    
    
    if (method[1] == 0L) {method <- 0L}
    else {method <- method}
    
    # method = 0 -> column-pivoted QR decomposition
    # method = 1 -> unpivoted QR decomposition
    # method = 2 -> LLT Cholesky
    # method = 3 -> LDLT Cholesky
    # method = 4 -> full pivoted QR decomposition
    # method = 5 -> Bidiagonal Divide and Conquer SVD
    link <- substitute(link)
    Call <- match.call()
    
    #loglikelihood of model
    loglik <- function(n, theta, mu, y, weights) {
        sum(weights * (
            lgamma(theta + y) - lgamma(theta) - lgamma(y + 1) + theta * log(theta) + y *
                log(mu + (y == 0)) - (theta + y) * log(theta + mu)
        ))
    }
    n <- NROW(y)
    
    #start fit with poisson regression
    
    fam0 <- do.call("poisson", list(link = link))
    
    fit0 <- fastglmPure(
        x = x,
        y = y,
        family = fam0,
        weights = weights,
        offset = offset,
        method = method
    )
    
    mu0 <- fit0$fitted.values
    
    theta0 <- MASS::theta.ml(
        y = y,
        mu = mu0,
        weights = weights,
        limit = maxit
    )
    
    init.theta <- theta0
    
    fam0 <- MASS::negative.binomial(theta = theta0, link = link)
    g <- fam0$linkfun
    
    eta0 <- g(mu0)
    
    iter <- 0
    
    d1 <- sqrt(2 * max(1, fit0$df.residual))
    d2 <- dev <- 1
    Lm <- loglik(n, theta0, mu0, y, weights)
    Lm0 <- Lm + 2 * d1
    
    # algorithm
    while ((abs(Lm0 - Lm) / d1 + abs(dev) / d2) > tol &&
           iter <= maxit) {
        fit <- fastglmPure(
            x = x,
            y = y,
            family = fam0,
            weights = weights,
            offset = offset,
            etastart = eta0,
            method = method
        )
        
        mu0 <- fit$fitted.values
        
        theta1 <- MASS::theta.ml(
            y = y,
            mu =  mu0,
            weights =  weights,
            limit = maxit
        )
        
        fam0 <- MASS::negative.binomial(theta = theta1, link = link)
        
        eta0 <- g(mu0)
        
        dev <- abs(theta1 - theta0)
        
        Lm0 <- Lm
        
        Lm <- loglik(n, theta1, mu0, y, weights)
        
        theta0 <- theta1
        
        iter <- iter + 1
    }
    
    class(fit) <- c("fastglm_nb", "fastglm")
    
    Call$init.theta <- signif(as.vector(init.theta), 10)
    Call$link <- link
    
    fit$call <- Call
    
    fit$theta <- theta0
    fit$SE.theta <- attr(theta0, "SE")
    fit$loglik <- loglik(n, fit$theta, fit$fitted, y, weights)
    fit$aic <- -2 * fit$loglik + 2 * fit$rank + 2
    fit$null.deviance <- fastglmPure(
        x = matrix(1, nrow = n),
        y = y,
        family = MASS::negative.binomial(theta = theta0),
        weights = weights,
        offset = offset
    )$deviance
    
    fit$df.null <- n - 1
    fit$model.matrix <- x
    fit$family <- fam0
    fit$dispersion <- 1
    fit
}

#### fastglm_nb methods ####

#' vcov method for fastglm_nb fitted objects
#'
#' @param object fastglm_nb fitted object
#' @param ... not used
#' @return the covariance matrix
#' @rdname vcov
#' @method vcov fastglm_nb
#' @export
vcov.fastglm_nb <- function(object, ...)
{
    object$XXinv * object$dispersion
}



#' summary method for fastglm_nb fitted objects
#'
#' @param object fastglm_nb fitted object
#' @param ... not used
#' @return a summary.fastglm_nb object
#' @rdname summary
#' @method summary fastglm_nb
#' @export
summary.fastglm_nb <- function(object, ...)
{
    dispersion <- 1
    p <- object$rank
    df.r <- object$df.residual
    aliased <- is.na(coef(object))  # used in print method
    
    if (p > 0)
    {
        coef   <- object$coefficients
        se     <- object$se
        zvalue <- coef / se
        
        dn <- c("Estimate", "Std. Error")
        
        pvalue <- 2 * pnorm(-abs(zvalue))
        coef.table <- cbind(coef, se, zvalue, pvalue)
        dimnames(coef.table) <- list(names(coef), c(dn, "z value", "Pr(>|z|)"))
        
        df.f <- length(aliased)
    } else
    {
        coef.table <- matrix(0, 0L, 4L)
        dimnames(coef.table) <-
            list(NULL, c("Estimate", "Std. Error", "z value", "Pr(>|Z|)"))
        covmat.unscaled <- covmat <- matrix(0, 0L, 0L)
        df.f <- length(aliased)
    }
    df.int <- if (object$intercept)
        1L
    else
        0L
    
    ## these need not all exist, e.g. na.action.
    keep <- match(
        c(
            "call",
            "terms",
            "family",
            "deviance",
            "aic",
            "contrasts",
            "df.residual",
            "null.deviance",
            "df.null",
            "iter",
            "na.action"
        ),
        names(object),
        0L
    )
    ans <- c(
        object[keep],
        list(
            deviance.resid = residuals(object, type = "deviance"),
            coefficients = coef.table,
            aliased = aliased,
            df = c(object$rank, df.r, df.f)
        )
    )
    #cov.unscaled = covmat.unscaled,
    #cov.scaled = covmat))
    
    class(ans) <- "summary.glm"
    return(ans)
}
