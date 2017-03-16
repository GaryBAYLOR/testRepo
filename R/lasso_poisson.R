## ======================================
## Poisson Regression with LASSO Penality
## ======================================

#' Standardize a vector
#'
#' This function standardizes a numeric vector to have mean zero and variacne one.
#'
#' @param x a numeric vector
#' @return A numeric vector of the same length as input vector.
#' @examples
#' x <- rnorm(100)
#' y <- standardize(x)
#' mean(y)
#' sum(y^2)
standardize <- function(x) {
	mu <- mean(x, na.rm = TRUE)
	sd <- sqrt(sum((x-mean(mu))^2)/length(x))
	(x - mu) / sd
}

#' Poisson regression with lasso regularity
#'
#' This function performs possion regression with lasso regularity.
#'
#' The coordinate descent is used as in the package \code{glmnet}.
#'
#' @param X a matrix for predictors
#' @param y a non-negative integer vector
#' @param lambda a positive numeric vector or scalar.
#' @param max_iter the maximum iteration in coordiate descent algorithm
#' @param tol the tolorance for the coordinate descent algorithm
#'
#' @return a matrix of class \code{dgCMatrix}.
#' @examples
#' x <- matrix(rnorm(1000), ncol = 5)
#' y <- rpois(200, 3)
#' fit <- lasso_poisson(x, y)
#'
lasso_poisson <- function(X, y, lambda, max_iter = 300L, tol = 1e-7) {
	p <- ncol(X)
	n <- nrow(X)

    X2 <- X
	for(j in 1:p) {
	    X2[, j] <- standardize(X2[, j])
     }

	lambda.max <- max(abs(t(X2) %*% y / n)[-1])

	if(missing(lambda)) {
		nlambda <- 100
	    lambda <- numeric(nlambda)
	    lambda[1] <- lambda.max
	    for(i in 2:nlambda) {
		    lambda[i] <- lambda[i - 1] / 1.097499
		    if(lambda[i] < 5e-4) break
	    }
	    lambda <- lambda[lambda > 0]
	}

	k <- length(lambda)
	if(k == 1) {
		ret <- lasso_poisson_single(X, y, lambda)
		ret <- c(ret, lambda)
		names(ret) <- c("intercept", paste0("V", 1:p), "lambda")
	} else {
		res <- matrix(0, ncol = p + 1, nrow = k)
	    res[1, ] <- lasso_poisson_singleC(X, y, lambda[1], max_iter = max_iter, tol = tol)
	    for(i in 2:k) {
		    res[i, ] <- lasso_poisson_singleC(X, y, lambda[i], res[i-1, ], max_iter = max_iter, tol = tol)
	    }
	    ret <- cbind(res, lambda)
	    colnames(ret) <- c("intercept", paste0("V", 1:p), "lambda")
	    ret <- as(ret, "dgCMatrix")
	}
	ret
}

#' Plotting Poisson regression with lasso regularity
#'
#' A function used to plot the Poisson regression with lasso regularity.
#' @param fit an object which is an output from the function \code{lasso_poisson}.
#' @examples
#' x <- matrix(rnorm(1000), ncol = 5)
#' y <- rpois(200, 3)
#' fit <- lasso_poisson(x, y)
#' plotx(fit)
#' ## the plot should be similar to the one in package glment
#' ## fit2 <- glmnet(x, y)
#' ## plot(fit2, "lambda")
plotx <- function(fit) {
	nc <- ncol(fit)
	lambda <- fit[, nc]
	rg <- range(fit[, -c(1, nc)])
	loglam <- log(lambda)
	plot(1, 1, xlim = range(loglam), ylim = rg, main = "Poisson LASSO",
	xlab = "log(lambda)", ylab = "Parameter Estimates")
	for(i in 2:(nc - 1)) {
		lines(loglam, fit[, i], col = i - 1)
	}
}








