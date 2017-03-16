# implementation of Poisson regression with lasso regularity

This R package is just a demo to show how to implement Poisson regression with lasso regularity. The algorithm implemented is the same as in the R package **glmnet** -- the coordinate descent algorithm. The core part is implemented in C++ via Rcpp. A plot function `plotx` is used to plot the result.

The following code can be used to compare the results of this package and `glmnet(..., "poisson")`
```{r}
library(glmnet)

set.seed(101)
X <- matrix(rnorm(1000, 4, 3), ncol = 10)
y <- rpois(100, 3)
fit1 <- lasso_poisson(X, y)
fit2 <- glmnet(X, y, "poisson")

par(mfrow = c(1, 2))
plotx(fit1)
plot(fit2, "lambda")
```
The two plots should look similar. 
