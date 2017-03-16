#include <Rcpp.h>
# include <iostream>
using namespace std;
using namespace Rcpp;

/////////////// BASIC FUNCTION DEFINITION ////////////////

// 3. calculate row mean and column mean
// [[Rcpp::export]]
NumericVector rcmean(NumericMatrix x, int margin) {
    if(margin == 1) {
        int nx = x.nrow();
        NumericVector res(nx);
        for(int i = 0; i < nx; i++) {
            res(i) = mean(x(i, _));
        }
        return res;
    } else {
        int ny = x.ncol();
        NumericVector res(ny);
        for(int j = 0; j < ny; j++) {
            res(j) = mean(x(_, j));
        }
        return res;
    }
}

// 4. calculate row sum and column sum
// [[Rcpp::export]]
NumericVector rcsum(NumericMatrix x, int margin) {
    if(margin == 1) {
        int nx = x.nrow();
        NumericVector res(nx);
        for(int i = 0; i < nx; i++) {
            res(i) = sum(x(i, _));
        }
        return res;
    } else {
        int ny = x.ncol();
        NumericVector res(ny);
        for(int j = 0; j < ny; j++) {
            res(j) = sum(x(_, j));
        }
        return res;
    }
}

// function to calculate standard deviation
// [[Rcpp::export]]
double stdC(NumericVector x) {
    return sqrt(sum(pow(x - mean(x), 2)) / x.size());
}
                     
// 5. function for standardization
// [[Rcpp::export]]
NumericVector standardizeC(NumericVector x) {
    double mu = mean(x);
    double sd = stdC(x);
    return (x - mu) / sd;
}

// rc standard deviation
// [[Rcpp::export]]
NumericVector rcstdC(NumericMatrix x, int margin) {
    if(margin == 1) {
        int nx = x.nrow();
        NumericVector res(nx);
        for(int i = 0; i < nx; i++) {
            res(i) = stdC(x(i, _));
        }
        return res;
    } else {
        int ny = x.ncol();
        NumericVector res(ny);
        for(int j = 0; j < ny; j++) {
            res(j) = stdC(x(_, j));
        }
        return res;
    }
}

// 6. function for soft-thresholding
// [[Rcpp::export]]
double SC(double z, double r) {
    if(z > 0 & r < abs(z)) {
        return z - r;
    } else if(z < 0 & r < abs(z)) {
        return z + r;
    } else {
        return 0.0;
    }
}


// function for n by p matrix times p by 1 matrix
// [[Rcpp::export]]
NumericVector mat_vec(NumericMatrix X, NumericVector y) {
    int n = X.nrow();
    NumericVector res(n);
    for(int i = 0; i < n; i++) {
        res(i) = sum(X(i, _) * y);
    }
    return res;
}

//lasso poisson for single lambda
// [[Rcpp::export]]
NumericVector lasso_poisson_singleC(NumericMatrix X, NumericVector y,double lambda = 0, NumericVector beta_ini = 0, int max_iter = 300, double tol = 1e-7) {
    NumericVector mu_x = rcmean(X, 2);
    NumericVector sd_x = rcstdC(X, 2);
    int n = X.nrow();
    int p = X.ncol() + 1;
    
    for(int j = 0; j < p - 1; j++) {
        X(_, j) = standardizeC(X(_, j));
    }
    
    NumericMatrix X1(n, p);
    X1(_, 0) = rep(1, n);
    for(int j = 1; j < p; j++) {
        X1(_, j) = X(_, j - 1);
    }
    double mu_y = mean(y);
    NumericVector beta(p), beta_new(p);
    if(beta_ini(0) == 0) {
        beta = rep(0.1, p);
    } else {
        beta = beta_ini;
    }
    NumericVector ti = mat_vec(X1, beta);
    NumericVector z = ti - 1 + y / exp(ti);
    NumericVector weight = exp(ti);
    weight = weight / sum(weight);
    
    NumericVector r;
    double Q, Q_new;
    int iter = 1;
    
    Q = sum(weight * pow(z - mat_vec(X1, beta), 2));
    while(iter < max_iter) {
        r = z - mat_vec(X1, beta);
        beta_new(0) = sum(weight * r + weight * beta(0));
        for(int j = 1; j < p; j++) {
            beta_new(j) = SC(sum(weight * X1(_, j) * r + weight * pow(X1(_, j), 2) * beta(j)), lambda / mu_y) / sum(weight * pow(X1(_, j), 2));
            beta(j) = beta_new(j);
        }
        beta(0) = beta_new(0);
        ti = mat_vec(X1, beta);
        z = ti - 1 + y / exp(ti);
        weight = exp(ti);
        weight = weight / sum(weight);
        Q_new = sum(weight * pow(z - mat_vec(X1, beta), 2));
        if(abs(Q_new - Q) < tol) break;
        Q = Q_new;
        iter += 1;
    }
    double intercept = beta(0);
    for(int i = 1; i < p; i++) {
        intercept += mu_x(i - 1) * beta(i) / sd_x(i - 1);
    }
    NumericVector ret(p);
    ret(0) = intercept;
    for(int i = 1; i < p; i++) {
        ret(i) = beta(i) / sd_x(i - 1);
    }
    return ret;
}




