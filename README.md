# Source code for SQR algorithm

## Instructions

It can be used to reproduce the simulation studies in the following paper:

**Ye Fan and Nan Lin**. *Sequential Quantile Regression for Stream Data by Least Squares.*

## Installation requirements 
```
Install the quantreg package before using the SQR function. 
```

## Code for testing PIQR
```
install.packages("quantreg")
library(quantreg)

###function for generating correlation matrix
gcov = function(p, rho){
  cov = diag(p)
  for(i in 1:p){
    for(j in 1:p){
      if(i<j) cov[i,j] = rho^{j-i}
      else cov[i,j] = cov[j,i]
    }
  }
  cov
}

###generate synthetic data
N = 100000
n = 500
M = N/n
p = 50
rhoX = 0.5
eps = 1e-06
tau = 0.75
set.seed(666)
X = matrix(rnorm(N*p), N, p)
cov = gcov(p, rhoX)
X = X%*%chol(cov)
Xinter = cbind(1, X)
beta = rep(c(-1, 1), p/2)
e = rnorm(N)
y = X%*%beta+e
beta_true = c(quantile(e, tau), beta)

###test SQR
result_SQR = SQR(Xinter, y, tau, n, eps) 
AE_SQR = sum(abs(result_SQR$beta-beta_true))
Time_SQR = result_SQR$Time

###output results
AE_SQR
Time_SQR
```
