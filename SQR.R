###SQR algorithm
SQR = function(X, Y, tau, n, eps){
  ptm = proc.time()  
  N = length(Y)
  M = N/n
  xi1 = (1-2*tau)/(tau*(1-tau))
  xi2 = sqrt(2/(tau*(1-tau)))
  bm = 2+xi1^2/xi2^2
  X1 = X[1:n,]
  y1 = y[1:n]
  result1 = rq.fit.fnb(X1, y1, tau, eps = eps)
  beta_BQR = result1$coefficients
  for(m in 2:M){
    Xm = X[(1+(m-1)*n):(m*n),]
    ym = y[(1+(m-1)*n):(m*n)]
    am = (ym-Xm%*%beta_BQR)^2/xi2^2
    vm = (sqrt(1+4*am*bm)-1)/(2*bm)
    c1m = xi1*vm
    ymnew = ym-c1m
    Xmnew = Xm
    betam = solve(t(Xmnew)%*%Xmnew)%*%(t(Xmnew)%*%ymnew)
    beta_BQR = betam/m+beta_BQR*(1-1/m)
  }
  beta_BQR[1] = quantile(ym-Xm[,-1]%*%beta_BQR[-1],tau)
  T_BQR = (proc.time() - ptm)[3]
  return(list(beta = beta_BQR, Time = T_BQR))
}

