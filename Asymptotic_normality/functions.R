
composite_likelihood_est = function(x,z,K,L,s) {
  n = length(z)
  R = 1 - abs(outer(x[c(1:K,(K+2):(K+L+1))],x[c(1:K,(K+2):(K+L+1))],'-'))^s
  r = 1 - abs( x[K+1] -x[c(1:K,(K+2):(K+L+1))]  )^s
  r = matrix(nrow=K+L,ncol=1,data=r)
  vsig2 = rep(1 - t(r)%*%solve(R)%*%r,n-K-L)
  vproj =  t(r)%*%solve(R)
  vpred = seq(from=0,to=0,length=n-K-L)
  for (k in 1:K) {
    vpred = vpred + vproj[K+1-k]*z[(K+1-k):(n-L-k)]
  }
  for (l in 1:L) {
    vpred = vpred + vproj[K+l]*z[(K+1+l):(n-L+l)]
  }
  verr = z[(K+1):(n-L)] - vpred
  mean( (verr^2) / vsig2)
}

finite_sample_variance = function(x,K,L,s) {
  n = length(x)
  R = 1 - abs(outer(x[1:(K+L+1)],x[1:(K+L+1)],'-'))^s
  iR = solve(R)
  sig2 = 1/(iR[K+1,K+1])
  Mloo =  diag(1/diag(iR))%*%iR
  vproj = matrix(nrow=1,ncol=K+L+1,data=Mloo[K+1,])
  vcov = seq(from=0,to=0,length=n-K-L-1)
  for (a in 1:(n-K-L-1)) {
    Rcross = 1 - abs(outer(x[1:(K+L+1)],x[(1+a):(K+L+1+a)],'-'))^s
    vcov[a] = vproj%*%Rcross%*%t(vproj)
  }
  (1/sig2^2)*(1/((n-K-L)^2))*(4*sum((vcov^2)*seq(from=(n-K-L-1),to=1,by=-1))+2*(n-K-L)*sig2^2)
}


asymptotic_variance = function(K,L,s,n) {
  B = matrix(nrow=K+L,ncol=K+L,data=0)
  for (i in 1:(K+L)) {
    for (j in 1:(K+L)) {
      B[i,j] = abs(i)^s + abs(j)^s - abs(i-j)^s
    }
  }
  b = solve(B)
  S = sum(b[K,]*(seq(from=1,to=K+L,by=1)^s))
  vt = seq(from=0,to=1,length=1000)
  int = mean((1-vt)*((1-vt^s)^2))
  (1/(n^(2*s)))*4*(1/(b[K,K]^2))*(S^4)*int
}

asymptotic_variance_K1L1 = function(s,n) {
  vt = seq(from=0,to=1,length=1000)
  int = mean((1-vt)*((1-vt^s)^2))
  (1/(n^(2*s)))*(1-2^(s-1))^4*(1/(1-2^(s-2))^2)*int
}


