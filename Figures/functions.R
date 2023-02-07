library(Rcpp)

fcov_matern = function(x1,x2,s,alpha) {
  nu = s/2
  Anu = pi / (gamma(nu)*sin(nu*pi)*2^(2*nu)*gamma(1+nu))
  M = Anu^(-1/(2*nu))*alpha*abs(outer(x1,x2,'-')) + 10^(-14)
  ((M^nu)/(2^(nu-1)*gamma(nu)))*besselK(M,nu)
}


fcov_tri = function(x1,x2,s) {
  1 - abs(outer(x1,x2,'-'))^s
}

composite_likelihood_est = function(x,z,K,L,fcov) {
  n = length(z)
  R = fcov(x[c(1:K,(K+2):(K+L+1))],x[c(1:K,(K+2):(K+L+1))])
  r = fcov(x[K+1],x[c(1:K,(K+2):(K+L+1))])
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

finite_sample_variance = function(x,K,L,fcov) {
  n = length(x)
  R = fcov(x[1:(K+L+1)],x[1:(K+L+1)])
  iR = solve(R)
  sig2 = 1/(iR[K+1,K+1])
  Mloo =  diag(1/diag(iR))%*%iR
  vproj = matrix(nrow=1,ncol=K+L+1,data=Mloo[K+1,])
  vcov = seq(from=0,to=0,length=n-K-L-1)
  for (a in 1:(n-K-L-1)) {
    Rcross = fcov(x[1:(K+L+1)],x[(1+a):(K+L+1+a)])
    vcov[a] = vproj%*%Rcross%*%t(vproj)
  }
  (1/sig2^2)*(1/((n-K-L)^2))*(4*sum((vcov^2)*seq(from=(n-K-L-1),to=1,by=-1))+2*(n-K-L)*sig2^2)
}

finite_sample_variance_2 = function(x,K,L,fcov) {
  n = length(x)
  R = fcov(x[1:(K+L+1)],x[1:(K+L+1)])
  iR = solve(R)
  sig2 = 1/(iR[K+1,K+1])
  Mloo =  diag(1/diag(iR))%*%iR
  vproj = matrix(nrow=1,ncol=K+L+1,data=Mloo[K+1,])
  vcov = seq(from=0,to=0,length=n-K-L-1)
  vr = fcov(x[1],x)
  for (a in 1:(n-K-L-1)) {
    Rcross = matrix(nrow=K+L+1,ncol=K+L+1)
    for (b in 1:(K+L+1)) {
      for (c in 1:(K+L+1)) {
          Rcross[b,c] = vr[1+abs(b-(a+c))]
        }
    }
    quad_form=0
    for (b in 1:(K+L+1))  {
      for (c in 1:(K+L+1)) {
        quad_form = quad_form + vproj[b]*vproj[c]*Rcross[b,c]
      }
    }
    vcov[a] = quad_form
  }
  (1/sig2^2)*(1/((n-K-L)^2))*(4*sum((vcov^2)*seq(from=(n-K-L-1),to=1,by=-1))+2*(n-K-L)*sig2^2)
}

finite_sample_variance_3 = function(x,K,L,fcov) {
  n = length(x)
  R = fcov(x[1:(K+L+1)],x[1:(K+L+1)])
  iR = solve(R)
  sig2 = 1/(iR[K+1,K+1])
  Mloo =  diag(1/diag(iR))%*%iR
  vproj = matrix(nrow=1,ncol=K+L+1,data=Mloo[K+1,])
  vr = fcov(x[1],x)
  Rcross = matrix(nrow=K+L+1,ncol=K+L+1,data=0)
  vcov = vcov_c(n,K,L,as.vector(vr),as.vector(vproj),Rcross)
  (1/sig2^2)*(1/((n-K-L)^2))*(4*sum((vcov^2)*seq(from=(n-K-L-1),to=1,by=-1))+2*(n-K-L)*sig2^2)
}

cppFunction('NumericVector vcov_c(int n, int K, int L, NumericVector vr, NumericVector vproj, NumericMatrix Rcross) {
  NumericVector out(n-K-L-1);
  for(int a = 0; a < (n-L-K-1); a++) {
    for(int b = 0; b < (K+L+1); b++) {
      for(int c = 0; c < (K+L+1); c++) {
        Rcross(b,c) = vr[fabs(b-(a+1+c))];
      }
    }
    double quad_form = 0;
    for(int b = 0; b < (K+L+1); b++) {
      for(int c = 0; c < (K+L+1); c++) {
        quad_form = quad_form + vproj[b]*vproj[c]*Rcross(b,c);
      }
    }
    out[a] = quad_form;
  }
  return out;
}')

finite_sample_variance_4 = function(n,K,L,fcov,v_batch_size) {
  xKLp1 = seq(from=0,by=1,length=K+L+1)/n
  R = fcov(xKLp1,xKLp1)
  iR = solve(R)
  sig2 = 1/(iR[K+1,K+1])
  Mloo =  diag(1/diag(iR))%*%iR
  vproj = matrix(nrow=1,ncol=K+L+1,data=Mloo[K+1,])
  Rcross = matrix(nrow=K+L+1,ncol=K+L+1,data=0)
  sum_total = 0
  for (i in 1:length(v_batch_size)) {
    batch_size=v_batch_size[i]
    if (i==1) {
      index_start=1
    }  else {
      index_start=sum(v_batch_size[1:(i-1)])+1 
    }
    xi = seq(from=(index_start-K-L),to=(index_start+batch_size-1+K+L),by=1)/n
    vri = fcov(0,xi)
    sum_total = sum_total + sum_c(n,batch_size,index_start,K,L,as.vector(vri),as.vector(vproj),Rcross)
  }
  (1/sig2^2)*(1/((n-K-L)^2))*(4*sum_total+2*(n-K-L)*sig2^2) 
}

cppFunction('double sum_c(int n, int batch_size, int index_start, int K, int L, NumericVector vr, NumericVector vproj, NumericMatrix Rcross) {
  NumericVector vcov(batch_size);
  for(int a = 0; a < batch_size; a++) {
    for(int b = 0; b < (K+L+1); b++) {
      for(int c = 0; c < (K+L+1); c++) {
        Rcross(b,c) = vr[fabs(index_start+a+c-b) - (index_start-K-L)];
      }
    }
    double quad_form = 0;
    for(int b = 0; b < (K+L+1); b++) {
      for(int c = 0; c < (K+L+1); c++) {
        quad_form = quad_form + vproj[b]*vproj[c]*Rcross(b,c);
      }
    }
    vcov[a] = quad_form;
  }
  double out = 0;
  for(int a = 0; a < batch_size; a++) {
    out = out + vcov[a]*vcov[a]*(n-K-L-index_start-a);
  }
  return out;
}')

finite_sample_variance_5 = function(n,K,L,fcov,v_batch_size) {
  Bn = matrix(nrow=K+L,ncol=K+L,data=0)
  for (i in 2:(K+L+1)) {
    for (j in 2:(K+L+1)) {
      Bn[i-1,j-1] = n^s*( fcov(i/n,j/n) - fcov(i/n,1/n)*fcov(j/n,1/n) )
    }
  }
  bn = solve(Bn)
  sum_total = 0
  for (i in 1:length(v_batch_size)) {
    batch_size=v_batch_size[i]
    if (i==1) {
      index_start=1
    }  else {
      index_start=sum(v_batch_size[1:(i-1)])+1 
    }
    xi1 = seq(from=1,to=K+L,by=1)/n
    vri1 = fcov(0,xi1)
    shift_min = max( 0 ,  index_start-K-L )
    shift_max = index_start+batch_size+K+L
    xi2 = seq(from=shift_min,to=shift_max,by=1)/n
    vri2 = fcov(0,xi2)
    sum_total = sum_total + sum_c_5(n^s,n,batch_size,index_start,K,L,as.vector(vri1),as.vector(vri2),bn,shift_min)
  }
  (1/((n-K-L)^2))*(4*(1/(bn[K,K])^2)*sum_total+2*(n-K-L)) 
}

cppFunction('double sum_c_5(double nps, int n, int batch_size, int index_start, int K, int L, NumericVector vr1, NumericVector vr2, NumericMatrix bn, int shift_min) {
  NumericVector vcov(batch_size);
  for(int a = 0; a < batch_size; a++) {
    double quad_form = 0;
    for(int k = 1; k < K+L+1; k++) {
      for(int l = 0; l < K+L+1; l++) {
        quad_form = quad_form + bn(K-1,k-1)*bn(K-1,l-1)*nps*(vr2[abs(index_start+a+l-k)-shift_min]-vr1[l-1]*vr2[abs(index_start+a-k)-shift_min]-vr1[k-1]*vr2[abs(index_start+a+l)-shift_min]+vr1[k-1]*vr1[l-1]*vr2[abs(index_start+a)-shift_min]);
      }
    }
    vcov[a] = quad_form;
  }
  double out = 0;
  for(int a = 0; a < batch_size; a++) {
    out = out + vcov[a]*vcov[a]*(n-L-K-index_start-a);
  }
  return out;
}')

finite_sample_variance_6 = function(n,K,L,fcov,v_batch_size) {
  Bn = matrix(nrow=K+L,ncol=K+L,data=0)
  for (i in 2:(K+L+1)) {
    for (j in 2:(K+L+1)) {
      Bn[i-1,j-1] = n^s*( fcov(i/n,j/n) - fcov(i/n,1/n)*fcov(j/n,1/n) )
    }
  }
  bn = solve(Bn)
  vcov=rep(-50,n-K-L-1)
  for (a in 1:(n-K-L-1)) {
    quad_form=0
    for (k in 1:(K+L)) {
      for (l in 1:(K+L)) {
        quad_form = quad_form + bn[K,k]*bn[K,l]*n^s*(fcov((a+l-k)/n,0)-fcov(l/n,0)*fcov((a-k)/n,0)
            -fcov(k/n,0)*fcov((a+l)/n,0) + fcov(k/n,0)*fcov(l/n,0)*fcov(a/n,0))
      }
    }
    vcov[a] = quad_form
  }
  sum_total = 0
  for (a in 1:(n-K-L-1)) {
    sum_total = sum_total+vcov[a]^2*(n-K-L-a)
  }
  (1/((n-K-L)^2))*(4*(1/(bn[K,K])^2)*sum_total+2*(n-K-L)) 
}

finite_sample_variance_7 = function(n,K,L,fcov,v_batch_size) {
  Bn = matrix(nrow=K+L,ncol=K+L,data=0)
  for (i in 2:(K+L+1)) {
    for (j in 2:(K+L+1)) {
      Bn[i-1,j-1] = n^s*( fcov(i/n,j/n) - fcov(i/n,1/n)*fcov(j/n,1/n) )
    }
  }
  bn = solve(Bn)
  sum_total = 0
  xi1 = seq(from=1,to=K+L,by=1)/n
  vri1 = fcov(0,xi1)
  shift_min = 0
  shift_max = n
  xi2 = seq(from=shift_min,to=shift_max,by=1)/n
  vri2 = fcov(0,xi2)
  sum_total =  sum_c_7(n^s,n,K,L,as.vector(vri1),as.vector(vri2),bn,shift_min)
  (1/((n-K-L)^2))*(4*(1/(bn[K,K])^2)*sum_total+2*(n-K-L)) 
}

cppFunction('double sum_c_7(double nps, int n, int K, int L, NumericVector vr1, NumericVector vr2, NumericMatrix bn, int shift_min) {
  NumericVector vcov(n-L-K-1);
  for(int a = 1; a < n-L-K; a++) {
    double quad_form = 0;
    for(int k = 1; k < K+L+1; k++) {
      for(int l = 1; l < K+L+1; l++) {
        quad_form = quad_form + bn(K-1,k-1)*bn(K-1,l-1)*nps*(vr2[abs(a+l-k)]-vr1[l-1]*vr2[abs(a-k)]-vr1[k-1]*vr2[abs(a+l)]+vr1[k-1]*vr1[l-1]*vr2[abs(a)]);
      }
    }
    vcov[a-1] = quad_form;
  }
  double out = 0;
  for(int a = 1; a < n-L-K; a++) {
    out = out + vcov[a-1]*vcov[a-1]*(n-L-K-a);
  }
  return out;
}')

asymptotic_variance = function(K,L,s,n,fcov) {
  B = matrix(nrow=K+L,ncol=K+L,data=0)
  for (i in 1:(K+L)) {
    for (j in 1:(K+L)) {
      B[i,j] = abs(i)^s + abs(j)^s - abs(i-j)^s
    }
  }
  b = solve(B)
  S = sum(b[K,]*(seq(from=1,to=K+L,by=1)^s))
  vt = seq(from=0,to=1,length=100000)
  int = mean((1-vt)*(fcov(0,vt)^2))
  (1/(n^(2*s)))*4*(1/(b[K,K]^2))*(S^4)*int
}

asymptotic_variance_K1L1 = function(s,n,fcov) {
  vt = seq(from=0,to=1,length=1000)
  int = mean((1-vt)*(fcov(0,vt)^2))
  (1/(n^(2*s)))*(1-2^(s-1))^4*(1/(1-2^(s-2))^2)*int
}


