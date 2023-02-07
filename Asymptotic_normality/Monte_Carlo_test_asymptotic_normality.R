source('functions.R')

K=2
L=1
vn=c(100,500,1000,2000)
s=0.25
nmc=10000
nbreaks=50

for (n in vn) {
  vest=seq(from=0,to=0,length=nmc)
  x = seq(from=0,to=1,length=n)
  R = 1-abs(outer(x,x,'-'))^s
  cR = chol(R)
  
  for (imc in 1:nmc) {
    if (imc/100 == floor(imc/100)) {
      cat("imc = ",imc,"\n")
    }
    z = t(cR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
    vest[imc] = composite_likelihood_est(x,z,K,L,s)
  }
  var_est = var(vest)
  mean_est = mean(vest)
  pdf(file = paste0("hist_n",n,".pdf"))
  hist(vest,breaks=nbreaks,freq=FALSE)
  xhist = seq(from=-4,to=4,length=1000)
  points(xhist,dnorm(x=xhist,mean=mean_est,sd=sqrt(var_est)),
         type="l",col="red")
  dev.off()
  cat("empirical variance ",var(vest),'\n')
  cat("finite_sample variance ",finite_sample_variance(x,K,L,s),'\n')
}










