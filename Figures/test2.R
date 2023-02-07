rm(list=ls(all=TRUE))
source('functions.R')

###############################
#reglages specifiques
###############################
name = "s0p3_matern"
s=0.3
K=2
L=3
n=10^7
alpha=1
fcov = function(x1,x2) {
  fcov_matern(x1,x2,s,alpha)
}
x = seq(from=0,length=n,by=1/n)
v_batch_size = c(rep(10^6,9),10^6-K-L-1)

start_time <- Sys.time()
var_finite_3 = finite_sample_variance_3(x,K,L,fcov)
cat("var_finite_3 = ", var_finite_3,"\n")
end_time <- Sys.time()
end_time - start_time
#
start_time <- Sys.time()
var_finite_4 = finite_sample_variance_4(n,K,L,fcov,v_batch_size)
cat("var_finite_4 = ", var_finite_4,"\n")
end_time <- Sys.time()
end_time - start_time
#
start_time <- Sys.time()
var_finite_5 = finite_sample_variance_5(n,K,L,fcov,v_batch_size)
cat("var_finite_4 = ", var_finite_4,"\n")
end_time <- Sys.time()
end_time - start_time

