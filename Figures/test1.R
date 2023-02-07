rm(list=ls(all=TRUE))
source('functions.R')

###############################
#reglages specifiques
###############################
name = "s0p3_matern"
s=0.3
K=2
L=3
n=30
alpha=1
fcov = function(x1,x2) {
  fcov_matern(x1,x2,s,alpha)
}
x = seq(from=0,length=n,by=1/n)
v_batch_size = c(10,10,10-K-L-1)  #total must be n-K-L-1 (number of elements in the sum)
#
start_time <- Sys.time()
var_finite_1 = finite_sample_variance(x,K,L,fcov)
cat("var_finite_1 = ", var_finite_1,"\n")
end_time <- Sys.time()
end_time - start_time
#
start_time <- Sys.time()
var_finite_2 = finite_sample_variance_2(x,K,L,fcov)
cat("var_finite_2 = ", var_finite_2,"\n")
end_time <- Sys.time()
end_time - start_time
#
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
var_finite_6 = finite_sample_variance_6(n,K,L,fcov,v_batch_size)
cat("var_finite_6 = ", var_finite_6,"\n")
end_time <- Sys.time()
end_time - start_time
#
start_time <- Sys.time()
var_finite_7 = finite_sample_variance_7(n,K,L,fcov,v_batch_size)
cat("var_finite_7 = ", var_finite_7,"\n")
end_time <- Sys.time()
end_time - start_time
#
start_time <- Sys.time()
var_finite_5 = finite_sample_variance_5(n,K,L,fcov,v_batch_size)
cat("var_finite_5 = ", var_finite_5,"\n")
end_time <- Sys.time()
end_time - start_time
