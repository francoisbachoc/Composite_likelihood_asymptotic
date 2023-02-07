rm(list=ls(all=TRUE))
source('functions.R')

###############################
#reglages specifiques
###############################
name = "s0p3_slepian_3"
s=0.3
vK=c(4)
vL=c(0)
vn=c(10^3,10^4,10^5,10^6,10^7,10^8,10^9)
l_batch_size = list()
l_batch_size[[1]] = c(rep(10^2,9),10^2-vK[1]-vL[1]-1)
l_batch_size[[2]] = c(rep(10^3,9),10^3-vK[1]-vL[1]-1)
l_batch_size[[3]] = c(rep(10^4,9),10^4-vK[1]-vL[1]-1)
l_batch_size[[4]] = c(rep(10^5,9),10^5-vK[1]-vL[1]-1)
l_batch_size[[5]] = c(rep(10^6,9),10^6-vK[1]-vL[1]-1)
l_batch_size[[6]] = c(rep(10^7,9),10^7-vK[1]-vL[1]-1)
l_batch_size[[7]] = c(rep(10^7,99),10^7-vK[1]-vL[1]-1)
fcov = function(x1,x2) {
  fcov_tri(x1,x2,s)
}

###############################
#calcul des variances
###############################
nKL = length(vK)
nn = length(vn)
m_var_finite = matrix(nrow=nKL,ncol=nn,data=-1)
m_var_ass = matrix(nrow=nKL,ncol=nn,data=-1)


for (i in 1:nKL) {
  for(j in 1:nn) {
    K = vK[i]
    L = vL[i]
    n = vn[j]
    v_batch_size = l_batch_size[[j]]
    cat("n= ",n,"K=",K,"L=",L,"\n")
    m_var_finite[i,j] = finite_sample_variance_5(n,K,L,fcov,v_batch_size)
    m_var_ass[i,j] = asymptotic_variance(K,L,s,n,fcov)
  }
}



###############################
#Sauvegarde des resultats
###############################
save(s,vK,vL,vn,fcov,m_var_finite,m_var_ass,file=paste0(name,".Rdata"))

