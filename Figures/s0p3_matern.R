source('functions.R')

###############################
#reglages specifiques
###############################
name = "s0p3_matern"
s=0.3
vK=c(2,2,4,4)
vL=c(0,2,0,4)
vn=c(10^3,10^4,10^5,10^6,10^7,10^8)
alpha=1
fcov = function(x1,x2) {
  fcov_matern(x1,x2,s,alpha)
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
    cat("n= ",n,"K=",K,"L=",L,"\n")
    x = seq(from=0,to=1,length=n)
    m_var_finite[i,j] = finite_sample_variance(x,K,L,fcov)
    m_var_ass[i,j] = asymptotic_variance(K,L,s,n,fcov)
  }
}



###############################
#Sauvegarde des resultats
###############################
save(s,vK,vL,vn,fcov,m_var_finite,m_var_ass,file=paste0(name,".Rdata"))

