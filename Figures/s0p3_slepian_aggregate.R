rm(list=ls(all=TRUE))
source('functions.R')

###############################
#reglages specifiques
###############################
name = "s0p3_slepian"
s=0.3
vK=c(2,2,4,4)
vL=c(0,2,0,4)
vn=c(10^3,10^4,10^5,10^6,10^7,10^8,10^9)
fcov = function(x1,x2) {
  fcov_tri(x1,x2,s)
}

###############################
#calcul des variances
###############################
nKL = length(vK)
nn = length(vn)
m_var_finite_tot = matrix(nrow=nKL,ncol=nn,data=-1)
m_var_ass_tot = matrix(nrow=nKL,ncol=nn,data=-1)


for (i in 1:nKL) {
  for(j in 1:nn) {
    load(paste0("s0p3_slepian_",i,".Rdata"))
    m_var_finite_tot[i,j] = m_var_finite[1,j]
    m_var_ass_tot[i,j] = m_var_ass[1,j]
  }
}
m_var_finite = m_var_finite_tot
m_var_ass =  m_var_ass_tot

###############################
#Sauvegarde des resultats
###############################
vK=c(2,2,4,4)
vL=c(0,2,0,4)
save(s,vK,vL,vn,fcov,m_var_finite,m_var_ass,file=paste0(name,".Rdata"))

