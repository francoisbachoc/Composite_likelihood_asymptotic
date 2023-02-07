rm(list=ls(all=TRUE))
vname = c("s0p15_slepian",
          "s0p3_slepian",
          "s0p3_matern",
          "s0p15_matern",
          "s0p2_matern")

ymin = 0.5
ymax=0
for (name in vname) {
  load(paste0(name,".Rdata"))
  ymax = max(ymax,(1.2*max(m_var_finite/m_var_ass)))
}
for (name in vname) {
  load(paste0(name,".Rdata"))
  if (name == "s0p15_matern") {
    m_var_finite = m_var_finite[,1:6]
    m_var_ass = m_var_ass[,1:6]
    vn = vn[1:6]
  }
  nKL = length(vK)
  pdf(file = paste0(name,".pdf"))
  plot(x=vn,y=m_var_finite[1,]/m_var_ass[1,],type="b",lty=1,col=2,
       ylim=c(ymin,ymax),xlab="n",ylab="ratio",log = "xy",lwd=2,
       cex.axis=2,cex.lab=2,xlim=c(10^3,10^9))
  abline(h=1,col=1)
  for ( i in 2:nKL) {
    points(x=vn,y=m_var_finite[i,]/m_var_ass[i,],type="b",lty=i,col=i+1,
           ylim=c(ymin,ymax),xlab="n",ylab="ratio",lwd=2,xlim=c(10^3,10^9))
  }
  legendS = rep("",nKL)
  for(i in 1:nKL) {
    legendS[i] = paste0("K ",vK[i],", L ",vL[i]) 
  }
  legend(x = "topleft",lty=1:nKL,col=2:(nKL+1),legend=legendS,cex=2,lwd=2)
  dev.off()
}

