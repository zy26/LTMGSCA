
### Gamma Normal from scimpute
calculate_weight <-function (x, paramt){
  pz1 = paramt[1] * dgamma(x, shape = paramt[2], rate = paramt[3])
  pz2 = (1 - paramt[1]) * dnorm(x, mean = paramt[4], sd = paramt[5])
  pz = pz1/(pz1 + pz2)
  pz[pz1 == 0] = 0
  return(cbind(pz, 1 - pz))
}

fn = function(alpha, target){
  log(alpha) - digamma(alpha) - target
}

dmix <-function (x, pars) {
  pars[1] * dgamma(x, shape = pars[2], rate = pars[3]) + 
    (1 - pars[1]) * dnorm(x, mean = pars[4], sd = pars[5])
}    


update_gmm_pars = function(x, wt){
  tp_s = sum(wt)
  tp_t = sum(wt * x)
  tp_u = sum(wt * log(x))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v, 
                      extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}


get_mix = function(xdata, point=log10(1.01)){
  inits = rep(0, 5)
  inits[1] = sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {inits[1] = 0.01}
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > point]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {inits[5] = 0}
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0
  
  while(eps > 0.5) {
    wt = calculate_weight(xdata, paramt)
    paramt[1] = sum(wt[, 1])/nrow(wt)
    paramt[4] = sum(wt[, 2] * xdata)/sum(wt[, 2])
    paramt[5] = sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))
    paramt[2:3] = update_gmm_pars(x=xdata, wt=wt[,1])
    
    loglik = sum(log10(dmix(xdata, paramt)))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100) 
      break
  }
  return(paramt)
}













pdf("GN_compare.pdf")
par(mfrow=c(3,5))
for (tg_gene in c("FAP","TIGIT","COL6A3")) {
  ccc<-log(File_data[tg_gene,])
  Zcut<-min(ccc[which(File_data[tg_gene,]!=0)])
  ccc[which(ccc<Zcut0)]<-Zcut0-2
  num_c<-File_LTMG[[1]][tg_gene,1]
  aaa<-rbind(File_LTMG[[2]][tg_gene,1:num_c],File_LTMG[[3]][tg_gene,1:num_c],File_LTMG[[4]][tg_gene,1:num_c])
  aaa<-t(aaa)
  
  y<-ccc
  if(sum(y<Zcut)>0)
  {
    zz<-min(y[which(y>Zcut)])-2
    y[which(y<=Zcut)]<-zz
  }
  mm<-min(y)
  MM<-max(y)
  diff_c<-MM-mm
  MM<-MM+diff_c*0.2
  mm<-mm-diff_c*0.2
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(y,breaks=30,plot=F)
  plot(h,col="lightblue",main=paste(tg_gene,"LTMG",sep = " "))
  n<-length(ccc)*(h$breaks[2]-h$breaks[1])
  for(i in 1:nrow(aaa))
  {
    z<-dnorm(x,aaa[i,2],aaa[i,3])*aaa[i,1]*n
    #abline(v=aaa[i,2],col=c(i+1),lwd=2)
    points(z~x,type="l",col=c(i+1),lwd=2)
  }
  
  y<-ccc
  y<-y[!y<Zcut]
  mm<-min(y)
  MM<-max(y)
  diff_c<-MM-mm
  MM<-MM+diff_c*0.2
  mm<-mm-diff_c*0.2
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(y,breaks=30,plot=F)
  plot(h,col="lightblue",main=paste(tg_gene,"LTMG",sep = " "))
  n<-length(ccc)*(h$breaks[2]-h$breaks[1])
  #aaa<-t(aaa)
  for(i in 1:nrow(aaa))
  {
    z<-dnorm(x,aaa[i,2],aaa[i,3])*aaa[i,1]*n
    #abline(v=aaa[i,2],col=c(i+1),lwd=2)
    points(z~x,type="l",col=c(i+1),lwd=2)
  }
  
  #plot_peak(log(File_data[tg_gene,]),t(aaa),Zcut=Zcut0,main0 = paste(tg_gene,"LTMG",sep = " "))
  
  ccc<-File_data[tg_gene,]
  test<-log10(ccc+1.01)
  para<-get_mix(test)
  y<-test
  mm<-min(test)
  MM<-max(test)
  diff_c<-MM-mm
  MM<-MM+diff_c*0.2
  mm<-mm-diff_c*0.2
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(test,breaks=30,plot=F)
  plot(h,col="lightblue",main=paste(tg_gene,"GN",sep = " "))
  n<-length(ccc)*(h$breaks[2]-h$breaks[1])
  i<-1
  z<-dgamma(x,para[2],para[3])*para[1]*n
  points(z~x,type="l",col=c(i+1),lwd=2)
  
  i<-2
  z<-dnorm(x,para[4],para[5])*(1-para[1])*n
  points(z~x,type="l",col=c(i+1),lwd=2)
  
  
  ccc<-File_data[tg_gene,]
  test<-log10(ccc[ccc>0]+1.01)
  y<-test
  mm<-min(test)
  MM<-max(test)
  diff_c<-MM-mm
  MM<-MM+diff_c*0.2
  mm<-mm-diff_c*0.2
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(test,breaks=30,plot=F)
  plot(h,col="lightblue",main=paste(tg_gene,"GN",sep = " "))
  n<-length(ccc)*(h$breaks[2]-h$breaks[1])
  i<-1
  z<-dgamma(x,para[2],para[3])*para[1]*n
  points(z~x,type="l",col=c(i+1),lwd=2)
  
  i<-2
  z<-dnorm(x,para[4],para[5])*(1-para[1])*n
  points(z~x,type="l",col=c(i+1),lwd=2)
  
  
  
  ccc<-log(File_data[tg_gene,])
  Zcut<-min(ccc[which(File_data[tg_gene,]!=0)])
  ccc[which(ccc<Zcut0)]<-Zcut0-2
  num_c<-File_LTMG[[1]][tg_gene,1]
  aaa<-rbind(File_ZIMG[[2]][tg_gene,1:num_c],File_ZIMG[[3]][tg_gene,1:num_c],File_ZIMG[[4]][tg_gene,1:num_c])
  aaa<-t(aaa)
  
  y<-ccc
  y<-y[!y<Zcut]
  mm<-min(y)
  MM<-max(y)
  diff_c<-MM-mm
  MM<-MM+diff_c*0.2
  mm<-mm-diff_c*0.2
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(y,breaks=30,plot=F)
  plot(h,col="lightblue",main=paste(tg_gene,"ZIMG",sep = " "))
  n<-length(y)*(h$breaks[2]-h$breaks[1])
  #aaa<-t(aaa)
  for(i in 1:nrow(aaa))
  {
    z<-dnorm(x,aaa[i,2],aaa[i,3])*aaa[i,1]*n
    #abline(v=aaa[i,2],col=c(i+1),lwd=2)
    points(z~x,type="l",col=c(i+1),lwd=2)
  }
  
  
  
}
dev.off()

