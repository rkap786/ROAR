imv.binary<-function(y, #outcomes
                     p1,#baseline
                     p2 #enhanced
) {
  ##
  ll<-function(x,p) {
    z<-log(p)*x+log(1-p)*(1-x)
    z<-sum(z)/length(x)
    exp(z)
  }    
  loglik1<-ll(y,p1)
  loglik2<-ll(y,p2)
  getcoins<-function(a) {
    f<-function(p,a) abs(p*log(p)+(1-p)*log(1-p)-log(a))
    nlminb(.5,f,lower=0.001,upper=.999,a=a)$par
  }
  c1<-getcoins(loglik1)
  c2<-getcoins(loglik2)
  ew<-function(p1,p0) (p1-p0)/p0
  imv<-ew(c2,c1)
  imv
}



getrespmatrix= function(data) {
  data = data |>
    dplyr::select(subj,item,correct) |>
    pivot_wider(values_from= correct,
                names_from = item)
  
  #names(data) = paste0("item",names(data))
  #names(data)[1] = "personID"
  data
}


getp<-function(m,ability='EAP',id,
               x #dataframe containing oos to merge in
) { 
  co<-coef(m)
  nms<-names(co)
  co<-do.call("rbind",co[-length(co)])
  item<-data.frame(item=nms[-length(nms)],co)
  ##
  th<-fscores(m,method=ability)
  stud<-data.frame(id=id,th=th[,1])
  ##
  x<-merge(x[x$oos==1,],stud)
  x<-merge(x,item)
  ##
  kk<-x$th*x$a+x$d
  kk<-exp(kk)
  x$p<-x$g+(x$u-x$g)*kk/(1+kk)
  x
}