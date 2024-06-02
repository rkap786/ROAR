
func.catSim <- function(resp=datasim, item.bank, item.bank.u=coef,method='ML', 
                        theta=th_u, stoplen=20,skipitems=3, pids){
  ## First three items are selected at random
  ## Following that, 3 items are selected and the one to present is selected at random
  ## ML/ MFI is used for estimation
  
  ni = ncol(resp)
  np = nrow(resp)
  list.thetas <- NULL
  list.se <- NULL
  list.pid <- NULL
  list.item <- NULL
  list.theta.final <- NULL
  list.trialNumBlock = NULL
  list.theta.se <- NULL
  list.th.u <- NULL
  list.th.u.se <- NULL
  
  # define cat
  #pids <- 1:nrow(resp)
  #resp <- resp %>% select(-pid)
  test <- list(method = 'ML', itemSelect = method, infoType = "Fisher")
  #test<-list(method='ML',itemSelect=method, randomesque=3, infoType = "Fisher") 
  #stop <- list(rule = 'length',thr = ni) #test stops after all items are given
  final <- list(method = 'ML')
  #next item is selected using MFI and ability is estimated using method. randomesque specifies no prefixed items at start
  #consider itemSelect="bOpt" for diff ability match
  #stop<-list(rule=c("length","precision"),thr=c(stoplen,0.3))
  stop<-list(rule="length",thr=stoplen)
  cl<-makeCluster(parallel::detectCores())
  registerDoSNOW(cl)
  for (i in 1:np){
    #print(i)
    pid <- pids[i]
    start <- list(nrItems=skipitems,theta = theta[i],startSelect="MFI", seed=123) #try startSelect bOpt
    #start <- list(theta = theta[i],startSelect="MFI", randomesque=5) #try startSelect bOpt
    #random_number <- sample(c(1, ncol(resp)), size = 1)
    #start <- list(nrItems = 10, theta = 0,randomesque=5, seed=1234)
    
    res <- randomCAT(itemBank = item.bank,
                     responses = as.numeric(resp[i,]),
                     start = start,
                     test = test,
                     final = final,
                     stop = stop #, cbControl = cbList
    )
    
    
    len=nrow(item.bank)
    list.thetas <- c(list.thetas,c(rep(NA,skipitems-1),res$thetaProv))
    list.pid <- c(list.pid, rep(pid, times = len))
    list.se <- c(list.se, c(rep(NA,skipitems-1),res$seProv))
    list.item <- c(list.item, res$testItems)
    list.trialNumBlock <- c(list.trialNumBlock,c(1:len))
    list.theta.final=c(list.theta.final,rep(res$thFinal, len))
    list.theta.se= c(list.theta.se, rep(res$seFinal,len))
    
    ## Get unidimensional overall theta
    start= list(fixItems=res$testItems,nrItems=len,theta = theta[i],startSelect="MFI")
    res <- randomCAT(itemBank = item.bank.u,
                     responses = as.numeric(resp[i,]),
                     start = start,
                     test = test,
                     final = final,
                     stop = stop #, cbControl = cbList
    )
    list.th.u=c(list.th.u, res$thFinal)
    list.th.u.se=c(list.th.u.se, res$seFinal)
    
  }
  stopCluster(cl)
  
  return(list(data.frame(pid = list.pid, 
                         trialNumTotal = list.trialNumBlock, 
                         item = list.item,
                         thetaEstimate = list.thetas, 
                         thetaSE = list.se,
                         thetaFinal=list.theta.final,
                         thetaSE.final= list.theta.se),
              data.frame(pid = pids,
                         thetaU.final =list.th.u,
                         thetaU.se.final=list.th.u.se,
                         test.length=rep(len,np)
              )
  ))
}


func.catSim.error <- function(resp=datasim, item.bank, item.bank.u=coef,method='ML', 
                        theta=th_u, stoplen=20,skipitems=3, pids){
  ## First three items are selected at random
  ## Following that, 3 items are selected and the one to present is selected at random
  ## ML/ MFI is used for estimation
  
  ni = ncol(resp)
  np = nrow(resp)
  list.thetas <- NULL
  list.se <- NULL
  list.pid <- NULL
  list.item <- NULL
  list.theta.final <- NULL
  list.trialNumBlock = NULL
  list.theta.se <- NULL
  list.th.u <- NULL
  list.th.u.se <- NULL
  
  ## Add error
  item.bank[,2]= item.bank[,2] + rnorm(ni,5,1)
  
  # define cat
  #pids <- 1:nrow(resp)
  #resp <- resp %>% select(-pid)
  test <- list(method = 'ML', itemSelect = method, infoType = "Fisher")
  #test<-list(method='ML',itemSelect=method, randomesque=3, infoType = "Fisher") 
  #stop <- list(rule = 'length',thr = ni) #test stops after all items are given
  final <- list(method = 'ML')
  #next item is selected using MFI and ability is estimated using method. randomesque specifies no prefixed items at start
  #consider itemSelect="bOpt" for diff ability match
  #stop<-list(rule=c("length","precision"),thr=c(stoplen,0.3))
  stop<-list(rule="length",thr=stoplen)
  cl<-makeCluster(parallel::detectCores())
  registerDoSNOW(cl)
  for (i in 1:np){
    #print(i)
    pid <- pids[i]
    start <- list(nrItems=skipitems,theta = theta[i],startSelect="MFI", seed=123) #try startSelect bOpt
    #start <- list(theta = theta[i],startSelect="MFI", randomesque=5) #try startSelect bOpt
    #random_number <- sample(c(1, ncol(resp)), size = 1)
    #start <- list(nrItems = 10, theta = 0,randomesque=5, seed=1234)
    
    res <- randomCAT(itemBank = item.bank,
                     responses = as.numeric(resp[i,]),
                     start = start,
                     test = test,
                     final = final,
                     stop = stop #, cbControl = cbList
    )
    
    
    len=nrow(item.bank)
    list.thetas <- c(list.thetas,c(rep(NA,skipitems-1),res$thetaProv))
    list.pid <- c(list.pid, rep(pid, times = len))
    list.se <- c(list.se, c(rep(NA,skipitems-1),res$seProv))
    list.item <- c(list.item, res$testItems)
    list.trialNumBlock <- c(list.trialNumBlock,c(1:len))
    list.theta.final=c(list.theta.final,rep(res$thFinal, len))
    list.theta.se= c(list.theta.se, rep(res$seFinal,len))
    
    ## Get unidimensional overall theta
    start= list(fixItems=res$testItems,nrItems=len,theta = theta[i],startSelect="MFI")
    res <- randomCAT(itemBank = item.bank.u,
                     responses = as.numeric(resp[i,]),
                     start = start,
                     test = test,
                     final = final,
                     stop = stop #, cbControl = cbList
    )
    list.th.u=c(list.th.u, res$thFinal)
    list.th.u.se=c(list.th.u.se, res$seFinal)
    
  }
  stopCluster(cl)
  
  return(list(data.frame(pid = list.pid, 
                         trialNumTotal = list.trialNumBlock, 
                         item = list.item,
                         thetaEstimate = list.thetas, 
                         thetaSE = list.se,
                         thetaFinal=list.theta.final,
                         thetaSE.final= list.theta.se),
              data.frame(pid = pids,
                         thetaU.final =list.th.u,
                         thetaU.se.final=list.th.u.se,
                         test.length=rep(len,np)
              )
  ))
}


metrics.mse= function(true, est) {
  sqrt(sum((est-true)^2)/length(true))
}

metrics.bias= function(true, est) {
  sum(est-true)/length(true)
}