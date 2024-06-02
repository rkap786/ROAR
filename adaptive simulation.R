setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
## Simulate data
## n items = 10, 20, 30, 40
## Generate response data, 
## items ordered either based on fisher information or random
library(dplyr)
library(catR)
library(tidyverse)
library(mirt)
library(parallel)
library(doSNOW)
setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
## SEM and Square error compared

simbf= function(n=100, nitems=15, ndim=3, sigma, aset="bf") {
  library(MASS)
  ### Map items to factors
  bound=round(nitems/ndim)
  item_mapping=c()
  
  for (i in 1:(ndim-1)) {
    item_mapping=c(item_mapping, rep(i,bound))
  }
  
  item_mapping=c(item_mapping, rep(ndim,(nitems-(bound*i))))
  print(table(item_mapping))
  
  ## Generate item discriminations
  if(aset=="bf"){
    a= matrix(data=NA, nrow=nitems, ncol=ndim+1)
    a[,1]= rep(1, nitems)
    for (i in 2:(ndim+1)) {
      a[,i]= ifelse(item_mapping==i-1,1,0)
    }}
  
  if(aset!="bf") {
    #a[,1]= rep(1, nitems)
    a= matrix(data=rep(1,nitems*(ndim+1)), nrow=nitems, ncol=ndim+1)
  }
  
  #diff= matrix(rnorm(nitems))
  diff= matrix(runif(nitems, min=0.6, max=0.9))
  #sigma=diag(ndim+1)
  ## 
  # sigma[2,3] =0.8
  # sigma[3,2] =0.8
  # sigma[3,4] =0.8
  # sigma[4,3] =0.8
  # sigma[4,2] =0.8
  # sigma[2,4] =0.8
  theta= mvrnorm(n, mu=rep(0,ndim+1), Sigma=sigma)
  
  datasim = mirt::simdata(N=n, a=a, d=diff, itemtype = 'dich', Theta=theta,
                    sigma= diag(ndim+1))
  
  return(list(datasim, theta))
}


  
# func.catSim <- function(resp=datasim, item.bank, item.bank.u=coef,method='ML', 
#                         theta=th_u, stoplen=20,skipitems=3){
#   ## First three items are selected at random
#   ## Following that, 3 items are selected and the one to present is selected at random
#   ## ML/ MFI is used for estimation
#   
#   ni = ncol(resp)
#   np = nrow(resp)
#   list.thetas <- NULL
#   list.se <- NULL
#   list.pid <- NULL
#   list.item <- NULL
#   list.theta.final <- NULL
#   list.trialNumBlock = NULL
#   list.theta.se <- NULL
#   list.th.u <- NULL
#   list.th.u.se <- NULL
#   
#   # define cat
#   pids <- 1:nrow(resp)
#   #resp <- resp %>% select(-pid)
#   test <- list(method = 'ML', itemSelect = method, infoType = "Fisher",randomesque=5)
#   #test<-list(method='ML',itemSelect=method, randomesque=3, infoType = "Fisher") 
#   #stop <- list(rule = 'length',thr = ni) #test stops after all items are given
#   final <- list(method = 'ML')
#   #next item is selected using MFI and ability is estimated using method. randomesque specifies no prefixed items at start
#   #consider itemSelect="bOpt" for diff ability match
#   #stop<-list(rule=c("length","precision"),thr=c(stoplen,0.3))
#   stop<-list(rule="length",thr=stoplen)
#   cl<-makeCluster(parallel::detectCores())
#   registerDoSNOW(cl)
#   for (i in 1:np){
#     pid <- pids[i]
#     start <- list(nrItems=skipitems,theta = theta[i],startSelect="MFI",
#                   seed=123) #try startSelect bOpt
#     #start <- list(theta = theta[i],startSelect="MFI", randomesque=5) #try startSelect bOpt
#     random_number <- sample(c(1, ncol(resp)), size = 1)
#     #start <- list(nrItems = 10, theta = 0,randomesque=5, seed=1234)
# 
#     res <- randomCAT(itemBank = item.bank,
#                      responses = as.numeric(resp[i,]),
#                      start = start,
#                      test = test,
#                      final = final,
#                      stop = stop #, cbControl = cbList
#                      )
#     
#     
#     len=nrow(item.bank)
#     list.thetas <- c(list.thetas,c(rep(NA,skipitems-1),res$thetaProv))
#     list.pid <- c(list.pid, rep(pid, times = len))
#     list.se <- c(list.se, c(rep(NA,skipitems-1),res$seProv))
#     list.item <- c(list.item, res$testItems)
#     list.trialNumBlock <- c(list.trialNumBlock,c(1:len))
#     list.theta.final=c(list.theta.final,rep(res$thFinal, len))
#     list.theta.se= c(list.theta.se, rep(res$seFinal,len))
#     
#     ## Get unidimensional overall theta
#     start= list(fixItems=res$testItems,nrItems=len,theta = theta[i],startSelect="MFI")
#     res <- randomCAT(itemBank = item.bank.u,
#                      responses = as.numeric(resp[i,]),
#                      start = start,
#                      test = test,
#                      final = final,
#                      stop = stop #, cbControl = cbList
#     )
#     list.th.u=c(list.th.u, res$thFinal)
#     list.th.u.se=c(list.th.u.se, res$seFinal)
#     
#   }
#   stopCluster(cl)
#   
#   return(list(data.frame(pid = list.pid, 
#                     trialNumTotal = list.trialNumBlock, 
#                     item = list.item,
#                     thetaEstimate = list.thetas, 
#                     thetaSE = list.se,
#                     thetaFinal=list.theta.final,
#                     thetaSE.final= list.theta.se),
#               data.frame(thetaU.final =list.th.u,
#                 thetaU.se.final=list.th.u.se
#                 )
#               ))
# }

metrics.mse= function(true, est) {
  sqrt(sum((est-true)^2)/length(true))
}

metrics.bias= function(true, est) {
  sum(est-true)/length(true)
}

### Simulated data
ndim=3
nitems=75
np=1000
data = simbf(n=np,nitems=nitems, ndim=ndim,sigma=diag(ndim+1))
datasim=data[[1]]
theta_true=data.frame(data[[2]])

## Map items to dimensions
pa = factanal(datasim, factors=3, rotation="varimax")
loads=data.frame(pa$loadings[1:nitems,])
maxcol= max.col(loads)

loads = ifelse(pa$loadings[1:nitems,]>0.3,1,0)

bound=round(nitems/ndim)

## This part is hardcoded for 3 dimensions
data.dim=list()
data.dim[[1]]= datasim[,(1:bound)]
data.dim[[2]]= datasim[,(bound+1):(bound*2)]
data.dim[[3]]= datasim[,(bound*2+1):nitems]


names(theta_true)= c("th.gen", paste0("th.dim",c(1:ndim)))
theta_true=cbind(pid=1:np,theta_true)
theta_true = theta_true |> 
  mutate(trueEst.dim1= th.gen+th.dim1,
         trueEst.dim2= th.gen+th.dim2,
         trueEst.dim3= th.gen+th.dim3
  ) 



#sigma2= matrix(c(1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1), nrow=4, ncol=4, byrow=T)

model_unidim= mirt(datasim,1,  itemtype = 'Rasch')
th_u= data.frame(fscores(model_unidim))
th_u=cbind(pid=1:np, th=th_u)
names(th_u) = c("pid","th")
th_u_se= fscores(model_unidim,full.scores.SE=TRUE)
coef = coef(model_unidim, simplify=T)$items
coef[,2]= -coef[,2]

###### Subscales
th_dims=c()
coef_dims=list()
for (i in 1:ndim) {
  data=data.dim[[i]]
  m_dim= mirt(data,1,  itemtype = 'Rasch')
  th= fscores(m_dim)
  th_dims=cbind(th_dims,th)
  coef_dims[[i]] = coef(m_dim, simplify=T)$items
  coef_dims[[i]][,2]= -coef_dims[[i]][,2]
  # th_dim_se= fscores(m_dim,full.scores.SE=TRUE)
  # se_unidim=c(se_unidim,unidim_subscale=as.numeric(empirical_rxx(th_dim_se)))
}

### Simulate as follows: Run separately for each dimension
### Dimension 1: Run CAT using unidimensional IRT for dim 1

# CAT
th_start= rep(0,np)
results.compare=list()
th_u_final=list()
len=nitems
for (i in 1:ndim) {
  
  results= func.catSim(resp=data.dim[[i]], item.bank= (coef_dims[[i]]),
                       item.bank.u=coef, method='MFI', theta=th_start, 
                       stoplen=len,skipitems=3, pids=1:np)
  
  results_random= func.catSim(resp=data.dim[[i]], item.bank= coef_dims[[i]], 
                              item.bank.u=coef, method='random', theta=rep(0,np),
                              stoplen=len ,skipitems=3, pids=1:np)
  
  results.compare[[i]]= bind_rows(list(results[[1]], results_random[[1]]),.id="type")
  results.compare[[i]]= results.compare[[i]] |> 
    mutate(type=if_else(type==1,"CAT","Random"))
  # results.compare[[i]]= results.compare[[i]] |> 
  #   mutate(stoplen= rep(len,2*np*ni))
  # 
  ## Add in true values
  results.compare[[i]]= left_join(results.compare[[i]],theta_true)
  ## Add in full estimated IRT dim theta values
  #results.compare[[i]]= left_join(results.compare[[i]],th_dims[,i])
  
  ## Results from unidimensional

  #results.compare[[i]]= left_join(results.compare[[i]],results1.th.u)
  
  #results.full= left_join(results.full,results1.th.u)
  
  ## Results from unidimensional
  th_u_final[[i]]=bind_rows(results[[2]],results_random[[2]],.id="type")
  th_u_final[[i]]= th_u_final[[i]] |> mutate(type=if_else(type==1,"CAT","Random"))
  th_u_final[[i]]= left_join(th_u_final[[i]],th_u)
  th_u_final[[i]]= left_join(th_u_final[[i]],theta_true)
  
  # Unidimensional model to give starting parameter for next dimension
  th_start= results[[2]]$thetaU.final
  
  # th_u_final=rbind(th_u_final,cbind(results[[2]],len=rep(len,np)))
  # 
  # th_u_final_r=rbind(th_u_final_r,cbind(results_random[[2]],len=rep(len,np)))
  #  results.compare[[i]]= bind_rows(list(results[[1]], results_random[[1]]),.id="type")
  # # results.compare[[i]]= results.compare[[i]] |> 
  #   mutate(stoplen= rep(len,2*np*ni))
  
}


#th_u_res= bind_rows(th_u_final,th_u_final_r, .id="type")
th_u_res= bind_rows(th_u_final, .id = "dim")

##################################################################################################
### Calculate metrics
### Hardcoded for 3 dimensions

quartiles1 <- quantile(results.compare[[1]]$trueEst.dim1, 
                      probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
quartiles2 <- quantile(results.compare[[2]]$trueEst.dim2, 
                       probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
quartiles3 <- quantile(results.compare[[3]]$trueEst.dim3, 
                       probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))


metrics.dim1= results.compare[[1]] |>
  mutate(theta.bin = cut(trueEst.dim1, 
                         breaks = quartiles1, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(thetaFinal, thetaEstimate), 
                   bias = Metrics::bias(thetaFinal, thetaEstimate)) 


metrics.dim2= results.compare[[2]] |>
  mutate(theta.bin = cut(trueEst.dim2, 
                         breaks = quartiles2, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(trueEst.dim2, thetaEstimate), 
                   bias = Metrics::bias(trueEst.dim2, thetaEstimate)) 


metrics.dim3= results.compare[[3]] |>
  mutate(theta.bin = cut(trueEst.dim3, 
                         breaks = quartiles3, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(trueEst.dim3, thetaEstimate), 
                   bias = Metrics::bias(trueEst.dim3, thetaEstimate)) 

### Plot estimated ability

results.compare[[1]] |>
  filter(pid==10) |>
  select(pid, type, trialNumTotal, thetaEstimate, trueEst.dim1) |> 
  ggplot(aes(x=trialNumTotal, y=thetaEstimate, group=type, color=type)) +
  geom_point() +
  geom_hline(aes(yintercept= unique(trueEst.dim1)))



file= list(metrics.dim1,metrics.dim2, metrics.dim3)
### Plot metrics
for (i in 1:3) {
  titlename= paste0("Reliability, Dimension",i)
  g1=file[[i]] |>
    ggplot(aes(x=trialNumTotal, y=reliability, group=type, color=type)) +
    geom_point() +
    labs(x = "Number of test items",
         y = "Mean squared error",
         title = titlename) + 
    scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
    theme( 
      legend.title = element_blank(),
      legend.position = "bottom")
  plotname=paste0("Plots/ROAR_simulation_d",i,"_N_rel.jpeg")
  ggsave(plotname, plot=g1)
  
  titlename= paste0("Bias, Dimension",i)
  g2=file[[i]] |>
    ggplot(aes(x=trialNumTotal, y=bias, group=type, color=type)) +
    geom_point() +
    labs(x = "Number of test items",
         y = "Bias",
         title = titlename) + 
    scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
    theme( 
      legend.title = element_blank(),
      plot.margin = margin(1,1,1,1, "cm"),
      panel.spacing = unit(1, "lines"), 
      legend.position = "bottom")
  plotname=paste0("Plots/ROAR_simulation_d",i,"_N_bias.jpeg")
  ggsave(plotname, plot=g2)
  
  titlename= paste0("MSE, Dimension",i)
  g3=file[[i]] |>
    ggplot(aes(x=trialNumTotal, y=mse, group=type, color=type)) +
    geom_point() +
    labs(x = "Number of test items",
         y = "MSE",
         title = titlename) + 
    scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
    theme( 
      legend.title = element_blank(),
      plot.margin = margin(1,1,1,1, "cm"),
      panel.spacing = unit(1, "lines"), 
      legend.position = "bottom")
  plotname=paste0("Plots/ROAR_simulation_d",i,"_N_mse.jpeg")
  ggsave(plotname, plot=g3)
  
}

##################################################################################################
##################################################################################################

