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

error.datasim= function(datasim, per=0.05) {
  ni=ncol(datasim)
  np= nrow(datasim)
  x= round(np*per,0)
  for (i in 1:ncol(datasim)) {
    rnd = sample.int(np,x)
    xnew = rbinom(x, 1, 0.5)
    datasim[rnd,i]=xnew
  }
  datasim
}


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

####### Add error 
datasim.original = datasim
datasim = error.datasim(datasim, per=0.05)


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

############ Generate noisy difficulty estimates




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

quartiles1 <- quantile(results.compare[[1]]$th.dim1, 
                      probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
quartiles_true <- quantile(results.compare[[2]]$trueEst.dim1, 
                       probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
# quartiles3 <- quantile(results.compare[[3]]$trueEst.dim3, 
#                        probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))


metrics.dim1= results.compare[[1]] |>
  mutate(theta.bin = cut(th.dim1, 
                         breaks = quartiles1, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(th.dim1, thetaEstimate), 
                   bias = Metrics::bias(th.dim1, thetaEstimate)) 


metrics.dim1_true= results.compare[[1]] |>
  mutate(theta.bin = cut(trueEst.dim1, 
                         breaks = quartiles_true, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(trueEst.dim1, thetaEstimate), 
                   bias = Metrics::bias(trueEst.dim1, thetaEstimate)) 


file= list(metrics.dim1,metrics.dim1_true)
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

