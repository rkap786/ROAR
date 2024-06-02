## Simulation
# simulate multidimensional theta
library(MASS)
library(mirt)
library(dplyr)
library(tidyverse)


simbf= function(n=100, nitems=15, ndim=3, sigma, aset="bf") {
#n=nstudents[1]
#nitems=nitems_list[2]
#ndim=3 # number of specific factors

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

diff= matrix(rnorm(nitems))
#sigma=diag(ndim+1)
## 
# sigma[2,3] =0.8
# sigma[3,2] =0.8
# sigma[3,4] =0.8
# sigma[4,3] =0.8
# sigma[4,2] =0.8
# sigma[2,4] =0.8
theta= mvrnorm(n, mu=rep(0,ndim+1), Sigma=sigma)

datasim = simdata(N=n, a=a, d=diff, itemtype = 'dich', Theta=theta,
             sigma= diag(ndim+1))


### Factor analysis
# pa = factanal(datasim, factors=3, rotation="varimax")
# print(pa, digits=2, cutoff=0.3, sort=TRUE)
# dim(pa$loadings)
# loads= round(pa$loadings[ 1:74,], 3)
# loads= data.frame(loads) |>
#   mutate(across(starts_with("Factor"), ~ifelse(.<0.3, NA, .)))
# loads$item= rownames(loads)

### 3 dimensional model as generated
model_unidim= mirt(datasim,1,  itemtype = 'Rasch')
m2_u= M2(model_unidim)



## Collect standard errors
se_unidim=c()
se_bf=c()
### Unidimensional model
th_u= fscores(model_unidim)
th_u_se= fscores(model_unidim,full.scores.SE=TRUE)
se_unidim=c(se_unidim,gen=as.numeric(empirical_rxx(th_u_se)))

cor(th_u, theta[,1])
p_u=data.frame(probtrace(model_unidim, th_u)) |>
  dplyr::select(contains("P.1")) |>
  mutate(id= 1:n) |>
  pivot_longer(
    cols=contains("P.1"),
    names_to = "item",
    values_to = "pr_u"
  ) |>
  mutate(item= str_remove(item, ".P.1"))

###### Subscales
th_dims=c()
for (i in 1:ndim) {
  m_dim= mirt(datasim[,(bound*(i-1)+1):(bound*i)],1,  itemtype = 'Rasch')
  th= fscores(m_dim)
  th_dims=cbind(th_dims,th)
  th_dim_se= fscores(m_dim,full.scores.SE=TRUE)
  se_unidim=c(se_unidim,unidim_subscale=as.numeric(empirical_rxx(th_dim_se)))
}


### Bifactor
specific= item_mapping

bfactormod <- bfactor(datasim, model=specific, 
                      quadpts = 9, TOL = 1e-3) #itemtype default is 2PL
m2_bf= M2(bfactormod, QMC=TRUE)
coef = coef(bfactormod, simplify=T)$items
th_bf= fscores(bfactormod,  QMC=TRUE)
th_bf_se= fscores(bfactormod,  QMC=TRUE,full.scores.SE=T)
se_bf=c() 
se_bf=c(se_bf,bf=as.numeric(empirical_rxx(th_bf_se)))
p_bf = data.frame(probtrace(bfactormod, th_bf)) |>
  dplyr::select(contains("P.1")) |>
  mutate(id=1:n) |>
  pivot_longer(
    cols=contains("P.1"),
    names_to = "item",
    values_to = "pr_bf"
  ) |>
  mutate(item= str_remove(item, ".P.1"))

probs = left_join(p_u,p_bf)


#all general discrim params divided by all general + specific
weight=c()
wt_general= sum(coef[,1])/ (sum(colSums(coef)[1:ndim+1]) ) 
weight=c(weight, wt_general)
for (i in 1:ndim) {
#wt_dim =sum(coef[((bound*i)-bound+1):(bound*i),i])/ (sum(colSums(coef)[1:ndim]) ) #only general FSM params
wt_dim =sum(coef[((bound*i)-bound+1):(bound*i),i+1])/ (sum(coef[((bound*i)-bound+1):(bound*i),i+1]) + sum(coef[((bound*i)-bound+1):(bound*i),1])) #only general FSM params
weight=c(weight, wt_dim)
}
th_bf_overall=weight[1]*th_bf[,1]   
# wt_dim1= sum(coef[1:25,2])/ (sum(coef[1:25,2]) +sum(coef[1:25,1]))
# wt_dim2= sum(coef[26:50,3])/ (sum(coef[26:50,3]) +sum(coef[26:50,1]))
# wt_dim3= sum(coef[51:75,4])/ (sum(coef[51:75,4]) +sum(coef[51:75,1]))



for (i in 2:ndim) {
  th_bf_overall= th_bf_overall + weight[i]*th_bf[,i]    
}
cor(th_u, rowSums(theta))
cor(th_u, theta[,1])
cor(th_bf_overall, theta[,1])
cor(th_bf[,1], theta[,1])
cor(th_bf_overall, rowSums(theta)) #bifactor score vs true theta
cor(th_u, th_bf_overall)
cor_uni=c()
cor_bf=c()
cor_bf_uni=c()

cor_bf =c(cor_bf, cor_bf_gen= cor(th_bf[,1], theta[,1]))
cor_uni =c(cor_uni, cor_uni_gen= cor(th_u, theta[,1]))
cor_bf_uni =c(cor_bf_uni, cor_gen= cor(th_u, th_bf[,1]))
## Dimesnions
metrics=c()
for (i in 1:ndim) {
  cor_uni=c(cor_uni,cor((theta[,(i+1)]), th_dims[,i]))
  metrics=c(metrics, paste0("cor_true_uni_dim",i)) 
  print(cor((theta[,(i+1)]), th_dims[,i]))
}
names(cor_uni)[2:(2+ndim-1)]= metrics

metrics=c()
for (i in 1:ndim) {
  cor_bf_uni=c(cor_bf_uni, cor((th_bf[,1]+th_bf[,i+1]), th_dims[,i]))
  metrics=c(metrics, paste0("cor_bf_uni_dim",i)) 
  print(cor((th_bf[,1]+th_bf[,i+1]), th_dims[,i]))
}
names(cor_bf_uni)[2:(2+ndim-1)] = metrics

metrics=c()
for (i in 1:ndim) {
  cor_bf=c(cor_bf, cor(th_bf[,i+1], theta[,i+1]))
  metrics=c(metrics, paste0("cor_true_bf_dim",i)) 
  print(cor(th_bf[,1]+th_bf[,i+1], (theta[,(i+1)]+ theta[,1])))
}
names(cor_bf)[2:(2+ndim-1)] = metrics

out= rbind(cbind(cor_uni,cor_bf), cbind(se_unidim, se_bf))

return(out)
}

results=data.frame()

ndim=3
sigma1=diag(ndim+1)
sigma2= matrix(c(1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1), nrow=4, ncol=4, byrow=T)
## Generate data with correlation between factors
# sigma2= diag(ndim+1)
# totaldim=ndim+1
# for (i in 1:totaldim) {
#   for (j in 1:totaldim) {
#     if(i!=j) {sigma2[i,j]= 0.5}
#   }
# }
# datasim2 = simdata(N=n, a=a, d=diff, itemtype = 'dich', Theta=theta,
#                    sigma= sigma2)


out=simbf(n=1000, nitems=75, ndim=3, sigma=sigma1,aset="bf")
out2=simbf(n=1000, nitems=75, ndim=3, sigma=sigma2,aset="bf")


results=bind_rows(results, 
              c(n=1000, nitems=75,ndim=3,dim_cor=0,out[[1]], out[[2]]))

results=bind_rows(results, 
                  c(n=200, nitems=10,ndim=3,dim_cor=0,out[[1]], out[[2]]))


results=bind_rows(results, 
                  c(
                    n=100, nitems=30,ndim=3,dim_cor=0,
                    simbf(n=100, nitems=75, ndim=3, sigma=diag(ndim+1), aset="bf"))
)

results=bind_rows(results, 
                  c(
                    n=1000, nitems=75,ndim=3,dim_cor=0,
                    simbf(n=1000, nitems=75, ndim=3, sigma=diag(ndim+1), aset="not bf"))
)

results=bind_rows(results, 
                  c(
                    n=100, nitems=30,ndim=3,dim_cor=0,
                    simbf(n=100, nitems=30, ndim=3, sigma=diag(ndim+1), aset="not bf"))
)



results=bind_rows(results, 
                  c(
                    n=100, nitems=30,ndim=3,dim_cor=0,
                    simbf(n=100, nitems=30, ndim=3, sigma=diag(ndim+1), aset="bf"))
)





## Simulate using mirt
### Item paramters