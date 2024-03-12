## Simulation
# simulate multidimensional theta
library(MASS)
library(mirt)
nstudents=c(1000,5000)
nitems_list= c(50,75,100)

n=nstudents[1]
nitems=nitems_list[2]
ndim=3 # number of specific factors
bound=round(nitems/ndim)
item_mapping=c()
for (i in 1:(ndim-1)) {
  item_mapping=c(item_mapping, rep(i,bound))
}
item_mapping=c(item_mapping, rep(ndim,(nitems-(bound*i))))
## Generate item discriminations
a= matrix(data=NA, nrow=nitems, ncol=ndim+1)
a[,1]= rep(1, nitems)
for (i in 2:(ndim+1)) {
  a[,i]= ifelse(item_mapping==i-1,1,0)
}


a[,2]= rep(1,75)
a[,3]= rep(1,75)


diff= matrix(rnorm(nitems))
sigma=diag(ndim+1)
## 
# sigma[2,3] =0.8
# sigma[3,2] =0.8
# sigma[3,4] =0.8
# sigma[4,3] =0.8
# sigma[4,2] =0.8
# sigma[2,4] =0.8
theta= mvrnorm(n, mu=rep(0,ndim+1), Sigma=sigma)

#prob = exp()

datasim = simdata(a=a, d=diff, N=n, itemtype = 'dich', Theta=theta,
             sigma= diag(ndim+1))

### Factor analysis
pa = factanal(datasim, factors=3, rotation="varimax")
print(pa, digits=2, cutoff=0.3, sort=TRUE)
dim(pa$loadings)
loads= round(pa$loadings[ 1:74,], 3)
loads= data.frame(loads) |>
  mutate(across(starts_with("Factor"), ~ifelse(.<0.3, NA, .)))
loads$item= rownames(loads)

### 3 dimensional model as generated
model_unidim= mirt(datasim,1,  itemtype = 'Rasch')
#m2_u= M2(model_unidim)


### Unidimensional model
th_u= fscores(model_unidim)
th_u_se= fscores(model_unidim,full.scores.SE=TRUE)
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
  #th_dim1_se= fscores(m_dim1,full.scores.SE=TRUE)
}


### Bifactor
specific= item_mapping

bfactormod <- bfactor(datasim, model=specific, 
                      quadpts = 9, TOL = 1e-3) #itemtype default is 2PL
coef = coef(bfactormod, simplify=T)$items
th_bf= fscores(bfactormod,  QMC=TRUE)
th_bf_se= fscores(bfactormod,  QMC=TRUE,full.scores.SE=T)
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

# wt_dim1= sum(coef[1:25,2])/ (sum(coef[1:25,2]) +sum(coef[1:25,1]))
# wt_dim2= sum(coef[26:50,3])/ (sum(coef[26:50,3]) +sum(coef[26:50,1]))
# wt_dim3= sum(coef[51:75,4])/ (sum(coef[51:75,4]) +sum(coef[51:75,1]))


th_bf_overall=weight[1]*th_bf[,1]   
for (i in 2:ndim) {
  th_bf_overall= th_bf_overall + weight[i]*th_bf[,i]    
}
cor(th_u, rowSums(theta))
cor(th_u, theta[,1])
cor(th_bf_overall, theta[,1])
cor(th_bf_overall, rowSums(theta)) #bifactor score vs true theta
cor(th_u, rowSums(theta)) #bifactor score vs true theta
cor(th_u, th_bf_overall)

## Dimesnions
for (i in 1:ndim) {
print(cor((theta[,(i+1)]+ theta[,1]), th_dims[,i]))
}

for (i in 1:ndim) {
  print(cor((th_bf[,1]+th_bf[,i+1]), th_dims[,i]))
}

for (i in 1:ndim) {
  print(cor(th_bf[,1]+th_bf[,i+1], (theta[,(i+1)]+ theta[,1])))
}







## Simulate using mirt
### Item paramters