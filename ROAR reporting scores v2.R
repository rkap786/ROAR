library(readr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(mirt)
library(ggmirt)
library(psych)
library(Matrix)

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

setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
df= read_csv("PA data/roar_pa_trialdata.csv")
df = df[,-1]
df = df |> mutate(removed = ifelse(is.na(removed), 0, removed)) |>
  rename(item= itemId)




getrespmatrix= function(data) {
  data = data |>
    dplyr::select(subj,item,correct) |>
    pivot_wider(values_from= correct,
                names_from = item)
  
  #names(data) = paste0("item",names(data))
  #names(data)[1] = "personID"
  data
}



table(df$block)
n_distinct(df$subj)
df_resp = df |>  getrespmatrix() |> drop_na()

## keep respondents who have answered all three subscales
resp_list= df_resp |> select(subj) |> unique() |> as.vector()
resp_list= resp_list[[1]]

df= df |> filter(subj %in% resp_list)


#Factor analysis
pa = factanal(df_resp[,-1], factors=3, rotation="varimax")
print(pa, digits=2, cutoff=0.3, sort=TRUE)
dim(pa$loadings)
loads= round(pa$loadings[ 1:74,], 3)
loads= data.frame(loads) |>
  mutate(across(starts_with("Factor"), ~ifelse(.<0.3, NA, .)))
loads$item= rownames(loads)
write_csv(loads, "PA/factor_analysis_results.csv")

## MIRT model
cov=matrix(c(1,0.3,0.3,0,1,0.3,0.3,0.3,1),nrow=3,ncol=3)
cov=tril(cov)
model= mirt.model('F1= 1-25,  
                      F2 = 26-50,
                      F3= 51-74,
                      COV=cov')
model_multidim = mirt(df_resp[,-1], 
                      model=model,
                      itemtype = 'Rasch', guess=0.3, method='QMCEM')
coef=coef(model_multidim,simplify=T,IRTpars = TRUE)$items
# eas = data.frame(
#   item=rownames(coef),
#   easy=coef[,4],
#   n= 1:nrow(eas))

summary(model_multidim)
m2_m= M2(model_multidim)
th_m= fscores(model_multidim, full.scores = T)

# Predicted probability from MIRT, where each column is an item
p_m=data.frame(probtrace(model_multidim, th_m)) |>
  dplyr::select(contains("P.1")) |>
  mutate(subj = df_resp$subj) |>
  pivot_longer(
    cols=contains("P.1"),
    names_to = "item",
    values_to = "pr_m"
  ) |>
  mutate(item= str_remove(item, ".P.1"))

## 1 factor IRT
model_unidim= mirt(df_resp[,-1],1,  itemtype = 'Rasch', guess=0.3)
m2_u= M2(model_unidim)
th_u= fscores(model_unidim)
p_u=data.frame(probtrace(model_unidim, th_u)) |>
  dplyr::select(contains("P.1")) |>
  mutate(subj = df_resp$subj) |>
  pivot_longer(
    cols=contains("P.1"),
    names_to = "item",
    values_to = "pr_u"
  ) |>
  mutate(item= str_remove(item, ".P.1"))

## IMV for multidimensional model, compared to mean correct
resp_mean = bind_cols(item=names(df_resp)[-1],
                  pr0=as.numeric(colMeans(df_resp[,-1])))
resp = df |> dplyr::select(subj,item,correct)
p_u = left_join(p_u,resp_mean)
p_u = left_join(p_u,resp)

probs = left_join(p_m,resp_mean)
probs = left_join(probs,p_u)

imv.binary(probs$correct,probs$pr0,probs$pr_u)
imv.binary(probs$correct,probs$pr0,probs$pr_m)
imv.binary(probs$correct,probs$pr_u,probs$pr_m)


(rbind("3F"=m2_m, "1F"=m2_u) |> t() |>round(4))[c(1,4,8,9),]

### Separate unidimensional models for each subscale

resp_del = df    |> filter(block=="DEL") |>  getrespmatrix()
resp_fsm = df |> filter(block=="FSM")  |>  getrespmatrix()
resp_lsm = df  |> filter(block=="LSM") |>  getrespmatrix()


#Estimate model - FSM
#resp_fsm = df_resp_FSM |> getrespmatrix()
m_fsm= mirt(resp_fsm[-1],1,  itemtype = 'Rasch', guess=0.3)
th_fsm= fscores(m_fsm)

#resp_del = df_resp_DEL |> dplyr::filter(item!="DEL_3") |> getrespmatrix()
m_del= mirt(resp_del[-1],1,  itemtype = 'Rasch', guess=0.3)
th_del= fscores(m_del)

#resp_lsm = df_resp_LSM |>  getrespmatrix()
m_lsm= mirt(resp_lsm[-1],1,  itemtype = 'Rasch', guess=0.3)
th_lsm= fscores(m_lsm)


  
################ Bifactor models
names(df_resp[-1])
loads
## First set is FSM, second set is LSM, third set is DEL
specific= c(rep(1,25),rep(2,25), rep(3,24))
guess=rep(0.3,74)

bfactormod_Rasch <- bfactor(df_resp[-1], model=specific, itemtype='Rasch',
                          quadpts = 9, TOL = 1e-3,
                          guess=guess)

bfactormod_2pl <- bfactor(df_resp[-1], model=specific, itemtype='2PL',
                      quadpts = 9, TOL = 1e-3,
                      guess=guess)

bfactormod <- bfactor(df_resp[-1], model=specific, 
                          quadpts = 9, TOL = 1e-3,
                          guess=guess)

## Compare loadings 

summary(bfactormod)
cbind(summary(bfactormod)[[1]], loads)
spec1_load=summary(bfactormod)[[1]][,2]
fa_2= loads[,2]
## Correlation between factor analysis domain and bifactor domain
cor(spec1_load, fa_2, use="pairwise.complete.obs")

coef_rasch= coef(bfactormod_Rasch, simplify=T)
coef_2PL = coef(bfactormod_2pl, simplify=T)$items
coef = coef(bfactormod, simplify=T)



th_bf= fscores(bfactormod,  QMC=TRUE)
head(th_bf)
th_bf_Rasch= fscores(bfactormod_Rasch,  QMC=TRUE)

summary(bfactormod)

## Compare ability scores
### Method 1: take general score as general score, and specific as subscale score
### Remember that specific factor 1 is FSM, specific factor 2 is LSM, 
######### specific factor 3 is DEL

cor(th_u,th_bf[,1])
cor(th_fsm,th_bf[,2])
cor(th_lsm,th_bf[,3])
cor(th_del,th_bf[,4])


### Method 2: take general score as general score + sum of all specific scores, and specific as subscale score
th_bf_overall= rowSums(th_bf)
cor(th_u,th_bf_overall) # about the same


th_fsm2= th_bf[,2]+th_bf[,1]
th_lsm2= th_bf[,3]+th_bf[,1]
th_del2= th_bf[,4]+th_bf[,1]
cor(th_fsm2,th_fsm)
cor(th_lsm2,th_lsm)
cor(th_del2,th_del)


### Method 3: weight scores using discrimination
## Overall score weights
wt_general= sum(coef_2PL[,1])/ (sum(colSums(coef_2PL)[1:4]) )
wt_fsm= sum(coef_2PL[1:25,1])/ (sum(colSums(coef_2PL)[1:4]) )
wt_lsm= sum(coef_2PL[26:50,1])/ (sum(colSums(coef_2PL)[1:4]) )
wt_del= sum(coef_2PL[51:74,1])/ (sum(colSums(coef_2PL)[1:4]) )

th_bf_overall2= wt_general*th_bf[,1] + wt_fsm*th_bf[,2] + wt_lsm*th_bf[,3] +
  wt_del*th_bf[,4]

cor(th_u,th_bf_overall2) # Pretty good 

## Specific domain weights
wt_fsm_sp= sum(coef_2PL[1:25,2])/ (sum(coef_2PL[1:25,2]) +sum(coef_2PL[1:25,1]))
wt_lsm_sp= sum(coef_2PL[26:50,3])/ (sum(coef_2PL[1:25,3]) +sum(coef_2PL[1:25,1]))
wt_del_sp= sum(coef_2PL[51:74,4])/ (sum(coef_2PL[1:25,4]) +sum(coef_2PL[1:25,1]))

th_fsm4= wt_fsm_sp*th_bf[,2]+wt_general*th_bf[,1]
th_lsm4= wt_lsm_sp*th_bf[,3]+wt_general*th_bf[,1]
th_del4= wt_del_sp*th_bf[,4]+wt_general*th_bf[,1]

cor(th_fsm4,th_fsm)
cor(th_lsm4,th_lsm)
cor(th_del4,th_del)



## Evaluate fit
#itemplot(bfactormod, 1, drop.zeros = TRUE)
#itemfit(bfactormod, QMC=TRUE)
M2(bfactormod, QMC=TRUE)
M2(bfactormod_Rasch, QMC=TRUE)
## Check fit stats
anova(model_unidim, bfactormod)
anova(bfactormod_2pl, bfactormod_Rasch) #very similar

#probs=probs |> dplyr::select(-pr_bf)

## Using 2PL bifactor

p_bf = data.frame(probtrace(bfactormod, th_bf)) |>
  dplyr::select(contains("P.1")) |>
  mutate(subj = df_resp$subj) |>
  pivot_longer(
    cols=contains("P.1"),
    names_to = "item",
    values_to = "pr_bf"
  ) |>
  mutate(item= str_remove(item, ".P.1"))

probs = left_join(probs,p_bf)

### Correlation between abilities
cor(probs$pr_u, probs$pr_bf)
cor(probs$pr_m, probs$pr_bf)


## Fit

imv.binary(probs$correct,probs$pr0,probs$pr_bf)
imv.binary(probs$correct,probs$pr_u,probs$pr_bf)
imv.binary(probs$correct,probs$pr_m,probs$pr_bf)

### Trying other loadings for bifactor
# 
# loads_max=data.frame(do.call(rbind,apply(loads,1, which.max)))
# loads_max$item= rownames(loads_max)
# specific2= data.frame(item=loads$item)
# specific2= left_join(specific2, loads_max, by="item")
# 
# ### Specific factors map onto Loadings 
# ### Does not work very well
# bfactormod2 <- bfactor(df_resp[-1], specific2[,2], guess=guess)
# M2(bfactormod2, QMC=TRUE)

# p_bf = data.frame(probtrace(bfactormod2, th_bf)) |>
#   dplyr::select(contains("P.1")) |>
#   mutate(subj = df_resp$subj) |>
#   pivot_longer(
#     cols=contains("P.1"),
#     names_to = "item",
#     values_to = "pr_bf2"
#   ) |>
#   mutate(item= str_remove(item, ".P.1"))
# 
# probs = left_join(probs,p_bf)

# imv.binary(probs$correct,probs$pr0,probs$pr_bf2)
# imv.binary(probs$correct,probs$pr_u,probs$pr_bf2)
# imv.binary(probs$correct,probs$pr_m,probs$pr_bf2)
# imv.binary(probs$correct,probs$pr_bf,probs$pr_bf2)
