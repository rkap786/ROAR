library(readr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(mirt)
library(ggmirt)
library(psych)
library(Matrix)


setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
source("PA/ROAR/helper fncs.R")
df= read_csv("PA data/roar_pa_trialdata.csv")
df = df[,-1] |> 
  mutate(removed = ifelse(is.na(removed), 0, removed)) |>
  rename(item= itemId) 

## keep respondents who have answered all three subscales
df_resp = df |>  getrespmatrix() |> drop_na()
n_distinct(df_resp$subj)
resp_list= df_resp |> select(subj) |> unique() |> as.vector()
resp_list= resp_list[[1]]

df= df |> filter(subj %in% resp_list)

## keep respondents who have answered all three subscales
resp_list= df_resp |> select(subj) |> unique() |> as.vector()
resp_list= resp_list[[1]]
df= df |> filter(subj %in% resp_list)

## Get datasets for each subscale
resp_del = df    |> filter(block=="DEL") |>  getrespmatrix()
resp_fsm = df |> filter(block=="FSM")  |>  getrespmatrix()
resp_lsm = df  |> filter(block=="LSM") |>  getrespmatrix()
n_distinct(resp_del$subj)
n_distinct(resp_fsm$subj)
n_distinct(resp_lsm$subj)

####################################
######### Unidimensional IRT

## 1 factor IRT
model_unidim= mirt(df_resp[,-1],1,  itemtype = 'Rasch', guess=0.3)
#m2_u= M2(model_unidim)
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

resp_mean = bind_cols(item=names(df_resp)[-1],
                      pr0=as.numeric(colMeans(df_resp[,-1])))
resp = df |> dplyr::select(subj,item,correct)
p_u = left_join(p_u,resp_mean)
p_u = left_join(p_u,resp)

#probs = left_join(p_m,resp_mean)
probs = left_join(resp_mean,p_u)

#imv.binary(probs$correct,probs$pr0,probs$pr_u)
### Unidimesnionalfor subscales
m_fsm= mirt(resp_fsm[-1],1,  itemtype = 'Rasch', guess=0.3)
th_fsm= fscores(m_fsm)

#resp_del = df_resp_DEL |> dplyr::filter(item!="DEL_3") |> getrespmatrix()
m_del= mirt(resp_del[-1],1,  itemtype = 'Rasch', guess=0.3)
th_del= fscores(m_del)

#resp_lsm = df_resp_LSM |>  getrespmatrix()
m_lsm= mirt(resp_lsm[-1],1,  itemtype = 'Rasch', guess=0.3)
th_lsm= fscores(m_lsm)

################## Bifactor
## First set is FSM, second set is LSM, third set is DEL
specific= c(rep(1,25),rep(2,25), rep(3,24))
guess=rep(0.3,74)


bfactormod <- bfactor(df_resp[-1], model=specific, 
                      quadpts = 9, TOL = 1e-3,
                      guess=guess)
coef = coef(bfactormod, simplify=T)$items
th_bf= fscores(bfactormod,  QMC=TRUE)
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


head(th_bf)
summary(bfactormod)

#Model fit - bifactor fits better
anova(model_unidim, bfactormod)
imv.binary(probs$correct,probs$pr_u,probs$pr_bf)




## Compare ability scores
### Method 1: take general score as general score, and specific as subscale score
### Remember that specific factor 1 is FSM, specific factor 2 is LSM, 
######### specific factor 3 is DEL

## Compare bifactor and unidimensional models
cor(probs$pr_u,probs$pr_bf) #0.87
cor(th_bf[,1],th_u) #0.92


cor(th_fsm,th_bf[,2]) #0.68
cor(th_lsm,th_bf[,3]) #0.86
cor(th_del,th_bf[,4]) #0.72


### Method 2: take general score as general score + sum of all specific scores, and specific as subscale score
th_bf_overall= rowSums(th_bf)
cor(th_u,th_bf_overall) # about the same, 0.94


th_fsm2= th_bf[,2]+th_bf[,1]
th_lsm2= th_bf[,3]+th_bf[,1]
th_del2= th_bf[,4]+th_bf[,1]
cor(th_fsm2,th_fsm) #0.973
cor(th_lsm2,th_lsm) #0.971
cor(th_del2,th_del) #0.961


### Method 3: weight scores using discrimination
## Overall score weights
wt_general= sum(coef_2PL[,1])/ (sum(colSums(coef_2PL)[1:4]) )
wt_fsm= sum(coef_2PL[1:25,1])/ (sum(colSums(coef_2PL)[1:4]) )
wt_lsm= sum(coef_2PL[26:50,1])/ (sum(colSums(coef_2PL)[1:4]) )
wt_del= sum(coef_2PL[51:74,1])/ (sum(colSums(coef_2PL)[1:4]) )

th_bf_overall2= wt_general*th_bf[,1] + wt_fsm*th_bf[,2] + wt_lsm*th_bf[,3] +
  wt_del*th_bf[,4]

cor(th_u,th_bf_overall2) # 0.98

## Specific domain weights
wt_fsm_sp= sum(coef_2PL[1:25,2])/ (sum(coef_2PL[1:25,2]) +sum(coef_2PL[1:25,1]))
wt_lsm_sp= sum(coef_2PL[26:50,3])/ (sum(coef_2PL[1:25,3]) +sum(coef_2PL[1:25,1]))
wt_del_sp= sum(coef_2PL[51:74,4])/ (sum(coef_2PL[1:25,4]) +sum(coef_2PL[1:25,1]))

th_fsm4= wt_fsm_sp*th_bf[,2]+wt_general*th_bf[,1]
th_lsm4= wt_lsm_sp*th_bf[,3]+wt_general*th_bf[,1]
th_del4= wt_del_sp*th_bf[,4]+wt_general*th_bf[,1]

cor(th_bf_overall2,th_u) # no change

cor(th_fsm4, th_fsm) #0.971
cor(th_lsm4,th_lsm) #0.98
cor(th_del4,th_del) # 0.96





