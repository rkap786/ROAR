library(readr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(mirt)
library(ggmirt)
library(psych)
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

#Factor analysis
pa = factanal(df_resp[,-1], factors=3, rotation="varimax")
print(pa, digits=2, cutoff=0.3, sort=TRUE)
dim(pa$loadings)
loads= round(pa$loadings[ 1:74,], 3)
loads= data.frame(loads) |>
  mutate(across(starts_with("Factor"), ~ifelse(.<0.3, NA, .)))

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
eas = data.frame(
  item=rownames(coef),
  easy=coef[,4],
  n= 1:nrow(eas))

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

(rbind("3F"=m2_m, "1F"=m2_u) |> t() |>round(3))[c(1,4,8,9),]
  

