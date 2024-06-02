## Simulate data
## n items = 10, 20, 30, 40
## Generate response data, 
## items ordered either based on fisher information or random

setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
source("PA/ROAR/cat.functions.R")
library(readxl)
library(dplyr)
library(catR)

## SEM and Square error compared


### ROAR overall theta estimates

getrespmatrix= function(data) {
  data = data |>
    dplyr::select(subj,item,correct) |>
    pivot_wider(values_from= correct,
                names_from = item)
  
  #names(data) = paste0("item",names(data))
  #names(data)[1] = "personID"
  data
}

ROAR_scores= read_csv("PA data/roar_pa_trialdata.csv")
ROAR_scores = ROAR_scores[,-1] |> 
  mutate(removed = ifelse(is.na(removed), 0, removed)) 

ROAR_subj_ids= unique(ROAR_scores$subj)

CTOPP_scores <- read_excel("PA data/CTOPP_scores.xlsx")
head(CTOPP_scores)
names(CTOPP_scores)
CTOPP_scores= CTOPP_scores |> dplyr::select(Subject, CTOPP_PA_raw) |> 
  dplyr::rename(subj=Subject) |> 
  filter(subj %in% ROAR_subj_ids) 

CTOPP_subj_ids= unique(CTOPP_scores$subj)
ROAR_scores= ROAR_scores |> filter(subj %in% CTOPP_subj_ids)

## Same unique respondent set, N=251
n_distinct(CTOPP_scores$subj)
n_distinct(ROAR_scores$subj)

##Unidim model
ROAR_scores=ROAR_scores |> dplyr::rename(item=itemId)
resp_roar = getrespmatrix(ROAR_scores )
model_unidim= mirt(resp_roar[,-1],1,  itemtype = 'Rasch', guess=0.3)
th_u_roar= fscores(model_unidim)
th_u_roar= cbind(th_u_roar,resp_roar[,1])
names(th_u_roar)=c("th","subj")
th_u_roar_se= fscores(model_unidim,full.scores.SE=TRUE)
empirical_rxx(as.matrix(tibble(F1 = th_u_roar_se[,1], SE_F1 = th_u_roar_se[,2])))

coef = coef(model_unidim, simplify=T)$items
coef[,2]= -coef[,2]
coef=data.frame(coef)
coef$block = rownames(coef)


item.bank.overall=list()
#del, fsm, lsm
item.bank.overall[['DEL']]=coef |> filter(str_detect(block,"DEL") ) |> dplyr::select(-block)
item.bank.overall[['FSM']]=coef |> filter(str_detect(block,"FSM") ) |> dplyr::select(-block)
item.bank.overall[['LSM']]=coef |> filter(str_detect(block,"LSM") ) |> dplyr::select(-block)


th_merged= th_u_roar
th_merged= merge(th_merged, CTOPP_scores, by="subj")

cor(th_merged$CTOPP_PA_raw, th_merged$th)
## lets use raw CTOPP, cor=0.73

## Sum score correlation
resp_roar_ss= rowSums(resp_roar[,-1], na.rm=T)
resp_roar_ss= cbind(ss=as.numeric(resp_roar_ss),subj=resp_roar[,1])
resp_roar_ss= merge(resp_roar_ss,CTOPP_scores, by="subj")
cor(resp_roar_ss$ss, resp_roar_ss$CTOPP_PA_raw)


### Separate unidimensional models for each subscale
resp_del = ROAR_scores |> filter(block=="DEL") |>  getrespmatrix()
resp_fsm = ROAR_scores |> filter(block=="FSM")  |>  getrespmatrix()
resp_lsm = ROAR_scores  |> filter(block=="LSM") |>  getrespmatrix()
subj.list= unique(resp_del$subj)

data.dim= list()
data.dim[['DEL']] = resp_del
data.dim[['FSM']] = resp_fsm
data.dim[['LSM']] = resp_lsm

## Item difficulties
d1=colMeans(resp_del[,-1], na.rm=T)
d2=colMeans(resp_lsm[,-1], na.rm=T)
d3=colMeans(resp_fsm[,-1], na.rm=T)
plot(density(d1))
lines(density(d2),col="red")
lines(density(d3),col="blue")


###### Subscales IRT calibration
th_dims=list()
coef_dims=list()
rel=c()
item.bank.overall=list()
for (i in c('DEL', 'FSM', 'LSM')) {
  data=data.dim[[i]]
  data = data |> filter(subj %in% subj.list)
  m_dim= mirt(data[,-1],1,  itemtype = 'Rasch')
  th= fscores(m_dim)
  th_se=fscores(m_dim, full.scores.SE=TRUE)
  rel=c(rel,empirical_rxx(as.matrix(tibble(F1 = th_se[,1], SE_F1 = th_se[,2]))))
  th_dims[[i]]=data.frame(th)
  coef_dims[[i]] = coef(m_dim, simplify=T)$items
  coef_dims[[i]][,2]= -coef_dims[[i]][,2]
  coef_dims[[i]]= data.frame(coef_dims[[i]])
  coef_dims[[i]]$block = rownames(coef_dims[[i]])
  item.bank.overall[[i]] = coef_dims[[i]]
  #del, fsm, lsm
  # th_dim_se= fscores(m_dim,full.scores.SE=TRUE)
  # se_unidim=c(se_unidim,unidim_subscale=as.numeric(empirical_rxx(th_dim_se)))
}
th_dims= dplyr::bind_cols(th_dims) 
names(th_dims)= c("th_del","th_fsm","th_lsm")
subj_final=data$subj
th_dims$pid=data$subj


### Simulate as follows: Run separately for each dimension
### Dimension 1: Run CAT using unidimensional IRT for dim 1
#item.bank.u=coef_dims[[1]]
np= n_distinct(subj_final)
th_start= rep(0,np)
results.compare=list()
th_u_final=list()


for (i in c('DEL', 'LSM', 'FSM')) {
  data=data.dim[[i]]
  data = data |> filter(subj %in% subj.list)
  item.bank= coef_dims[[i]]
  coef= item.bank.overall[[i]]
  df = data |> filter(subj %in% subj.list) |> select(-subj)
  ni= ncol(df)
  results= func.catSim(resp=df, item.bank= item.bank,
                       item.bank.u=coef, method='MFI', theta=th_start, 
                       stoplen=ni,skipitems=1, pids=pids=data$subj)
  ## First element in list returns CAT results; 
  ### second is final unidim theta estimates for overall ability
  
  results_random= func.catSim(resp=df, item.bank= item.bank,
                       item.bank.u=coef, method='random', theta=th_start, 
                       stoplen=ni,skipitems=1,pid=pids=data$subj)
  
  results.compare[[i]]= bind_rows(list(results[[1]], results_random[[1]]),.id="type")
  results.compare[[i]]= results.compare[[i]] |> 
    mutate(type=if_else(type==1,"CAT","Random"))
  # results.compare[[i]]= results.compare[[i]] |> 
  #   mutate(stoplen= rep(len,2*np*ni))
  
  ## Add in true values
  results.compare[[i]]= left_join(results.compare[[i]],th_dims[,c(4,i)])
  
  
  ## Results from unidimensional
  th_u_final[[i]]=bind_rows(results[[2]],results_random[[2]],.id="type")
  th_u_final[[i]]= th_u_final[[i]] |> mutate(type=if_else(type==1,"CAT","Random"))
  th_u_final[[i]]= left_join(th_u_final[[i]],th_u_roar)
  
  # Unidimensional model to give starting parameter for next dimension
  th_start= results[[2]]$thetaU.final
  
}

#### Overall theta reporting


#th_u_final= do.call("rbind",th_u_final)
th_u_df= bind_rows(th_u_final, .id = "dim")
th_dim_merged=bind_rows(results.compare, .id="dim") |>
  group_by(dim) |>
  mutate(maxtrial= max(trialNumTotal)) |>
  ungroup() |>
  filter(type=="CAT", trialNumTotal==maxtrial) |>
  dplyr::select(dim, pid, thetaFinal)  |> 
  rename(subj=pid) |>
  pivot_wider(names_from = dim, 
              values_from = thetaFinal)

th_u_df |> filter(type=="CAT") |> dplyr::select(th,thetaU.final) |> cor()
th_u_df_filter = th_u_df |> filter(type=="CAT", dim==3) |> rename(subj=pid)
th_merged= left_join(th_merged,th_u_df_filter)
th_merged= left_join(th_merged,th_dim_merged)
th_merged$th_total=th_merged$thetaU.final +th_merged$DEL + th_merged$LSM + th_merged$FSM


th_merged  |> dplyr::select(th,CTOPP_PA_raw) |> cor()
th_merged  |> dplyr::select(thetaU.final,CTOPP_PA_raw) |> filter(!is.na(thetaU.final)) |> cor()
th_merged  |> dplyr::select(thetaU.final,th) |> filter(!is.na(thetaU.final)) |> cor() #Full unidimensional vs CAT theta
th_merged  |> dplyr::select(th_total,th) |> filter(!is.na(th_total)) |> cor()
th_merged  |> dplyr::select(th_total,CTOPP_PA_raw) |> filter(!is.na(th_total)) |> cor()

##################################################################################################
### Calculate metrics
### Hardcoded for 3 dimensions

quartiles1 <- quantile(results.compare[[1]]$th1,
                       probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
quartiles2 <- quantile(results.compare[[2]]$th2, 
                       probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
quartiles3 <- quantile(results.compare[[3]]$th3, 
                       probs = c(0, 0.20, 0.40, 0.60, 1))

#data.dim= list(resp_del,resp_lsm,resp_fsm)
metrics.dim1= results.compare[[1]] |>
  mutate(theta.bin = cut(th1, 
                         breaks = quartiles1, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(th1, thetaEstimate), 
                   bias = Metrics::bias(th1, thetaEstimate)) 


metrics.dim2= results.compare[[2]] |>
  mutate(theta.bin = cut(th2, 
                         breaks = quartiles2, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(th2, thetaEstimate), 
                   bias = Metrics::bias(th2, thetaEstimate)) 


metrics.dim3= results.compare[[3]] |>
  mutate(theta.bin = cut(th3, 
                         breaks = quartiles3, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(th3, thetaEstimate), 
                   bias = Metrics::bias(th3, thetaEstimate)) 

### Plot estimated ability

results.compare[[1]] |>
  filter(pid==946) |>
  select(pid, type, trialNumTotal, thetaEstimate, th1) |> 
  ggplot(aes(x=trialNumTotal, y=thetaEstimate, group=type, color=type)) +
  geom_point() +
  geom_hline(aes(yintercept= unique(th1)))




### Plot metrics
metrics.dim1 |>
  ggplot(aes(x=trialNumTotal, y=reliability, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Reliability",
       title = "Reliability, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")


metrics.dim1 |>
  ggplot(aes(x=trialNumTotal, y=bias, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Bias",
       title = "Bias, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")

metrics.dim1 |>
  ggplot(aes(x=trialNumTotal, y=mse, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "MSE",
       title = "MSE, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")

##################################################################################################
##################################################################################################

# Real data


metrics.dim2 |>
  ggplot(aes(x=trialNumTotal, y=reliability, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Reliability",
       title = "Reliability, LSM") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")


metrics.dim2 |>
  ggplot(aes(x=trialNumTotal, y=bias, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Bias",
       title = "Bias, LSM") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")

metrics.dim2 |>
  ggplot(aes(x=trialNumTotal, y=mse, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "MSE",
       title = "MSE, LSM") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")


##################################################################################################
##################################################################################################

# Real data


metrics.dim3 |>
  ggplot(aes(x=trialNumTotal, y=reliability, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Mean squared error",
       title = "Reliability, FSM") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")


metrics.dim3 |>
  ggplot(aes(x=trialNumTotal, y=bias, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Bias",
       title = "Bias, FSM") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")

metrics.dim3 |>
  ggplot(aes(x=trialNumTotal, y=mse, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "MSE",
       title = "MSE, FSM") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")


#############################################################################
#############################################################################
#############################################################################

### Overall ability

for (l in c(5,10,15,20,24)) {
  
}




th_u_final= bind_rows(results[[2]],results_random[[2]], .id="type")
th_u_final= left_join(th_u_final, CTOPP_scores)
th_u_final= th_u_final |> mutate(type=if_else(type==1, "CAT", "Random"))

th_u_final |> filter(type=="CAT") |> dplyr::select(thetaU.final,CTOPP_PA_raw) |> cor()
th_u_final |> filter(type=="Random") |> dplyr::select(thetaU.final,CTOPP_PA_raw) |> cor()



#item.bank.u=coef_dims[[1]]
np= n_distinct(subj_final)
th_start= rep(0,np)
results.compare=list()
th_u_final=data.frame()
th_u_final_r=data.frame()
results.compare= list()


for (i in 1:3) {
  data=data.dim[[i]]
  item.bank= coef_dims[[i]]
  item.bank[,2]= -item.bank[,2]
  df = data |> dplyr::filter(subj %in% subj.list) |> select(-subj)
  ni= ncol(df)
  for (stop in c(5,10,15,20,ni)) {
   results= func.catSim(resp=df, item.bank= item.bank,
                       item.bank.u=item.bank.overall[[i]], method='MFI', 
                       theta=th_start, 
                       stoplen=stop,skipitems=1, pid=subj_final)
  
  results_random= func.catSim(resp=df, item.bank= item.bank,
                              item.bank.u=item.bank.u[[i]], method='random', 
                              theta=th_start, 
                              stoplen=stop,skipitems=1,pid=subj_final)
  
  results.compare[[i]]= bind_rows(list(results[[1]], results_random[[1]]),.id="type")
  results.compare[[i]]= results.compare[[i]] |>
    mutate(type=if_else(type==1,"CAT","Random"))
  results.compare[[i]]= results.compare[[i]] |> filter(trialNumTotal==stop) |> 
    distinct() |>
    mutate(stoplen= rep(stop,2*np), dim=rep(i, 2*np))
  
  ## Add in true values
  
  ## Results from unidimensional
  results1.th.u= cbind(results[[2]], test.length=rep(stop,np), dim=i)
  results1.th.u.r= cbind(results_random[[2]], test.length=rep(stop,np),dim=i)
  #results.compare[[i]]= left_join(results.compare[[i]],th_dims[,c(4,i)])
  
  #results.full= left_join(results.full,results1.th.u)
  
  # Unidimensional model to give starting parameter for next dimension
  th_start= results[[2]]$thetaU.final
  
  th_u_final=rbind(th_u_final,results1.th.u)
  th_u_final_r=rbind(th_u_final_r,results1.th.u.r)
  }
}

th_uni_dim= bind_rows(th_u_final,th_u_final_r, .id="type")
th_uni_dim= th_uni_dim |> mutate(type=ifelse(type==1, "CAT", "Random"))
th_uni_dim=  left_join(th_uni_dim, th_merged)
th_uni_dim=  left_join(th_uni_dim, th_dims)

results_plot= th_uni_dim |> filter(dim==3) 

quartiles <- quantile(results_plot$th, 
                      probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))
#data.dim= list(resp_del,resp_lsm,resp_fsm)
metrics.dim1= results_plot |>
  mutate(theta.bin = cut(th, 
                         breaks = quartiles, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(test.length,type) |>
  dplyr::summarise(sem = mean(thetaU.se.final), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaU.final, SE_F1 = thetaU.se.final))),
                   mse = Metrics::mse(th, thetaU.se.final), 
                   bias = Metrics::bias(th, thetaU.se.final)) 


metrics.dim1 |>
  ggplot(aes(x=test.length, y=mse, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "MSE",
       title = "MSE, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")


metrics.dim1 |>
  ggplot(aes(x=test.length, y=reliability, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Reliability",
       title = "MSE, Overall") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")

# Correlations
th_uni_dim |> filter(dim==3,test.length=20) |> cor(CTOPP_PA_raw,thetaU.final)
