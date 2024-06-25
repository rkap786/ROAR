
setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
source("PA/ROAR/cat.functions.R")

library(readxl)
library(dplyr)
library(catR)
library(Metrics)
library(readr)

############################################################

## Import data
ROAR_scores= read_csv("PA data/roar_pa_trialdata.csv")
ROAR_scores = ROAR_scores[,-1] |> 
  mutate(removed = ifelse(is.na(removed), 0, removed)) 

ROAR_subj_ids= unique(ROAR_scores$subj)

CTOPP_scores <- read_excel("PA data/CTOPP_scores.xlsx")

### Filter to get the same respondent set for ROAR and CTOPP
CTOPP_scores= CTOPP_scores |> dplyr::select(Subject, CTOPP_PA_raw) |> 
  dplyr::rename(subj=Subject) |> 
  filter(subj %in% ROAR_subj_ids) 

CTOPP_subj_ids= unique(CTOPP_scores$subj)
ROAR_scores= ROAR_scores |> filter(subj %in% CTOPP_subj_ids)

## Same unique respondent set, N=251
n_distinct(CTOPP_scores$subj)
n_distinct(ROAR_scores$subj)


### ROAR overall theta estimates
##Unidim overll model
ROAR_scores=ROAR_scores |> dplyr::rename(item=itemId)
resp_roar = getrespmatrix(ROAR_scores )
model_unidim= mirt(resp_roar[,-1],1,  itemtype = 'Rasch', guess=0.3)
th_u_roar= fscores(model_unidim,method = "EAP",theta_lim = c(-4, 4))
th_u_roar= cbind(th_u_roar,resp_roar[,1])
names(th_u_roar)=c("th","subj")

## Overall item paramters
coef = coef(model_unidim, simplify=T)$items
coef[,2]= -coef[,2]
coef=data.frame(coef)
coef$itemId = rownames(coef)
#### Save item parameter file
write_csv(coef,"PA/Data/itemparam_overall_unidim.csv")


item.bank.overall=list()
#del, fsm, lsm
item.bank.overall[['DEL']]=coef |> filter(str_detect(itemId,"DEL") ) |> dplyr::select(-itemId)
item.bank.overall[['FSM']]=coef |> filter(str_detect(itemId,"FSM") ) |> dplyr::select(-itemId)
item.bank.overall[['LSM']]=coef |> filter(str_detect(itemId,"LSM") ) |> dplyr::select(-itemId)


th_merged= merge(th_u_roar, CTOPP_scores, by="subj")
cor(th_merged$CTOPP_PA_raw, th_merged$th)
## lets use raw CTOPP, cor=0.73



# results all in one place
result.merged= list()
for (subtask in c('DEL', 'LSM', 'FSM')) {
  
data = ROAR_scores |> filter(block==subtask) |>  filter(subj %in% subj.list) |>
  getrespmatrix()
subj.list= unique(resp_del$subj)
df = data |> select(-subj)

## Unidim model overall
m_dim= mirt(df,1,  itemtype = 'Rasch', guess= 0.3)
th= fscores(m_dim, method = "EAP",theta_lim = c(-4, 4))
th_se=fscores(m_dim, full.scores.SE=TRUE,method = "EAP",theta_lim = c(-4, 4))
rel=empirical_rxx(as.matrix(tibble(F1 = th_se[,1], SE_F1 = th_se[,2])))
th_dims=data.frame(th)
coef_dims = coef(m_dim, simplify=T)$items
coef_dims[,2]= -coef_dims[,2] #convert easiness to difficulty

coef_output= data.frame(coef_dims) |>
  mutate(itemId = rownames(coef_dims))
filename= paste0("PA/Data/irtparam_", subtask, ".csv")
write_csv(coef_output,filename)

names(th_dims)= c("th1")
subj_final=data$subj
th_dims$subj=data$subj

### Adaptive simulation
## Set up paramters
np= n_distinct(subj_final)
th_start= rep(0,np)
results.compare=data.frame()
th_u_final=data.frame()
item.bank= coef_dims
coef= item.bank.overall[[subtask]]
dim(coef)

ni= ncol(df)
results= func.catSim(resp=df, item.bank= item.bank,
                     item.bank.u=coef, method='MFI', theta=th_start, 
                     stoplen=ni,skipitems=1, pids=data$subj)
## This outputs a list with two elements
## First element in list returns CAT results; 
### Second element is final unidim theta estimates for overall ability

results_random= func.catSim(resp=df, item.bank= coef,
                            item.bank.u=coef, method='random', theta=th_start, 
                            stoplen=ni,skipitems=1,pids=data$subj)

results.compare= bind_rows(list(results[[1]], results_random[[1]]),.id="type")
results.compare= results.compare |> 
  mutate(type=if_else(type==1,"CAT","Random")) |>
  rename(subj=pid)

## Add in true values
results.compare= left_join(results.compare,th_dims)
result.merged[[subtask]]=results.compare

## This is more useful after stopping rule is applied
final.theta.unidim= bind_rows(list(results[[2]], results_random[[2]]),.id="type")
final.theta.unidim= final.theta.unidim |> 
  mutate(type=if_else(type==1,"CAT","Random")) |>
  rename(subj=pid)

}



##################################################################################################
### Calculate metrics
### Hardcoded for 3 dimensions

### First deletion
for (subtask in c('DEL', 'LSM', 'FSM')) {
  
results.compare= result.merged[[subtask]]
# quartiles1 <- quantile(results.compare$th1,
#                        probs = c(0, 0.3, 0.5, 0.7, 1))

metrics.dim1= results.compare  |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(th1, thetaEstimate), 
                   bias = Metrics::bias(th1, thetaEstimate)) 


### Plot metrics
plot.rel= metrics.dim1 |>
  ggplot(aes(x=trialNumTotal, y=reliability, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Reliability",
       title = "Reliability, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  ylim(0,1) +
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")
plotname= paste0("Plots/ROAR_plot.rel",".",subtask,".jpeg")
ggsave(plotname, plot=plot.rel)

plot.bias = metrics.dim1 |>
  ggplot(aes(x=trialNumTotal, y=bias, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "Bias",
       title = "Bias, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  ylim(-0.1,0.25) +
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")

plotname= paste0("Plots/ROAR_plot.bias",".",subtask,".jpeg")
ggsave(plotname, plot=plot.bias)

plot.mse = metrics.dim1 |>
  ggplot(aes(x=trialNumTotal, y=mse, group=type, color=type)) +
  geom_point() +
  labs(x = "Number of test items",
       y = "MSE",
       title = "MSE, Deletion") + 
  scale_color_manual(values=c("#8F993E", "#E05A1D")) + 
  ylim(0,2.5) +
  theme( 
    legend.title = element_blank(),
    plot.margin = margin(1,1,1,1, "cm"),
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom")
plotname= paste0("Plots/ROAR_plot.mse",".",subtask,".jpeg")
ggsave(plotname, plot=plot.mse)


}
##################################################################################################
##################################################################################################
