
setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
library(readxl)
library(dplyr)
library(catR)
library(Metrics)
###############Functions we will need
getrespmatrix= function(data) {
  data = data |>
    dplyr::select(subj,item,correct) |>
    pivot_wider(values_from= correct,
                names_from = item)
  
  #names(data) = paste0("item",names(data))
  #names(data)[1] = "personID"
  data
}


func.catSim <- function(resp=datasim, item.bank, item.bank.u=coef,method='ML', 
                        theta=th_u, stoplen=20,skipitems=3, pids){
  
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
  test <- list(method = 'ML', itemSelect = method, infoType = "Fisher")
  #test<-list(method='ML',itemSelect=method, randomesque=3, infoType = "Fisher") 
  final <- list(method = 'ML')
  stop<-list(rule="length",thr=stoplen)
  cl<-makeCluster(parallel::detectCores())
  registerDoSNOW(cl)
  for (i in 1:np){
    pid <- pids[i]
    start <- list(nrItems=skipitems,theta = theta[i],startSelect="MFI") 
    #start <- list(theta = theta[i],startSelect="MFI", randomesque=5) 
    #random_number <- sample(c(1, ncol(resp)), size = 1)
    #start <- list(nrItems = 10, theta = 0,randomesque=5, seed=1234)
    
    res <- randomCAT(itemBank = item.bank,
                     responses = as.numeric(resp[i,]),
                     start = start,
                     test = test,
                     final = final,
                     stop = stop 
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
th_u_roar= fscores(model_unidim)
th_u_roar= cbind(th_u_roar,resp_roar[,1])
names(th_u_roar)=c("th","subj")

## Overall item paramters
coef = coef(model_unidim, simplify=T)$items
coef[,2]= -coef[,2]
coef=data.frame(coef)
coef$block = rownames(coef)

item.bank.overall=list()
#del, fsm, lsm
item.bank.overall[['DEL']]=coef |> filter(str_detect(block,"DEL") ) |> dplyr::select(-block)


th_merged= merge(th_u_roar, CTOPP_scores, by="subj")
cor(th_merged$CTOPP_PA_raw, th_merged$th)
## lets use raw CTOPP, cor=0.73

############# Lets focus on deletion subtask
data = ROAR_scores |> filter(block=="DEL") |>  filter(subj %in% subj.list) |>
  getrespmatrix()
df = data |> select(-subj)

## Unidim model overall
m_dim= mirt(df,1,  itemtype = 'Rasch')
th= fscores(m_dim)
th_se=fscores(m_dim, full.scores.SE=TRUE)
rel=empirical_rxx(as.matrix(tibble(F1 = th_se[,1], SE_F1 = th_se[,2])))
th_dims=data.frame(th)
coef_dims = coef(m_dim, simplify=T)$items
coef_dims[,2]= -coef_dims[,2] #convert easiness to difficulty
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
coef= item.bank.overall[['DEL']]
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


##################################################################################################
### Calculate metrics
### Hardcoded for 3 dimensions

quartiles1 <- quantile(th_dims$th1,
                       probs = c(0, 0.20, 0.40, 0.60, 0.80, 1))

metrics.dim1= results.compare |>
  mutate(theta.bin = cut(th1, 
                         breaks = quartiles1, include.lowest = TRUE, 
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))) |>
  group_by(trialNumTotal,type) |>
  dplyr::summarise(sem = mean(thetaSE), 
                   reliability = empirical_rxx(as.matrix(tibble(F1 = thetaEstimate, SE_F1 = thetaSE))),
                   mse = Metrics::mse(th1, thetaEstimate), 
                   bias = Metrics::bias(th1, thetaEstimate)) 





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


