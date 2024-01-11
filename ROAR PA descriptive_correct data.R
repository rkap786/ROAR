library(readr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(mirt)
library(ggmirt)
library(here)
setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR")
here::here()
source(here("ROAR/PA/ROAR","imv.binary.R"))
df= read_csv("PA/roar_pa_trialdata.csv")
df = df[,-1]
df = df |> mutate(removed = ifelse(is.na(removed), 0, removed)) |>
  rename(item= itemId)

table(df$block)
n_distinct(df$subj)
df_resp = df |> dplyr::select(subj, block, correct, rt)
df_resp_DEL = df |> filter(block=="DEL")
df_resp_FSM = df |> filter(block=="FSM")
df_resp_LSM = df |> filter(block=="LSM")

n_distinct(df_resp_DEL$item) # 24 items
n_distinct(df_resp_FSM$item) # 25 items
n_distinct(df_resp_LSM$item) # 25 items
n_distinct(df_resp_DEL$subj) # 269
n_distinct(df_resp_FSM$subj) # 272
n_distinct(df_resp_LSM$subj) # 272

## Check item ID and item order
df= df |> mutate(item2 = strsplit(item, "_")) |> unnest(item2)  |> 
  filter(!item2 %in% c("FSM", "DEL", "LSM")) 



table(df_resp_DEL$item,df_resp_DEL$item_order)
  
 
getrespmatrix= function(data) {
 data = data |>
   dplyr::select(subj,item,correct) |>
   pivot_wider(values_from= correct,
               names_from = item)
 
 #names(data) = paste0("item",names(data))
 #names(data)[1] = "personID"
 data
}


itemfitPlot = function(m1, data_col, fitstat="infit", fitvalue=c("infit", "outfit"), 
                       metric=c(0.5,1,1.5), limit=c(0,2), title) {
  #metric=c(0.5,1,1.5)  
  # limit=c(0,2)  
  fit = itemfit(m1, fitstat)
  fit = fit %>%
    dplyr::select(item, !!fitvalue) 
  fit=  left_join(fit, data_col)
  
  fit |>
    gather(key, value, -item, -removed) |>
    ggplot(aes(x = item, y = value, color = factor(removed))) +
    geom_point(size = 3,  shape = shape) +
    geom_line() +
    geom_hline(yintercept = metric[1], color = "darkgrey", linetype = "dashed") +
    geom_hline(yintercept = metric[2], color = "darkgrey") +
    geom_hline(yintercept = metric[3], color = "darkgrey", linetype = "dashed") +
    scale_y_continuous(breaks = metric, limits = limit) +
    facet_grid(~key) +
    coord_flip() +
    theme_minimal() +
    labs(y = "", x = "", 
         caption = "Note: Items with values within dashed lines are considered to be productive for measurement",
         title=title
    ) + theme(legend.position = "bottom")
  
}

# df2= read_csv("pa_block_clean.csv")
# head(df2)
# names(df2)


table_summary= df |> group_by(item_order, block, removed) |> 
  summarise(mean_correct=mean(correct), 
            mean_rt= mean(rt)/1000) 


### Correctness vs 
data.frame(table_summary) |>
  ggplot(aes(x=factor(item_order, levels=sort(as.numeric(unique(df$item_order)))), 
             y=mean_correct, group=factor(removed), color=factor(removed))) + 
  geom_point() + facet_wrap(~ block, ncol=1) + ylim(0.5,1) + theme(legend.position = "bottom") +
  labs(x="Item Order")
  
data.frame(table_summary) |>
  ggplot(aes(x=factor(item_order, levels=sort(as.numeric(unique(df$item_order)))), 
             y=mean_rt, group=factor(removed), color=factor(removed))) + 
  geom_point() + facet_wrap(~ block, ncol=1) + ylim(2,13) + theme(legend.position = "bottom") +
  labs(x="Item Order")



########


#Estimate model - FSM
resp_fsm = df_resp_FSM |> getrespmatrix()
fsm_remove_item = df_resp_FSM |> dplyr::select(item, removed) |> distinct() 
cronbach.alpha(resp_fsm[,-1], standardized = T, CI=T)
#colMeans(resp_del[,-1], na.rm=T)
m1= mirt(resp_fsm[-1],1,  itemtype = 'Rasch', guess=0.3)
M2(m1)
m1_coef = as.data.frame(coef(m1,simplify=TRUE, IRTpars = TRUE))
m1_coef$item = row.names(m1_coef)
m1_coef = merge(m1_coef, fsm_remove_item, by="item")

m1_fsm_time = df_resp_FSM |> group_by(item) |> summarise(mean_rt=mean(rt, na.rm=T)/1000)
m1_coef_time= merge(m1_coef, m1_fsm_time, by="item")
m1_coef_time |> ggplot(aes(items.b, mean_rt)) + geom_point() +
  geom_smooth(method = "lm", se = T)
  
M2(m1, na.rm=T) #model fit is good: M2 is not significant, RMSEA<0.05, CFI~1
# del_f1 <- fscores(m1)
# del_m1.coef <- as.data.frame(coef(m1,simplify=TRUE, IRTpars = TRUE))
m1_itemfit = itemfit(m1, c('X2', 'S_X2', 'infit'))
plot(m1, type='infotrace')


#itemfit(m1, S_X2.plot=10, na.rm=T)
#itemfit(m1, group.bins=15, empirical.plot = 17, method="ML")
# fit = itemfit(m1, c('S_X2'))
# fit = fit |>
#   dplyr::select(item, infit, outfit)
# 
# fit=  left_join(fit, del_resp_order)
itemfitPlot(m1, data_col=fsm_remove_item,
            fitstat="infit", fitvalue=c("infit", "outfit"), 
         metric=c(0.5,1,1.5), limit=c(0,2), title="FSM: Item Infit and Outfit Statistics")

itemfitPlot(m1, data_col=fsm_remove_item,
            fitstat="infit", fitvalue=c("z.infit", "z.outfit"), 
            metric=c(2,0, -2), limit=c(-3,3), title="Item Infit and Outfit Statistics (Standardized)")


itemfitPlot(m1, data_col=fsm_remove_item,
            fitstat="S_X2", fitvalue="RMSEA.S_X2", 
            metric=c(-0.06,0,0.06), limit=c(-0.1,0.1), title="RMSEA S_X2")

itemfitPlot(m1, data_col=fsm_remove_item,
            fitstat="X2", fitvalue="RMSEA.X2", 
            metric=c(-0.06,0,0.06), limit=c(-0.15,0.15), title="RMSEA X2")


personfitPlot(m1)
personfit(m1) %>%
  summarize(infit.outside = prop.table(table(z.infit > 1.96 | z.infit < -1.96)),
            outfit.outside = prop.table(table(z.outfit > 1.96 | z.outfit < -1.96)))

plot(m1, type = 'trace')
plot(m1, type = 'trace', facet_items=F)
plot(m1, type = 'trace', which.items = c(1,3, 9,10,16), facet_items=F)

itempersonMap(m1)
#itemInfoPlot(m1)
#tracePlot(m1, facet = F, legend = T) + scale_color_brewer(palette = "Set3")
plot(m1, type='infoSE')
#itemInfoPlot(m1, facet = T)
#testInfoPlot(m1, adj_factor = 2)


itemfitPlot(m1)


#Estimate model - DELETION
resp_del = df_resp_DEL |> dplyr::filter(item!="DEL_3") |> getrespmatrix()
cronbach.alpha(resp_del[,-1], standardized = T, CI=T)

del_resp_order = df_resp_DEL |> dplyr::select(item, removed) |> distinct() 
resp_del_2 = df_resp_DEL |>  getrespmatrix()
#colMeans(resp_del[,-1], na.rm=T)
m1= mirt(resp_del[-1],1,  itemtype = 'Rasch', guess=0.3)
m1_coef = as.data.frame(coef(m1,simplify=TRUE, IRTpars = TRUE))
m1_coef$item = row.names(m1_coef)
m1_coef = merge(m1_coef, del_resp_order, by="item")



M2(m1, na.rm=T) #model fit is good: M2 is not significant, RMSEA<0.05, CFI~1
# del_f1 <- fscores(m1)
# del_m1.coef <- as.data.frame(coef(m1,simplify=TRUE, IRTpars = TRUE))
m1_itemfit = itemfit(m1, c('X2', 'S_X2', 'infit'), na.rm=T)
plot(m1, type='infotrace')


#itemfit(m1, S_X2.plot=10, na.rm=T)
#itemfit(m1, group.bins=15, empirical.plot = 17, method="ML")


itemfitPlot(m1, data_col=del_resp_order,
            fitstat="infit", fitvalue=c("infit", "outfit"), 
            metric=c(0.5,1,1.5), limit=c(0,2), title="LSM: Item Infit and Outfit Statistics")

itemfitPlot(m1, data_col=del_resp_order,
            fitstat="infit", fitvalue=c("z.infit", "z.outfit"), 
            metric=c(2,0, -2), limit=c(-3,3), title="Item Infit and Outfit Statistics (Standardized)")


itemfitPlot(m1, data_col=del_resp_order,
            fitstat="S_X2", fitvalue="RMSEA.S_X2", 
            metric=c(-0.06,0,0.06), limit=c(-0.1,0.1), title="RMSEA S_X2")

itemfitPlot(m1, data_col=del_resp_order,
            fitstat="X2", fitvalue="RMSEA.X2", 
            metric=c(-0.06,0,0.06), limit=c(-0.15,0.15), title="RMSEA X2")


itempersonMap(m1)
#itemInfoPlot(m1)
plot(m1, type = 'trace')
plot(m1, type = 'trace', which.items = c(1,4,11, 12, 16), facet_items=F)
#tracePlot(m1, facet = F, legend = T) + scale_color_brewer(palette = "Set3")
plot(m1, type='infoSE')
itemInfoPlot(m1)
#testInfoPlot(m1, adj_factor = 2)


itemfitPlot(m1)
personfitPlot(m1)
personfit(m1) %>%
  summarize(infit.outside = prop.table(table(z.infit > 1.96 | z.infit < -1.96)),
            outfit.outside = prop.table(table(z.outfit > 1.96 | z.outfit < -1.96)))



########## Last sound matching

resp_lsm = df_resp_LSM |>  getrespmatrix()
summary(resp_lsm)
cronbach.alpha(resp_lsm[,-1], standardized = T, CI=T)

lsm_remove_item = df_resp_LSM |> dplyr::select(item, removed) |> distinct() 

#colMeans(resp_del[,-1], na.rm=T)
m1= mirt(resp_lsm[-1],1,  itemtype = 'Rasch', guess=0.3)
m1_coef = as.data.frame(coef(m1,simplify=TRUE, IRTpars = TRUE))
m1_coef$item = row.names(m1_coef)
m1_coef = merge(m1_coef, lsm_remove_item, by="item")



M2(m1) #model fit is good: M2 is not significant, RMSEA<0.05, CFI~1
# del_f1 <- fscores(m1)
# del_m1.coef <- as.data.frame(coef(m1,simplify=TRUE, IRTpars = TRUE))
m1_itemfit = itemfit(m1, c('X2', 'S_X2', 'infit'))
plot(m1, type='infotrace')


#itemfit(m1, S_X2.plot=10, na.rm=T)
#itemfit(m1, group.bins=15, empirical.plot = 17, method="ML")

itemfitPlot(m1, data_col=lsm_remove_item,
            fitstat="infit", fitvalue=c("infit", "outfit"), 
            metric=c(0.5,1,1.5), limit=c(0,2), title="LSM: Item Infit and Outfit Statistics")

itemfitPlot(m1, data_col=lsm_remove_item,
            fitstat="infit", fitvalue=c("z.infit", "z.outfit"), 
            metric=c(2,0, -2), limit=c(-3,3), title="Item Infit and Outfit Statistics (Standardized)")


itemfitPlot(m1, data_col=lsm_remove_item,
            fitstat="S_X2", fitvalue="RMSEA.S_X2", 
            metric=c(-0.06,0,0.06), limit=c(-0.1,0.1), title="RMSEA S_X2")

itemfitPlot(m1, data_col=lsm_remove_item,
            fitstat="X2", fitvalue="RMSEA.X2", 
            metric=c(-0.06,0,0.06), limit=c(-0.15,0.15), title="RMSEA X2")




itempersonMap(m1)
itemInfoPlot(m1)
plot(m1, type = 'trace')
plot(m1, type = 'trace', which.items = c(1,4,9, 23,25), facet_items=F)
#tracePlot(m1, facet = F, legend = T) + scale_color_brewer(palette = "Set3")
plot(m1, type='infoSE')
#itemInfoPlot(m1, facet = T)
#testInfoPlot(m1, adj_factor = 2)


personfitPlot(m1)
personfit(m1) %>%
  summarize(infit.outside = prop.table(table(z.infit > 1.96 | z.infit < -1.96)),
            outfit.outside = prop.table(table(z.outfit > 1.96 | z.outfit < -1.96)))




