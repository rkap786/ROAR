setwd("/Users/radhika/Library/CloudStorage/GoogleDrive-rkap786@stanford.edu/My Drive/0. Projects - Stanford/ROAR/")
source("PA/ROAR/cat.functions.R")

library(readxl)
library(dplyr)
library(catR)
library(Metrics)

############################################################

## Import data
#ROAR_scores= read_csv("PA data/roar_pa_trialdata.csv")
ROAR_items= read_csv("PA data/test (1).csv")
ROAR_items= ROAR_items |>
  mutate(ItemId = paste0(trial_type, "_", original_name)) |>
  rename(block=trial_type)
names(ROAR_items)

param.unidim= read_csv("PA/Data/itemparam_overall_unidim.csv")
param.undim.del=read_csv("PA/Data/irtparam_DEL.csv")
param.undim.lsm=read_csv("PA/Data/irtparam_LSM.csv")
param.undim.fsm=read_csv("PA/Data/irtparam_FSM.csv")

names(param.unidim) = c("a.overall", "easiness.overall", "guess.overall", 
                        "u.overall", "ItemId")
names(param.undim.del)= c("a.subtask", "easiness.subtask", "guess.subtask", "u.subtask", "ItemId")
names(param.undim.lsm)= c("a.subtask", "easiness.subtask", "guess.subtask", "u.subtask", "ItemId")
names(param.undim.fsm)= c("a.subtask", "easiness.subtask", "guess.subtask", "u.subtask", "ItemId")

param.undim.subtask= bind_rows(param.undim.del,param.undim.lsm,param.undim.fsm)

ROAR_items= left_join(ROAR_items, param.unidim)
ROAR_items= left_join(ROAR_items, param.undim.subtask)
write_csv(ROAR_items,"PA/Data/items_with_param.csv")
