################################################################################
#                                                                              #
# SEER SEMISYNTH PREPROCESSING                                                 #
#                                                                              #
# Contains additional preprocessing for the semi-synthetic analysis of the     #
# SEER data                                                                    #
#                                                                              #
################################################################################

remove(list=ls()) 

#path to output pre-processed ds
out.path <- "..."

load("breast04.RData")

library(data.table)
library(forcats)

breast04semis <- as.data.table(breast04)

#recode marstat (marital status) to include 3 categories only (single, married, other)
breast04semis$marstatrec <- fct_collapse(breast04semis$marstat, "3" = c("3","4","5"))
#recode racrecy (race) to include 3 categories only (white, black, other)
breast04semis$racrecyrec <- fct_collapse(breast04semis$racrecy, "3" = c("3","4","9"))
#recode radiatn (radiation therapy) to include 2 categories only (none, any kind of radiation therapy)
breast04semis$radiatnrec <- fct_collapse(breast04semis$radiatn, "0" = c("0","7"), "1" = c("1","2","3","4","5"))

#drop unbalanced varaibles (radsurgrec (order of radiation therapy/surgery), adjtm6valuerec (adjusted T/M variable)) & old variables
breast04semis <- breast04semis[,c("marstat","racrecy","radiatn","radsurgrec","adjtm6valuerec") := NULL]

#save to .RData
save(breast04semis, list = "breast04semis", file = file.path(out.path,"breast04semis.RData"))

#clean workspace, remove breast04
rm(breast04)
gc()
