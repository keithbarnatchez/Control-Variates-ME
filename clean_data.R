# code for proccessing NHANES data
rm(list=ls())
library(haven)
library(tidyverse)

demo <- read_xpt('data/raw/DEMO_H.XPT')
cot <- read_xpt('data/raw/COT_H.XPT')
smoke <- read_xpt('data/raw/SMQRTU_H.XPT')
smoker <- read_xpt('data/raw/SMQ_H.XPT')

# First, categorize smoking status
# non-smoker if < 100 in lifetime or >= 100 but not current smoker
# smoker if otherwise
smoker <- smoker %>% mutate(non_smoker = ifelse( (SMQ020==2 | SMQ040==3),
                                                 1,0)) %>% 
  mutate(non_smoker = replace(non_smoker,
                              is.na(SMQ020),
                              NA)) %>% 
  select(c(SEQN,non_smoker,SMQ020, SMQ040)) %>% 
  mutate(smoker_report=1-non_smoker)

# Next, get lab status
cot <- cot %>% mutate(smoker_cot = ifelse(LBXCOT>=100,1,0)) %>%
  select(c(SEQN,smoker_cot,LBXCOT))

smoking_merged <- inner_join(smoker, cot, by='SEQN') 

smoking_inc <- smoking_merged %>% filter(smoker_report == 0, smoker_cot==1)
