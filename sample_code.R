#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(lubridate)
library(mgcv)
library(gamm4)
library(fields)

#=============================================================================
#Load data
#=============================================================================

master_wSUMMED_BNF= read.csv("Nact_master_datafile.csv")

#Dependent variables: NPP
#Dependent variables: BNF 
#Ratio of algae/dias
#48 tanks in 12 metacommunities: 
#Only used no dispersal
#8 Dates, measured
#3 community types: h20 water column, peri(phyton), sed(iment)
#adj_eth_area - measure of 
# bluegreen:yellow_substances are measures from water column samples
#chl_a_fluoro  = sum of algaes
#total_adj_eth_area = within a tank, sediment BNF+water BNF + peri BNF
#bluegreen is a diazatroph 
#Does temperature response of N fixation effect the response of NPP, given that when 
#N is low, NF rates determine the rate of N flow into the sytem. 
#MT doesn't take into account absolute quantities of nutrients. 
#You can say photosynthesis is limiting respiration, because of Carbon. 
#But don't have a piece in metabolic theory for limiting nutrient supply. 


algae_npp1 <- lmer(NPP ~ mean_chla + (1 | tank), data = master_wSUMMED_BNF) 
alg_gam = gam(NPP ~ mean_chla + te(tank, date, bs="re"), data = master_wSUMMED_BNF )
#I think this is the better way: 
alg_gam = gam(NPP ~ mean_chla + s(date, bs="re")+s(tank, bs="re"), data = master_wSUMMED_BNF )

#=============================================================================
#Define the growth rate functions for the LG population model(s)
#=============================================================================
#The model for the resident with only the resident species
LG_res = formula (Ndiff_res ~ lambda/(1+alpha_ii*N_res) + 0.5)
m1=nls(LG_res, data = mydata1, start=list(lambda = 1.1, alpha_ii=0.5) )

