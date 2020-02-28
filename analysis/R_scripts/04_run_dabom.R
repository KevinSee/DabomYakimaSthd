# Author: Kevin See
# Purpose: prep and run DABOM
# Created: 2/27/20
# Last Modified: 2/27/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)

spp = "Steelhead"
#-----------------------------------------------------------------
# load configuration and site_df data

#-----------------------------------------------------------------
# Load required DABOM data
#-----------------------------------------------------------------
# set year
yr = 2019

load(paste0('analysis/data/derived_data/PITcleanr/PRO_', spp, '_', yr, '.rda'))

proc_ch <- proc_list$ProcCapHist %>%
  mutate(UserProcStatus = AutoProcStatus) %>%
  filter(UserProcStatus)


# file path to the default and initial model
basic_modNm = 'analysis/model_files/PRO_DABOM.txt'

writeDABOM_PRO(file_name = basic_modNm)

#------------------------------------------------------------------------------
# Alter default model code for species and year of
# interest; sets prior for some detection node efficiencies at 0 or 100%
# based on actual tag detection data; 0% if no tags were seen
#------------------------------------------------------------------------------

# filepath for specific JAGS model code for species and year
mod_path = paste0('analysis/model_files/PRO_', spp, '_', yr, '.txt')

#writes species and year specific jags code
fixNoFishNodes(basic_modNm,
               mod_path,
               proc_ch,
               proc_list$NodeOrder)

#------------------------------------------------------------------------------
# Create capture history matrices for each main branch to be used in
# the JAGS data list
#------------------------------------------------------------------------------
dabom_list = createDABOMcapHist(proc_ch,
                                proc_list$NodeOrder,
                                split_matrices = T)


# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues_PRO(dabom_list)

#Create all the input data for the JAGS model
jags_data = createJAGSinputs_PRO(dabom_list)


#------------------------------------------------------------------------------
# Tell JAGS which parameters in the model that it should save.
# the fnc is hard coded and needs to be updated if there are changes!
#------------------------------------------------------------------------------

jags_params = setSavedParams(model_file = mod_path)
