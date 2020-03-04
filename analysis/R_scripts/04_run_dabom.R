# Author: Kevin See
# Purpose: prep and run DABOM
# Created: 2/27/20
# Last Modified: 3/4/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)
library(jagsUI)

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

# writes species and year specific jags code
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

# add fish origin
# call them all wild for now
dabom_list$fishOrigin = rep(1, nrow(dabom_list[[1]]))

# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues_PRO(dabom_list,
                                mod_path,
                                parent_child)

# Create all the input data for the JAGS model
jags_data = createJAGSinputs_PRO(dabom_list,
                                 mod_path,
                                 parent_child)


#------------------------------------------------------------------------------
# Tell JAGS which parameters in the model that it should save.
# the fnc is hard coded and needs to be updated if there are changes!
#------------------------------------------------------------------------------

jags_params = setSavedParams(model_file = mod_path)


#------------------------------------------------------------------------------
# Run the model

# Recommended MCMC parameters are:
# * `n.chains`: 4
# * `n.iter`: 5,000
# * `n.burnin`: 2,500
# * `n.thin`: 10
# ( 4*(5000-2500) ) / 10 = 1000 samples

set.seed(12)
dabom_mod <- jags.basic(data = jags_data,
                        inits = init_fnc,
                        parameters.to.save = jags_params,
                        model.file = mod_path,
                        n.chains = 4,
                        n.iter = 5000,
                        n.burnin = 2500,
                        n.thin = 10,
                        # n.chains = 1,
                        # n.iter = 2,
                        # n.burnin = 1,
                        DIC = T)

# save some stuff
proc_list[["proc_ch"]] <- proc_ch

save(dabom_mod, dabom_list, proc_list,
     file = paste0("analysis/data/derived_data/model_fits/PRO_DABOM_", spp, '_', yr,'.rda'))



#------------------------------------------------------------------------------
# diagnostics
library(mcmcr)

# pull out mcmc.list object
my_mod = dabom_mod

#---------------------------------------
# using mcmcr
anyNA(my_mod)
my_mcmcr = as.mcmcr(my_mod)

# get Rhat statistics for all parameters
rhat_df = rhat(my_mcmcr,
               by = 'parameter',
               as_df = T) %>%
  mutate(type = if_else(grepl('_p$', parameter),
                        'Detection',
                        if_else(grepl('^p_pop', parameter) |
                                  grepl('^phi', parameter),
                                'Movement',
                                'Other')))

# plot histogram of Rhat statistics
rhat_df %>%
  ggplot(aes(x = rhat)) +
  geom_histogram(fill = 'blue',
                 bins = 40) +
  facet_wrap(~ type,
             scales = 'free')

# which parameters have converged and which haven't?
convg_df = converged(my_mcmcr,
                     by = 'parameter',
                     as_df = T)

janitor::tabyl(convg_df,
               converged)

# look at parameters that have not converged
convg_df %>%
  filter(!converged) %>%
  left_join(rhat_df)
