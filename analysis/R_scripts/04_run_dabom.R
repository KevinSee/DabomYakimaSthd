# Author: Kevin See
# Purpose: prep and run DABOM
# Created: 2/27/20
# Last Modified: 4/10/23
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)
library(rjags)
library(magrittr)
library(lubridate)
library(here)

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data/site_config.rda'))

#-----------------------------------------------------------------
# Load required DABOM data
#-----------------------------------------------------------------
# set year
yr = 2020

# for(yr in 2011:2019) {
#   cat(paste("Working on", yr, "\n\n"))

# load processed detection histories
load(here('analysis/data/derived_data/PITcleanr',
          paste0('PRO_Steelhead_', yr, '.rda')))

# filter to keep only the observations you want to keep
filter_obs = prepped_ch %>%
  mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                 auto_keep_obs,
                                 user_keep_obs)) %>%
  filter(user_keep_obs)


# determine origin of each fish
fish_origin = bio_df %>%
  filter(tag_code %in% unique(filter_obs$tag_code)) %>%
  mutate(origin = if_else(spp_code == "wsth",
                          "W",
                          "H")) %>%
  select(tag_code, origin) %>%
  distinct()

# fix all origins to wild
fish_origin <- fish_origin %>%
  mutate(origin = "W")

# file path to the default and initial model
basic_modNm = here('analysis/model_files', "PRO_DABOM.txt")

writeDABOM(file_name = basic_modNm,
           parent_child = parent_child,
           configuration = configuration,
           time_varying = F)

#------------------------------------------------------------------------------
# Alter default model code for species and year of
# interest; sets prior for some detection node efficiencies at 0 or 100%
# based on actual tag detection data; 0% if no tags were seen
#------------------------------------------------------------------------------

# filepath for specific JAGS model code for species and year
mod_path = here('analysis/model_files',
                paste0('PRO_Steelhead_', yr, '.txt'))

# writes species and year specific jags code
fixNoFishNodes(init_file = basic_modNm,
               file_name = mod_path,
               filter_ch = filter_obs,
               parent_child = parent_child,
               configuration = configuration,
               fish_origin = fish_origin)

#------------------------------------------------------------------------------
# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues(filter_obs,
                            parent_child,
                            configuration)

# Create all the input data for the JAGS model
jags_data = createJAGSinputs(filter_ch = filter_obs,
                             parent_child = parent_child,
                             configuration = configuration,
                             fish_origin = fish_origin)

# Tell JAGS which parameters in the model that it should save.
jags_params = setSavedParams(model_file = mod_path,
                             time_varying = F)


# run the model
set.seed(5)
jags = jags.model(mod_path,
                  data = jags_data,
                  inits = init_fnc,
                  # n.chains = 1,
                  # n.adapt = 5)
                  n.chains = 4,
                  n.adapt = 5000)


#--------------------------------------
# test the MCMC outcome and summary functions
dabom_mod = coda.samples(jags,
                         jags_params,
                         # n.iter = 10)
                         n.iter = 5000,
                         thin = 10)


save(dabom_mod, jags_data, filter_obs, bio_df,
     file = here("analysis/data/derived_data/model_fits",
                 paste0('PRO_DABOM_Steelhead_', yr,'.rda')))

# rm(dabom_mod, jags_data, filter_obs)
# }

#------------------------------------------------------------------------------
# diagnostics
#------------------------------------------------------------------------------
# load model run
load(here("analysis/data/derived_data/model_fits",
          paste0('PRO_DABOM_Steelhead_', yr,'.rda')))

# using mcmcr package
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
  full_join(esr(my_mcmcr,
                by = 'parameter',
                as_df = T)) %>%
  mutate(type = if_else(grepl('_p$', parameter),
                        'Detection',
                        if_else(grepl('^psi', parameter) |
                                  grepl('^phi', parameter),
                                'Movement',
                                'Other')))

# plot histogram of Rhat statistics
rhat_df %>%
  ggplot(aes(x = rhat)) +
  geom_histogram(fill = 'blue',
                 # bins = 40) +
                 binwidth = 0.001) +
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
  # filter(!converged) %>%
  left_join(rhat_df) %>%
  arrange(esr)

#---------------------------------------
# using postpack
library(postpack)

# what parameters were tracked?
get_params(my_mod,
           type = 'base_only')

# some summary statistics
post_summ(my_mod,
          '_p$') %>%
  t() %>%
  as_tibble(rownames = 'param')

post_summ(my_mod,
          '^phi') %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(mean > 0)

post_summ(my_mod,
          '^psi') %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(mean > 0)




param_chk = c('psi_ROZ')
param_chk = convg_df %>%
  filter(!converged) %>%
  pull(parameter)

post_summ(my_mod,
          param_chk) %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  mutate(cv = sd / mean) %>%
  arrange(desc(cv))


diag_plots(post = my_mod,
           p = param_chk,
           save = F,
           file = here('outgoing/PRO_diagnostics.pdf'))

# calculate Brooks-Gelman-Rubin Potential Scale Reduction Factor (Rhat)
# if ratio is close to 1, the chains have converged to the same distribution
# <1.10 is generally considered converged
post_summ(my_mod,
          # '_p$',
          get_params(my_mod,
                     type = 'base_only'),
          neff = T, # effective sample size
          Rhat = T)[c("Rhat", "neff"),] %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(!is.na(Rhat)) %>%
  arrange(neff)

# find and remove params where Rhat == "NaN"
all_params = get_params(my_mod,
                        type = 'base_only')


post_summ_nas = post_summ(my_mod,
                          # '_p$',
                          all_params,
                          neff = T, # effective sample size
                          Rhat = T)[c("Rhat", "neff"),] %>%
  t() %>%
  as.data.frame() %>%
  as_tibble(rownames = 'param') %>%
  filter(Rhat == "NaN") %>%
  pull(param)

param_chk = get_params(my_mod,
                       type = 'base_only')[grep('_p$', get_params(my_mod, type = 'base_only'))]
param_chk = param_chk[!param_chk %in% post_summ_nas]

# diagnostic plots for remaining params
diag_plots(post = my_mod,
           p = param_chk,
           ext_device = T)
