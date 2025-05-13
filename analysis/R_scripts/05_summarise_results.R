# Author: Kevin See
# Purpose: summarise DABOM results
# Created: 3/4/20
# Last Modified: 4/10/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(DABOM)
library(PITcleanr)
library(tidyverse)
library(jagsUI)
library(STADEM)
library(writexl)
library(moments)
library(coda)
library(here)

#-----------------------------------------------------------------
# set species
spp = "Steelhead"
# set year
yr = 2020

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data/site_config.rda'))

# load JAGS MCMC results
load(here("analysis/data/derived_data/model_fits",
          paste0("PRO_DABOM_", spp, '_', yr,'.rda')))

# summarise detection probabilities
detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   filter_ch = filter_obs)

# which sites had detection probabilities fixed at 0% or 100%
detect_summ %>%
  filter(sd == 0)

# look at all the other sites
detect_summ %>%
  filter(sd > 0) %>%
  arrange(desc(sd))

# compile all movement probabilities, and multiply them appropriately
trans_df = compileTransProbs_PRO(dabom_mod,
                                 parent_child) %>%
  mutate(origin = recode(origin,
                         "2" = "H",
                         "1" = "W"))

# summarize transition probabilities
trans_summ = trans_df %>%
  group_by(origin, param) %>%
  summarise(mean = mean(value),
            median = median(value),
            mode = estMode(value),
            sd = sd(value),
            skew = moments::skewness(value),
            kurtosis = moments::kurtosis(value),
            lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2],
            .groups = "drop") %>%
  mutate(across(c(mean, median, mode, sd, matches('CI$')),
                ~ if_else(. < 0, 0, .)))

#-----------------------------------------------------------------
# total escapement past Prosser
# window count
tot_win_cnt = getWindowCounts(dam = 'PRO',
                              spp = 'Steelhead',
                              start_date = paste0(yr-1, '0701'),
                              end_date = paste0(yr, '0630')) %>%
  summarize(
    across(
      win_cnt,
      ~ sum(.,
            na.rm = T)
    )
  ) %>%
  pull(win_cnt)

# total counts from Yakima Nation (use these)
yak_cnts = tibble(year = c(2012:2014,
                           2019,
                           2020),
                  tot_win_cnt = c(6359,
                                  4787,
                                  4143,
                                  1132,
                                  1657))

if(yr %in% yak_cnts$year) {
  tot_win_cnt = yak_cnts %>%
    filter(year == yr) %>%
    pull(tot_win_cnt)
}

# translate movement estimates to escapement
escape_post <- trans_df %>%
  filter(origin == "W") %>%
  mutate(tot_escp = tot_win_cnt,
         escp = value * tot_escp)

# summary statistics
escape_summ = escape_post %>%
  group_by(location = param) %>%
  summarise(mean = mean(escp),
            median = median(escp),
            mode = estMode(escp),
            sd = sd(escp),
            skew = moments::skewness(escp),
            kurtosis = moments::kurtosis(escp),
            lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2]) %>%
  mutate(across(c(mean, median, mode, sd, matches('CI$')),
                ~ if_else(. < 0, 0, .))) %>%
  ungroup()

# generate population level estimates
pop_summ = escape_post %>%
  select(-value) %>%
  pivot_wider(names_from = param,
              values_from = escp) %>%
  mutate(Status = SAT,
         Toppenish = TOP,
         Naches = LNR + AH1,
         `Upper Yakima` = LWC + ROZ,
         Sunnyside_bb = SUN_bb) %>%
  select(chain:tot_escp, Status:Sunnyside_bb) %>%
  pivot_longer(cols = Status:Sunnyside_bb,
               names_to = "pop",
               values_to = "escp") %>%
  mutate(pop = factor(pop,
                      levels = c('Status',
                                 'Toppenish',
                                 'Naches',
                                 'Upper Yakima',
                                 'Sunnyside_bb'))) %>%
  group_by(pop) %>%
  summarise(mean = mean(escp),
            median = median(escp),
            mode = estMode(escp),
            sd = sd(escp),
            skew = moments::skewness(escp),
            kurtosis = moments::kurtosis(escp),
            lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2]) %>%
  mutate_at(vars(mean, median, mode, sd, matches('CI$')),
            list(~ if_else(. < 0, 0, .))) %>%
  ungroup()



# compare with dam counts at Roza dam
roz_win_cnt = getWindowCounts(dam = 'ROZ',
                              spp = 'Steelhead',
                              start_date = paste0(yr-1, '0701'),
                              end_date = paste0(yr, '0630')) %>%
  summarize(
    across(
      win_cnt,
      ~ sum(.,
            na.rm = T)
    )
  ) %>%
  pull(win_cnt)
roz_win_cnt

escape_summ %>%
  filter(location == 'ROZ')


#-----------------------------------------------------------------
# summarize tag data
tag_summ = summarizeTagData(filter_obs,
                            bio_df)

#-----------------------------------------------------------------
# save some of these objects
save(tag_summ,
     bio_df,
     trans_df,
     trans_summ,
     tot_win_cnt,
     escape_post,
     escape_summ,
     detect_summ,
     configuration,
     flowlines,
     parent_child,
     sites_sf,
     file = here("analysis/data/derived_data/estimates",
                 paste0("PRO_Sthd_DABOM_", yr, ".rda")))


#-----------------------------------------------------------------
# write results to an Excel file
save_list = list('Population Escapement' = pop_summ %>%
                   select(-skew, -kurtosis) %>%
                   mutate(across(-pop,
                                 ~ round(.,
                                         digits = 1))),
                 'All Escapement' = escape_summ %>%
                   select(-skew, -kurtosis) %>%
                   mutate(across(-location,
                                 ~ round(.,
                                         digits = 1))),
                 'Detection' = detect_summ %>%
                   mutate(across(-node,
                                 ~ round(.,
                                         digits = 3))),
                 "Tag Summary" = tag_summ)

write_xlsx(x = save_list,
           path = here('outgoing/estimates',
                       paste0('PRO_est_Steelhead_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx')))
