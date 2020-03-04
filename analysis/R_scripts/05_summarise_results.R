# Author: Kevin See
# Purpose: summarise DABOM results
# Created: 3/4/20
# Last Modified: 3/4/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(DABOM)
library(tidyverse)
library(jagsUI)
library(STADEM)
library(WriteXLS)

#-----------------------------------------------------------------
# set species
spp = "Steelhead"
# set year
yr = 2019

#-----------------------------------------------------------------
# load JAGS MCMC results
load(paste0("analysis/data/derived_data/model_fits/PRO_DABOM_", spp, '_', yr,'.rda'))

# summarise detection probabilities
detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   capHist_proc = proc_list$proc_ch) %>%
  filter(!is.na(mean))

# which sites had detection probabilities fixed at 0% or 100%
detect_summ %>%
  filter(sd == 0)

# look at all the other sites
detect_summ %>%
  filter(sd > 0) %>%
  arrange(desc(sd))

# compile all movement probabilities, and multiply them appropriately
trans_df = compileTransProbs_PRO(dabom_mod,
                                 parent_child)

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
            upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2]) %>%
  mutate_at(vars(mean, median, mode, sd, matches('CI$')),
            list(~ if_else(. < 0, 0, .))) %>%
  ungroup()


#-----------------------------------------------------------------
# total escapement past Prosser
# window count
tot_win_cnt = getWindowCounts(dam = 'PRO',
                              spp = 'Steelhead',
                              start_date = paste0(yr-1, '0701'),
                              end_date = paste0(yr, '0630')) %>%
  summarise_at(vars(win_cnt),
               list(sum)) %>%
  pull(win_cnt)

# Chris Frederiksen says total count was 1132
tot_win_cnt = 1132


# translate movement estimates to escapement
escape_summ = trans_df %>%
  filter(origin == 1) %>%
  mutate(tot_escp = tot_win_cnt,
         escp = value * tot_escp) %>%
  group_by(location = param) %>%
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

# generate population level estimates
pop_summ = trans_df %>%
  filter(origin == 1) %>%
  mutate(tot_escp = tot_win_cnt,
         escp = value * tot_escp) %>%
  select(-value) %>%
  spread(param, escp) %>%
  mutate(Status = SAT,
         Toppenish = TOP,
         Naches = LNR + AH1,
         `Upper Yakima` = LWC + ROZ,
         Sunnyside_bb = SUN_bb) %>%
  select(chain:tot_escp, Status:Sunnyside_bb) %>%
  gather(pop, escp, Status:Sunnyside_bb) %>%
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
  summarise_at(vars(win_cnt),
               list(sum)) %>%
  pull(win_cnt)

escape_summ %>%
  filter(location == 'ROZ')

#-----------------------------------------------------------------
# write results to an Excel file
save_list = list('Population Escapement' = pop_summ %>%
                   select(-skew, -kurtosis) %>%
                   mutate_at(vars(-pop),
                             list(round),
                             digits = 1),
                 'All Escapement' = escape_summ %>%
                   select(-skew, -kurtosis) %>%
                   mutate_at(vars(-location),
                             list(round),
                             digits = 1),
                 'Detection' = detect_summ %>%
                   mutate_at(vars(-Node),
                             list(round),
                             digits = 3))

WriteXLS(x = save_list,
         ExcelFileName = paste0('outgoing/estimates/PRO_est_', spp, '_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx'),
         AdjWidth = T,
         AutoFilter = F,
         BoldHeaderRow = T,
         FreezeRow = 1)
