# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 2/21/20
# Last Modified: 4/10/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(lubridate)
library(janitor)
library(readxl)
library(magrittr)
library(here)

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data',
          'site_config.rda'))

# which spawn year are we dealing with?
yr = 2022

# for(yr in c(2012:2014, 2019, 2020)) {

#-----------------------------------------------------------------
# start date is July 1 of previous year
start_date = paste0(yr - 1, '0701')
# when is the last possible observation date?
# we chose one entire year from the start date
max_obs_date = as.character(ymd(start_date) + years(1) - days(1))

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
ptagis_file = here("analysis/data/raw_data/PTAGIS",
                   paste0("PTAGIS_", yr, ".csv"))

# recode the PTAGIS observations of double tagged fish so that the tag code matches the TagID (not TagOther)
ptagis_obs = readCTH(ptagis_file)

# any orphaned or disowned tags?
qcTagHistory(ptagis_obs, T)

# compress and process those observations with PITcleanr
prepped_ch = PITcleanr::prepWrapper(ptagis_file = ptagis_obs,
                                    configuration = configuration,
                                    parent_child = parent_child %>%
                                      addParentChildNodes(configuration = configuration),
                                    min_obs_date = start_date,
                                    max_obs_date = max_obs_date,
                                    ignore_event_vs_release = F,
                                    add_tag_detects = T,
                                    save_file = F,
                                    file_name = here('outgoing/PITcleanr',
                                                     paste0('PRO_Steelhead_', yr, '.xlsx')))



# load and filter biological data by spawn year
bio_df = read_rds(here('analysis/data/derived_data/Bio_Data_2012_2020.rds')) %>%
  filter(spawn_year == yr)

# save some stuff
save(parent_child, configuration, start_date, bio_df, prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('PRO_Steelhead_', yr, '.rda')))

# }

#-------------------------------------------
# NEXT STEPS
#-------------------------------------------
# open that Excel file, and filter on the column user_keep_obs, looking for blanks. Fill in each row with TRUE or FALSE, depending on whether that observation should be kept or not. The column auto_keep_obs provides a suggestion, but the biologist's best expert judgment should be used based on detection dates, detection locations before and after, etc.

load(here('analysis/data/derived_data/PITcleanr',
          paste0('PRO_Steelhead_', yr, '.rda')))

# read in the updated version of the PITcleanr output Excel file
yn_df = read_excel(here('analysis/data/derived_data/YakimaNation',
                          paste0('PRO_Steelhead_', yr, '.xlsx')))

filter_obs = yn_df %>%
  mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                 auto_keep_obs,
                                 user_keep_obs)) %>%
  filter(user_keep_obs)

# construct all valid paths
all_paths = buildPaths(addParentChildNodes(parent_child,
                                           configuration))

tag_path = summarizeTagData(filter_obs,
                            bio_df) %>%
  select(tag_code, spawn_node) %>%
  distinct() %>%
  left_join(all_paths,
            by = c('spawn_node' = 'end_loc')) %>%
  rename(tag_path = path)

# check if any user definied keep_obs lead to invalid paths
error_tags = filter_obs %>%
  left_join(tag_path) %>%
  rowwise() %>%
  mutate(node_in_path = str_detect(tag_path, node)) %>%
  ungroup() %>%
  filter(!node_in_path) %>%
  select(tag_code) %>%
  distinct()

nrow(error_tags)
error_tags %>%
  left_join(yn_df) %>%
  as.data.frame()


prepped_ch %>%
  select(-user_keep_obs) %>%
  left_join(yn_df %>%
              select(tag_code:max_det,
                     user_keep_obs)) %>%
  tabyl(auto_keep_obs,
        user_keep_obs)

prepped_ch = prepped_ch %>%
  select(-user_keep_obs) %>%
  left_join(yn_df %>%
              select(tag_code:max_det,
                     user_keep_obs))

save(parent_child, configuration, start_date, bio_df, prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('UC_Steelhead_', yr, '.rda')))



#-----------------------------------------------------------------
# tag summaries
#-----------------------------------------------------------------
# use auto_keep_obs for the moment
tag_summ = summarizeTagData(prepped_ch %>%
                              mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                                             auto_keep_obs,
                                                             user_keep_obs)),
                            bio_df)

# any duplicated tags?
sum(duplicated(tag_summ$tag_code))
# tag_summ %>%
#   filter(tag_code %in% tag_code[duplicated(tag_code)]) %>%
#   as.data.frame()

# where are tags assigned?
janitor::tabyl(tag_summ, spawn_node) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# preliminary estimate of node efficiency
node_eff = prepped_ch %>%
  mutate(user_keep_obs = auto_keep_obs) %>%
  filter(user_keep_obs) %>%
  estNodeEff(node_order = buildNodeOrder(addParentChildNodes(parent_child, configuration)))

node_eff %>%
  filter(tags_at_node > 0,
         eff_est < 1)


#-----------------------------------------------------------------
# examine some of the output
#-----------------------------------------------------------------
# which tags have "strange" capture histories?
prepped_ch %>%
  summarise(n_tags = n_distinct(tag_code),
            n_weird = n_distinct(tag_code[direction == "unknown"]),
            n_fix = n_distinct(tag_code[is.na(user_keep_obs)]),
            prop_weird = n_weird / n_tags,
            prop_fix = n_fix / n_tags)


# look at which branch each tag was assigned to for spawning
brnch_df <- addParentChildNodes(parent_child,
                                configuration) %>%
  buildNodeOrder() %>%
  mutate(branch_nm = if_else(node == "PRO",
                             "Start",
                             if_else(str_detect(path, "TOP"),
                                     "Toppenish",
                                     if_else(str_detect(path, "SUN") & str_detect(path, "ROZ", negate = T),
                                             "Naches",
                                             if_else(str_detect(path, "ROZ"),
                                                     "UpperYakima",
                                                     if_else(str_detect(path, "SAT"),
                                                             "Satus",
                                                             "Downstream")))))) %>%
  mutate(across(branch_nm,
                as.factor))

tag_summ %<>%
  left_join(brnch_df %>%
              select(spawn_node = node,
                     branch_nm))

# how many tags in each branch?
tag_summ %>%
  janitor::tabyl(branch_nm) %>%
  janitor::adorn_pct_formatting() %>%
  arrange(desc(n))
