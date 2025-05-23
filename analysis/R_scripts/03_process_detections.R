# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 2/21/20
# Last Modified: 5/16/25
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
yr = 2024

# for(yr in c(2012:2022)) {
  cat(paste("Working on", yr, "\n"))

#-----------------------------------------------------------------
# start date is July 1 of previous year
start_date = paste0(yr - 1, '0701')
# when is the last possible observation date?
# we chose one entire year from the start date
max_obs_date = as.character(ymd(start_date) + years(1) - days(1))

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
ptagis_file = here("analysis/data/raw_data/PTAGIS",
                   # paste0("PRO_Sthd_", yr, "_CTH.csv"))
                   paste0("PTAGIS_", yr, ".csv"))

# recode the PTAGIS observations of double tagged fish so that the tag code matches the TagID (not TagOther)
ptagis_obs = readCTH(ptagis_file)

# any orphaned or disowned tags?
qcTagHistory(ptagis_obs,
             ignore_event_vs_release = T)

# compress and process those observations with PITcleanr
prepped_ch = PITcleanr::prepWrapper(cth_file = ptagis_obs,
                                    file_type = "PTAGIS",
                                    configuration = configuration,
                                    parent_child = parent_child %>%
                                      addParentChildNodes(configuration = configuration),
                                    min_obs_date = start_date,
                                    max_obs_date = max_obs_date,
                                    ignore_event_vs_release = F,
                                    filter_orphan_disown_tags = FALSE,
                                    add_tag_detects = T,
                                    save_file = T,
                                    file_name = here('outgoing/PITcleanr',
                                                     paste0('PRO_Steelhead_', yr, '.xlsx')))



# load and filter biological data by spawn year
# bio_df = read_rds(here('analysis/data/derived_data/Bio_Data_2012_2020.rds')) %>%
# bio_df = read_rds(here('analysis/data/derived_data/Bio_Data_2012_2022.rds')) |>
bio_df = read_rds(here('analysis/data/derived_data/Bio_Data_2012_2024.rds')) |>
  filter(spawn_year == yr)


# save some stuff
save(parent_child, configuration, start_date, bio_df, prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('PRO_Steelhead_', yr, '.rda')))

rm(prepped_ch,
   bio_df,
   start_date,
   ptagis_obs,
   ptagis_file)

# }



#-------------------------------------------
# NEXT STEPS
#-------------------------------------------
# open that Excel file, and filter on the column user_keep_obs, looking for blanks. Fill in each row with TRUE or FALSE, depending on whether that observation should be kept or not. The column auto_keep_obs provides a suggestion, but the biologist's best expert judgment should be used based on detection dates, detection locations before and after, etc.

load(here('analysis/data/derived_data/PITcleanr',
          paste0('PRO_Steelhead_', yr, '.rda')))

# how many tags to look at and "fix"
prepped_ch |>
  summarize(n_tags = n_distinct(tag_code),
            n_fix = n_distinct(tag_code[is.na(user_keep_obs)]),
            n_weird = n_distinct(tag_code[direction == "unknown"]),
            perc_fix = n_fix / n_tags,
            n_weird = n_weird / n_tags)

# an example
prepped_ch |>
  filter(is.na(user_keep_obs)) |>
  select(tag_code) |>
  distinct() |>
  slice_sample(n = 1) |>
  left_join(prepped_ch)



# read in the updated version of the PITcleanr output Excel file
yn_df <-
  read_excel(here('analysis',
                  "data",
                  "raw_data",
                  "YakimaNation",
                  # paste0('PRO_Steelhead_', yr, '.xlsx'))) |>
                  paste0('PRO_Steelhead_', yr, '_KS.xlsx'))) |>
  mutate(across(c(duration,
                  travel_time),
                as.numeric),
         across(c(duration,
                  travel_time),
                ~ as.difftime(., units = "mins")))

identical(dim(prepped_ch),
          dim(yn_df))


if(!"node_order" %in% names(yn_df)) {
  no <- parent_child |>
    addParentChildNodes(configuration) |>
    buildNodeOrder()

  yn_df <-
    yn_df |>
    left_join(no)
}


filter_obs <-
  yn_df %>%
  # prepped_ch |>
  mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                 auto_keep_obs,
                                 user_keep_obs)) %>%
  filter(user_keep_obs)

# construct all valid paths
all_paths = buildPaths(addParentChildNodes(parent_child,
                                           configuration))

tag_path <-
  estimateFinalLoc(filter_obs) |>
  select(tag_code,
         final_node) %>%
  distinct() %>%
  left_join(all_paths,
            by = join_by(final_node == end_loc)) %>%
  rename(tag_path = path)

# check if any user defined keep_obs lead to invalid paths
error_tags <-
  filter_obs %>%
  left_join(tag_path,
            by = join_by(tag_code)) %>%
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

prepped_ch <-
  prepped_ch %>%
  select(-user_keep_obs) %>%
  left_join(yn_df %>%
              select(tag_code:max_det,
                     user_keep_obs))

save(parent_child, configuration, start_date, bio_df, prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('PRO_Steelhead_', yr, '.rda')))



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
janitor::tabyl(tag_summ, final_node) %>%
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
  mutate(branch_nm = case_when(node == "PRO" ~ "Start",
                               str_detect(path, "TOP") ~ "Toppenish",
                               str_detect(path, "SUN") &
                                 str_detect(path, "ROZ", negate = T) ~ "Naches",
                               str_detect(path, "ROZ") ~ "UpperYakima",
                               str_detect(path, "SAT") ~ "Satus",
                               .default = "Downstream")) |>
  mutate(across(branch_nm,
                as.factor))

tag_summ <-
  tag_summ |>
  left_join(brnch_df %>%
              select(final_node = node,
                     branch_nm),
            by = join_by(final_node))

# how many tags in each branch?
tag_summ %>%
  janitor::tabyl(branch_nm) %>%
  janitor::adorn_pct_formatting() %>%
  arrange(desc(n))
