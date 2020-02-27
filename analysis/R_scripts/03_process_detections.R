# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 2/21/20
# Last Modified: 2/27/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
# library(WriteXLS)

#-----------------------------------------------------------------
# load configuration and site_df data
load('analysis/data/derived_data/site_config.rda')

# which spawn year are we dealing with?
yr = 2019
# start date is July 1 of the previous year
start_date = paste0(yr - 1, '0701')

# build parent-child table
parent_child = createParentChildDf(site_df,
                                   configuration,
                                   startDate = start_date)

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
observations = read_csv(paste0('analysis/data/raw_data/PTAGIS/PTAGIS_', yr-1, '_', str_sub(as.character(yr), 3, 4), '.csv'))

# process those observations with PITcleanr, using Yakima-specific function
proc_list = processCapHist_PRO(start_date,
                               configuration = configuration,
                               parent_child = parent_child,
                               observations = observations,
                               # for earlier years, may want to use the code below
                               last_obs_date = format(lubridate::ymd(start_date) + lubridate::years(1), "%Y%m%d"),
                               site_df = site_df,
                               save_file = T,
                               file_name = paste0('outgoing/PITcleanr/PRO_Steelhead_', yr, '.xlsx'))

# save some stuff
save(yr, start_date, parent_child, proc_list,
     file = paste0('analysis/data/derived_data/PITcleanr/PRO_Steelhead_', yr, '.rda'))

#-----------------------------------------------------------------
# tag summaries
#-----------------------------------------------------------------
bio_df = read_rds('analysis/data/derived_data/Bio_2018_19.rds') %>%
  filter(TagID %in% proc_list$ProcCapHist$TagID)

# Fix UserProcStatus, and summarise tag data
tag_summ = proc_list$ProcCapHist %>%
  filter(AutoProcStatus) %>%
  mutate(UserProcStatus = AutoProcStatus) %>%
  summariseTagData(trap_data = bio_df)

# any duplicated tags?
tag_summ %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  as.data.frame()

# where are tags assigned?
janitor::tabyl(tag_summ, AssignSpawnSite) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# which branch are tags assigned to?
tag_summ %>%
  mutate(Branch = fct_explicit_na(Group,
                                 'PRO_bb')) %>%
  janitor::tabyl(Branch) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# preliminary estimate of node efficiency
node_eff = proc_list$ProcCapHist %>%
  filter(AutoProcStatus) %>%
  mutate(UserProcStatus = AutoProcStatus) %>%
  estNodeEff(node_order = proc_list$NodeOrder)

node_eff %>%
  filter(tagsAtNode > 0,
         detEff < 1)

node_eff %>%
  xtabs(~ (!is.na(detEff)) + (detEff_SE > 0), .)

node_eff %>%
  filter(!is.na(detEff),
         detEff_SE > 0)

node_eff %>%
  filter(grepl('^TOP', Node))

#-----------------------------------------------------------------
# examine some of the output
#-----------------------------------------------------------------
proc_ch = proc_list$ProcCapHist

# which tags have "strange" capture histories?
weird_tags = proc_ch %>%
  filter(UserProcStatus == '') %>%
  pull(TagID) %>%
  unique()

length(weird_tags)
length(weird_tags) / n_distinct(proc_ch$TagID)

proc_ch %>%
  filter(TagID %in% weird_tags) %>%
  select(TagID:AutoProcStatus) %>%
  as.data.frame()

proc_ch %>%
  filter(TagID == weird_tags[[1]])

configuration %>%
  filter(SiteID == 'BCC') %>%
  as.data.frame()

tag_id = "3DD.00773715CE"
observations %>%
  filter(`Tag Code` == tag_id) %>%
  as.data.frame()

proc_ch %>%
  filter(TagID == tag_id) %>%
  as.data.frame()
