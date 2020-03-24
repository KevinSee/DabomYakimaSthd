# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 2/21/20
# Last Modified: 3/24/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(readxl)
library(magrittr)

#-----------------------------------------------------------------
# load configuration and site_df data
load('analysis/data/derived_data/site_config.rda')

# which spawn year are we dealing with?
yr = 2012

# start date is July 1 of the previous year
start_date = paste0(yr - 1, '0701')

# build parent-child table
parent_child = createParentChildDf(site_df,
                                   configuration,
                                   startDate = start_date)

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
# observations = read_csv(paste0('analysis/data/raw_data/PTAGIS/PTAGIS_', yr-1, '_', str_sub(as.character(yr), 3, 4), '.csv'))
observations = read_csv(paste0('analysis/data/raw_data/PTAGIS/PTAGIS_', yr, '.csv'))

# deal with some double tagged fish
if(yr %in% 2012:2014) {
  bio_df = read_rds('analysis/data/derived_data/Bio_2012_14.rds')
  dbl_tag = bio_df %>%
    filter(Year == yr,
           !is.na(AltTag))

  observations %<>%
    left_join(dbl_tag %>%
                select(`Tag Code` = AltTag,
                       TagID)) %>%
    # filter(!is.na(TagID))
    mutate(`Tag Code` = if_else(!is.na(TagID),
                                TagID,
                                `Tag Code`)) %>%
    select(-TagID)
}


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
# Read in the file returned by the Yakima
#-----------------------------------------------------------------
load(paste0('analysis/data/derived_data/PITcleanr/PRO_Steelhead_', yr, '.rda'))

proc_ch = read_excel(paste0('analysis/data/raw_data/YakimaNation/PRO_Steelhead_', yr, '.xlsx')) %>%
  mutate_at(vars(AutoProcStatus:ValidPath),
            list(as.logical)) %>%
  filter(UserProcStatus)

proc_list$ProcCapHist = proc_ch

# re-save some stuff
save(yr, start_date, parent_child, proc_list,
     file = paste0('analysis/data/derived_data/PITcleanr/PRO_Steelhead_', yr, '.rda'))


#-----------------------------------------------------------------
# tag summaries
#-----------------------------------------------------------------
bio_df = read_rds('analysis/data/derived_data/Bio_2018_19.rds')

bio_df = read_rds('analysis/data/derived_data/Bio_2012_14.rds') %>%
  filter(Year == yr)

# Fix UserProcStatus, and summarise tag data
tag_summ = proc_list$ProcCapHist %>%
  filter(AutoProcStatus) %>%
  mutate(UserProcStatus = AutoProcStatus) %>%
  summariseTagData(trap_data = bio_df)

# any duplicated tags?
tag_summ %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  as.data.frame()

# where did hatchery fish go?
tag_summ %>%
  filter(SppCode == 'hsth')

# save as csv file to send out
tag_summ %>%
  select(-BranchNum, -Group) %>%
  write_csv(paste0('outgoing/tag_summary/PRO_Steelhead_TagSummary_', yr, '.csv'))


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


tag_summ %>%
  mutate(Population = Group,
         Population = if_else(AssignSpawnSite == 'LWC' | grepl('ROZ', TagPath),
                              'Upper Yakima',
                              Population),
         Population = if_else(AssignSpawnSite %in% c('LNR', 'AH1'),
                              'Naches',
                              Population),
         Population = if_else(is.na(Group),
                              'PRO_bb',
                              Population)) %>%
  filter(Population %in% c('Status', 'Naches', 'Toppenish', 'Upper Yakima')) %>%
  ggplot(aes(x = PassDate,
             color = Population,
             fill = Population)) +
  # geom_density(alpha = 0.2) +
  geom_histogram() +
  facet_wrap(~ Population) +
  theme_bw() +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1')

# ggsave('outgoing/figures/RunTiming_2019.pdf',
#        width = 6,
#        height = 6)
