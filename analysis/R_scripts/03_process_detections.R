# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 2/21/20
# Last Modified: 2/21/20
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
startDate = paste0(yr - 1, '0701')

# build parent-child table
parent_child = createParentChildDf(site_df,
                                   configuration,
                                   startDate = startDate)

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
observations = read_csv(paste0('analysis/data/raw_data/PTAGIS/PTAGIS_', yr-1, '_', str_sub(as.character(yr), 3, 4), '.csv'))

#-----------------------------------------------------------------
# need a new processCapHist for Yakima. The code below is the begining of that function, which will be moved to the PITcleanr package

tagging_site = c('PRO')

# construct valid paths
cat('Constructing valid paths.\n')
valid_paths = getValidPaths(parent_child)

# create dataframe describing node order
cat('Creating node order')
node_order = createNodeOrder(valid_paths = valid_paths,
                             configuration = configuration,
                             site_df = site_df,
                             step_num = 2)

# this could be used to create different matrices to feed to DABOM
node_order = node_order %>%
  select(-Group) %>%
  left_join(tibble(BranchNum = 1:8,
                   Group = c(rep("Downstream", 5),
                             'Status',
                             'Toppenish',
                             'Sunnyside')))

# pull out tag ID and trap date at Prosser
cat('Getting trap date.\n')
valid_tag_df = observations %>%
  filter(`Event Site Code Value` %in% tagging_site) %>%
  mutate_at(vars(`Event Date Time Value`, `Event Release Date Time Value`),
            list(lubridate::mdy_hms)) %>%
  mutate(ObsDate = if_else(!is.na(`Event Release Date Time Value`) &
                             is.na(`Antenna ID`),
                           `Event Release Date Time Value`,
                           `Event Date Time Value`)) %>%
  filter(ObsDate >= lubridate::ymd(startDate)) %>%
  group_by(TagID = `Tag Code`) %>%
  summarise(TrapDate = min(lubridate::floor_date(ObsDate,
                                                 unit = 'days'))) %>%
  ungroup()

# translate in nodes and simplify consecutive hits on the same node
cat('Assigning nodes.\n')
valid_obs = assignNodes(valid_tag_df,
                        observations,
                        configuration,
                        parent_child,
                        truncate = T)


# drop detecions of BelowJD1 if other detections exist
# this eliminates detections of kelts moving down through the mainstem Columbia
onlyBelowJD1_tags = anti_join(valid_obs,
                              valid_obs %>%
                                group_by(TagID) %>%
                                filter(!Node %in% c('BelowJD1', tagging_site)) %>%
                                ungroup() %>%
                                select(TagID) %>%
                                distinct()) %>%
  pull(TagID) %>%
  unique()

valid_obs = valid_obs %>%
  filter(TagID %in% onlyBelowJD1_tags |
           (!TagID %in% onlyBelowJD1_tags & Node != 'BelowJD1'))

# check if any tags have been dropped incorrectly along the way
if(n_distinct(valid_obs$TagID) != n_distinct(observations$`Tag Code`)) {
  cat('Error: some tags being dropped')
  # warning('Error: some tags being dropped')
}


cat('Processing assigned nodes\n')
save_df = writeCapHistOutput(valid_obs,
                             valid_paths,
                             node_order,
                             last_obs_date = NULL,
                             save_file = F)

# which tags have "strange" capture histories?
weird_tags = save_df %>%
  filter(UserProcStatus == '') %>%
  pull(TagID) %>%
  unique()

length(weird_tags)
length(weird_tags) / n_distinct(save_df$TagID)

save_df %>%
  filter(TagID %in% weird_tags) %>%
  select(TagID:AutoProcStatus) %>%
  as.data.frame()

save_df %>%
  filter(TagID == weird_tags$TagID[[1]])

configuration %>%
  filter(SiteID == 'BCC') %>%
  as.data.frame()

tag_id = "3DD.00773715CE"
observations %>%
  filter(`Tag Code` == tag_id) %>%
  as.data.frame()

save_df %>%
  filter(TagID == tag_id) %>%
  as.data.frame()
