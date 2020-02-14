# Author: Kevin See
# Purpose: Develop configuration file for DABOM
# Created: 2/14/20
# Last Modified: 2/14/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(magrittr)

#-----------------------------------------------------------------
# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  mutate(Node = if_else(SiteID %in% c('PRO'),
                        'PRO',
                        Node),
         Node = if_else(SiteID %in% c("NFTEAN", "TEANAR", "TEANM", "TEANWF"),
                       "TEAN",
                       Node),
         Node = if_else(SiteID == 'ROZ',
                        if_else(AntennaID %in% c('01', '02', '03'),
                                Node,
                                as.character(NA)),
                        Node),
         Node = if_else(SiteID == 'TAN' & ConfigID %in% c(120, 130),
                        "TANB0",
                        Node),
         Node = if_else(SiteID %in% c('MC1', 'MC2', 'MCJ', 'MCN'),
                       'MCN',
                       Node),
         Node = if_else(SiteID == 'ICH',
                       'ICHB0',
                       Node),
         Node = if_else(grepl('522\\.', RKM) & RKMTotal > 538,
                       'ICHA0',
                       Node),
         Node = if_else(SiteID == 'JD1',
                       'JD1B0',
                       Node),
         Node = if_else(SiteID %in% c('30M', 'BR0', 'JDM', 'SJ1', 'SJ2', 'MJ1'),
                       'JD1A0',
                       Node),
         Node = if_else(SiteID != 'JD1' & as.integer(stringr::str_split(RKM, '\\.', simplify = T)[,1]) < 351,
                       'BelowJD1',
                       Node),
         Node = if_else(SiteID == 'PRA',
                        'PRAB0',
                        Node),
         Node = if_else(SiteID != 'PRA' & as.integer(stringr::str_split(RKM, '\\.', simplify = T)[,1]) >= 639,
                        'PRAA0',
                        Node))


# Node network for DABOM
site_df = writePRONodeNetwork()

# Save file.
save(configuration, site_df, file = 'analysis/data/derived_data/site_config.rda')


#-----------------------------------------------------------------
# which sites are in site_df, but not in the PTAGIAS configuration file?
site_df %>%
  filter(!(SiteID %in% configuration$SiteID |
             SiteID %in% configuration$Node))

configuration %>%
# org_config %>%
  filter(grepl("PRA", SiteID)) %>%
  # filter(ConfigID == 130) %>%
  # as.data.frame()
  select(SiteID, ConfigID, SiteName, Node, AntennaID, AntennaGroup, SiteDescription)


#-----------------------------------------------------------------
yr = 2019
# start date is July 1 of the previous year
startDate = paste0(yr - 1, '0701')

# build parent-child table
parent_child = createParentChildDf(site_df,
                                   configuration,
                                   startDate = startDate)

parent_child %>%
  filter(ParentNode == 'PRO')
parent_child %>%
  filter(grepl('ROZ', ChildNode))

parent_child %>%
  filter(grepl('539\\.', RKM))
