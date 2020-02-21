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
