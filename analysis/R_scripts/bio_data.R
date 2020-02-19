# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 2/19/20
# Last Modified: 2/19/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(readxl)
library(janitor)

#-----------------------------------------------------------------
# read in biological data from trap
bio_df = read_excel('analysis/data/raw_data/YakimaNation/Denil 2018_19.xlsx') %>%
  filter(!is.na(LadCode))

# pull out PIT tag numbers
tags = bio_df %>%
  filter(SppCode == 'wsth') %>%
  filter(!is.na(PitTag)) %>%
  select(PitTag)

# save tags to upload to PTAGIS
write_delim(tags,
            path = 'analysis/data/raw_data/tag_lists/Tags_2018_19.txt',
            delim = '\n',
            col_names = F)
