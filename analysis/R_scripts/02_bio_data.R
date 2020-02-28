# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 2/19/20
# Last Modified: 2/27/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(janitor)

#-----------------------------------------------------------------
# read in biological data from trap
bio_df = read_excel('analysis/data/raw_data/YakimaNation/Denil 2018_19.xlsx') %>%
  rename(TagID = PitTag) %>%
  mutate_at(vars(PassTime),
            list(as.numeric)) %>%
  filter(!is.na(LadCode)) %>%
  filter(!is.na(TagID)) %>%
  # fix one tag code
  mutate(TagID = if_else(TagID == "389.1C2E70563A",
                         "3D9.1C2D70563A",
                         TagID)) %>%
  # filter out duplicate tags by keeping the first record
  arrange(PassDate, TagID) %>%
  group_by(TagID) %>%
  slice(1) %>%
  ungroup()

# pull out PIT tag numbers
tags = bio_df %>%
  # filter(SppCode == 'wsth') %>%
  select(TagID)

# save tags to upload to PTAGIS
write_delim(tags,
            path = 'analysis/data/raw_data/tag_lists/Tags_2018_19.txt',
            delim = '\n',
            col_names = F)

# save biological data for later
write_rds(bio_df,
          path = 'analysis/data/derived_data/Bio_2018_19.rds')
