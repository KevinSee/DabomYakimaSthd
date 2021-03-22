# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 2/19/20
# Last Modified: 3/22/21
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(janitor)
library(magrittr)


#-----------------------------------------------------------------
# for spawn year 2020
# read in data from 2019-2020
bio_df = read_excel('analysis/data/raw_data/YakimaNation/2019_2020_sthd_denil_PITtag.xlsx') %>%
  rename(TagID = PitTag) %>%
  mutate(across(PassTime,
                as.numeric)) %>%
  filter(!is.na(LadCode)) %>%
  filter(!is.na(TagID)) %>%
  # # fix one tag code
  # mutate(TagID = if_else(TagID == "389.1C2E70563A",
  #                        "3D9.1C2D70563A",
  #                        TagID)) %>%
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
            path = 'analysis/data/raw_data/tag_lists/Tags_2020.txt',
            delim = '\n',
            col_names = F)

# save biological data for later
write_rds(bio_df,
          path = 'analysis/data/derived_data/Bio_2020.rds')



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

#-----------------------------------------------------------------
# for older years
bio_df = 2012:2014 %>%
  as.list %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(yr) {

           df = read_excel(paste0('analysis/data/raw_data/YakimaNation/Prosser_steelhead_', yr, '.xlsx'))

           if("WDFW Age" %in% names(df)) {
             df %<>%
               rename(age = `WDFW Age`)
           }

           if("ProsserDate...8" %in% names(df)) {
             df %<>%
               rename(ProsserDate = `ProsserDate...8`)
           }

           df %>%
             mutate(TagID = if_else(!is.na(pittag),
                                    pittag,
                                    jvpittag)) %>%
             rename(SppCode = SpeciesCode) %>%
             mutate(alt_tag = if_else(!is.na(pittag) & !is.na(jvpittag),
                                      jvpittag,
                                      as.character(NA))) %>%
             select(PassDate = ProsserDate,
                    TagID,
                    alt_tag,
                    SppCode,
                    Status = status,
                    age,
                    sex,
                    forklgth,
                    mehlgth,
                    weight,
                    AdiposeClip,
                    comments) %>%
             clean_names(case = 'upper_camel') %>%
             rename(TagID = TagId)
         })

tabyl(bio_df, Year)
xtabs(~ Year + is.na(AltTag), bio_df)

bio_df %>%
  filter(Year != 2012,
         !is.na(AltTag))

# pull out PIT tag numbers
tags = bio_df %>%
  split(list(.$Year)) %>%
  map(.f = function(x) {
    x %>%
      select(TagID) %>%
      bind_rows(x %>%
                  filter(!is.na(AltTag)) %>%
                  select(TagID = AltTag))
  })

# save tags to upload to PTAGIS
for(yr in names(tags)) {
  write_delim(tags[[yr]],
              path = paste0('analysis/data/raw_data/tag_lists/Tags_', yr, '.txt'),
              delim = '\n',
              col_names = F)
}

# save biological data for later
write_rds(bio_df,
          path = 'analysis/data/derived_data/Bio_2012_14.rds')
