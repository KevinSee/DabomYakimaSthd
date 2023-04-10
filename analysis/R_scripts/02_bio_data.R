# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 2/19/20
# Last Modified: 4/10/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(here)


#-----------------------------------------------------------------
# for spawn year 2020
# read in data from 2019-2020
bio_2020 = read_excel(here('analysis/data/raw_data',
                         'YakimaNation/2019_2020_sthd_denil_PITtag.xlsx'),
                    2,
                    skip = 1,
                    col_names = 'tag_code') %>%
  distinct() %>%
  left_join(read_excel(here('analysis/data/raw_data/YakimaNation',
                            '2019_2020_sthd_denil_PITtag.xlsx')) %>%
              clean_names() %>%
              mutate(tag_code = if_else(!is.na(pit_tag),
                                        pit_tag,
                                        jv_pit_tag)),
            by = join_by(tag_code),
            multiple = "all") %>%
  mutate(pass_date_time = ymd_hms(paste(str_split(pass_date, " ", simplify = T)[,1],
                                        str_split(pass_time, " ", simplify = T)[,2]))) %>%
  relocate(pass_date_time,
           .after = "pass_time") %>%
  select(-pass_date,
         -pass_time) %>%
  # fix a few tag codes
  mutate(across(tag_code,
                ~ str_replace(.,
                              "\\.\\.",
                              "\\."))) %>%
  left_join(read_csv(here("analysis/data/raw_data/YakimaNation",
                          "Tags_Not_In_PTAGIS_cf_corrections.csv")) %>%
              select(tag_code = TagID,
                     new_tag_code = `Corrected PITtag ID`)) %>%
  mutate(tag_code = if_else(!is.na(new_tag_code),
                            new_tag_code,
                            tag_code)) %>%
  select(-new_tag_code) %>%
  # filter out duplicate tags by keeping the first record
  arrange(pass_date_time, tag_code) %>%
  group_by(tag_code) %>%
  slice(1) %>%
  ungroup() %>%
  add_column(spawn_year = 2020,
             .before = 0)

# for spawn year 2019
# read in biological data from trap
bio_2019 = read_excel(here('analysis/data/raw_data/YakimaNation',
                         'Denil 2018_19.xlsx')) %>%
  clean_names() %>%
  rename(tag_code = pit_tag) %>%
  mutate(across(pass_time,
                as.numeric)) %>%
  filter(!is.na(lad_code)) %>%
  filter(!is.na(tag_code)) %>%
  # fix one tag code
  mutate(tag_code = if_else(tag_code == "389.1C2E70563A",
                            "3D9.1C2D70563A",
                            tag_code)) %>%
  # filter out duplicate tags by keeping the first record
  arrange(pass_date, tag_code) %>%
  group_by(tag_code) %>%
  slice(1) %>%
  ungroup() %>%
  add_column(spawn_year = 2019,
             .before = 0) %>%
  select(any_of(names(bio_2020)),
         everything())

# for older years
bio_old = 2012:2014 %>%
  as.list %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
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
             select(pass_date = ProsserDate,
                    tag_code = TagID,
                    alt_tag,
                    spp_code = SppCode,
                    status,
                    age,
                    sex,
                    forklgth,
                    mehlgth,
                    weight,
                    ad_clip = AdiposeClip,
                    comments) %>%
             arrange(pass_date,
                     tag_code) %>%
             group_by(tag_code) %>%
             slice(1) %>%
             ungroup()
         }) %>%
  mutate(
    across(
      spawn_year,
      as.numeric))

# put all years biological data together
bio_df <- bio_2020 %>%
  bind_rows(bio_2019) %>%
  bind_rows(bio_old) %>%
  arrange(spawn_year,
          pass_date,
          tag_code)


#-----------------------------------------------------------------
# for tag lists
#-----------------------------------------------------------------
# put bounds around years
min_yr = min(bio_df$spawn_year)
max_yr = max(bio_df$spawn_year)


# pull out PIT tag numbers
tag_list = bio_df %>%
  split(list(.$spawn_year)) %>%
  map(.f = function(x) {
    x %>%
      pivot_longer(cols = starts_with("tag"),
                   names_to = "source",
                   values_to = "tag_code") %>%
      filter(!is.na(tag_code)) %>%
      select(tag_code)
  })

# save tags to upload to PTAGIS

# for(yr in names(tag_list)) {
#   write_delim(tag_list[[yr]],
#               file = here('analysis/data/raw_data/tag_lists',
#                           paste0('Prosser_Tags_', yr, '.txt')),
#               delim = '\n',
#               col_names = F)
# }

# just write the latest year
write_delim(tag_list[[as.character(max_yr)]],
            file = here('analysis/data/raw_data/tag_lists',
                        paste0('Prosser_Tags_', max_yr, '.txt')),
            delim = '\n',
            col_names = F)

# save biological data for later
write_rds(bio_df,
          file = here('analysis/data/derived_data',
                      paste0('Bio_Data_', min_yr, '_', max_yr, '.rds')))
