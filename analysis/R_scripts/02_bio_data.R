# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 2/19/20
# Last Modified: 4/10/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(here)

#-----------------------------------------------------------------

# read data directly from trap database exported from Yakima Nation
trap_data <-
  read_excel(here("analysis/data/raw_data",
                  "YakimaNation",
                  "prossertrap_sth.xlsx"),
             1,
             col_types = c(rep("guess", 8),
                           rep("text", 8),
                           rep("numeric", 4))) |>
  clean_names() |>
  mutate(trap_date =
           case_when(!is.na(pass_time) ~ pass_date +
                       hours(hour(pass_time)) +
                       minutes(minute(pass_time)) +
                       seconds(second(pass_time)),
                     is.na(pass_time) ~ pass_date,
                     .default = NA_Date_),
         pit_tag = case_when(!is.na(tagged) ~ tagged,
                              is.na(tagged) &
                                !is.na(jv_pit_tag) ~ jv_pit_tag,
                              .default = NA_character_),
         across(pit_tag,
                ~ case_when(str_detect(., "\\.", negate = T) ~ paste0(str_sub(., 1, 3),
                                                                      ".",
                                                                      str_sub(., 4)),
                            .default = .))) |>
  mutate(spawn_year = case_when(month(trap_date) >= 7 ~ year(trap_date) + 1,
                                month(trap_date) < 7 ~ year(trap_date),
                                .default = NA_real_)) |>
  relocate(trap_date,
           spawn_year,
           .before = "pass_date") |>
  select(-starts_with("pass_"),
         -sy) |>
  mutate(across(age,
                ~ case_when(nchar(.) > 7 ~ as.character(round(as.numeric(.), 1)),
                            .default = .)),
         across(origin,
                ~ case_when(spp_code == "hsth" ~ "H",
                            spp_code == "wsth" ~ "W",
                            .default = NA_character_)))
trap_data <-
  trap_data |>
  filter(spawn_year < year(today()))

# fix a few tag codes
# tags with the wrong number of characters in tag code
# trap_data |>
#   mutate(n_chr = nchar(pit_tag)) |>
#   filter(n_chr != 14) |>
#   select(trap_date,
#          pit_tag)
#
# trap_data |>
#   select(spawn_year,
#          trap_date,
#          pit_tag) |>
#   mutate(n_chr = nchar(pit_tag),
#          tag_code1 = str_sub(pit_tag, 1, 8),
#          tag_code2 = str_sub(pit_tag, 9)) |>
#   filter(#spawn_year == 2019,
#          # str_detect(tag_code1, "3D9.1C2"),
#          # str_detect(tag_code1, "D$", negate = T))
#          n_chr < 14)

trap_data <-
  trap_data |>
  mutate(n_chr = nchar(pit_tag)) |>
  mutate(across(pit_tag,
                ~ case_when(n_chr == 14 ~ .,
                            spawn_year == 2019 &
                              n_chr == 13 &
                              str_detect(pit_tag, "^3D9.1C2") ~ paste0(str_sub(., 1, 8),
                                                                        "D",
                                                                        str_sub(., 9)),
                            spawn_year == 2019 &
                              n_chr == 13 &
                              str_detect(pit_tag, "C2D75CCE7$") ~ paste0(str_sub(., 1, 4),
                                                                       "1",
                                                                       str_sub(., 5)),
                            .default = .))) |>
  select(-n_chr)

# pull out PIT tag numbers
tag_codes <-
  trap_data |>
  select(spawn_year,
         contains("tag")) |>
  nest(tags = contains("tag"),
       .by = spawn_year) |>
  mutate(tag_list = map(tags,
                        .f = function(x) {
                          x |>
                            pivot_longer(cols = contains("tag"),
                                         names_to = "source",
                                         values_to = "tag_code") |>
                            filter(!is.na(tag_code)) |>
                            select(tag_code) |>
                            distinct()
                          })) |>
  select(-tags)

# save tags to upload to PTAGIS
for(yr in tag_codes$spawn_year) {
  tag_codes |>
    filter(spawn_year == yr) |>
    pull(tag_list) |>
    magrittr::extract2(1) |>
    write_delim(file = here('analysis',
                            "data",
                            "raw_data",
                            "tag_lists",
                            paste0('Prosser_Sthd_Tags_', yr, '.txt')),
                delim = '\n',
                col_names = F)
}

bio_data <-
  trap_data |>
  select(spawn_year,
         trap_date,
         spp_code,
         pit_tag,
         origin,
         age,
         sex,
         forklgth,
         weight,
         status)

# save biological data for later
yr_range = range(bio_data$spawn_year)
write_rds(bio_data,
          file = here('analysis/data/derived_data',
                      paste0('Bio_Data_', yr_range[1], '_', yr_range[2], '.rds')))

#-----------------------------------------------------------------
# tag_list_folder <-
#   here('analysis',
#        "data",
#        "raw_data",
#        "tag_lists")
#
# old_tags <-
#   tibble(file_nm = list.files(tag_list_folder)) |>
#   filter(str_detect(file_nm, "Sthd", negate = T)) |>
#   mutate(spawn_year = str_extract(file_nm, "[:digit:]+"),
#          across(spawn_year,
#                 as.numeric)) |>
#   mutate(tags = map(file_nm,
#                     .f = function(x) {
#                       read_csv(str_glue("{tag_list_folder}/", x),
#                                col_names = "tag_code",
#                                show_col_types = FALSE) |>
#                         mutate(across(tag_code,
#                                       ~ str_remove(., "^'")))
#                     },
#                     .progress = T)) |>
#   unnest(tags) |>
#   select(-file_nm)
#
# old_tags |>
#   anti_join(tag_codes |>
#               unnest(tag_list))
#
# tag_codes |>
#   filter(spawn_year %in% old_tags$spawn_year) |>
#   unnest(tag_list) |>
#   anti_join(old_tags) |>
#   tabyl(spawn_year)

# tags with the wrong number of characters in tag code
trap_data |>
  mutate(n_chr = nchar(pit_tag)) |>
  filter(n_chr != 14) |>
  select(trap_date,
         tagged,
         jv_pit_tag,
         pit_tag)

# any duplicated tags in a single spawn year?
trap_data |>
  unite(col = "tag_yr",
        pit_tag,
        spawn_year,
        remove = F) |>
  filter(tag_yr %in% tag_yr[duplicated(tag_yr)]) |>
  arrange(pit_tag,
          trap_date) |>
  group_by(pit_tag) |>
  mutate(n_trap_dates = n_distinct(trap_date)) |>
  ungroup() |>
  filter(n_trap_dates == 1) |>
  pull(pit_tag)

tagging_df |>
  anti_join(trap_data |>
               select(spawn_year,
                      pit_tag) |>
               distinct()) |>
  filter(spawn_year == 2014)
  tabyl(spawn_year)

trap_data |>
  select(spawn_year,
         pit_tag) |>
  distinct() |>
  anti_join(tagging_df) |>
  mutate(n_chr = nchar(pit_tag)) |>
  filter(n_chr == 14)

trap_data |>
  select(spawn_year,
         trap_date,
         pit_tag,
         forklgth) |>
  distinct() |>
  # correct any lengths that were mistakenly entered as cm, convert to mm
  mutate(across(forklgth,
                ~ if_else(. < 100,
                          . * 10,
                          .))) |>
  inner_join(tagging_df,
             by = join_by(spawn_year,
                          pit_tag)) |>
  filter(forklgth != length) |>
  select(spawn_year,
         pit_tag,
         trap_date,
         event_date,
         release_date,
         trap_date,
         forklgth,
         length)



#-----------------------------------------------------------------
sthd_tags <-
  read_csv(here("analysis/data/raw_data",
                # "tagging_recapture",
                "Prosser_Mark_Recapture_2012_2024.csv"),
           show_col_types = F) |>
  clean_names() |>
  mutate(across(contains("_date_mm"),
                mdy),
         across(contains("date_time"),
                mdy_hms)) |>
  mutate(spawn_year = case_when(!is.na(event_release_date_mmddyyyy) &
                                  month(event_release_date_mmddyyyy) < 7 ~ year(event_release_date_mmddyyyy),
                                !is.na(event_release_date_mmddyyyy) &
                                  month(event_release_date_mmddyyyy) >= 7 ~ year(event_release_date_mmddyyyy) + 1,
                                is.na(event_release_date_mmddyyyy) &
                                  month(event_date_mmddyyyy) < 7 ~ year(event_date_mmddyyyy),
                                is.na(event_release_date_mmddyyyy) &
                                  month(event_date_mmddyyyy) >= 7 ~ year(event_date_mmddyyyy) + 1,
                                .default = NA_real_)) |>
  # grab only tags from the adult ladder
  filter(event_capture_method_code == "LADDER")

# remove tags that are spawning in a future year
sthd_tags <-
  sthd_tags |>
  filter(spawn_year < year(today()))

# pull out MRR data about all PIT tags from those MRR files
all_tags <-
  sthd_tags |>
  select(spawn_year,
         event_file_name) |>
  distinct() |>
  arrange(spawn_year,
          event_file_name) |>
  mutate(tag_file = map(event_file_name,
                        .f = function(x) {
                          out <-
                            tryCatch(queryMRRDataFile(x),
                                     error =
                                       function(cond) {
                                         message(paste("Error with file", x))
                                         message("Error message:")
                                         message(cond)
                                         return(NULL)
                                       },
                                     warning =
                                       function(cond) {
                                         message(paste("Warning with file", x))
                                         message("Warning message:")
                                         message(cond)
                                         return(NULL)
                                       })
                          if("spawn_year" %in% names(out)) {
                            out <- out |>
                              select(-spawn_year)
                          }

                          return(out)
                        },
                        .progress = T)) |>
  unnest(tag_file) |>
  mutate(year = if_else(month(event_date) < 7,
                        year(event_date),
                        year(event_date) + 1)) |>
  relocate(year,
           .before = 0) |>
  distinct() |>
  arrange(year,
          event_date) |>
  # # pull adult out steelhead tags
  # filter(life_stage == "Adult",
  #        str_detect(species_run_rear_type, "^3"),
  #        # ignore tags identified as rainbow trout
  #        str_detect(species_run_rear_type, "^30", negate = T),
  #        # ignore strange tag numbers
  #        str_detect(pit_tag, "\\.\\.\\.", negate = T)) |>
  # fix a few metrics
  # mutate(across(poh,
  #               as.numeric)) |>
  # assign sex based on conditional comments
mutate(sex = case_when(str_detect(conditional_comments, "FE") ~ "Female",
                       str_detect(conditional_comments, "MA") ~ "Male",
                       .default = NA_character_)) |>
  # determine CWT and ad-clip status from conditional comments
  mutate(cwt = if_else(str_detect(conditional_comments, "CP") |
                         str_detect(conditional_comments, "CW"),
                       T, F),
         ad_clip = case_when(str_detect(conditional_comments, "AD") ~ T,
                             str_detect(conditional_comments, "AI") ~ F,
                             .default = FALSE)) |>
  # correct any lengths that were mistakenly entered as cm, convert to mm
  mutate(across(length,
                ~ if_else(. < 100,
                          . * 10,
                          .)))

if(!"second_pit_tag" %in% names(all_tags)) {
  all_tags <-
    all_tags |>
    add_column(second_pit_tag = NA_character_,
               .after = "pit_tag")
}

# are all tags from first file accounted for?
sthd_tags$tag_code[!sthd_tags$tag_code %in% all_tags$pit_tag]
sthd_tags$tag_code[!(sthd_tags$tag_code %in% all_tags$pit_tag |
                       sthd_tags$tag_code %in% na.omit(all_tags$second_pit_tag))]


# what tags are not in the first file, but appear to be adult steelhead?
extra_tags <-
  all_tags |>
  filter(spawn_year %in% unique(sthd_tags$spawn_year)) |>
  filter(!(pit_tag %in% sthd_tags$tag_code |
             second_pit_tag %in% sthd_tags$tag_code)) |>
  filter(life_stage == "Adult",
         str_detect(species_run_rear_type, "^3"),
         # ignore tags identified as rainbow trout
         str_detect(species_run_rear_type, "^30", negate = T),
         # ignore strange tag numbers
         str_detect(pit_tag, "\\.\\.\\.", negate = T))

# pull out data for tags from first file (sthd_tags)
# and a few additional tags from extra_tags
tagging_df <-
  all_tags |>
  # filter(pit_tag %in% sthd_tags$tag_code |
  #          second_pit_tag %in% sthd_tags$tag_code |
  #          pit_tag %in% unique(extra_tags$pit_tag[extra_tags$event_type != "Recovery"])) |>
  # filter(pit_tag %in% sthd_tags2$tag_code) |>
  # put a few filters on
  filter(#life_stage == "Adult",                                # filter for adults
         str_detect(species_run_rear_type, "^3"),              # filter for O. mykiss
         str_detect(species_run_rear_type, "^30", negate = T), # ignore tags identified as rainbow trout
         str_detect(conditional_comments, "KL", negate = T),   # filter out kelts
         str_detect(pit_tag, "\\.\\.\\.", negate = T))         # ignore strange tag numbers

tagging_df |>
  count(spawn_year,
        name = "tagging_n") |>
  left_join(sthd_tags2 |>
              count(spawn_year,
                    name = "sthd_tags2_n")) |>
  mutate(diff = sthd_tags2_n - tagging_n)

tagging_df |>
  filter(spawn_year <= 2021) |>
  anti_join(sthd_tags2 |>
              select(-spawn_year),
            by = join_by(pit_tag == tag_code)) |>
  filter(spawn_year == 2013)



sthd_tags2 |>
  anti_join(tagging_df |>
              select(-spawn_year),
            by = join_by(tag_code == pit_tag)) |>
  select(tag_code) |>
  left_join(all_tags,
            by = join_by(tag_code == pit_tag))

# pull out certain columns to use going forward
bio_df <-
  tagging_df |>
  # determine origin
  mutate(origin = str_extract(species_run_rear_type, "[:alpha:]$")) |>
  # round off the hours/minutes/seconds
  mutate(across(ends_with("_date"),
                ~ floor_date(., unit = "days"))) |>
  # filter out any recovery mortalities
  filter(event_type %in% c("Mark",
                           "Recapture")) |>
  select(spawn_year,
         event_file_name,
         pit_tag,
         second_pit_tag,
         event_date,
         release_date,
         event_type,
         species_run_rear_type,
         origin,
         sex,
         cwt,
         ad_clip,
         length,
         # poh,
         conditional_comments,
         text_comments) |>
  mutate(across(release_date,
                ~ if_else(is.na(.),
                          event_date,
                          .)))



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

# add 2022 data
bio_2022 <- read_excel(here('analysis/data/raw_data/YakimaNation',
                            'Prosser sampled steelhead 2021_22 run.xlsx')) %>%
  clean_names() %>%
  rename(tag_code = pit_tag) %>%
  mutate(across(pass_time,
                as.numeric)) %>%
  filter(!is.na(lad_code)) %>%
  filter(!is.na(tag_code)) %>%
  mutate(pass_date_time = pass_date) %>%
  arrange(pass_date, tag_code) %>%
  group_by(tag_code) %>%
  slice(1) %>%
  ungroup() %>%
  add_column(spawn_year = 2022,
             .before = 0)



# put all years biological data together
bio_df <- bio_2020 %>%
  bind_rows(bio_2019) %>%
  bind_rows(bio_old) %>%
  bind_rows(bio_2022) %>%
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
