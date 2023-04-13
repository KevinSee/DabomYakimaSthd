# Author: Kevin See
# Purpose: Develop configuration file for DABOM
# Created: 2/14/20
# Last Modified: 4/10/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(sf)
library(here)

#-----------------------------------------------------------------
# set starting point
root_site = "PRO"

#-----------------------------------------------------------------
# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  mutate(node = if_else(site_code %in% c('PRO'),
                        'PRO',
                        node),
         node = if_else(site_code %in% c("NFTEAN", "TEANAR", "TEANM", "TEANWF"),
                       "LMTA0",
                       node),
         node = if_else(site_code == 'ROZ',
                        if_else(antenna_id %in% c('01', '02', '03'),
                                node,
                                as.character(NA)),
                        node),
         node = if_else(site_code %in% c('MC1', 'MC2', 'MCJ', 'MCN'),
                       'MCN',
                       node),
         node = if_else(site_code == 'ICH',
                       'ICHB0',
                       node),
         node = if_else(grepl('522\\.', rkm) & rkm_total > 538,
                       'ICHA0',
                       node),
         node = if_else(site_code == 'JD1',
                       'JD1B0',
                       node),
         node = if_else(site_code %in% c('30M', 'BR0', 'JDM', 'SJ1', 'SJ2', 'MJ1'),
                       'JD1A0',
                       node),
         node = if_else(site_code != 'JD1' & as.integer(stringr::str_split(rkm, '\\.', simplify = T)[,1]) < 351,
                       'JDA',
                       node),
         node = if_else(site_code == 'PRA',
                        'PRAB0',
                        node),
         node = if_else(site_code != 'PRA' & as.integer(stringr::str_split(rkm, '\\.', simplify = T)[,1]) >= 639,
                        'PRAA0',
                        node),
         latitude = if_else(site_code == "SWK" & is.na(latitude),
                            unique(latitude[site_code == "SWAUKC"]),
                            latitude),
         longitude = if_else(site_code == "SWK" & is.na(longitude),
                            unique(longitude[site_code == "SWAUKC"]),
                            longitude))


# Node network for DABOM

# get spatial object of sites used in model
sites_sf = writeOldNetworks()$Prosser %>%
  mutate(
    across(
      c(SiteID, Step2),
      ~ recode(.,
               "BelowJD1" = "JDA")),
    path = str_replace(path, "BelowJD1", "JDA")) %>%
  rename(site_code = SiteID) %>%
  left_join(configuration %>%
              group_by(site_code) %>%
              filter(config_id == max(config_id)) %>%
              ungroup(),
            multiple = "all") %>%
  select(site_code,
         site_name,
         site_type = site_type_name,
         type = site_type,
         rkm,
         site_description = site_description,
         latitude, longitude) %>%
  distinct() %>%
  filter(!is.na(latitude)) %>%
  st_as_sf(coords = c("longitude",
                      "latitude"),
           crs = 4326) %>%
  st_transform(crs = 5070)

# tmp <- writeOldNetworks()$Prosser
# pro_sites <- tmp$SiteID
# miss_sites <- pro_sites[!pro_sites %in% sites_sf$site_code]
# tmp |>
#   filter(SiteID %in% miss_sites)


#-----------------------------------------------------------------
# download the NHDPlus v2 flowlines
# do you want flowlines downstream of root site? Set to TRUE if you have downstream sites
dwn_flw = T
nhd_list = queryFlowlines(sites_sf = sites_sf,
                          root_site_code = root_site,
                          min_strm_order = 2,
                          dwnstrm_sites = dwn_flw,
                          dwn_min_stream_order_diff = 4)

# compile the upstream and downstream flowlines
flowlines = nhd_list$flowlines
if(dwn_flw) {
  flowlines %<>%
    rbind(nhd_list$dwn_flowlines)
}

#-----------------------------------------------------------------
# plot the flowlines and the sites
ggplot() +
  geom_sf(data = flowlines,
          aes(color = as.factor(StreamOrde),
              size = StreamOrde)) +
  scale_color_viridis_d(direction = -1,
                        option = "D",
                        name = "Stream\nOrder",
                        end = 0.8) +
  scale_size_continuous(range = c(0.2, 1.2),
                        guide = 'none') +
  geom_sf(data = nhd_list$basin,
          fill = NA,
          lwd = 2) +
  geom_sf(data = sites_sf,
          size = 4,
          color = "black") +
  geom_sf_label(data = sites_sf,
                aes(label = site_code)) +
  geom_sf_label(data = sites_sf %>%
                  filter(site_code == root_site),
                aes(label = site_code),
                color = "red") +
  theme_bw() +
  theme(axis.title = element_blank())


#-----------------------------------------------------------------
# build parent child table
parent_child = sites_sf |>
  buildParentChild(flowlines,
                   # rm_na_parent = T,
                   add_rkm = F) |>
  editParentChild(fix_list =
                    list(c("JDA", 'MCN', "PRO"),
                         c("JDA", 'JD1', "PRO"),
                         c(NA, "JDA", "PRO"),
                         c("MCN", 'PRA', "PRO"),
                         c("MCN", "ICH", "PRO")),
                  switch_parent_child =
                    list(c("MCN", "PRO")))

# add RKMs from configuration file (since we had to fix at least one from PTAGIS)
parent_child %<>%
  left_join(configuration %>%
              select(parent = site_code,
                     parent_rkm = rkm) %>%
              distinct(),
            by = "parent") %>%
  left_join(configuration %>%
              select(child = site_code,
                     child_rkm = rkm) %>%
              distinct(),
            by = "child") %>%
  distinct()


#-----------------------------------------------------------------
# Save some configuration stuff
save(configuration,
     sites_sf,
     flowlines,
     parent_child,
     file = here('analysis/data/derived_data/site_config.rda'))


#-----------------------------------------------------------------
# pull out configuration info about all sites in the model
pro_sites <- configuration %>%
  filter(site_code %in% sites_sf$site_code) %>%
  select(node) %>%
  distinct() %>%
  left_join(configuration %>%
              select(site_code,
                     node,
                     site_name,
                     site_type,
                     site_description) %>%
              distinct(),
            multiple = "all") %>%
  select(site_code,
         node,
         site_name,
         site_type,
         site_description) %>%
  distinct() %>%
  arrange(node, site_code)

write_csv(pro_sites,
          file = here("analysis/data/derived_data",
                      "PRO_DABOM_sites_nodes.csv"))

# save flowlines
st_write(flowlines,
         dsn = here("analysis/data/derived_data",
                    "PRO_flowlines.gpkg"))

# save parent child table
parent_child %>%
  write_csv(file = here("analysis/data/derived_data",
                        "PRO_DABOM_ParentChild.csv"))

# save relevant parts of configuration file
configuration %>%
  filter(site_code %in% sites_sf$site_code) %>%
  write_csv(file = here("analysis/data/derived_data",
                        "PRO_DABOM_Configuration.csv"))


#-----------------------------------------------------------------
# Build network diagram
# simple
pc_graph = plotNodes(parent_child,
                     layout = "tree")
pc_graph

# add nodes
pc_nodes_graph = parent_child %>%
  addParentChildNodes(configuration) %>%
  plotNodes()
pc_nodes_graph


#-----------------------------------------------------------------
# A fancier version
#-----------------------------------------------------------------
library(ggraph)

# assign branch names to each site
node_order = buildNodeOrder(parent_child) %>%
  mutate(branch_nm = if_else(node == "PRO",
                             "Start",
                             if_else(str_detect(path, "TOP"),
                                     "Toppenish",
                                     if_else(str_detect(path, "SUN") & str_detect(path, "ROZ", negate = T),
                                             "Naches",
                                             if_else(str_detect(path, "ROZ"),
                                                     "UpperYakima",
                                                     if_else(str_detect(path, "SAT"),
                                                             "Satus",
                                                             "Downstream")))))) %>%
  mutate(
    across(
      branch_nm,
      ~ factor(.,
               levels = c("Start",
                          "Satus",
                          "Toppenish",
                          "Naches",
                          "UpperYakima",
                          "Downstream"))
    )
  )

# construct nodes and edges of a graph
nodes = buildNodeGraph(parent_child) %>%
  as_tibble() %>%
  left_join(node_order %>%
              select(label = node,
                     branch_nm))
edges = parent_child %>%
  left_join(nodes, by = c('parent' = 'label')) %>%
  rename(from = index) %>%
  left_join(nodes, by = c('child' = 'label')) %>%
  rename(to = index) %>%
  select(from, to)

node_graph = tidygraph::tbl_graph(nodes = nodes,
                                  edges = edges)

node_p = node_graph %>%
  ggraph(layout = "tree",
         circular = F,
         flip.y = F) +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm'),
                               type = "closed"),
                 end_cap = circle(4, 'mm')) +
  geom_node_point(size = 7,
                  aes(color = branch_nm)) +
  theme_graph(base_family = 'Times') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette = "Set2",
                     na.value = "black") +
  geom_node_label(aes(label = label),
                  size = 1.5,
                  label.padding = unit(0.1, 'lines'),
                  label.size = 0.1)

node_p

# save as pdf
ggsave(here("analysis/figures",
            "Prosser_DABOM_sites.pdf"),
       node_p,
       width = 9,
       height = 6)

# save as png
ggsave(here("analysis/figures",
            "Prosser_DABOM_sites.png"),
       node_p,
       width = 9,
       height = 6)
