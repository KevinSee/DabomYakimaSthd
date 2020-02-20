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
# which sites are in site_df, but not in the PTAGIS configuration file?
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

parent_child %>%
  filter(nodeOrder > 1) %>%
  select(ChildNode) %>%
  mutate_at(vars(ChildNode),
            list(~ paste0(., "_p ~ dbeta(1,1)"))) %>%
  write_delim(path = '/Users/seek/Desktop/PRO_sites.txt',
              delim = '\n')


parent_child %>%
  filter(nodeOrder == 2) %>%
  select(ChildNode) %>%
  mutate_at(vars(ChildNode),
            list(~ str_remove(., "B0$")))

library(DABOM)
setBranchNums(parent_child)
valid_paths = getValidPaths(parent_child,
                            root_site = "PRO")
node_order = createNodeOrder(valid_paths,
                             configuration,
                             site_df,
                             step_num = 2) %>%
  select(-Group) %>%
  left_join(tibble(BranchNum = 1:8,
                   Group = c(rep("Downstream", 5),
                             'Status',
                             'Toppenish',
                             'Sunnyside')))
node_order %>%
  filter(NodeOrder == 2)

createNodeList(node_order)



branch_nums  = setBranchNums(parent_child) %>%
  unlist() %>%
  enframe(name = 'var',
          value = "n_branch") %>%
  mutate(SiteID = str_remove(var, 'n_pops_')) %>%
  select(SiteID, n_branch)


site = "ROZ"

n_branch_site = branch_nums %>%
  filter(SiteID == site) %>%
  pull(n_branch)


node_order %>%
  filter(grepl(site, Path)) %>%
  select(NodeSite) %>%
  distinct() %>%
  slice(1:n_branch_site) %>%
  filter(NodeSite != site) %>%
  mutate(nm = paste0('past_', NodeSite)) %>%
  pull(nm) %>%
  c(., paste0(site, "_bb"))

parent_child %>%
  mutate(ParentSite = str_remove(str_remove(ParentNode, "B0$"), "A0$"),
         ChildSite = str_remove(str_remove(ChildNode, "B0$"), "A0$"))

parent_child %>%
  mutate_at(vars(ParentNode, ChildNode),
            list(~ str_remove(str_remove(., "B0$"), "A0$"))) %>%
  rename(ParentSite = ParentNode,
         ChildSite = ChildNode) %>%
  filter(ParentSite != ChildSite) %>%
  group_by(ParentSite) %>%
  summarise(n_childs = n_distinct(ChildSite),
            child_sites = list(ChildSite)) %>%
  # filter(n_childs > 1) %>%
  mutate(vec_nm = map(child_sites,
                      .f = function(x) {
                        paste0('past_', x)
                      })) %>%
  mutate(vec_nm = map2(vec_nm,
                       ParentSite,
                       .f = function(x, y) {
                         if(length(x) > 1) {
                           c(x, paste0(y, "_bb"))
                         } else {
                           x
                         }
                       })) %>%
  mutate(param_nm = paste0("p_pop_", ParentSite)) %>%
  unnest(vec_nm) %>%
  group_by(param_nm) %>%
  mutate(vec_pos = 1:n()) %>%
  ungroup() %>%
  mutate(param_nm = if_else(n_childs == 1,
                            str_replace(vec_nm, 'past', 'phi'),
                            param_nm)) %>%
  mutate(param_nm = if_else(n_childs > 1,
                            paste0(param_nm, "[", vec_pos, "]"),
                            param_nm)) %>%
  select(param_nm, vec_nm)
