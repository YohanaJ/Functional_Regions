rm(list = ls())
gc()

library(paletteer)
library(doParallel)
library(arrow)
library(RColorBrewer)
library(dynamicTreeCut)
library(foreach)
library(doMC)
library(dggridR)
library(tidyverse)



# specify the resolution of the map
myres <- 7

# are we doing all traits together or each separately?
do_each_trait <- TRUE

# read in the grids
ll_grids <- read_csv("grids_dggridR_lat_longs_forested.csv")

# construct the grid for the given resolution 
dggs <- dgconstruct(res = myres)

# create the suffix for writing files
suffix <- paste0("non_natives", myres, "")

# read in the community data
comm <- arrow::read_feather("comm_equal_area.feather") %>%
  select(grid_id, accepted_bin, paste0("new_grid_id_", myres)) %>%
  setNames(c("grid_id", "accepted_bin", "new_grid_id")) %>%
  distinct()

# get the lat long info
ll <- ll_grids %>% select(paste0("new_grid_id_", myres), paste0("Latitude_", myres), paste0("Longitude_", myres)) %>% 
  distinct() %>% setNames(c("new_grid_id", "new_lat", "new_lon")) %>% distinct()

# read in the trait data
dt_mat <- read_csv("traits18.csv") %>% select(-1) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'accepted_bin')

#  clean up trait matrix
dt_mat_sub <- dt_mat %>%
  filter(rownames(.) %in% comm$accepted_bin)

# trim the community data to only species in the trait data
comm_coor <- comm %>%
  filter(accepted_bin%in%rownames(dt_mat_sub)) %>% 
  group_by(new_grid_id) %>% mutate(n = length(unique(accepted_bin))) %>% ungroup %>%
  select(-n)

# summary info
spp <- sort(unique(comm_coor$accepted_bin))
grids <- sort(unique(comm_coor$new_grid_id))
(n <- length(grids))

# log/scale the data
dt_stand <- scale(log(dt_mat_sub))


# create the overall distance matrix -- this take A LOT of memory
dst <- as.matrix(dist(dt_stand, method = "euclidean"))

# cycle through each trait and get the distance matrix
for(i in 1:ncol(dt_stand)){
  print(i)
  mydst <- as.matrix(dist(dt_stand[,i], method = "euclidean"))
  arrow::write_feather(mydst %>% data.frame() %>% as_tibble(), paste0("data/julia/traits/distance_mat_", gsub("\\.", "_", tolower(colnames(dt_stand)[i])), "_", suffix, ".arrow"))
}

# make sure no missing
sum(is.na(dst))

# inspect the objects
dim(dst)
length(grids)

# create the lookup tables for the julia script to match species/grids to integers
spp_look <- tibble(col = 1:ncol(dst), accepted_bin = colnames(dst))
grid_look <- comm_coor %>% select(new_grid_id) %>% arrange(new_grid_id) %>% distinct() %>% mutate(grid_num = 1:nrow(.))

###########################################################
## write the output
arrow::write_feather(dst %>% data.frame() %>% as_tibble(), paste0("data/julia/distance_mat_", suffix, ".arrow"))
arrow::write_feather(comm_coor %>% left_join(spp_look)  %>% left_join(grid_look) %>% select(grid_num, col) %>% mutate(abundance = 1), 
                     paste0("data/julia/community_mat_", suffix, ".arrow"))
write_csv(spp_look, paste0("data/julia/spp_lookup_", suffix, ".csv"))
write_csv(ll, paste0("data/julia/lat_long_", suffix, ".csv"))
arrow::write_feather(grid_look, paste0("data/julia/grid_lookup_",suffix, ".arrow"))

