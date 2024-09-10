# Clear all objects from the workspace 
rm(list = ls())
gc()

# Load required libraries
library(paletteer)      # For additional color palettes
library(doParallel)     # For parallel processing
library(arrow)          # For reading and writing parquet and feather files

library(RColorBrewer)   # For color palettes
library(dggridR)        # For creating equal-area grids

library(sf)             # For handling spatial data
library(tidyverse)      # Collection of R packages for data manipulation and visualization
library(dplyr)          # Part of tidyverse, for data manipulation
library(cluster)        # For clustering algorithms





# Set the working directory
setwd("D:/")


# Combine the data into one dataframe
distib <- read_feather("equal_area/distance_matrix_7.arrow")


# Read latitude and longitude data
ll <- read_csv("lat_long.csv")

# Ensure that the column names align with the data
sum(distib$col != as.numeric(names(distib)[-1]))

# Convert the dataframe to a matrix, excluding the 'col' column
dmat <- as.matrix(distib %>% select(-col))

# Convert the matrix to a distance matrix object
distmat <- as.dist(dmat)



# Test different clustering algorithms using cophenetic correlation coefficients
# This evaluates how well the clustering preserves the pairwise distances
var.list <- c("single", "complete", "median", "mcquitty", "average", "centroid", "ward.D2")
co <- do.call(c, lapply(1:length(var.list), function(j) {
  hc <- hclust(distmat, method = var.list[j])
  return(cor(distmat, cophenetic(hc)))
}))
names(co) <- var.list
co

# Select clustering algorithms with cophenetic correlation coefficients greater than 0.7
var.list <- c("complete", "mcquitty", "average", "centroid", "ward.D2")
total_variance <- sum((dmat[upper.tri(dmat)] - mean(dmat[upper.tri(dmat)]))^2)

# Create a sequence of minimum grid sizes
size_seq <- seq(5, 500, by = 2)

# Register parallel backend to use multiple processors
registerDoParallel(8)

results_df <- data.frame()

# Loop through each clustering method and evaluate
for(method in var.list) {
  print(paste("Processing method:", method))
  
  # Compute hierarchical clustering
  hc <- hclust(distmat, method = method)
  
  # Perform analysis for each grid size
  method_results <- foreach(k = 1:length(size_seq), .combine = rbind, .inorder = FALSE) %dopar% {
    
    # Implement clustering, pruning unstable groups
    res <- dbscan::extractFOSC(hc, minPts = size_seq[k])
    
    # Get the clusters and remove the zero cluster
    clusters <- res$cluster
    uclusts <- unique(clusters)
    uclusts <- uclusts[uclusts != 0]
    
    # Calculate within-cluster sum of squares for each cluster
    within_cluster_variance <- sapply(uclusts, function(i) {
      cluster_points <- which(clusters == i)
      mydmat <- dmat[cluster_points, cluster_points]
      sum((mydmat[upper.tri(mydmat)] - mean(mydmat[upper.tri(mydmat)]))^2)
    })
    
    # Output the results
    data.frame(
      minsize = size_seq[k], 
      method = method,
      nclusts = length(uclusts), 
      missing = length(which(clusters == 0)) / length(clusters), 
      val = 1 - (sum(within_cluster_variance) / total_variance)
    )
  }
  
  # Combine results for all methods
  results_df <- rbind(results_df, method_results)
}

# Convert results to wide format (optional)
results_wide <- reshape(results_df, timevar = "method", idvar = "minsize", direction = "wide")

# Plot variance explained by number of clusters
ggplot(results_df, aes(x = nclusts , y = val, group = method, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "Cluster Analysis by Method",
       x = "Number of Clusters" ,
       y = "Variance Explained") +
  theme_minimal() +
  xlim(0, 50) + ylim(0.7, 1)

# Plot missing clusters by number of clusters
ggplot(results_df, aes(x = nclusts , y = missing, group = method, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "Cluster Analysis by Method",
       x = "Number of Clusters" ,
       y = "Missing Clusters") +
  theme_minimal() +
  xlim(0, 50) + guides(fill="none", color="none")

# Perform clustering with a specific minPts with the best method (considering variance explained and missing clusters)
hc <- hclust(distmat, method = "ward.D2")

# Perform clustering analysis using extractFOSC with minPts set to 355, the minimum number of points required to 
# explain 85% of the variance. 
res1 <- dbscan::extractFOSC(hc, minPts = 355, prune_unstable = TRUE)  # 85% variance explained, functional realms
#res1 <- dbscan::extractFOSC(hc, minPts = 161,  prune_unstable = T) ### 90 variance explained, functional biomes
#res1 <- dbscan::extractFOSC(hc, minPts =13,  prune_unstable = T)##### 99 avriance explained, functional ecoregions

hc <- res1$hc
clusters <- res1$cluster

# Get the number of clusters
(nclust <- length(unique(res1$cl)))

# Add the clusters to the latitude/longitude data
clusts <- tibble(new_grid_id = as.numeric(attr(distmat, "Labels")), grps = as.factor(as.numeric(res1$cluster))) %>% 
  left_join(ll %>% dplyr::select(new_grid_id, new_lat, new_lon) %>% distinct()) %>% 
  rename(Latitude = new_lat, Longitude = new_lon)

# Define colors for functional realms
mycols <- c("#0107C8","#3CB44B","#E6191B", "#FFE119", "#911EB4")

# Create the wrapped polygon grid
myres <- 7
dggs <- dgconstruct(res = myres)
grid <- dgcellstogrid(dggs, clusts$new_grid_id)
grid <- merge(grid, clusts, by.x = "seqnum", by.y = "new_grid_id")
wrapped_grid <- st_wrap_dateline(grid, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
countries <- map_data("world")

# Plot functional realms
ggplot() +
  geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = wrapped_grid, aes(fill = grps, color = grps)) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw() + guides(fill = "none", color = "none")






#################### Silhouette Analysis

# This section is used to evaluate the quality of clusters by calculating the silhouette width,
# which measures how similar an object is to its own cluster compared to other clusters.

############################################ Run this part if we have missing clusters

# Convert the distance matrix to a standard matrix format
dmat1 <- as.matrix(dmat)

# Set row and column names for the matrix to ensure they match
rownames(dmat1) <- as.character(colnames(dmat1))
colnames(dmat1) <- as.character(colnames(dmat1))

# Extract the clusters and labels from the hierarchical clustering result
clusters <- res1$cluster
labels <- res1$hc$labels

# Combine clusters and labels into a named vector
res <- setNames(clusters, labels)


# Calculate silhouette widths, which measure the separation distance between the resulting clusters
sil_widths <- silhouette(res, distmat)

# Plot the silhouette widths for visual assessment
windows() 
plot(sil_widths)

# Convert silhouette results to a data frame for further analysis
sil_widths <- as.data.frame(sil_widths)

# Assign the silhouette width and neighbor information to the grid data
grid$sil_width <- sil_widths$sil_width
grid$neighbor <- sil_widths$neighbor

# Recode neighbor clusters and handle missing silhouette width data
grid$neighbor <- as.factor(ifelse(grid$grps == 0, 0, grid$neighbor))
grid$sil_width <- ifelse(grid$grps == 0, NA, grid$sil_width)
nclust <- length(unique(grid$neighbor))

# Wrap the grid around the dateline for proper global plotting
wrapped_grid1 <- st_wrap_dateline(grid, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
str(wrapped_grid1)

############### Plot Silhouette Score

# Plot the silhouette width for each cluster, showing how well-separated each cluster is
ggplot() +
  geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = wrapped_grid1, aes(fill = sil_width)) +
  scale_fill_viridis_c(option = "C", guide = guide_legend(title = "Silhouette Width")) +
  theme_bw()

############# Plot Neighbor Clusters

# Plot the neighbor clusters to visualize the relationship between adjacent clusters
ggplot() +
  geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill = NA, color = "black", linewidth = 0.1) +
  geom_sf(data = wrapped_grid1, aes(fill = neighbor, color = neighbor)) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  theme_bw()



####################################################### Plot Dendrogram from Clusters

# Load libraries for dendrogram visualization
library(dendextend)
library(dbscan)
library(ggdendro)

# Convert hierarchical clustering result to a dendrogram object
dend <- as.dendrogram(res1$hc)

# Map cluster labels to specific colors for better visual differentiation
label_colors <- rep("black", length(res1$cluster))  # Set all initial colors to black
label_colors[res1$cluster == 1] <- "blue"
label_colors[res1$cluster == 5] <- "green"
label_colors[res1$cluster == 6] <- "red"
label_colors[res1$cluster == 9] <- "yellow"
label_colors[res1$cluster == 10] <- "purple"

# Assign the colors to the dendrogram's leaves (labels)
labels_colors(dend) <- label_colors[order.dendrogram(dend)]

# Further adjust the dendrogram for visualization
dend <- dend %>%
  set("labels_cex", 0.5) %>%  # Reduce the size of the labels
  set("leaves_pch", 19)       # Change the style of the leaves (points)

# Plot the dendrogram with the customized settings
plot(dend, main = "Dendrogram with Colored Clusters", ylab = "Height")


########################### Network Diagram (neighbor)

# Load library for creating alluvial diagrams (flow diagrams)
library(ggalluvial)

# Prepare the grid data for plotting
table <- grid %>%
  distinct() %>%
  mutate(grps = as.character(grps),
         neighbor = as.character(neighbor))

# Create a mapping for cluster numbers to a specific range
unique_clusters <- unique(c(table$grps, table$neighbor))
cluster_mapping <- setNames(seq_along(unique_clusters), unique_clusters)

# Apply the cluster mapping to the grid data
table <- table %>%
  mutate(grps_num = factor(cluster_mapping[grps], levels = 1:length(unique_clusters)),
         nghbr_num = factor(cluster_mapping[neighbor], levels = 1:length(unique_clusters)))

# Calculate the frequency of each cluster-neighbor relationship
cluster_relationships <- table %>%
  count(grps_num, nghbr_num) %>%
  ungroup()

# Map the cluster numbers to descriptive names for better readability
name_mapping <- c('1' = "Cool temperate", 
                  '2' = "Boreal", 
                  '3' = "Warm temperate", 
                  '4' = "Neotropical", 
                  '5' = "Paleotropical")

# Apply the name mapping to the cluster relationships data
cluster_relationships <- cluster_relationships %>%
  mutate(grps_name = name_mapping[as.character(grps_num)],
         nghbr_name = name_mapping[as.character(nghbr_num)])

# Further organize the cluster relationships for plotting
cluster_relationships <- cluster_relationships %>%
  mutate(grps_name = factor(name_mapping[as.character(grps_num)], 
                            levels = c("Boreal", "Warm temperate", "Cool temperate", "Neotropical", "Paleotropical")),
         nghbr_name = factor(name_mapping[as.character(nghbr_num)], 
                             levels = c("Boreal", "Warm temperate", "Cool temperate", "Neotropical", "Paleotropical")))

# Create and plot the alluvial diagram to visualize the relationships between clusters and their neighbors
ggplot(cluster_relationships,
       aes(axis1 = grps_name, axis2 = nghbr_name, y = n)) +
  geom_alluvium(aes(fill = grps_name), width = 0.1) +
  geom_stratum(width = 0.4, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Cluster", "Neighbor"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c('#0000FF','darkorchid', '#FFFF00','#00FF00','red')) +
  theme_minimal() +
  labs(title = "Cluster Relationships", x = "", y = "Count")
