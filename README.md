# Delineating Global Tree Functional Regions: A Trait-Based Approach
## Project Overview
This repository contains the code and datasets used in the project Delineating Global Tree Functional Regions: A Trait-Based Approach. The project aims to create a comprehensive functional regionalization of global forests by analyzing 18 key functional traits.
## Getting Started
###  Prerequisites
To run the analysis, you will need to have the following software installed:
- **Julia v1.10.3** (for functional distance calculation). Instructions to download and install Julia can be found here: https://julialang.org/downloads/
- **R 4.3.2** (for clustering analysis and figure generation). Instructions to download and install R can be found here: https://posit.co/download/rstudio-desktop/

### Dataset
The dataset used in this project is provided in the data directory. It includes all the necessary information for running the analysis, such as trait data and species distribution
### Functional Trait Imputation
The information on how functional traits were imputed is available [here](https://www.nature.com/articles/s41467-022-30888-2).
### Data Collection Information
The methodology and data on species distribution are available [here](https://onlinelibrary.wiley.com/doi/10.1111/geb.13877)


### Running the Analysis
Functional Distance Calculation (Julia)
Navigate to the `functional_distance` directory.  
Run the Julia script to calculate the functional distances.

**Note:** For our data, the script may take approximately **2 hours** to complete on a typical desktop computer, depending on the system's performance.


### Clustering and Visualization (R)
Navigate to the clustering_and_visualization directory.
Run the R script to perform clustering and generate figures.
