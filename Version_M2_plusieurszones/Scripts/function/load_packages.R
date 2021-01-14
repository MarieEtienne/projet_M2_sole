############################
## Load and install packages
############################
#' @title ipak
#' 
#' @param pkg  vector of package to install or load
#' 

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){
    install.packages(new.pkg)} # ,repos = "http://cran.us.r-project.org" # dependencies = TRUE,type = "source"
  sapply(pkg, require, character.only = TRUE)
}



####### Package list ########
# careful VAST, TMB, INLA VMStools are more difficult to install : several packages must be installed by hand
packages <- c("chron",
              "classInt",
              "cluster",
              "cowplot",
              "curl",
              "data.table",
              # "DATRAS",
              "devtools",
              "doBy",
              "dplyr",
              "DBI",
              "FactoMineR",
              # "FLCore",
              # "fmsb",
              "geosphere",
              "ggforce",
              # "ggfortify",
              "ggmap",
              "ggplot2",
              "ggspatial",
              # "ggplotFL",
              "githubinstall",
              "grid",
              "gridExtra",
              "gstat",
              "gt",
              "gtools",
              "icesDatras",
              "INLA",
              "kableExtra",
              "latticeExtra",
              # "leaflet",
              # "lwgeom",
              "maps",
              "mapdata",
              "maptools",
              "mgcv",
              "nngeo",
              "paletteer",
              "PBSmapping",
              "plotly",
              "RandomFields",
              "raster",
              "RColorBrewer",
              "rnaturalearth",
              "rnaturalearthdata",
              "RNetCDF",
              "RPostgreSQL",
              "rworldmap",
              "sf",
              "sp",
              "spatstat",
              "spacetime",
              "stringr",
              "tidyr",
              # "tidyverse",
              "tmap",
              "TMB",
              # "tools",
              "viridis",
              #"vmstools",
              "worms",
              "ggpubr")

ipak(packages)

