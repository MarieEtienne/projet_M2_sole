#######################
## Load commercial data
#######################

## Load tacsatEflalo
setwd(folder_phd_codes)
filter_with_Eflalo <- F
source("r/real_data/load_data/extract_commercial_data.R")

## Load ObsMer
HH_obsmer <- read.csv(paste0(ObsMer_file,"/HH_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
HL_obsmer <- read.csv(file = paste0(ObsMer_file,"/HL_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
SL_obsmer <- read.csv(file = paste0(ObsMer_file,"/SL_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
TR_obsmer <- read.csv(file = paste0(ObsMer_file,"/TR_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)

## Load scientific data
if(species_to_plot == "Solea_solea") load("data/tidy/fitting_intermediate/bob_cs/Solea_solea/OTB/with_size/mature_survey_data.RData")
if(species_to_plot == "Merluccius_merluccius") load("data/tidy/fitting_intermediate/bob_cs/Merluccius_merluccius/OTB/with_size/mature_survey_data.RData")

setwd(folder_project)

## 'VMS x logbooks'
#------------------
source("Scripts/source/shape_VMS_logbooks.R")

## ObsMer data
#-------------
source("Scripts/source/shape_ObsMer.R")

## Scientific data
#-----------------
source("Scripts/source/shape_scientific.R")
