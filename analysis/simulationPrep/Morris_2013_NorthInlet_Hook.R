# Hook script for Morris et al. 2013 North Inlet Data.
# This is a hook for a digitized figure

## 1. Workspace Prep ############

# Load necessary packages
library(tidyverse)
library(lubridate)
library(fuzzyjoin)

# Load internal functions
# Internal Functions Needed
source("analysis/supportScripts/qa_functions.R")
source("analysis/supportScripts/getTidalDatums.R")

# Load data
morris_2013 <- read_csv("data/Morris_2013/intermediate/Morris_2013_WPD_output_hand_edit.csv") 

##### 2. Calculate Tidal Datums ###########

# North Inlet data are means of collections made during 2005–2010
# Oyster Landing is nearest with NAVD88 referene and complete data over time period. 
oyster_landing_datums <- getDatumsForGauge(station_id = 8662245,
                  startDate = "2005-01-01 00:00", 
                  endDate = "2010-12-31 23:59",
                  graph = T)

oyster_landing_datums_wide <- oyster_landing_datums[[1]] %>% 
  select(Datum, meters) %>% 
  spread(key = Datum, value = meters)

oyster_inundation <- oyster_landing_datums[[2]]

##### 3. Rename and Add Attributes ######

morris_curated <- morris_2013 %>% 
  mutate(study_id = "Morris_et_al_2013",
         site_id = "North Inlet",
         plot_id = NA,
         longitude = -79.166846,
         latitude = 33.328097,
         elevation =  round(elevation_round / 100, 2),
         datum = "NAVD88",
         elevation_method = NA) %>% 
  bind_cols(oyster_landing_datums_wide) %>% 
  mutate(tidal_datum_source = "NOAA gauge 8662245 2005 to 2010.",
         year = 2010,
         is_peak_biomass = "yes",
         field_or_manipulation_code = "marsh organ",
         genus = "Spartina",
         species = "alterniflora",
         harvest_or_allometry = "harvest",
         plot_area = NA) %>% 
  rename(AGB_total = mean,
         AGB_total_se = se) %>%  
  mutate(AGB_total_n = 6*5) # 6 replicates and 5 years

  
morris_output <-  order_columns(morris_curated, category="plant_plot") 

# Attributes all match
test_colnames(list(plant_plot = morris_output))

# Required and conditional attributes all present.
test_requirements(list(plant_plot = morris_output))

# Variables all match
test_variables(list(plant_plot = morris_output))        

write_csv(morris_output, "data/Morris_2013/derivative/Morris_et_al_2013_plant_plot.csv")

##### 5. Data Vis ######
morris_vis <- morris_output %>% 
  mutate(zStar = (elevation-MSL)/(MLHW-MSL))

ggplot(data = morris_vis, aes(x=zStar, y=AGB_total)) +
  geom_point() +
  geom_segment(aes(xend = zStar, y = AGB_total+AGB_total_se, yend = AGB_total-AGB_total_se))


