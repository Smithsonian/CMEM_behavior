## Figure summarising rCMEM behavior. 

require(rCMEM)
require(tidyverse)
require(gridExtra)

all_sites_and_variables_wide <- read_csv("data/parameterSet/all_sites_and_variables_wide.csv")

View(all_sites_and_variables_wide)

charleston <- all_sites_and_variables_wide %>% 
  filter(Name == "Charleston, SC")

# Should equal orthometric elevation 
initElvMin <- rCMEM::convertZStarToZ(zStar = -0.470,
                                     meanSeaLevel = charleston$meanSeaLevel[1],
                                     meanHighWater = charleston$meanSeaLevel[1] +
                                       (charleston$meanHighWaterDatum[1]-charleston$meanSeaLevelDatum[1])
                                     )

# Generic input will run from mean sea level to highest tide + sea-level rise 
transectoutput <- runCohortMemTransect(startYear = 2000, 
                                       relSeaLevelRiseInit=charleston$relSeaLevelRiseInit[1], 
                                       relSeaLevelRiseTotal=charleston$relSeaLevelRiseTotal[1],
                                       initElv=charleston$initElv[1], 
                                       meanSeaLevel=charleston$meanSeaLevel[1],
                                       meanSeaLevelDatum = charleston$meanSeaLevelDatum[1],
                                       meanHighWaterDatum=charleston$meanHighWaterDatum[1], 
                                       meanHighHighWaterDatum=charleston$meanHighHighWaterDatum[1],
                                       meanHighHighWaterSpringDatum=charleston$meanHighHighWaterSpringDatum[1], 
                                       suspendedSediment=3e-05, 
                                       lunarNodalAmp=charleston$lunarNodalAmp[1], 
                                       lunarNodalPhase = charleston$lunarNodalPhase[1],
                                       bMax = 866 * 0.0001, # square meters to square cm
                                       zVegMin = -0.470,
                                       zVegMax = 2.08,
                                       zVegPeak = 0.831,
                                       plantElevationType = "dimensionless",
                                       rootToShoot=2,
                                       rootTurnover=0.5, rootDepthMax=30, 
                                       omDecayRate=0.5,
                                       recalcitrantFrac=0.2,
                                       shape = "linear",
                                       captureRate = charleston$captureRate[1],
                                       elvIntervals=8,
                                       initElvMin = initElvMin,
                                       overrideAnalyticSolution = T)

scenarioTransectGraph <- transectoutput[[1]] %>% 
  filter(surfaceElevation != initElv) %>% 
  mutate(above_below_msl = ifelse(surfaceElevation >= meanSeaLevel, "above MSL", "below MSL"))

# Which members to visualize in Cohort Graph?

# initElv == 27
# year %in% c(2015, 2065, 2114)
cohortSubset <- transectoutput[[2]] %>% 
  filter(initElv == median(unique(transectoutput[[2]]$initElv)),
         year %in% c(2001, 2050, 2099)) %>% 
  mutate(labs1 = recode(as.character(year),
                        "2001" = "B", 
                        "2050"="C",
                        "2099"="D"),
         labs2 = paste(labs1, as.character(year), sep=". Cohort depth series, "))

scenarioSubset <- transectoutput[[1]] %>% 
  filter(initElv == median(unique(transectoutput[[2]]$initElv)),
         year %in% c(2001, 2050, 2099)) %>% 
  mutate(labs1 = recode(as.character(year),
                        "2001" = "B", 
                        "2050"="C",
                        "2099"="D"),
         labs2 = paste(labs1, as.character(year), sep=". Cohort depth series, "),
         above_below_msl = ifelse(surfaceElevation >= meanSeaLevel, "above MSL", "below MSL")) 
  

scenarioTransect <- ggplot(data=scenarioTransectGraph, aes(x=year, y=surfaceElevation, color=above_below_msl)) +
  geom_line(aes(group=as.character(initElv)), size = 1) +
  # geom_point(aes(shape=above_below_msl)) +
  # scale_shape_manual(values=c(24, 25)) +
  geom_point(data = scenarioSubset, shape = 21, size = 2, color = "red", fill = "white") +
  geom_label(data = scenarioSubset, aes(label = labs1), nudge_x = -3, color = "red", fill = "white") +
  theme_minimal() +
  ylab("Surface Elevation (cm NAVD88)") +
  scale_color_manual(values=c("black", "grey")) +
  theme(legend.title = element_blank(), 
        legend.position = "bottom") +
  ggtitle("A. CohortMEM Behavior Over Range of Initial Elevations")

(scenarioTransect)


# First reshape the mass cohorts so that they're in long form
mass_cohorts <- cohortSubset %>%
  dplyr::select(-cumCohortVol, -respired_OM) %>% 
  # dplyr::filter(complete.cases(.)) %>% 
  dplyr::group_by(year, initElv) %>% 
  dplyr::mutate(cohortIndex = length(age):1) %>% 
  ungroup() %>% 
  tidyr::gather(key = "mass_pool", value = "mass_fraction", 
                -age, -year, -layer_top, -layer_bottom, -cohortIndex, -initElv, -labs1, -labs2, -elvMin, -elvMax) %>%
  dplyr::group_by(year, age, cohortIndex, initElv) %>%
  dplyr::mutate(mass_pool = str_replace(mass_pool, "_", " "),
                mass_pool = factor(mass_pool, 
                                   levels=c("mineral",
                                            "root mass",
                                            "fast OM",
                                            "slow OM"
                                            ))) %>%
  dplyr::arrange(year, age, mass_pool) %>%
  dplyr::mutate(max_mass = cumsum(mass_fraction),
                min_mass = ifelse(mass_pool==first(mass_pool),0,lag(max_mass)),
                mass_pool = as.character(mass_pool)) %>%
  # Join mass cohorts with scenario table to convert depths to referenced elevations
  dplyr::ungroup() %>%
  dplyr::left_join(transectoutput[[1]]) %>%
  filter(initElv == median(unique(transectoutput[[2]]$initElv)),
         year %in% c(2001, 2050, 2099)) %>% 
  dplyr::mutate(layer_top=surfaceElevation-layer_top, 
                layer_bottom=surfaceElevation-layer_bottom) 

# For each unique core. 
years <- c(2001, 2050, 2099)
for (i in 1:length(years)) {
  coreCandidate_temp <- filter(cohortSubset, year == years[i])
  labs1 <- first(coreCandidate_temp$labs1)
  labs2 <- first(coreCandidate_temp$labs2)
  initElv <- first(coreCandidate_temp$initElv)
  
  core_temp <- simulateSoilCore(cohortSubset, coreYear = years[i]) %>% 
    mutate(labs1 = labs1,
           labs2 = labs2,
           initElv = initElv,
           year = years[i])
  
  if (i == 1) {
    all_cores <- core_temp
  } else {
    all_cores <- bind_rows(all_cores, core_temp)
  }
}

core_graph <- all_cores %>%
  dplyr::select(-age, -input_yrs, -om_fraction, -dry_bulk_density, -oc_fraction) %>% 
  # dplyr::filter(complete.cases(.)) %>% 
  dplyr::group_by(year, initElv) %>% 
  dplyr::mutate(cohortIndex = n():1) %>% 
  ungroup() %>% 
  tidyr::gather(key = "mass_pool", value = "mass_fraction", 
                -year, -layer_top, -layer_bottom, -cohortIndex, -initElv, -labs1, -labs2) %>%
  dplyr::group_by(year, cohortIndex, initElv, layer_top, layer_bottom) %>%
  dplyr::mutate(mass_pool = str_replace(mass_pool, "_", " "),
                mass_pool = factor(mass_pool, 
                                   levels=c("mineral",
                                            "root mass",
                                            "fast OM",
                                            "slow OM"
                                   ))) %>%
  dplyr::arrange(year, layer_bottom, mass_pool) %>%
  dplyr::mutate(mass_fraction = mass_fraction/sum(mass_fraction)) %>% 
  dplyr::mutate(max_mass = cumsum(mass_fraction),
                min_mass = ifelse(mass_pool==first(mass_pool),0,lag(max_mass)),
                mass_pool = as.character(mass_pool)) %>%
  # Join mass cohorts with scenario table to convert depths to referenced elevations
  dplyr::ungroup() %>%
  dplyr::left_join(transectoutput[[1]]) %>%
  filter(initElv == median(unique(transectoutput[[2]]$initElv)),
         year %in% c(2001, 2050, 2099)) %>% 
  dplyr::mutate(layer_top=surfaceElevation-layer_top, 
                layer_bottom=surfaceElevation-layer_bottom,
                labs2 = str_replace_all(labs2, "Cohort d", "D")) 

# Reshape the scenario table
tides <- scenarioSubset %>%
  # Track any elevation threholds in the animation speciefied.
  # meanSeaLevel and meanHighWater are the defaults
  dplyr::select(year, meanSeaLevel, meanHighWater, labs2) %>%
  tidyr::gather(value="WaterLevel", key="datum", -year, -labs2) %>%
  # dplyr::rename(year=years) %>%
  dplyr::arrange(year) %>%
  dplyr::filter(complete.cases(.)) %>% 
  mutate(datum = recode(datum, 
                             "meanSeaLevel"= "MSL",
                             "meanHighWater"="MHW")) %>% 
  filter(! (datum == "MSL" & year %in% c(2001, 2050))) %>% 
  mutate(labs2 = str_replace_all(labs2, "Cohort d", "D"))

# get rid of any NA values.               

chPalette = c("#56B4E9", "#999999", "#E69F00", "#009E73")

# gganimate stuff
graph_cores <- ggplot2::ggplot(data = core_graph, 
                                      aes(xmin = min_mass, xmax = max_mass, 
                                          ymin = layer_top, ymax = layer_bottom
                                      )) +
  ggplot2::geom_rect(aes(fill = mass_pool), color = rgb(0,0,0, alpha = 0.1)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values=chPalette) +
  ggplot2::geom_hline(data=tides, aes(yintercept=WaterLevel, lty=datum), color="blue") +
  ggplot2::ylab("Depth (cm NAVD88)") +
  ggplot2::xlab("Fractional Mass") +
  facet_wrap(.~labs2) +
  # scale_x_log10() +
  theme(legend.position = "right",
        legend.title = element_blank())

(graph_cores)

grid.arrange(scenarioTransect, graph_cores, nrow = 2)
cMemFig <- arrangeGrob(scenarioTransect, graph_cores)  
ggsave("figs/cMemBehaviorFig_core.pdf", width=7.25, height=7.25, dpi=300, cMemFig)
ggsave("figs/cMemBehaviorFig_core.jpg", width=7.25, height=7.25, dpi=300, cMemFig)
