# Figure for plotting RCMEM behavior of geographic gradient

require(rCMEM)
require(tidyverse)
require(gridExtra)
require(patchwork)

all_sites_and_variables_wide <- read_csv("data/parameterSet/all_sites_and_variables_wide.csv")

gauge_df <- data.frame(Name = c("Seattle, WA", "San Francisco, CA", "Port Isabel, TX",
                                "Pensacola, FL", "Charleston, SC", "Annapolis, MD",
                                "Portland, ME"),
                       station_id = c(9447130, 9414290, 8779770,
                                      8729840, 8665530, 8575512,
                                      8418150),
                       Latitude = c(47+36.1/60, 37+48.4/60, 26+3.7/60,
                                    30+24.3/60, 32+46.8/60, 38+59.0/60,
                                    43+39.5/60),
                       Longitude = c(-(122+2.04/60), -(122+28.0/60), -(97+12.9/60),
                                     -(87+12.7/60), -(79+55.4/60), -(76+28.9/60),
                                     -(70+14.7/60)))

for (i in 1:nrow(all_sites_and_variables_wide)) {
  mem_output_temp <- runCohortMem(startYear = 2000,
                                  relSeaLevelRiseInit = all_sites_and_variables_wide$relSeaLevelRiseInit[i],
                                  relSeaLevelRiseTotal = all_sites_and_variables_wide$relSeaLevelRiseTotal[i],
                                  initElv = all_sites_and_variables_wide$initElv[i],
                                  meanSeaLevel = all_sites_and_variables_wide$meanSeaLevel[i],
                                  meanSeaLevelDatum = all_sites_and_variables_wide$meanSeaLevelDatum[i],
                                  meanHighWaterDatum = all_sites_and_variables_wide$meanHighWaterDatum[i],
                                  meanHighHighWaterDatum = all_sites_and_variables_wide$meanHighHighWaterDatum[i],
                                  meanHighHighWaterSpringDatum = all_sites_and_variables_wide$meanHighHighWaterSpringDatum[i],
                                  suspendedSediment = 30 * 1e-06, # mg / l to grams / cc (1g / 1000 mg)  (1l / 1000 cc) 
                                  lunarNodalAmp = all_sites_and_variables_wide$lunarNodalAmp[i],
                                  lunarNodalPhase = all_sites_and_variables_wide$lunarNodalPhase[i],
                                  nFloods = all_sites_and_variables_wide$nFloods[i],
                                  bMax = 866 * 0.0001, # square meters to square cm
                                  zVegMin = -0.470,
                                  zVegMax = 2.08,
                                  zVegPeak = 0.831,
                                  plantElevationType = "dimensionless",
                                  overrideAnalyticSolution = T,
                                  rootToShoot=2,
                                  rootTurnover=0.5, rootDepthMax=30, 
                                  omDecayRate=0.5, # Lowered OM decay rate from 0.8 to 0.5 so that we could see fast pool OM in figure
                                  recalcitrantFrac=0.2,
                                  shape = "linear",
                                  captureRate = all_sites_and_variables_wide$captureRate[i]
  )
  
  annualTimeStepsTemp <- mem_output_temp$annualTimeSteps %>% 
    mutate(station_id = all_sites_and_variables_wide$station_id[i],
           Name = all_sites_and_variables_wide$Name[i])
  
  # plot(annualTimeStepsTemp$year, annualTimeStepsTemp$surfaceElevation - annualTimeStepsTemp$meanSeaLevel)
  # 
  # plot(annualTimeStepsTemp$year, annualTimeStepsTemp$cSequestration*10000)
  
  if (i == 1) {
    allTimeSteps <- annualTimeStepsTemp
  } else {
    allTimeSteps <- bind_rows(allTimeSteps, annualTimeStepsTemp)
  }
  
  # Store cohort information for Annapolis
  if(i == 2){
    mem_output_temp$cohorts %>% filter(year == 2099) -> cohorts_Annapolis
  }
  # Store cohort information for Seattle
  if(i == 7){
    mem_output_temp$cohorts %>% filter(year == 2099) -> cohorts_Seattle
  }
    
}

geo_order <- gauge_df %>% 
  arrange(-Latitude) %>% 
  mutate(Name2 = Name) %>% 
  separate(Name2, c("City", "State"), ", ")

plot_these <- allTimeSteps %>% 
  mutate(`Elevation - MSL (cm)` = surfaceElevation - meanSeaLevel,
         `Carbon flux (g m^{-2} yr^{-1})` = cFlux * 10000) %>% 
  select(station_id, Name, year, `Elevation - MSL (cm)`, `Carbon flux (g m^{-2} yr^{-1})`) %>% 
  gather(key = "variable", value = "value", -station_id, -Name, -year) %>% 
  left_join(gauge_df) %>% 
  separate(Name, c("City", "State"), ", ") %>% 
  mutate(City = factor(City, geo_order$City)) %>% 
  mutate(variable = str_replace_all(variable, " ", "~"))

plot_these_alt <- allTimeSteps %>% 
  mutate(`Z^{"*"}` = (surfaceElevation - meanSeaLevel)/(meanHighWater-meanSeaLevel),
         `Carbon flux (g m^{-2} yr^{-1})` = cFlux * 10000) %>% 
  select(station_id, Name, year, `Z^{"*"}`, `Carbon flux (g m^{-2} yr^{-1})`) %>% 
  gather(key = "variable", value = "value", -station_id, -Name, -year) %>% 
  left_join(gauge_df) %>% 
  separate(Name, c("City", "State"), ", ") %>% 
  mutate(City = factor(City, geo_order$City)) %>% 
  mutate(variable = str_replace_all(variable, " ", "~"))

plot_these_alt2 <- allTimeSteps %>%
  dplyr::mutate(`Relative Tidal Elevation` = (surfaceElevation - meanSeaLevel)/(meanHighWater-meanSeaLevel),
                `Carbon flux (g m^{-2} yr^{-1})` = cFlux * 10000,
                `MHW-MSL` = meanHighWater-meanSeaLevel) %>% 
  dplyr::group_by(Name) %>% 
  dplyr::mutate(min_mhw = min(`MHW-MSL`, na.rm = T),
                max_mhw = max(`MHW-MSL`, na.rm = T ),
                mean_mhw = mean(`MHW-MSL`, na.rm = T),
                init_msl = first(meanSeaLevel)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(`Nodal Modulation (cm)` = (`MHW-MSL`-mean_mhw),
                `Sea Level Rise (cm)` = meanSeaLevel - init_msl) %>% 
  select(station_id, Name, year,  `Sea Level Rise (cm)`, `Relative Tidal Elevation`, `Carbon flux (g m^{-2} yr^{-1})`, `Nodal Modulation (cm)`) %>% 
  gather(key = "variable", value = "value", -station_id, -Name, -year) %>% 
  left_join(gauge_df) %>% 
  separate(Name, c("City", "State"), ", ") %>% 
  mutate(City = factor(City, geo_order$City)) %>% 
  mutate(variable = str_replace_all(variable, " ", "~"),
         value = ifelse((variable == "Nodal~Modulation~(percent)") & (City %in% c("Annapolis",
                                                                                  "San Francisco")),
                        0,
                        value),
         variable2 = factor(variable, c("Sea~Level~Rise~(cm)", "Nodal~Modulation~(cm)", "Relative~Tidal~Elevation",
                                        "Carbon~flux~(g~m^{-2}~yr^{-1})")))

library(RColorBrewer)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# devtools::install_github("https://github.com/yutannihilation/ggsflabel")
library(ggsflabel)

library(plyr)

my_label_parsed <- function (variable, value) {
  if (variable == "City") {
    return(as.character(value))
  } else {
    llply(as.character(value), function(x) parse(text = x))    
  }
}

forecasts <- ggplot(plot_these, aes(x = year, y=value, color = City, group = City)) +
  geom_hline(yintercept = 0) +
  facet_grid(variable~City, scale = "free_y", 
             labeller = my_label_parsed) +
  geom_line() +
  geom_point(aes(shape = City)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(palette = "Paired") 

(forecasts)

forecasts_alt <- ggplot(plot_these_alt, aes(x = year, y=value, color = City,
                                            shape = City)) +
  geom_hline(yintercept = 0) +
  facet_grid(variable~City, scale = "free_y", 
             labeller = my_label_parsed) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(palette = "Paired") 

plot_these_alt2_points <- plot_these_alt2 %>% 
  filter(year %in% seq(min(plot_these_alt2$year),
                       max(plot_these_alt2$year),
                       by = 3))

forecasts_alt2 <- ggplot(plot_these_alt2, aes(x = year, y=value, color = City)) +
  geom_hline(yintercept = 0) +
  facet_wrap(variable2~., scale = "free_y", 
             labeller = my_label_parsed) +
  geom_line() +
  geom_point(data = plot_these_alt2_points, aes(shape = City)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.title = element_blank()) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(palette = "Paired")  +
  theme(legend.position = "right")

(forecasts_alt2)

# MLV edits to just include relative tidal elevation and carbon flux in plot

# Relative tidal elevation plot
rel_tide <- ggplot(plot_these_alt2 %>% filter(variable == "Relative~Tidal~Elevation"),
                         aes(x = year, y=value, color = City)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point(data = plot_these_alt2_points %>% filter(variable == "Relative~Tidal~Elevation"),
             aes(shape = City)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.title = element_blank()) +
  xlab("Year") +
  ylab("Relative Tidal Elevation") +
  scale_color_brewer(palette = "Paired")  +
  theme(legend.position = "bottom")

# Carbon flux plot
carbon_flux <- ggplot(plot_these_alt2 %>% filter(variable == "Carbon~flux~(g~m^{-2}~yr^{-1})"),
                      aes(x = year, y=value, color = City)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point(data = plot_these_alt2_points %>% filter(variable == "Carbon~flux~(g~m^{-2}~yr^{-1})"),
             aes(shape = City)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.title = element_blank()) +
  xlab("Year") +
  ylab(expression(paste("Carbon Flux (g ", m^{-2}, "y", r^{-1}, ")"))) +
  scale_color_brewer(palette = "Paired")  +
  theme(legend.position = "bottom")

sites <- gauge_df %>% 
  dplyr::rename(lat = Latitude,
                lon = Longitude) %>% 
  select(Name, lat, lon) %>% 
  mutate(Name = factor(Name, geo_order$Name))

map.sf <- ne_countries(scale = 'medium', type = 'map_units',
                       returnclass = 'sf')

map.na.sf <- map.sf[map.sf$continent == 'North America',]

map.na.sf.aea <- st_transform(map.na.sf, 
                              crs = "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0")

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

site_data_summary_plot_sf <- st_as_sf(sites,
                                      coords = c("lon", "lat"),
                                      crs = projcrs)

site_data_summary_plot_aea <- st_transform(site_data_summary_plot_sf,
                                           crs = "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0")

b <- st_bbox(site_data_summary_plot_aea)

mapFig <- ggplot(data = site_data_summary_plot_aea) +
  geom_sf(data=map.na.sf.aea, color="black", size=0.1, fill="darkgrey") +
  geom_sf(aes(fill=Name), color="black",
          show.legend = T, pch = 21, size = 3) +
  geom_sf_label_repel(aes(label = Name, fill = Name),
                      color = 'black',
                      # box.padding = NA,
                      point.padding = NA,
                      show.legend = FALSE,
                      size = 3, 
                      force = 50) +
  scale_fill_brewer(palette = "Paired") +
  coord_sf(xlim = c(b["xmin"]-200000, b["xmax"]+80000), ylim = c(b["ymin"]-0.5, b["ymax"]+0.5),
           crs="+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0") +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  guides(color=guide_legend(ncol=2)) 

# Bring together plots
metrics <- (rel_tide + carbon_flux) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
png("figs/Geography of MEM simple.png", height = 8.5, width = 8.5, res = 300, units = "in")
mapFig / metrics + plot_layout(heights = c(3,2)) + plot_annotation(tag_levels = "A")
dev.off()
# gridExtra::grid.arrange(mapFig, forecasts, nrow=2, ncol=1,
#                         heights = c(2,2))
# 
# g <- gridExtra::arrangeGrob(mapFig, forecasts, nrow=2, ncol=1,
#                             heights = c(2,2))
# ggsave(file="figs/Geography of MEM.pdf", g,
#        width = 7.25,
#        height = 7.25) #saves g
# 
# ggsave(file="figs/Geography of MEM.jpg", g,
#        width = 7.25,
#        height = 7.25) #saves g
# 
# g2 <- gridExtra::arrangeGrob(mapFig, forecasts_alt, nrow=2, ncol=1,
#                              heights = c(2,2))
# 
# (g2)
# 
# ggsave(file="figs/Geography of MEM alt.pdf", g2,
#        width = 7.25,
#        height = 7.25) #saves g
# 
# ggsave(file="figs/Geography of MEM alt.jpg", g2,
#        width = 7.25,
#        height = 7.25) #saves g
# 
# g3 <- gridExtra::arrangeGrob(mapFig, forecasts_alt2, nrow=2, ncol=1,
#                              heights = c(2,2))
# 
# (g3)
# 
# ggsave(file="figs/Geography of MEM alt2.pdf", g3,
#        width = 7.25,
#        height = 7.25) #saves g
# 
# ggsave(file="figs/Geography of MEM alt2.jpg", g3,
#        width = 7.25,
#        height = 7.25) #saves g
size_set <- 4
stroke_set <- 0.7
## Create plot that shows differences in cohort dynamics across two locations ####

# Fast organic matter
fast_om <- tibble(site = rep(c("Annapolis", "Seattle"),
                             c(nrow(cohorts_Annapolis), nrow(cohorts_Seattle))),
                  fast_om = c(cohorts_Annapolis$fast_OM, cohorts_Seattle$fast_OM),
                  layer_top = c(cohorts_Annapolis$layer_top, cohorts_Seattle$layer_top))

fast_om %>% 
  ggplot(aes(fast_om, layer_top, color = site, shape = site)) +
  geom_point(size = size_set, stroke = stroke_set) +
  theme_bw(base_size = 14) +
  scale_y_reverse(lim = c(50,0)) +
  scale_color_manual(values = c("#B2DF8A", "#A6CEE3")) +
  scale_shape_manual(values = c(0,1)) +
  xlab("Fast Organic Matter (g)") +
  ylab("") +
  theme(axis.text.y = element_blank()) -> fast_om_plot

# Slow organic matter
slow_om <- tibble(site = rep(c("Annapolis", "Seattle"),
                             c(nrow(cohorts_Annapolis), nrow(cohorts_Seattle))),
                  slow_om = c(cohorts_Annapolis$slow_OM, cohorts_Seattle$slow_OM),
                  layer_top = c(cohorts_Annapolis$layer_top, cohorts_Seattle$layer_top))

slow_om %>% 
  ggplot(aes(slow_om, layer_top, color = site, shape = site)) +
  geom_point(size = size_set, stroke = stroke_set) +
  theme_bw(base_size = 14) +
  scale_y_reverse(lim = c(50,0)) +
  scale_color_manual(values = c("#B2DF8A", "#A6CEE3")) +
  scale_shape_manual(values = c(0,1)) +
  xlab("Slow Organic Matter (g)") +
  ylab("") +
  theme(axis.text.y = element_blank()) -> slow_om_plot

# Root mass
root_mass <- tibble(site = rep(c("Annapolis", "Seattle"),
                             c(nrow(cohorts_Annapolis), nrow(cohorts_Seattle))),
                    root_mass = c(cohorts_Annapolis$root_mass, cohorts_Seattle$root_mass),
                  layer_top = c(cohorts_Annapolis$layer_top, cohorts_Seattle$layer_top))

root_mass %>% 
  ggplot(aes(root_mass, layer_top, color = site, shape = site)) +
  geom_point(size = size_set, stroke = stroke_set) +
  theme_bw(base_size = 14) +
  scale_y_reverse(lim = c(50,0)) +
  scale_x_continuous(breaks = c(0, 0.001, 0.002)) +
  scale_color_manual(values = c("#B2DF8A", "#A6CEE3")) +
  scale_shape_manual(values = c(0,1)) +
  xlab("Belowground Biomass (g)") +
  ylab("Depth Below Marsh Surface (cm)")-> root_plot

# Mineral mass
mineral <- tibble(site = rep(c("Annapolis", "Seattle"),
                             c(nrow(cohorts_Annapolis), nrow(cohorts_Seattle))),
                  mineral = c(cohorts_Annapolis$mineral, cohorts_Seattle$mineral),
                  layer_top = c(cohorts_Annapolis$layer_top, cohorts_Seattle$layer_top))

mineral %>% 
  ggplot(aes(mineral, layer_top, color = site, shape = site)) +
  geom_point(size = size_set, stroke = stroke_set) +
  theme_bw(base_size = 14) +
  scale_y_reverse(lim = c(50,0)) +
  scale_color_manual(values = c("#B2DF8A", "#A6CEE3")) +
  scale_shape_manual(values = c(0,1)) +
  xlab("Mineral Mass (g)") +
  ylab("") +
  theme(axis.text.y = element_blank()) -> mineral_plot

root_plot + fast_om_plot + slow_om_plot + mineral_plot +
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(tag_levels = "A")&
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 14)) -> cohort_plot_combined

png("figs/geography_cohorts.png", height = 4.8, width = 14.8, res = 300, units = "in")
cohort_plot_combined
dev.off()
