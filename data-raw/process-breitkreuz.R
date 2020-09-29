## process Breitkreuz et al d18O, salinity and temperature data for amplitudes and means
#fl <- system.file("/extdata/breitkreuz.tbl.RData", package = "ecusdata")
load(fl)

GetAmp <- function(x) {diff(range(x, na.rm = TRUE))}

library(tidyverse)

ecustools::d18OcFromd18OwTemp(1.3, 5)

breitkreuz.tbl.2 <- breitkreuz.tbl %>%
  mutate(d18Oc = ecustools::d18OcFromd18OwTemp(d18O, potential.temperature))


breitkreuz.amp <- breitkreuz.tbl.2 %>%
  filter(depth >= -670) %>%
  rename(p.T = potential.temperature) %>%
  filter(complete.cases(d18O, salinity, p.T)) %>%
  select(-month) %>%
  group_by(longitude, latitude, depth) %>%
  summarise_all(funs(amp = GetAmp)) %>%
  ungroup()

breitkreuz.depth.tbl <- breitkreuz.depth.tbl %>%
  mutate(depth.range = depth.upr - depth.lwr)

breitkreuz.amp <- breitkreuz.amp %>%
  left_join(., breitkreuz.depth.tbl)

## Add mean temperature to breitkreuz.amp
breitkreuz.mean <- breitkreuz.tbl.2 %>%
  filter(depth >= -670) %>%
  rename(p.T = potential.temperature) %>%
  filter(complete.cases(d18O, salinity, p.T)) %>%
  select(-month) %>%
  group_by(longitude, latitude, depth) %>%
  summarise_all(funs(mean = mean)) %>%
  ungroup()

breitkreuz.amp <- left_join(breitkreuz.amp, breitkreuz.mean)

#usethis::use_data(breitkreuz.amp, breitkreuz.depth.tbl, overwrite = TRUE, internal = FALSE)



# 
# 
# breitkreuz.amp %>%
#   filter(longitude == 0.5) %>%
#   ggplot(aes(x = latitude, y = p.T_amp, colour = factor(depth))) +
#   geom_point()
# 
# 
# breitkreuz.tbl.2 %>%
#   filter(longitude == 10.5, latitude == 35.5, depth >= -670) %>%
#   ggplot(aes(x = potential.temperature, y = d18Oc, colour = factor(depth))) +
#   geom_point()
# 
# 
# breitkreuz.amp %>%
#   filter(longitude == 0.5) %>%
#   ggplot(aes(x = d18Oc_amp, y = p.T_amp, colour = factor(depth))) +
#   geom_point() +
#   facet_wrap(~depth, scales = "free") +
#   geom_abline(intercept = 0, slope = 4.8)
# 
# 
# breitkreuz.mean %>%
#   filter(longitude == 0.5) %>%
#   ggplot(aes(x = latitude, y = p.T_mean, colour = depth, group = depth)) +
#   geom_line()
# 
# 
# PlotWorld <- function(){
#   world <- map_data("world")
#   worldmap <- ggplot(world, aes(x = long, y = lat, group = group)) +
#     geom_polygon(fill = "Grey") +
#     coord_quickmap() +
#     scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 60))+
#     scale_y_continuous(limits = c(-90, 90), breaks = c(-45, 0, 45)) +
#     labs(x = "Longitude", y = "Latitude") +
#     theme_bw()
#   worldmap
# }
# 
# 
# PlotWorld() +
#   geom_tile(data = filter(breitkreuz.amp, depth == -25), aes(x = longitude, y = latitude,
#                             group = d18O_amp, fill = d18O_amp)) +
#   scale_fill_viridis_c(option = "inferno") +
#   expand_limits(fill = 0)
# 
# PlotWorld() +
#   geom_tile(data = filter(breitkreuz.amp, depth == -25),
#             aes(x = longitude, y = latitude,
#                 group = d18Oc_amp, fill = d18Oc_amp)) +
#   scale_fill_viridis_c(option = "inferno") +
#   expand_limits(fill = 0)
# 
# PlotWorld() +
#   geom_tile(data = filter(breitkreuz.amp, depth == -25),
#             aes(x = longitude, y = latitude,
#                 group = d18O_amp, fill = p.T_amp)) +
#   scale_fill_viridis_c(option = "inferno") +
#   expand_limits(fill = 0)
# 
# 
# PlotWorld() +
#   geom_tile(data = breitkreuz.amp, aes(x = longitude, y = latitude,
#                                                   group = p.T_amp,
#                                                   fill = p.T_amp)) +
#   scale_fill_viridis_c(option = "inferno") +
#   expand_limits(fill = 0) +
#   facet_wrap(~depth, labeller = label_both)
