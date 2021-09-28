## Extract foram abundances from PLAFOM for Breitkreuz locations -----

# 21.09.2021
# Andrew Dolman

# Reference:
#   Fraile, I., M. Schulz, S. Mulitza, and M. Kucera. “Predicting the Global Distribution of Planktonic Foraminifera Using a Dynamic Ecosystem Model.” Biogeosciences 5, no. 3 (2008): 891–911.

library(ncdf4)
library(dplyr)
library(tidyr)
library(ecustools)

# Open Fraile data -----------------
filename <- system.file("extdata/BioScalar_12feb07.nc", package = "ecusdata")
nc1 <- nc_open(filename, readunlim = FALSE)

dput(names(nc1$var))

# c("no3", "amm", "iron", "silica", "sphytoC", "sphytoN", "sphytoFe",
#   "sphytoChl", "lphytoC", "lphytoN", "lphytoFe", "lphytoChl", "lphytoSi",
#   "zooC", "zooN", "zooFe", "ldetrC", "ldetrN", "ldetrFe", "ldetrSi",
#   "sdetrC", "sdetrN", "sdetrFe", "po4", "sphytoP", "lphytoP", "zooP",
#   "ldetrP", "sdetrP", "diazC", "diazN", "diazFe", "diazChl", "diazP",
#   "sphytoCaco3", "ldetrCaco3", "pachydermaSG", "pachydermaD", "pachydermaSW",
#   "ruber", "bulloides", "sacculifer")

lats <- ncvar_get(nc1, varid = "Latitude")
lons <- ncvar_get(nc1, varid = "Longitude")
lons <- ifelse(lons > 180, lons - 360, lons)

nc_close(nc1)


# Make a list of variables that can be indexed ---------
fraile.vars <- list(bulloides, sacculifer, ruber, pachydermaD, pachydermaSG)
names(fraile.vars) <- c("bulloides", "sacculifer", "ruber", "pachydermaD", "pachydermaSG")


## Find closest coordinates in Fraile data
req.coords <- psem:::breitkreuz.coords

available.coords <- crossing(lons, lats) %>% 
  distinct()

system.time(
  {dist <- geosphere::distm(req.coords[, 1:2],
                            available.coords[, 1:2],
                            fun = geosphere::distHaversine)/1000
  
  nearest.coords <- available.coords[apply(dist, 1, which.min), ]}
)

bulloides <- ncvar_get(nc1, varid = "bulloides")
sacculifer <- ncvar_get(nc1, varid = "sacculifer")
ruber <- ncvar_get(nc1, varid = "ruber")
pachydermaD <- ncvar_get(nc1, varid = "pachydermaD")
pachydermaSG <- ncvar_get(nc1, varid = "pachydermaSG")

# Make a list of variables that can be indexed ---------
fraile.vars <- list(bulloides, sacculifer, ruber, pachydermaD, pachydermaSG)
names(fraile.vars) <- c("bulloides", "sacculifer", "ruber", "pachydermaD", "pachydermaSG")


breitkreuz.plafom.abundance <- req.coords %>% 
  mutate(lon.idx = sapply(nearest.coords$lons, function(x) which(lons == x)),
         lat.idx = sapply(nearest.coords$lats, function(x) which(lats == x))
  ) %>% 
  crossing(., taxon = names(fraile.vars)) %>% 
  group_by(longitude, latitude, taxon) %>% 
  summarise(
    month = 1:12,
    plafom.ab = fraile.vars[[taxon]][lon.idx, lat.idx, ]
  )

#usethis::use_data(breitkreuz.plafom.abundance, internal = TRUE)

fl <- system.file("/extdata/breitkreuz.tbl.RData", package = "ecusdata")
load(fl)

GetAmp <- function(x) {diff(range(x, na.rm = TRUE))}

library(tidyverse)

#ecustools::d18OcFromd18OwTemp(1.3, 5)

breitkreuz.tbl.2 <- breitkreuz.tbl %>%
  mutate(d18Oc = ecustools::d18OcFromd18OwTemp(d18O, potential.temperature)) %>% 
  mutate(cell = paste(longitude, latitude)) %>% 
  arrange(cell, depth, month)


breitkreuz.var <- breitkreuz.tbl.2 %>%
  #filter(depth >= -670) %>%
  rename(p.T = potential.temperature) %>%
  filter(complete.cases(d18O, salinity, p.T)) %>%
  select(-month) %>%
  group_by(longitude, latitude, depth) %>%
  summarise_all(funs(var = var)) %>%
  ungroup()


breitkreuz.plafom.weighted.var <- breitkreuz.tbl.2 %>% 
  #slice(1:1e03) %>% 
  rename(p.T = potential.temperature) %>%
  left_join(., breitkreuz.plafom.abundance) %>% 
  filter(complete.cases(plafom.ab)) %>%
  #ungroup() %>% 
  #arrange(taxon, cell, depth, month) %>% 
  group_by(cell, longitude, latitude, depth, taxon) %>%
  summarise(across(c(p.T, d18O, salinity, d18Oc),
                   ~WeightedVar(.x, plafom.ab), .names = "{.col}_wtd.var")) %>%
  ungroup()


#save(breitkreuz.plafom.abundance, file = "breitkreuz.plafom.abundance.rdata")
#save(breitkreuz.plafom.weighted.var, file = "breitkreuz.plafom.weighted.var.rdata")


tmp <- breitkreuz.plafom.weighted.var %>% 
  #filter(taxon == "pachydermaD") %>% 
  #rename(p.T_wtd.var = p.T_var.wtd) %>% 
  left_join(., breitkreuz.var) %>% 
  filter(depth >= -25, 
         complete.cases(p.T_wtd.var) 
         ) 
tmp %>% 
  ggplot(aes(x = p.T_var, y = p.T_wtd.var)) +
  geom_point() +
  geom_abline(intercept = 0, slope = seq(0.5, 1.1, 0.1), colour = "red") +
  facet_wrap(~taxon)

tmp %>% 
  ggplot(aes(x = latitude, y = sqrt(p.T_var))) +
  geom_point(alpha = 0.025) +
  geom_point(aes(y = sqrt(p.T_wtd.var)), colour = "Red", alpha = 0.025) +
  facet_wrap(~taxon)

tmp %>% 
  ggplot(aes(x = latitude, y = sqrt(p.T_wtd.var) / sqrt(p.T_var))) +
  geom_point(alpha = 0.025) +
  facet_wrap(~taxon)



a <- PlotWorld() +
  geom_tile(data = tmp, aes(x = longitude, y = latitude,
                            group = loc, fill = sqrt(p.T_var))) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, 6)) +
  expand_limits(fill = 0)+
  facet_wrap(~taxon)

b <- PlotWorld() +
  geom_tile(data = tmp, aes(x = longitude, y = latitude,
                            group = loc, fill = sqrt(p.T_wtd.var))) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, 6)) +
  expand_limits(fill = 0)+
  facet_wrap(~taxon)

egg::ggarrange(a, b, ncol = 1)



PlotWorld() +
  geom_tile(data = tmp, aes(x = longitude, y = latitude,
                            group = loc, fill = ((p.T_wtd.var/p.T_var)))) +
  #scale_fill_viridis_c(option = "inferno") +
  scale_fill_gradient2(midpoint = 0, trans = "log10") +
  expand_limits(fill = 0)


PlotWorld() +
  geom_tile(data = tmp, aes(x = longitude, y = latitude,
                            group = loc, fill = sqrt(p.T_wtd.var) - sqrt(p.T_var))) +
  #scale_fill_viridis_c(option = "inferno") +
  scale_fill_gradient2(midpoint = 0) +
  expand_limits(fill = 0)


left_join(psem:::breitkreuz.amp, breitkreuz.var) %>% 
  filter(depth >= -525) %>% 
  mutate(amp.var = psem:::VarSine(p.T_amp)) %>% 
  ggplot(aes(x = amp.var, y = p.T_var)) +
  geom_point() +
  geom_abline(intercept = 0, slope = c(0.9, 1.1, 1), colour = "red") +
  facet_wrap(~depth, scales = "free")







