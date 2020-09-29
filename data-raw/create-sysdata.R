# create sysdata

source("data-raw/create-lookup-tables-and-palettes.R")
source("data-raw/import-empirical-spectra.R")
source("data-raw/process-breitkreuz.R")

usethis::use_data(error.components, spec.colrs, linetype.scl, spec.linetypes, proxy.types,
parameter.description, parameter.source, parameter.symbol, parameter.units,
spectable,
breitkreuz.amp, breitkreuz.depth.tbl, internal = TRUE, overwrite = TRUE)

