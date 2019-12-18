library(tidyverse)
# Generate standard palettes etc.
error.components <- c("Climate", "Reference.climate", "Bioturbation",
                      "Seasonal.bias", "Orbital.mod.seas.bias",
                      "Seasonal.bias.unc.",
                      "Orbital.mod.seas.unc.",
                      "Calibration.unc.", "Aliasing.seasonal",
                      "Aliasing.stochastic", "Individual.variation",
                      "Meas.error", "Total.error")


# colours from RColorBrewer::brewer.pal(9, "Set1")
spec.colrs <-
  c(
    Climate = "#E41A1C",
    Reference.climate = "#E41A1C",
    Bioturbation = "#A65628",
    #Seasonal.bias = "#FFFF33", #Yellow,
    Seasonal.bias = "#FF7F00", # Orange
    Orbital.mod.seas.bias = "#FF7F00", # Orange
    Seasonal.bias.unc. = "Blue",
    Orbital.mod.seas.unc. = "Lightblue",
    Calibration.unc. = "#F781BF",
    Aliasing.seasonal = "#4DAF4A",
    Aliasing.stochastic = "#984EA3",
    Individual.variation = "#999999",
    Meas.error = "#377EB8",
    Total.error = "Black"
  )


linetype.scl  <- c(3, 5, rep(1, 10), 5)
names(linetype.scl) <- names(spec.colrs)

spec.linetypes <- linetype.scl

usethis::use_data(error.components, spec.colrs, linetype.scl, spec.linetypes, internal = FALSE, overwrite = TRUE)


proxy.types <- c("T_deg_Mg_Ca", "T_deg_Uk37", "Mg_Ca", "Uk37", "d18O")
usethis::use_data(proxy.types, overwrite = TRUE, internal = TRUE)



# table of parameters
parameter.symbol <-
  tribble(
    ~parameter, ~latex, ~label,
    "delta_t", "\\(\\Delta{t}\\)", "Delta[t]",
    #"mu_delta_t", "\\(\\overline\\Delta{t}\\)", "bar(Delta[t])",
    "tau_r", "\\(\\tau_{r}\\)", "tau[r]",
    "T", "\\(T\\)", "T",
    "tau_b", "\\(\\tau_{b}\\)", "tau[b]",
    "tau_s", "\\(\\tau_{s}\\)", "tau[s]",
    "tau_p", "\\(\\tau_{p}\\)", "tau[p]",
    "phi_c", "\\(\\langle \\phi_c \\rangle\\)",  "<phi[c]>",
    "delta_phi_c", "\\(\\Delta \\phi_c\\)",  "Delta~phi[c]",
    "N", "\\(N\\)", "italic(N)",
    #"seas.amp", "\\(A\\)", "italic(A)",
    "sig.sq_c", "\\(\\sigma_c^2\\)",  "sigma[c]^2",
    "sig.sq_a", "\\(\\sigma_a^2\\)",  "sigma[a]^2",
    "phi_a", "\\(\\phi_a\\)",  "phi[a]",
    "sigma.meas", "\\(\\sigma_{meas}\\)",  "sigma[meas]",
    "sigma.ind", "\\(\\sigma_{ind}\\)",  "sigma[ind]",
    "sigma.cal", "\\(\\sigma_{cal}\\)",  "sigma[cal]"
  )

parameter.units <-
  tribble(
    ~parameter, ~units,
    "delta_t", "years",
    #"mu_delta_t", "years",
    "tau_r", "years",
    "T", "years",
    "tau_b", "years",
    "tau_s", "years",
    "tau_p", "years (<= 1)",
    "phi_c", "\\(-\\pi \\leq \\langle \\phi_c \\rangle \\leq \\pi; 0 = \\text{midsummer}\\)",
    "delta_phi_c", "\\(0 - 2\\pi\\)",
    "N", "No. signal carriers",
    #"seas.amp", "Proxy units",
    "sig.sq_a", "?",
    "sig.sq_c", "\\(\\text{proxy units}^2\\)",
    "phi_a", "",
    "sigma.meas", "proxy units",
    "sigma.ind", "proxy units",
    "sigma.cal", "proxy units"
  )

parameter.description <-
  tribble(
    ~parameter, ~description,
    "delta_t", "The sampling frequency of the proxy record [years]",
    # "mu_delta_t", "Mean sampling frequency of an irregular timeseries [years]",
    "tau_r", "Interpreted timescale of the proxy timeseries [years]",
    "T", "Total length of the proxy record [years]",
    "tau_b", "Age heterogeneity of signal carriers due to bioturbation [years]",
    "tau_s", "Thickness of a sediment slice from which signal carrier are extracted [years]",
    "tau_p", "Proportion of the year during which signal carriers are created",
    "phi_c", "Expected phase of the signal carrier production period relative to the seasonal cycle [-pi, pi].",
    "delta_phi_c", "Uncertainty in the phase of the signal carrier production [0, 2pi].",
    "N", "No. of signal carriers per proxy measurement",
    #"seas.amp", "Peak-to-peak amplitude of the seasonal cycle [proxy units]",
    "sig.sq_c", "Variance of the seasonal cycle [\\(\\text{proxy units}^2\\)]",
    "sig.sq_a", "Variance of the orbital modulation of the seasonal cycle amplitude",
    "phi_a", "Phase of the proxy record in relation to the orbital solar radiation cycle",
    "sigma.meas", "Measurement error [proxy units]",
    "sigma.ind", "Inter-individual variation [proxy units]",
    "sigma.cal", "Calibration error [proxy units]"
  )

parameter.source <-
  tribble(
    ~parameter, ~source,
    "delta_t", "Approximated by the mean sampling frequency of an irregular timeseries",
    #"mu_delta_t", "The mean age difference between time-points in the proxy record",
    "tau_r", "Equal to \\(\\Delta{t}\\) unless explicitly estimated",
    "T", "Odd multiple of \\(\\Delta{t}\\) closest to the length of proxy record",
    "tau_b", "The bioturbation depth (estimated) divided by the sedimentation rate, or age-heterogeneity estimated from replicated radiocarbon dates",
    "tau_s", "Sediment slice thickness divided by the sedimentation rate",
    "tau_p", "Sedimentation trap data or predictions from a planktonic foraminifera model such as PLAFOM 2.0, FORAMCLIM, or FAME",
    "phi_c", "Sediment-trap data or predictions from a planktonic foraminifera model such as PLAFOM 2.0, FORAMCLIM, or FAME",
    "delta_phi_c", "",
    "N", "No. of signal carriers per proxy measurement",
    #"seas.amp", "Modern climatology, reanalysis data, e.g. HadSST",
    "sig.sq_c", "Calculated from the modern climatological amplitude of the seasonal cycle estimated from instrumental data e.g. HadSST or reanalysis data",
    "sig.sq_a", "Inferred from orbital variation in incoming solar radiation",
    "phi_a", "Frequency of orbital cycle being modelled, e.g. procession 1/23 kyr",
    "sigma.meas", "Reproducibility of measurements on real world material",
    "sigma.ind", "Individual foraminifera studies",
    "sigma.cal", "Standard error of the intercept term of a calibration regression model"
  )

usethis::use_data(parameter.description, parameter.source, parameter.symbol, parameter.units, overwrite = TRUE)


