test_that("spec pipe line works", {
  
  err.cache <- structure(
    list(
      smoothed.resolution = c("100", "100", "100", "100","100", "100", "100",
                              "100", "100", "100"),
                              component = c("Aliasing.seasonal", "Aliasing.stochastic",
                                            "Bioturbation", "Calibration.unc.",
                                            "Individual.variation","Meas.error",
                                            "Reference.climate", "Seasonal.bias",
                                            "Seasonal.bias.unc.",  "Total.error"),
                              f.zero = c(0.000392894339637473, 0.0251216652936918, 0, 0.3, 
                                         0.0273861278752583, 0.026, NA, 1.97722503894702, 0, 
                                         2.00036899272913),
                              inc.f.zero = c(0.00394853924559903, 0.252469611597032, 
                                             0.525258573452768,  0.3, 0.275227178890458, 
                                             0.261296766149143, 0.927987688415085,
                                             1.98063144982691, 0, 2.12051408973654),
                              exl.f.zero = c(0.00392894339637473, 0.251216652936918,
                                             0.525258573452768, 0, 0.273861278752583, 0.26, 
                                             NA, 0.116112382648807, 0, 0.703636196978993)),
                         class = c("proxy.error", "data.frame"), row.names = c(NA, -10L))
  
  spec.pars <- GetSpecPars("Mg_Ca", tau_p = 1 / 12, phi_c = 0, seas.amp = 4, T = 100 * 101)
  
  spec.obj <- do.call(ProxyErrorSpectrum, spec.pars)
  
  var.obj <- IntegrateErrorSpectra(spec.obj)
  
  
  err <- GetProxyError(var.obj, timescale = 100)

  expect_equal(err, err.cache)
 
   })


