#' Create list of required parameters for ProxyErrorSpectrum
#' @param proxy.type "Mg_Ca", "T_deg_Mg_Ca", "T_deg_Uk37", "Uk37" or "d18O"
#' @param ... Arguments to ProxyErrorSpectrum, e.g. delta_t = 1000
#' @description Returns a list of arguments to \code{\link{ProxyErrorSpectrum}}
#'  where arguments in the call replace the defaults. Default values are returned
#'  for all arguments not in the call.
#' @return A list of arguments to \code{\link{ProxyErrorSpectrum}}
#' @export
#' @examples
#' \dontrun{
#' GetSpecPars()
#' }
#' GetSpecPars(proxy.type = "Mg_Ca")
#' GetSpecPars(proxy.type = "Uk37")
#' do.call(ProxyErrorSpectrum, GetSpecPars("Mg_Ca"))
GetSpecPars <- function(proxy.type, ...) {

  proxy.type <- match.arg(proxy.type, choices = proxy.types)

  args <- list(...)

  legal.parnames <- c(
    "delta_t", "tau_r", "tau_b", "tau_s", "T", "tau_p", "N",
    "n.k", "clim.spec.fun", "clim.spec.fun.args", "seas.amp", "sig.sq_a",
    "sig.sq_c", "nu_a", "nu_c", "phi_a", "phi_c", "delta_phi_c",
    "sigma.meas", "sigma.ind", "sigma.cal"
  )

  bad.names.i <- which((names(args) %in% legal.parnames) == FALSE)

  if (length(bad.names.i) > 0)
    stop(
      paste0(names(args)[bad.names.i]), " is not a valid parameter name for ProxyErrorSpectrum. Valid names include: ",
      paste(legal.parnames, collapse = ", ")
    )

  if ("sig.sq_c" %in% names(args) & "seas.amp" %in% names(args)) {
    stop("Specifiy only 1 of either sig.sq_c or seas.amp, the other will be calculated.")
  }

  if ("delta_phi_c" %in% names(args)) {
    stopifnot(args$delta_phi_c <= 2 * pi & args$delta_phi_c >= 0)
  }


  default.mg.ca.pars <- list(
    delta_t = 100,
    tau_r = 100,
    tau_b = 200,
    tau_s = 20,
    T = 1e04 + 100,
    tau_p = 4 / 12,
    N = 30,
    n.k = 15,
    clim.spec.fun = "ClimPowerFunction",
    clim.spec.fun.args = list(a = 0.1, b = 1),
    seas.amp = 4,
    sig.sq_a = 0.004,
    sig.sq_c = 2, phi_c = 0, delta_phi_c = 0,
    nu_a = 1 / 23e03, nu_c = 1 / 1,
    phi_a = pi / 2, sigma.meas = 0.26, sigma.ind = 1.5,
    sigma.cal = 0.3
  )

  default.Uk37.pars <- list(
    delta_t = 100,
    tau_r = 100,
    tau_b = 200,
    tau_s = 20,
    T = 1e04 + 100,
    tau_p = 8 / 12,
    N = Inf,
    n.k = 15,
    clim.spec.fun = "ClimPowerFunction",
    clim.spec.fun.args = list(a = 0.1, b = 1),
    seas.amp = 4,
    sig.sq_a = 0.004,
    sig.sq_c = 2, phi_c = 0, delta_phi_c = 0,
    nu_a = 1 / 23e03, nu_c = 1 / 1,
    phi_a = pi / 2, sigma.meas = 0.23, sigma.ind = 0,
    sigma.cal = 0.5
  )


  default.d18O.pars <- list(
    delta_t = 100,
    tau_r = 100,
    tau_b = 200,
    tau_s = 20,
    T = 1e04 + 100,
    tau_p = 4 / 12,
    N = 10,
    n.k = 15,
    clim.spec.fun = "ClimPowerFunction",
    clim.spec.fun.args = list(a = 0.1 / 25, b = 1),
    seas.amp = 0.8,
    sig.sq_a = 0.004,
    sig.sq_c = 0.08, phi_c = 0, delta_phi_c = 0,
    nu_a = 1 / 23e03, nu_c = 1 / 1,
    phi_a = pi / 2, sigma.meas = 0.08, sigma.ind = 0.3,
    sigma.cal = 0
  )

  default.pars <- switch(proxy.type,
    "T_deg_Mg_Ca" = default.mg.ca.pars,
    "T_deg_Uk37" = default.Uk37.pars,
    "Mg_Ca" = default.mg.ca.pars,
    "Uk37" = default.Uk37.pars,
    "d18O" = default.d18O.pars
  )

  # print(args)
  default.pars[names(args)] <- args

  if ("seas.amp" %in% names(args)) {
    default.pars[["sig.sq_c"]] <- VarSine(args$seas.amp)
  } else if ("sig.sq_c" %in% names(args)) {
    default.pars[["seas.amp"]] <- 2 * sqrt(2 * args$sig.sq_c)
  }


  return(default.pars)
}


#' Get the discrete frequencies of the power spectrum of a regularly sampled
#'  timeseries.
#'
#' @inheritParams ProxyErrorSpectrum
#' @return a vector of discrete frequencies
#' @export
#'
#' @examples
#' GetNu(1e04, 100)
GetNu <- function(T, delta_t) {
  max.abs.m <- floor(T / (2 * delta_t))
  m <- seq(-max.abs.m, max.abs.m)
  nu_m <- m / T
  return(nu_m)
}



#' Power Spectral Density of the Error in a Sediment Archived Proxy Timeseries
#' @param nu frequency
#' @param delta_t sampling frequency of the sediment core / climate timeseries
#' @param nu_a 1/tau_a = frequency of the orbital variation, e.g. for precession 1/23e03 yrs
#' @param tau_s sediment slice thickness in years (layer.width / sedimentation
#'   rate)
#' @param tau_b timescale of bioturbation (bioturbation depth / sedimentation
#'   rate) (L/sr)
#' @param tau_r width of moving average filter that represents the "interpreted"
#'   resolution of the timeseries
#' @param T length of a finite time series, e.g. 1e04
#' @param nu_c 1/1 (yr) = frequency of the seasonal cycle
#' @param tau_p length of proxy carrier "growth season" (< 12 months)
#' @param N number of signal carriers (e.g. individual foraminifera)
#' @param sig.sq_a variance of precessionary amplitude modulation
#' @param sig.sq_c variance of the seasonal cycle
#' @param clim.spec.fun a function to return climate power spectral density as a
#'   function of frequency, nu
#' @param clim.spec.fun.args arguments to the named climate power spectrum
#'   function
#' @param phi_a phase of precessionary amplitude modulation relative to centre
#'   of finite timeseries of length T
#' @param phi_c phase of carrier growth season relative to the maximum of the
#'   seasonal cycle
#' @param delta_phi_c Uncertainty in the phase of the signal carrier production.
#'   A value between 0 and 2Pi.
#' @param n.nu.prime number of discrete frequencies at which to evaluate the PSD
#'   of the error
#' @param sigma.meas the standard deviation of the per sample measurement error
#' @param sigma.ind the standard deviation of error between individuals (e.g.
#'   foraminifera) e.g. due to "vital effects" or calcification depth
#' @param sigma.cal the 1-sigma (standard deviation) of the calibration error in
#'  the same units as the proxy
#' @param n.k the number of aliasers used when estimating the error spectrum of 
#' the stochastic climate signal. Defaults to 15, does not normally need to be 
#' adjusted.
#'   
#' @return a dataframe of frequencies and spectral power
#' @export
#'
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", tau_p = 4 / 12, delta_phi_c = pi)
#' spec.obj <- do.call(ProxyErrorSpectrum, spec.pars)
#' PlotSpecError(spec.obj)
ProxyErrorSpectrum <- function(nu = NULL, tau_s, tau_b, tau_p, tau_r, T, delta_t,
                               N, n.k,
                               clim.spec.fun, clim.spec.fun.args,
                               sig.sq_a, sig.sq_c,
                               nu_a, nu_c,
                               phi_a,
                               phi_c, delta_phi_c,
                               sigma.meas, sigma.ind, sigma.cal,
                               n.nu.prime = 1000, ...) {
  spec.pars <- c(as.list(environment()), list(...))

  n.dt <- round(T / delta_t)

  if (n.dt %% 2 == 0) {
    n.dt <- n.dt + 1
    T <- n.dt * delta_t
    warning(paste0("Rounding T to ", T, " so that T is an odd multiple of delta_t"))
  } else {
    T <- n.dt * delta_t
    warning(paste0("Rounding T to ", T, " so that T is an odd integer multiple of delta_t"))
  }

  spec.pars$T <- T

  if (is.null(nu)) nu <- GetNu(T = T, delta_t = delta_t)

  # function for spectrum of climate
  ClimSpec <- match.fun(clim.spec.fun)
  ClimSpec.args <- append(list(nu), clim.spec.fun.args)

  Climate <- do.call(ClimSpec, ClimSpec.args)$spec

  # Bioturbation error (infinite.N.part) and Aliasing.stochastic (finite.N.part)
  S_E <- S_E(
    nu = nu, n.k = n.k,
    clim.spec.fun = clim.spec.fun,
    clim.spec.fun.args = clim.spec.fun.args,
    tau_s = tau_s, tau_b = tau_b, tau_r = tau_r,
    delta_t = delta_t, N = N, n.nu.prime = n.nu.prime
  )


  Orbital <- OrbitalError(
    tau_s = tau_s, tau_b = tau_b, tau_p = tau_p,
    delta_t = delta_t, T = T,
    sig.sq_c = sig.sq_c, phi_c = phi_c, delta_phi_c = delta_phi_c,
    nu_c = nu_c,
    sig.sq_a = sig.sq_a, phi_a = phi_a, nu_a = nu_a,
    N = N
  )


  # Measurement error and individual variation.
  white.noise <- WhiteNoise(
    sigma.meas = sigma.meas, sigma.ind = sigma.ind,
    N = N, delta_t = delta_t, T = T
  )


  # Get the sign of the seasonal bias
  Sign.SB <- sign(Orbital$time$B)
  if ((diff(range(Sign.SB)) < 2) == FALSE) stop("Seasonal bias changes sign.")
  Sign.SB <- median(Sign.SB)


  out <- data.frame(
    nu = nu,
    Climate = Climate,
    Reference.climate = Climate * sinc(nu * tau_r),
    Bioturbation = S_E$inf.N.part, # Bioturbation smoothing
    Aliasing.stochastic = S_E$finite.N.part, # redistributed smoothed climate variance
    Aliasing.seasonal = Orbital$freq$SU4, # Aliasing.seasonal
    Seasonal.bias = Orbital$freq$SB, # Seasonal.bias
    Seasonal.bias.unc. = Orbital$freq$SU3, # Seasonal.bias.unc.
    Meas.error = white.noise$E_meas,
    Individual.variation = white.noise$E_ind
  )

  # Calibration uncertainty as a simple constant.
  out$Calibration.unc. <- 0
  out$Calibration.unc.[out$nu == 0] <- sigma.cal^2 / abs(diff(nu)[1])


  out <- within(out, {
    Total.error = rowSums(cbind(
      Bioturbation,
      Aliasing.stochastic, Aliasing.seasonal,
      Seasonal.bias, Seasonal.bias.unc.,
      Meas.error, Individual.variation,
      Calibration.unc.
    ),
    na.rm = TRUE
    )
  })

  out[out$nu == 0, c("Bioturbation")] <- 0

  class(out) <- c("tbl_df", "tbl", class(out))

  out <- list(spec.pars = spec.pars, proxy.error.spec = out, Sign.SB = Sign.SB)
  class(out) <- c("proxy.error.spec", class(out))

  return(out)
}


#' Plot method for error spec
#'
#' @param x proxy.error.spec object
#' @param ... other parameters passed to PlotSpecError
plot.proxy.error.spec <- function(x, ...){
  PlotSpecError(x, ...)
}


#' Get the orbital error components
#'
#' @inheritParams ProxyErrorSpectrum
#'
#' @return A list with the frequency domain and time domain error components
#' @export
#'
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", phi_c = pi, delta_phi_c = pi, sig.sq_a = 0.1, tau_p = 1/12)
#' do.call(OrbitalError, spec.pars[names(spec.pars) %in% names(formals(OrbitalError))])
OrbitalError <- function(nu = NULL,
                         tau_s, tau_b, tau_p, delta_t, T,
                         sig.sq_c, phi_c, delta_phi_c, nu_c,
                         sig.sq_a, phi_a, nu_a, N){

  if (is.null(nu))  nu <- GetNu(T = T, delta_t = delta_t)

  fax <- nu
  ### time scale parameters (all the same units - in the following examples: years)
  taus <- tau_s
  taub <- tau_b
  taup <- tau_p

  Deltat <- delta_t
  T <- T     ## T must be an odd multiple of Deltat

  ### amplitudes and phases of the seasonal cycle and the orbital modulation 
  ### (nus in 1/time, here 1/years; phases in radians)

  sigc <- sqrt(sig.sq_c)
  phic.ev <- phi_c
  Deltaphic <- delta_phi_c
  nuc <- nu_c
  siga <- sqrt(sig.sq_a)
  phia <- phi_a
  nua <- nu_a


  nh <- (T/Deltat - 1) / 2     ## defined 4 lines after Eq. (94)

  tax <- seq(from=-(nh * Deltat), to=+(nh * Deltat), by=Deltat)   ## time axis

  #fax <- seq(-1/2+1/(2*T/Deltat), 1/2-1/(2*T/Deltat), length.out=T/Deltat) / Deltat     ## frequency axis

  phib1 <- 2*pi*nua*taub - atan(2*pi*nua*taub)     ## defined by Eq. (77)
  phib2 <- 4*pi*nua*taub - atan(4*pi*nua*taub)     ## defined by Eq. (A17)

  Mb1 <- 1 / sqrt(1 + (2*pi*nua*taub)^2)     ## defined by Eq. (78)
  Mb2 <- 1 / sqrt(1 + (4*pi*nua*taub)^2)     ## defined by Eq. (A18)

  snc.cp  <- sinc(nuc*taup)     ## some frequently used variables to make the lines shorter below
  snc.2cp <- sinc(2*nuc*taup)
  snc.as  <- sinc(nua*taus)
  snc.2as <- sinc(2*nua*taus)
  snc.D   <- sinc(Deltaphic / (2*pi))
  snc.2D  <- sinc(2*Deltaphic / (2*pi))

  cos.c <- cos(phic.ev)     ## ...some more...
  cos.2c <- cos(2*phic.ev)

  gamma1 <- cos.c * snc.D * snc.cp     ## defined by Eq. (88)
  gamma2 <- 1 - snc.D^2 + cos.2c * (snc.2D - snc.D^2)     ## defined by Eq. (90)

  orb1   <- cos(2*pi*nua*tax +   phia+phib1)     ## the time dependent cosines that appear in Eq. (A12); this first one appears also in Eq. (80)
  orb2d  <- cos(4*pi*nua*tax + 2*phia+phib2)
  orb2dd <- cos(4*pi*nua*tax + 2*phia+phib1)

  And <- 1 + siga*sqrt(2) * Mb1 * snc.as * orb1     ## defined by Eq. (80)

  VB0   <- sigc^2 * (1 - snc.cp^2 + cos.2c * snc.2D * (snc.2cp - snc.cp^2)) + sigc^2*siga^2 * (1 - Mb1^2 * snc.as^2 * snc.cp^2 + cos.2c * snc.2D * (snc.2cp - Mb1^2 * snc.as^2 * snc.cp^2))     ## the VB terms, defined by Eqs. (A13-A16)
  VB1   <- sigc^2*siga*sqrt(8) * (Mb1 * snc.as * (1 - snc.cp^2) + cos.2c * snc.2D * Mb1 * snc.as * (snc.2cp - snc.cp^2))
  VB2d  <- sigc^2*siga^2 * (Mb2 * snc.2as * (1 + cos.2c * snc.2D * snc.2cp))
  VB2dd <- -sigc^2*siga^2 * (Mb1^2 * snc.as^2 * snc.cp^2 * (1 + cos.2c * snc.2D))

  B <- sigc*sqrt(2) * gamma1 * And     ## the bias, defined by Eq. (87)

  U3.sq <- sigc^2 * gamma2 * And^2     ## third squared uncertainty component, defined by Eq. (89)

  var.Bnj.1   <- VB1   * orb1     ## preparing for the next step...
  var.Bnj.2d  <- VB2d  * orb2d
  var.Bnj.2dd <- VB2dd * orb2dd

  var.Bnj <- VB0 + var.Bnj.1 + var.Bnj.2d + var.Bnj.2dd     ## defined by Eq. (A12)

  var.Bnj.mean <- VB0 + var.Bnj.1[nh+1] * sinc(nua*T) + var.Bnj.2d[nh+1] * sinc(2*nua*T) + var.Bnj.2dd[nh+1] * sinc(2*nua*T)     ## defined by Eq. (A19)

  U4.sq <- var.Bnj / N     ## fourth squared uncertainty component, defined by Eq. (92)

  dm <- rep(0, T/Deltat)     ## the single-argument Kronecker delta, defined one line after Eq. (99)
  dm[nh+1] <- 1

  xi.p <- asinc(fax + nua, T, Deltat)     ## defined by Eq. (99)
  xi.m <- asinc(fax - nua, T, Deltat)

  Sc <- dm     ## defined by Eq. (97), left
  Sca <- dm * siga*sqrt(2) * Mb1 * snc.as * 2 * xi.p[nh+1] * cos(phia+phib1)     ## defined by Eq. (97), right
  Sa <- (siga^2/2) * Mb1^2 * snc.as^2 * (xi.p^2 + xi.m^2 + 2 * xi.p * xi.m * cos(2 * (phia+phib1)))     ## defined by Eq. (98)
  S0 <- T * (Sc + Sca + Sa)     ## defined by Eq. (96)

  SB <- 2 * sigc^2 * gamma1^2 * S0     ## defined by Eq. (100)

  SU3 <- sigc^2 * gamma2 * S0     ## defined by Eq. (101)

  SU4 <- rep(Deltat/N * var.Bnj.mean, T/Deltat)     ## defined by Eq. (102)


  out.f <- data.frame(nu = nu,
                      SB = SB,
                      SU3 = SU3,
                      SU4 = SU4)

  out.t <- data.frame(tax, B, U3.sq, U4.sq)

  return(list(freq = out.f, time = out.t))

}


#' Normalized sinc function
#'
#' @param x numeric
#' @param normalized logical, return the normalized or non-normalized sinc
#'  function, defaults to TRUE
#' @description Return the normalized or non-normalized sinc function. Returns 1
#' for x == 0 
#' sinc(x) = sin(pi * x) / (pi * x)
#' sinc(x) = sin(x) / (x)
#'
#' @export
#'
#' @examples
#' x <- seq(-10, 10, length.out = 1000)
#' plot(x, sinc(x), type = "l")
#' lines(x, sinc(x)^2, col = "Red")
#' abline(h = 0)
#'
#' x <- seq(-10, 10, length.out = 1000)
#' plot(x, sinc(x, normalized = FALSE), type = "l")
#' lines(x, sinc(x, normalized = FALSE)^2, col = "Red")
#' abline(h = 0)
sinc <- function(x, normalized = TRUE) {
  if (normalized) {
    y <- sin(pi * x) / (pi * x)
  } else if (normalized == FALSE) {
    y <- sin(x) / (x)
  }
  y[x == 0] <- 1
  return(y)
}

#' Aliased sinc function
#'
#' @param x numeric
#' @inheritParams ProxyErrorSpectrum
#'
#' @export
#'
#' @examples
#' plot(asinc(seq(-1 / 2, 1 / 2, length.out = 1001), T = 10, delta_t = 1), type = "l")
asinc <- function(x, T, delta_t) {
  ## the aliased sinc function, defined 3 lines before Eq. (95)

  n.dt <- round(T / delta_t)
  # if (n.dt %% 2 != 1) stop("T must be an odd multiple of delta_t")

  s <- sin(pi * x * T) / (sin(pi * x * delta_t) * T / delta_t)
  s[x == 0] <- 1
  return(s)
}


#' Power function for climate
#'
#' @param nu frequency
#' @param a scaling of the power function
#' @param b exponent of power function
#'
#' @return vector of power
#' @export
#'
#' @examples
#' ClimPowerFunction(1/10, 0.125, 1)
ClimPowerFunction <- function(nu, a, b){
  list(nu = nu, spec = a * abs(nu)^-b)
}

#' Variance of a Sine wave
#'
#' @param full.amp Peak to peak amplitude of a sine wave
#'
#' @return numeric, the variance of a sine wave
#' @export
#'
#' @examples
#' x <- sin(seq(0, 2 * pi, length.out = 1000))
#' plot(x, type = "l")
#' var(x)
#' VarSine(2)
#' x5 <- x * 5
#' var(x5)
#' VarSine(diff(range(x5)))
VarSine <- function(full.amp) {
  (full.amp / 2)^2 / 2
}

#' Get error components for stochastic climate signal
#'
#' @param plot.integral Plot the PSD of the error of the bioturbated climate
#' timeseries. This PSD is numerically integrated to give the variance for the
#' aliasing of climate variation.

#' @inheritParams ProxyErrorSpectrum
#'
#' @return list of frequencies and error spectrum components
#' @export
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", phi_c = pi, delta_phi_c = pi,
#'  sig.sq_a = 0.1, tau_p = 1/12)
#' spec.pars$nu <- GetNu(spec.pars$T, spec.pars$delta_t)
#' do.call(S_E, spec.pars[names(spec.pars) %in% names(formals(S_E))])
#' @seealso \code{\link{ProxyErrorSpectrum}}
S_E <- function(nu, tau_s, tau_b, tau_r, T, delta_t, N, n.k,
                clim.spec.fun, clim.spec.fun.args, n.nu.prime = 1000,
                plot.integral = FALSE){

  # Nyquist frequency = 1 / (2*delta_t)
  nu_star <- 1/(2*delta_t)

  # function for spectrum of climate
  S_x <- match.fun(clim.spec.fun)

  # index of aliasers
  k <- -n.k:n.k

  # matrix of nu_k, rows are distinct nu, column corresponding nu_k
  nu_k <- vapply(k, function(k) nu + 2*nu_star*k,
                 FUN.VALUE = vector(mode = "numeric", length = length(nu)))

  S_x.args <- append(list(nu_k), clim.spec.fun.args)
  S_x_nu_k <- do.call(S_x, S_x.args)$spec

  InfNPart <- function(nu_k, tau_b, tau_s, tau_r, S_x_nu_k) {
    stopifnot(dim(nu_k) == dim(S_x_nu_k))
    i <- 0 + 1i
    abs(sinc(nu_k * tau_s) *
          (tau_b^-1 / (tau_b^-1 + i * 2 * pi * nu_k)) *
          exp(i * 2 * pi * nu_k * tau_b) - sinc(nu_k * tau_r))^2 * S_x_nu_k
  }

  FiniteNPart <- function(N, T, nu_star, tau_s, tau_b){

    nu.prime.min <- 1/(1e03*tau_b)
    # highest frequency is 12 (i.e. monthly climate variance)
    nu.prime <- exp(seq(log(nu.prime.min), log(12/1), length.out = n.nu.prime))

    d.nu.prime <- diff(nu.prime)
    d.nu.prime <- c(d.nu.prime[1], d.nu.prime)

    S_x.args <- append(list(nu.prime), clim.spec.fun.args)
    S_x_nu.prime <- do.call(S_x, S_x.args)$spec

    if (is.list(S_x_nu.prime)){
      S_x_nu.prime <- S_x_nu.prime$spec
    }

    to.integrate <- (1 - sinc(nu.prime*tau_s)^2 *
                       (tau_b^-2 / (tau_b^-2 + (2 * pi * nu.prime)^2))) *
      S_x_nu.prime * 2
    # extra factor of 2 because I'm only using the positive frequencies

    if (plot.integral) {plot(nu.prime, to.integrate, type = "l", log = "x")}

    1/N * 1/(2*nu_star) * sum(to.integrate*d.nu.prime)
  }


  inf.N.part.M <- InfNPart(nu_k, tau_b, tau_s, tau_r, S_x_nu_k)
  inf.N.part <- rowSums(inf.N.part.M)

  finite.N.part <- FiniteNPart(N, T, nu_star, tau_s, tau_b)

  out <- list(nu = nu,
              inf.N.part = inf.N.part,
              finite.N.part = finite.N.part)
  return(out)
}


#' Get spectrum of white noise error components
#'
#' @inheritParams ProxyErrorSpectrum
#' @keywords internal
WhiteNoise <- function(sigma.meas, sigma.ind, N, delta_t, T) {
  nu_m <- GetNu(T = T, delta_t = delta_t)

  c_a <- tail(nu_m, 1) - head(nu_m, 1)

  # Measurement error spectrum -------
  E_meas <- sigma.meas^2 / c_a
  E_ind <- (sigma.ind^2 / N) / c_a
  return(list(nu = nu_m, E_meas = E_meas, E_ind = E_ind))
}


#' @title Order the stages/components of the error
#' @param vec vector of names of error stages/components
#' @param rev logical, reverse the order of the error stages/components
#' @return An ordered factor
OrderStages <- function(vec, rev = FALSE) {
  fac.levels <- c(
    "Climate", "Reference.climate",
    "Seasonal.bias", "Seasonal.bias.unc.", "Calibration.unc.",
    "Bioturbation",
    "Orbital.mod.seas.bias", "Orbital.mod.seas.unc.",
    "Aliasing.seasonal", "Aliasing.stochastic",
    "Individual.variation", "Meas.error",
    "Total.error"
  )
  if (rev) {
    factor(vec,
      levels = rev(fac.levels),
      ordered = T
    )
  } else {
    factor(vec, levels = fac.levels, ordered = T)
  }
}

# Process analytical error spectra --------

#' Integrate Error Spectra
#'
#' Integrates error spectra to get the variance(s) and optionally filters to
#' simulate the case where a running mean has been passed over a proxy record.
#'
#' @param pes Object of class proxy.error.spec, e.g. output from \link{ProxyErrorSpectrum}
#' @param method sinc filter is the only currently supported method
#' @param tau_smooth temporal resolution to which the proxy timeseries is
#' smoothed.  If NULL, variance is returned at all odd multiples of delta_t
#' @param thin.tau_smooth If TRUE and the vector of tau_smooth values is longer
#'  than 1000, these are thinned to 1000 values equally spaced on a log-scale
#' 
#' @return timescale dependent variance object (dataframe)
#' @export
#'
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", tau_p = 1 / 12, phi_c = 0, seas.amp = 4, T = 100 * 101)
#' spec.obj <- do.call(ProxyErrorSpectrum, spec.pars)
#' PlotSpecError(spec.obj, show.low.power.panel = FALSE)
#' var.obj <- IntegrateErrorSpectra(spec.obj)
#' PlotTSDVariance(var.obj)
#'
#' var.obj <- IntegrateErrorSpectra(spec.obj, tau_smooth = 1000)
#' 
IntegrateErrorSpectra <- function(pes, method = "sincfilter",
                                  tau_smooth = NULL, thin.tau_smooth = TRUE) {
  if (is.proxy.error.spec(pes) == FALSE)
    stop("pes must be a proxy.error.spec object")

  spec.pars <- pes[["spec.pars"]]
  pes <- pes[["proxy.error.spec"]]

  method <- match.arg(method, choices = c("sincfilter"))

  if (method == "sincfilter") {
    df <- pes
    d.nu <- df$nu[2] - df$nu[1]

    # Fix climate power at nu == 0
    df$Climate[df$nu == 0] <- NA # df$Climate[df$nu == min(df$nu[df$nu > 0])]
    df$Reference.climate[df$nu == 0] <- NA # df$Climate[df$nu == min(df$nu[df$nu > 0])]

    df.0 <- df[df$nu == 0, ]
    var.0 <- df.0 * d.nu
    var.0$smoothed.resolution <- Inf

    # Nyquist frequency
    nu_star <- 1 / (2 * spec.pars$delta_t)

    # Return variance for a specific tau_smooth
    if (is.null(tau_smooth) == FALSE) {
      flt <- asinc(df$nu, T = tau_smooth, delta_t = spec.pars$delta_t)^2
      df.m <- d.nu * df * flt

      var.df <- colSums(df.m, na.rm = TRUE)
      var.df <- c("smoothed.resolution" = tau_smooth, var.df)

      var.df <- var.df[setdiff(names(var.df), c("nu"))]
      var.0 <- var.0[, setdiff(names(var.df), c("nu"))]
      var.df <- rbind(var.0, var.df)


      # Return variance for all odd multiples of delta_t
    } else {
      T <- 1 / min(df$nu[df$nu > 0])
      
      
      tau_smooth <- seq(spec.pars$delta_t, T, 2 * spec.pars$delta_t)
      tau_smooth <- head(tau_smooth, -1)
      
      # thin tau_smooth
      n.tau_smooth <- length(tau_smooth)
      if (thin.tau_smooth == TRUE & n.tau_smooth > 1000) {
        warning("Thinning frequencies for performance reasons.")
        tau_smooth <- tau_smooth[unique(round(exp(seq(log(1), log(n.tau_smooth), length.out = 1000))))]
        }
      
      flt.lst <- lapply(tau_smooth, function(s) {
        asinc(df$nu, T = s, delta_t = spec.pars$delta_t)^2
      })
      
      df2 <- as.matrix(df) * d.nu
      df.m.lst <- lapply(flt.lst, function(flt) df2 * flt)
      
      var.df.lst <- lapply(df.m.lst, function(df.m) {
        colSums(df.m, na.rm = TRUE)
      })
      
      var.df <- bind_rows(var.df.lst)
      var.df$smoothed.resolution = tau_smooth

      var.df <- var.df[, setdiff(names(var.df), c("nu"))]
      var.0 <- var.0[, setdiff(names(var.0), c("nu"))]
      var.df <- rbind(var.0, var.df)
    }

    var.df <- list(spec.pars = spec.pars, proxy.error.var = var.df)
    class(var.df) <- c("proxy.error.var", class(var.df))

    return(var.df)
  }
}


#' Get error from proxy.error.var object
#'
#' @param var.obj timescale dependent variance object from \link{IntegrateErrorSpectra}
#' @param timescale the temporal resolution at which the proxy values are being interpreted
#' @param exclude character vector of error components to exclude from the total
#' @param include.f.zero include or exclude power at nu = 0
#' @param format format of the output
#' @return 1 row dataframe with 1 SD error components
#' @export
#'
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", tau_p = 1 / 12, phi_c = 0, seas.amp = 4, T = 100 * 101)
#' spec.obj <- do.call(ProxyErrorSpectrum, spec.pars)
#' PlotSpecError(spec.obj)
#' var.obj <- IntegrateErrorSpectra(spec.obj)
#' PlotTSDVariance(var.obj)
#' GetProxyError(var.obj, timescale = 100)
GetProxyError <- function(var.obj, timescale, exclude = NULL,
                          include.f.zero = TRUE, format = "long"){

  if (is.vector(var.obj)){
    stop("Method for vector var.obj is deprecated")
  }

  spec.pars <- var.obj[["spec.pars"]]
  var.obj <- var.obj[["proxy.error.var"]]
  
  if (timescale %in% var.obj$smoothed.resolution == FALSE) 
    stop("Error not available at this timescale / frequency")

  if (include.f.zero == FALSE){
    f.zero <- var.obj[is.infinite(var.obj$smoothed.resolution),]
    f.zero$smoothed.resolution <- 0
    var.obj <- var.obj - as.list(f.zero)
  }

  if (format == "long"){
    error <- var.obj %>%
      select(-Climate) %>%
      filter(smoothed.resolution %in% c(timescale,  "Inf")) %>%
      gather(component, value, -smoothed.resolution) %>%
      spread(smoothed.resolution, value) %>%
      gather(smoothed.resolution, value, -"Inf", -component) %>%
      rename(inc.f.zero = value,
             f.zero = `Inf`) %>%
      mutate(exl.f.zero = inc.f.zero - f.zero) %>%
      mutate_if(is.numeric, sqrt) %>%
      select(smoothed.resolution, component, everything())
  }else{
    error <- var.obj %>%
      filter(smoothed.resolution == timescale) %>%
      select(-smoothed.resolution, -Climate) %>%
      gather(component, value, -Total.error) %>%
      filter(component %in% exclude == FALSE) %>%
      mutate(Total.inc.error = sum(value, na.rm = TRUE)) %>%
      spread(component, value) %>%
      mutate_all(sqrt) %>%
      select(Total.inc.error, everything())
  }

  class(error) <- c("proxy.error", class(error))

  return(error)
}



#' Confidence interval of a correlation coefficient
#'
#' @param rho correlation coefficient
#' @param n number of data points used to estimate correlation coefficient
#' @param alpha confidence level
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' rhoCI(rho = seq(0.1, 0.9, 0.1), n = c(100))
rhoCI <- function(rho, n, alpha = 0.05){

  z1 <- qnorm(1-alpha/2, 0, 1)

  se <- sqrt(1/(n-3))

  z <- 0.5*log((1+rho)/(1-rho))

  upr <- z + se*z1
  lwr <- z - se*z1

  ci <- tanh(cbind(lwr, upr))
  out <- cbind(alpha, n, rho, ci, width = upr - lwr)
  out <- data.frame(out)
  return(out)
}

# Expected correlation -----

#' Expected correlation between replicate proxy records and between a proxy record and the true climate
#' @param pes Object of class proxy.error.spec, e.g. output from \link{ProxyErrorSpectrum}
#' @param spec.pars Parameters of the proxy error spectrum, these are taken
#' from proxy.error.spec if it is a proxy.error.spec object. Can be passed here to
#' allow calculation on compatible none proxy.error.spec objects.
#'
#' @return a data.frame / tibble
#' @export
#'
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", tau_b = 1000 * 10 / 2, sigma.meas = 1)
#' spec.obj <- do.call(ProxyErrorSpectrum, spec.pars)
#'
#' exp.corr <- ExpectedCorrelation(spec.obj)
#' plot(rho~smoothed.resolution, data = exp.corr, type = "l", log = "x")
ExpectedCorrelation <- function(pes, spec.pars = NULL) {
  if (is.proxy.error.spec(pes) & is.null(spec.pars) == FALSE) {
    warning("spec.pars are being overridden by those contained in the proxy.error.spec object")
  }

  if (is.proxy.error.spec(pes)) {
    spec.pars <- pes[["spec.pars"]]
    pes.full <- pes
    pes <- pes[["proxy.error.spec"]]
  } else {
    stop("Not an object of class 'proxy.error.spec'.")
  }

  pes[pes == Inf] <- 0

  # Fourier transform of bioturbation and sediment slice thickness filter
  TF.f_bs <- function(nu, tau_b, tau_s) {
    i <- 0 + 1i
    sinc(nu * tau_s) * (1 + i * 2 * pi * nu * tau_b)^-1 * exp(i * 2 * pi * nu * tau_b)
  }

  filt <- abs(TF.f_bs(pes$nu, spec.pars$tau_b, spec.pars$tau_s))^2


  pes.full[["proxy.error.spec"]]$Climate.biot <- pes$Climate * filt

  var.obj <- IntegrateErrorSpectra(pes.full)

  rho_rep <- with(as.list(var.obj[["proxy.error.var"]]), {

    var.noise = colSums(rbind(Aliasing.stochastic, Aliasing.seasonal,
                              Meas.error, Individual.variation),
                        na.rm = TRUE)

    Climate.biot / (Climate.biot + var.noise)

  })

  rho_clim <- with(as.list(var.obj[["proxy.error.var"]]), {
    var.noise = colSums(rbind(Bioturbation,
                              Aliasing.stochastic, Aliasing.seasonal,
                              Meas.error, Individual.variation),
                        na.rm = TRUE)

    sqrt(Climate.biot / (Climate.biot + var.noise))
  })


  if (is.vector(var.obj[["proxy.error.var"]])) {
    rho_rep
  } else if (is.data.frame(var.obj[["proxy.error.var"]])) {
    d1 <- data.frame(
      correlation.type = "proxy-proxy",
      smoothed.resolution = var.obj[["proxy.error.var"]]$smoothed.resolution,
      #n = spec.pars$T / var.obj[["proxy.error.var"]]$smoothed.resolution,
      rho = rho_rep,
      stringsAsFactors = FALSE)

    d2 <- data.frame(
      correlation.type = "proxy-clim",
      smoothed.resolution = var.obj[["proxy.error.var"]]$smoothed.resolution,
      #n = spec.pars$T / var.obj[["proxy.error.var"]]$smoothed.resolution,
      rho = rho_clim,
      stringsAsFactors = FALSE)

    out <- rbind(d1, d2)

    #ci <- rhoCI(out$rho, out$n)

    #out <- cbind(out, ci[,c("lwr", "upr")])


    return(out)
  }
}



# Plotting functions ---------


#' @title Transform axis to reversed log scale
#' @description Custom axis transformation function
#' @param base base for the log transformation, Default: exp(1)
#' @details Copied here from ecustools to avoid dependency
#' @seealso
#'  \code{\link[scales]{trans_new}},\code{\link[scales]{breaks_log}}
#' @rdname reverselog_trans
#' @source https://stackoverflow.com/a/11054781/1198125
#' @keywords internal
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}


#' Plot a proxy error spectrum
#'
#' @param pes Object of class proxy.error.spec, e.g. output from ProxyErrorSpectrum
#' @param show.low.power.panel Show components with very low power in an
#' additional panel (otherwise they are excluded)
#'
#' @return a ggplot object
#' @import dplyr ggplot2
#' @export
#' @examples
#' spec.pars <- psem::GetSpecPars("Mg_Ca", T = 1e04)
#' spec.obj <- do.call(psem::ProxyErrorSpectrum, spec.pars)
#' PlotSpecError(spec.obj)
PlotSpecError <- function(pes, show.low.power.panel = FALSE) {
  if (is.logical(show.low.power.panel) == FALSE) stop("show.low.power.panel should be TRUE or FALSE")

  if (is.proxy.error.spec(pes)) {
    spec.pars <- pes[["spec.pars"]]
    pes <- pes[["proxy.error.spec"]]
  } else {
    warning("Not an object of class 'proxy.error.spec', attempting to plot anyway.")
  }

  df <- dplyr::filter(pes, nu >= 0)
  n.row <- nrow(df)

  df <- tidyr::gather(df, component, spec, -nu)

  min.pos.nu <- min(df$nu[df$nu != 0])

  # fake.nu <- median(df$nu)
  df$nu[df$nu == 0] <- min.pos.nu / 2

  df$ax.grp <- ifelse(df$nu < min.pos.nu, "nu == 0", "nu > 0")

  # If number of discrete frequencies is very large, interpolate to 1000 to
  # keep plotting fast
  if (n.row > 1000) {
    prec <- diff(range(log(df$nu), na.rm = TRUE)) / 1000
    df <- df %>%
      group_by(component) %>%
      mutate(nu.r = round(log(nu) / prec) * prec) %>%
      group_by(component, nu.r) %>%
      # quick and dirty option of taking first value rather than averaging or interpolating
      slice(1L)
  }

  # Midpoint on log scale of total.error, for setting plotting scales
  tot.midpoint <- exp(diff(range(log(pes$Total.error), na.rm = TRUE)))

  df$spec[is.infinite(df$spec)] <- NA

  df <- df %>%
    group_by(component, ax.grp) %>%
    summarise(
      max.spec = max(spec, na.rm = TRUE),
      min.spec = min(spec, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(max.spec.grp = ifelse(max.spec < median(min.spec) / 10, "low", "high")) %>%
    ungroup() %>%
    left_join(df, .) %>%
    group_by(component) %>%
    #filter(spec > max(spec, na.rm = TRUE) / 100000) %>%
    filter(complete.cases(spec),
           spec > 0) %>%
    ungroup() %>%
    mutate(component = OrderStages(component))

  base_breaks <- function(x, n = 3) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, nint = n)
  }

  brks <- c(min.pos.nu / 2, base_breaks(df$nu[df$nu >= min.pos.nu]))
  lbls <- brks
  lbls[lbls == min.pos.nu / 2] <- 0

  if (show.low.power.panel == FALSE) df <- subset(df, max.spec.grp == "high")
  
  p <- ggplot(data = df, aes(x = nu, y = spec, colour = component, linetype = component)) +
    geom_line() +
    scale_x_continuous(trans = "log10", breaks = brks, labels = lbls) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual("", values = spec.colrs) +
    scale_linetype_manual("", values = spec.linetypes) +
    theme_bw() +
    annotation_logticks(sides = "l") +
    labs(x = "Frequency [1/years]", y = "Power Spectral Density") +
    facet_grid(max.spec.grp ~ ax.grp, scales = "free", space = "free_x",
               labeller = label_parsed) +
    theme(strip.background = element_blank(), strip.text.y = element_blank())

  df.0 <- df[df$nu < min.pos.nu & df$component %in% c(
    "Seasonal.bias",
    "Seasonal.bias.unc.",
    "Calibration.unc."
  ), ]
  if (sum(df.0$spec, na.rm = TRUE) > 0) {
    p <- p + geom_linerange(
      data = df.0, aes(ymin = 0, ymax = spec),
      lwd = 2, show.legend = FALSE, position = position_dodge(width = 0.2)
    )
  }

  # Expand limits only for nu == 0 facet
  l <- expand_limits(x = c(min.pos.nu / 3, 0.75 * min.pos.nu))
  l$data <- data.frame(x = c(min.pos.nu / 3, 0.75 * min.pos.nu), ax.grp = c("nu == 0"))
  p <- p + l

  # Add logticks only for nu > 0
  a <- annotation_logticks(sides = 'b')
  a$data <- data.frame(x = NA, ax.grp = c("nu > 0"))
  p <- p + a

  p
}



#' Plot Time-scale Dependent Variance Object.
#'
#' @param var.obj Object of class proxy.error.var, e.g. output from \link{IntegrateErrorSpectra}
#' @param components Which error components to include or exclude
#' @param include Include or exclude components named in \code{components}
#' @param include.constant.errors Include the constant (nu = 0) error components
#' @param units Are the units degC or d18O, used for axis labelling
#' @param thin.var.obj If TRUE and the variance object contains more than 1000 frequencies
#' these are thinned to 1000 frequencies equally spaced on a log-scale
#'
#' @return ggplot object
#' @import dplyr
#'
#' @export
#'
#' @examples
#' spec.pars <- GetSpecPars("Mg_Ca", delta_phi_c = pi)
#' spec.obj <- do.call(ProxyErrorSpectrum, spec.pars)
#'
#' var.obj <- IntegrateErrorSpectra(spec.obj)
#' PlotTSDVariance(var.obj, units = "degC")
#' PlotTSDVariance(var.obj, units = "d18O")
#' PlotTSDVariance(var.obj, include = FALSE, c("Seasonal.bias", "Calibration.unc."))
#' PlotTSDVariance(var.obj, include.constant.error = TRUE)
PlotTSDVariance <- function(var.obj,
                            components = c("Bioturbation", "Seasonal.bias",
                                           #"Orbital.mod.seas.bias",
                                           "Seasonal.bias.unc.",
                                           #"Orbital.mod.seas.unc.",
                                           "Calibration.unc.",
                                           "Aliasing.seasonal", "Aliasing.stochastic",
                                           "Individual.variation",
                                           "Meas.error"),
                            include = TRUE,
                            include.constant.errors = FALSE,
                            units = c("degC", "d18O"),
                            thin.var.obj = TRUE){

  if ("proxy.error.var" %in% class(var.obj)){
    spec.pars <- var.obj[["spec.pars"]]
    var.obj <-  var.obj[["proxy.error.var"]]
  }else{
    warning("Not an object of class 'proxy.error.var', attempting to plot anyway.")
  }

  n.row <- nrow(var.obj)

  cmpnts <- c("Bioturbation", "Aliasing.stochastic", "Aliasing.seasonal",
              "Seasonal.bias", "Seasonal.bias.unc.", "Meas.error", "Individual.variation",
              "Calibration.unc.", "Total.error")

  if (include.constant.errors == FALSE){
    var.obj[,][cmpnts] <- var.obj[,][cmpnts] - as.list(var.obj[var.obj$smoothed.resolution==Inf,][cmpnts])
    }

  var.obj[var.obj == Inf] <- NA

  all.cats <- c("Bioturbation", "Seasonal.bias",
                #"Orbital.mod.seas.bias",
                "Seasonal.bias.unc.",
                #"Orbital.mod.seas.unc.",
                "Calibration.unc.",
                "Aliasing.seasonal", "Aliasing.stochastic",
                "Individual.variation",
                "Meas.error")

  components <- match.arg(components,
                          choices = all.cats,
                          several.ok = TRUE)

  exclude <- if (include) {
    all.cats[all.cats %in% components == FALSE]} else if (include == FALSE){
      all.cats[all.cats %in% components == TRUE]
      }

  exclude <- c(exclude, "Total.error")

  cats <- all.cats[all.cats %in% exclude == FALSE]


  var.obj <- var.obj %>%
    #dplyr::select(-Climate, -exclude) %>%
    tidyr::gather(component, var, cats)


  # if variance object is very large, thin out to 1k points evenly spread
  # in log frequency space
  if (n.row > 1000 & thin.var.obj == TRUE){
    message("Thinning variance object for plotting")
    prec <- diff(range(log(1/var.obj$smoothed.resolution), na.rm = TRUE)) / 1000
    var.obj <- var.obj  %>%
      group_by(component) %>%
      mutate(smoothed.resolution.r = round(log(1/smoothed.resolution) / prec) * prec) %>%
      group_by(component, smoothed.resolution.r) %>%
      slice(1L)
  }

  var.obj <- var.obj %>%
    ungroup() %>%
    mutate(component = OrderStages(component, rev = TRUE))

  units <- match.arg(units, choices = c("degC", "d18O"))

  ylab <- switch(units,
         degC = expression("Error variance [\u00B0C"^2 ~"]"),
         d18O = expression("Error variance [\u2030 d"^18*"O"^2 ~"]"))



  var.obj %>%
    ggplot(aes(x = smoothed.resolution, y = var,
               fill = component, colour = component, group = component)) +
    geom_area(alpha = 0.75) +
    scale_x_continuous(trans = reverselog_trans(10)) +
    scale_color_manual("", values = spec.colrs) +
    scale_fill_manual("", values = spec.colrs) +
    labs(y = ylab, x = "Timescale [years]") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}



# Methods and classes -------

is.proxy.error.spec <- function(x) inherits(x, "proxy.error.spec")
is.proxy.error.var <- function(x) inherits(x, "proxy.error.var")



