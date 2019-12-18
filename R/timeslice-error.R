# To do
# 1. deal with freq zero - need to be included in the integral of the spectrum
# 2. Aliasing of the sinc parts of the second term


#' Sort a power spectrum into the correct order for fft(., inverse = T)
#'
#' @param spec.obj proxy.error.spec object, list (or dataframe) of frequencies (nu)
#'  and power (spec)
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' nu <- seq(-1/2, 1/2, 1/10)
#' dat <- data.frame(nu = nu, spec = 0.1*nu^1)
#' SortForFFT(dat)
SortForFFT <- function(spec.obj){

  # make sure nu is sorted in normal way
  spec.obj <- spec.obj[order(spec.obj$nu),]

  # use indexing to get in FFT order (0, pos.freq, neg.freq)
  n.row.spec <- nrow(spec.obj)
  nu.0.i <- ceiling(n.row.spec/2)
  idx <- c(nu.0.i:n.row.spec, 1:(nu.0.i-1))

  spec.obj[idx,]
}


#' Estimate the error on the difference between two time-slices
#'
#' @param tau_1 Length of first timeslice
#' @param tau_2 Length of second timeslice
#' @param delta_ts Time interval between centres of the two timeslices
#' @inheritParams IntegrateErrorSpectra
#' @inheritParams ProxyErrorSpectrum
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' spec.pars <- psem::GetSpecPars("Mg_Ca", T = 1e04+100, delta_t = 100,
#' tau_r = 100, sig.sq_a = 0.1,
#' seas.amp = 6, N = 10,
#' tau_p = 4/12,
#' phi_c = 0,
#' phi_a = pi/2, sigma.cal = 0.3)
#' spec.obj <- do.call(psem::ProxyErrorSpectrum, spec.pars)
#'
#' ErrDiff2TimeSlices(spec.obj, 1100, 1100, 13*spec.pars$delta_t)
ErrDiff2TimeSlices <- function(pes, tau_1, tau_2, delta_ts){

  if (is.proxy.error.spec(pes)){
    pes.full <- pes
    spec.pars <- pes[["spec.pars"]]
    pes <- pes[["proxy.error.spec"]]
    delta_t <- spec.pars$delta_t
  }else{
    stop("Not an object of class 'proxy.error.spec'.")
  }


  delta_ts <- abs(delta_ts)

  if ((delta_ts%%delta_t == 0) == FALSE)
    stop("distance between time-slices is not an exact multiple of delta_t")

  d_nu <- mean(abs(diff(pes$nu)))
  pes <- SortForFFT(pes)

  pes$total <- with(pes,
                         rowSums(cbind(Bioturbation,
                                       Aliasing.stochastic, Aliasing.seasonal,
                                       Seasonal.bias.unc., Seasonal.bias,
                                       #Orbital.mod.seas.bias, Orbital.mod.seas.unc.,
                                       Meas.error, Individual.variation,
                                       Calibration.unc.),
                                 na.rm = TRUE))

  # Variance of timeslices
  var.tau.1 <- IntegrateErrorSpectra(pes.full,
                                     tau_smooth = tau_1)[["proxy.error.var"]][2,]
  var.tau.2 <- IntegrateErrorSpectra(pes.full,
                                     tau_smooth = tau_2)[["proxy.error.var"]][2,]

  var.tau <- bind_rows(var.tau.1 = var.tau.1, var.tau.2 = var.tau.2, .id = "tau") %>%
    select(-smoothed.resolution) %>%
    gather(., component, var, -tau) %>%
    spread(., tau, var) %>%
    mutate(var.tau.12 = var.tau.1 + var.tau.2)

  f1 <- asinc(pes$nu, tau_1, delta_t)
  f2 <- asinc(pes$nu, tau_2, delta_t)


  pes.long <- pes %>%
    gather(component, spec, -nu) %>%
    group_by(component) %>%
    mutate(to.fft = spec * f1 * f2)

  d_ts.i <- 1 + delta_ts / delta_t

  fft.obj <- pes.long %>%
    do({
      data.frame(
        cov = Re(fft(.$to.fft, inverse = TRUE)[d_ts.i]) * d_nu
        )
    })

  out <- left_join(var.tau, fft.obj) %>%
    mutate(var.diff = var.tau.12 - 2*cov,
           sigma.diff = sqrt(var.diff)) %>%
    filter(component != "Climate")


  class(out) <- c("proxy.timeslice.error", class(out))
  return(out)

}


# Methods
print.proxy.timeslice.error <- function(x){
  i <- sapply(x, is.numeric)
  x[, i] <- round(x[, i], 3)
  print.data.frame(x)
}
