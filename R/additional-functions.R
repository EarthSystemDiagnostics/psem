#' Regularise and smooth an irregular timeseries.
#'
#' @param df A dataframe with columns of ages and values
#' @param time.var A quoted string giving the name of the age column
#' @param value.var A quoted string giving the name of the value column
#' @param tau_smooth Width of the boxcar / uniform smoothing window
#' @param d_t Difference between time steps of new regular time series
#'
#' @export
#'
#' @examples
#' set.seed(20191205)
#' n <- 100
#' dat <- data.frame(time = sort(runif(n, 1, 1e04)),
#'  value = arima.sim(list(ar = 0.9), n = n)
#'  )
#' plot(dat, type = "b")
#' dat.smoothed <- SmoothIrregular(dat, "time", "value", tau_smooth = 500)
#' dat.smoothed.2 <- SmoothIrregular(dat, "time", "value", tau_smooth = 500,
#'  d_t = 500)
#' lines(mean.value~age, data = dat.smoothed, col = "Red")
#' lines(mean.value~age, data = dat.smoothed.2, col = "Blue")
#' head(dat.smoothed)
SmoothIrregular <- function(df, time.var, value.var, tau_smooth,
                            d_t = NULL){
  
  if (is.null(d_t)) {
    d_t <- mean(diff(df[[time.var]]))
  }
  
  strt.time <- d_t * round((min(df[[time.var]]) - d_t/2)/d_t)
  end.time <-  d_t * round((max(df[[time.var]]) + d_t/2)/d_t)
  
  age <- seq(strt.time, end.time, by = d_t)
  
  lst <- lapply(age, function(a){
    df.sub <- subset(df, df[[time.var]] > a - tau_smooth / 2 & df[[time.var]] < a + tau_smooth / 2)
    list(mean.value = mean(df.sub[[value.var]]), n = nrow(df.sub),
         mean.d_t = tau_smooth / nrow(df.sub), tau_smooth = tau_smooth)
  })
  
  df.2 <- do.call(rbind.data.frame, lst)
  cbind(age, df.2)
}


#' @title Return the fraction of the orbital variations in the seasonal
#' amplitude relative to the mean seasonal amplitude
#' @param latitude Latitude in degN
#' @param maxTimeKYear Maximum time in kyr BP which should be analysed
#' @param minTimeKYear Minimum time in kyr BP which should be analysed
#' @param bPlot logical, plot the max, min and mean daily insolation
#' @importFrom orbitalforcing DailyInsolation
#' @return fraction of orbital vs. seasonal amplitude
#' @author Thomas Laepple
#' @export
#' @examples
#' RelativeAmplitudeModulation(34)
#' \dontrun{
#' library(dplyr)
#' library(tidyr)
#' df <- crossing(lat = seq(-90, 90, by = 10), maxT = c(23, 100, 1000)) %>%
#' group_by(lat, maxT) %>%
#'   mutate(sig_a = RelativeAmplitudeModulation(lat, maxTimeKYear = maxT, minTimeKYear = 0)$sig_a) %>%
#'   ungroup() %>%
#'   mutate(A_a = sig_a * sqrt(2))
#'
#'
#' df %>%
#'   ggplot(aes(x = lat, y = A_a)) +
#'   geom_line(aes(colour = factor(maxT)))
#'
#'
#' df %>%
#'   ggplot(aes(x = lat, y = sig_a^2)) +
#'   geom_line(aes(colour = factor(maxT)))
#'
#'
#' df %>%
#'   ggplot(aes(x = lat, y = sig_a)) +
#'   geom_line(aes(colour = factor(maxT)))
#' }
RelativeAmplitudeModulation <- function(latitude, maxTimeKYear = 100,
                                        minTimeKYear = 1, bPlot = FALSE) {

  kYear <- minTimeKYear:maxTimeKYear

  rangeInsolation <- matrix(NA, length(kYear), 2)
  meanInsolation <- vector()
  for (i in 1:length(kYear)) {
    temp <- orbitalforcing::DailyInsolation(kYear[i], latitude,
                                            day = 1:365)$Fsw
    rangeInsolation[i, ] <- range(temp)
    meanInsolation[i] <- mean(temp)

  }

  if (bPlot) {
    plot(rangeInsolation[, 1],
         xlab = "Time (kyrBP)", ylab = "Envelope of daily insolation (W/m^2)",
         main = paste("Lat:", latitude, "degN", sep = ""),
         type = "l", lwd = 2, ylim = range(0, rangeInsolation))
    lines(rangeInsolation[, 2], lwd = 2)
    lines(meanInsolation, lwd = 2, col = "red")
  }

  amplitude <- apply(rangeInsolation, 1, diff)

  out <- list(
    mean.amp = mean(amplitude), max.amp = max(amplitude), min.amp = min(amplitude),
    relOrbitalAmplitude = diff(range(amplitude))/mean(amplitude),
    relOrbitalMean = diff(range(amplitude))/mean(meanInsolation)
  )

  out <- c(out, sig_a = ((out$max.amp - out$mean.amp)/out$mean.amp) / sqrt(2))
  out <- c(out, sig.sq_a = out$sig_a^2)

  return(out)
}




#' Get amplitude of seasonal cycle from location
#'
#' @param longitude Longitude
#' @param latitude Latitude
#' @param depth.upr Upper end of depth range, as negative metres below sealevel
#' @param depth.lwr Lower end of depth range, as negative metres below sealevel
#' @param proxy.type Proxy type, decC or d18O
#'
#' @return A list
#' @export
#'
#' @examples
#' AmpFromLocation(120.5, 0.5, 0, -120, proxy.type = "d18O")
AmpFromLocation <- function(longitude, latitude, depth.upr, depth.lwr,
                            proxy.type=c("degC", proxy.types)){

  proxy.type <- match.arg(proxy.type, choices = c("degC", proxy.types))

  breit.lons <- unique(breitkreuz.amp$longitude)
  nearest.lon <- breit.lons[which.min(abs(longitude - breit.lons))]

  breit.lats <- unique(breitkreuz.amp$latitude)
  nearest.lat <- breit.lats[which.min(abs(latitude - breit.lats))]

  if (longitude != nearest.lon | latitude != nearest.lat){
    message(paste0("Returning for closest available coordinates: longitude = ",
                   nearest.lon, ", latitude = ", nearest.lat))
  }

  depth.sub <- breitkreuz.amp[breitkreuz.amp$longitude == nearest.lon, ]
  depth.sub <- depth.sub[depth.sub$latitude == nearest.lat, ]
  depth.sub <- depth.sub[depth.sub$depth.upr <= depth.upr, ]
  depth.sub <- depth.sub[depth.sub$depth.lwr >= depth.lwr, ]

  depth.sub$depth.wt <- with(depth.sub, depth.range / sum(depth.range))

  if (proxy.type == "d18O"){
    mean.amp <- sum(depth.sub$d18Oc_amp * depth.sub$depth.wt)
    mean.d18Oc <- sum(depth.sub$d18Oc_mean * depth.sub$depth.wt)
    sig.sq_c <- VarSine(mean.amp)
    out <- list(mean.amp=mean.amp, sig.sq_c=sig.sq_c, mean.p.T = NA,
                mean.d18Oc = mean.d18Oc)
  }else{
    mean.amp <- sum(depth.sub$p.T_amp * depth.sub$depth.wt)
    mean.p.T <- sum(depth.sub$p.T_mean * depth.sub$depth.wt)
    sig.sq_c <- VarSine(mean.amp)
    out <- list(mean.amp=mean.amp, sig.sq_c=sig.sq_c, mean.p.T = mean.p.T,
                mean.d18Oc = NA)
  }

  return(out)
}




