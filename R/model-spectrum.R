#' @title Create a model SST spectrum based on the observed spectra and a
#'   powerlaw spectrum
#' @param spectable lookup table with spectra by latitude
#'   list(spec[iLat,iFreq],lat[iLat],freq[iFreq])
#' @param latitude latitude in degN
#' @param freq desired frequency vector
#' @param beta desired scaling
#' @param variable climate variable for which to return the power spectrum,
#'   currently temperature or d18O
#' @param freq.match.lower lower frequency boundary at which the variance is
#'   matched and the spectrum is switched from instrumental to powerlaw (default
#'   0.03 as lowest reliable estimate from the instrumental record given the
#'   detrending and the multitaper bias)
#' @param freq.match.upper higher frequency boundary at which the variance is
#'   matched (default 0.1 = decadal)
#' @param bPlot TRUE = provide diagnostic plot
#' @return model spectrum list(spec,freq)
#' @author Thomas Laepple
#' @export
#' @examples
#' ModelSpectrum(latitude = 20,beta=1, bPlot=TRUE)

#' ModelSpectrum(freq = GetNu(T = 10000, delta_t = 100), latitude = 20, beta=1, bPlot=FALSE)
#' ModelSpectrum(freq = c(0.1, GetNu(T = 10000, delta_t = 100), 0.1), latitude = 20, beta=1, bPlot=FALSE)
#' ModelSpectrum(freq = c(0.6, 1), latitude = 20, beta=1, bPlot=FALSE)

ModelSpectrum <- function(freq = NULL, latitude, spectable = psem:::spectable, beta = 1,
                          variable = c("temperature", "d18O", "T_deg_Mg_Ca", "T_deg_Uk37"),
                          freq.match.lower = 0.03,
                          freq.match.upper = 0.1, bPlot = FALSE){

  variable <- match.arg(variable)

  req.freq <- freq
  freq <- freq[freq>=0]
  default.freq <- seq(from = 0.001, to = 0.5, by = 0.001)

  if (is.null(freq)) {freq = default.freq} else {
    freq = sort(c(freq, default.freq))
  }

  #abs.freq <- abs(freq)

  # Search the closest exisiting latitude in the lookup table.
  # For latitudes > 75N/S, 75N/S is used as there is not enough
  # data further north/south

  #spectable <- spectable

  iLat <- which.min(abs(latitude - spectable$lat))

  fast.spec <- approx(spectable$freq, spectable$spec[, iLat], freq)$y
  slow.spec <- 1/(freq^beta)

  # Search the indices in the frequency vector for the
  # overlapping period used to rescale
  index.match.lower <- which.min(abs(freq.match.lower - freq))
  index.match.higher <- which.min(abs(freq.match.upper - freq))

  index.match.range <- index.match.lower:index.match.higher

  ## adapt the range if there is no spectral estimate at this
  ## frequency (only applies to high southern latitudes)
  index.match.range <- index.match.range[!is.na(fast.spec[index.match.range])]


  slow.spec.scl <- (sum(fast.spec[index.match.lower:index.match.higher],
                        na.rm = TRUE)/sum(slow.spec[index.match.lower:index.match.higher]))

  slow.spec <- slow.spec * slow.spec.scl

  spec <- c(slow.spec[1:index.match.lower],
            fast.spec[-1 * (1:index.match.lower)])

  spec[is.nan(spec)] <- Inf
  spec[freq > 0.5] <- 0

  spec[is.infinite(spec)] <- NA


  if (variable == "d18O"){
    spec <- spec / (4.8^2)
  }

  if (is.null(req.freq)){
    composite <- list(freq = freq, spec = spec)
  }else{
    req.freq.i <- match(abs(req.freq), freq)
    composite <- list(freq = freq[req.freq.i], spec = spec[req.freq.i])
    #composite <- list(freq = freq, spec = spec)
    }


  if (bPlot) {

    slow.fast.spec <- c(slow.spec, fast.spec)
    slow.fast.spec <- slow.fast.spec[is.finite(slow.fast.spec)]

    plot(freq, slow.spec, type = "l", log = "xy", ylim = range(slow.fast.spec,
      fast.spec, na.rm = TRUE))
    lines(freq, fast.spec, lwd = 2, col = "red")
    lines(composite$freq, composite$spec, lwd = 3, col = "Green")
  }

  
  if (variable == "d18O") slow.spec.scl <- slow.spec.scl / (4.8^2)
  
  composite$alpha = slow.spec.scl
  composite$beta = beta
  composite$freq.match.lower = freq.match.lower

  return(composite)
}


# Code for the figures

# quartz(width=10,height=6) par(mfcol=c(1,2))
# plot(1:10,type='n',xlim=c(1/1000,0.5),ylim=c(0.01,100),log='xy',xlab='f
# (1/yr)',ylab='PSD',main='model spec, beta=1, NH')
# cLat<-seq(from=0,to=70,by=5) for (iLat in 1:length(cLat)) {
# r<-ModelSpectrum(spectable,cLat[iLat],beta=1)
# lines(r$freq,r$spec,type='l',col=rbow(15)[iLat],lwd=2) }
# legend('bottomleft',col=rbow(15)[1:7],lwd=2,title='Latitude
# +-10deg',paste(cLat[1:7]),bty='n')
# legend('bottom',col=rbow(15)[8:15],lwd=2,paste(cLat[8:15]),bty='n')
# plot(1:10,type='n',xlim=c(1/1000,0.5),ylim=c(0.01,100),log='xy',xlab='f
# (1/yr)',ylab='PSD',main='model spec, beta=1, SH')
# cLat<-seq(from=0,to=-70,by=-5) for (iLat in 1:length(cLat))
# { r<-ModelSpectrum(spectable,cLat[iLat],beta=1)
# lines(r$freq,r$spec,type='l',col=rbow(15)[iLat],lwd=2) }
# legend('bottomleft',col=rbow(15)[1:7],lwd=2,title='Latitude
# +-10deg',paste(cLat[1:7]),bty='n')
# legend('bottom',col=rbow(15)[8:15],lwd=2,paste(cLat[8:15]),bty='n')


#' Plot a spliced empirical / power-law climate spectrum
#'
#' @param model.spec Output from ModelSpectrum
#'
#' @export
#'
#' @examples
#' clim.spec <- ModelSpectrum(
#'   freq = NULL,
#'   latitude = 20,
#'   variable = "temperature",
#'    beta = 1)
#'
#' PlotModelSpectrum(clim.spec)

PlotModelSpectrum <- function(model.spec){

  f.lwr <- model.spec$freq.match.lower

  # get min pos freq to construct axis
  min.pos.freq <- min(model.spec$freq[model.spec$freq > 0])

  max.exp <- round(log10(1/min.pos.freq))
  exps <- 10^(1:max.exp)
  brks <- 1/sort(c(2, exps, 5*exps))
  
  arrow.pos.y <- exp(log(max(model.spec$spec)) - 
    (0.1 * (log(max(model.spec$spec)) - log(min(model.spec$spec)))))

  model.spec[c("freq", "spec")] %>%
    as.data.frame() %>%
    ggplot(aes(x = freq, y = spec)) +
    geom_line() +
    geom_segment(aes(x = exp(log(f.lwr) - 1), xend = exp(log(f.lwr) + 1),
                     y = arrow.pos.y,
                     yend = arrow.pos.y),
                 colour = "Black",
                 arrow = arrow(
                   length = unit(0.015, "npc"), ends = "both",
                   type = "closed"
                 )
    ) +
    annotate("text", x = exp(log(f.lwr) - 1), y = max(model.spec$spec), label = "Theoretical") +
    annotate("text", x = exp(log(f.lwr) + 1), y = max(model.spec$spec), label = "Empirical") +
    scale_x_continuous(expression("Frequency ["*year^-1*"]"),
                       trans = "log10",
                       breaks = brks,
                       labels = function(n) format(n, drop0trailing = T)
    ) +
    scale_y_continuous(expression("Power spectral density ["*K^2~"year]"),
                       trans = "log10") +
    geom_vline(xintercept = f.lwr, linetype = 4) +
    annotation_logticks() +
    theme_bw() +
    theme(panel.grid = element_blank())
}


