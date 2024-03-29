% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additional-functions.R
\name{RelativeAmplitudeModulation}
\alias{RelativeAmplitudeModulation}
\title{Return the fraction of the orbital variations in the seasonal
amplitude relative to the mean seasonal amplitude}
\usage{
RelativeAmplitudeModulation(
  latitude,
  maxTimeKYear = 100,
  minTimeKYear = 1,
  bPlot = FALSE
)
}
\arguments{
\item{latitude}{Latitude in degN}

\item{maxTimeKYear}{Maximum time in kyr BP which should be analysed}

\item{minTimeKYear}{Minimum time in kyr BP which should be analysed}

\item{bPlot}{logical, plot the max, min and mean daily insolation}
}
\value{
fraction of orbital vs. seasonal amplitude
}
\description{
Return the fraction of the orbital variations in the seasonal
amplitude relative to the mean seasonal amplitude
}
\examples{
RelativeAmplitudeModulation(34)
\dontrun{
library(dplyr)
library(tidyr)
df <- crossing(lat = seq(-90, 90, by = 10), maxT = c(23, 100, 1000)) \%>\%
group_by(lat, maxT) \%>\%
  mutate(sig_a = RelativeAmplitudeModulation(lat, maxTimeKYear = maxT, minTimeKYear = 0)$sig_a) \%>\%
  ungroup() \%>\%
  mutate(A_a = sig_a * sqrt(2))


df \%>\%
  ggplot(aes(x = lat, y = A_a)) +
  geom_line(aes(colour = factor(maxT)))


df \%>\%
  ggplot(aes(x = lat, y = sig_a^2)) +
  geom_line(aes(colour = factor(maxT)))


df \%>\%
  ggplot(aes(x = lat, y = sig_a)) +
  geom_line(aes(colour = factor(maxT)))
}
}
\author{
Thomas Laepple
}
