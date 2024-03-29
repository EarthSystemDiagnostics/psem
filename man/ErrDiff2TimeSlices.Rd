% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timeslice-error.R
\name{ErrDiff2TimeSlices}
\alias{ErrDiff2TimeSlices}
\title{Estimate the error on the difference between two time-slices}
\usage{
ErrDiff2TimeSlices(pes, tau_1, tau_2, delta_ts)
}
\arguments{
\item{pes}{Object of class proxy.error.spec, e.g. output from \link{ProxyErrorSpectrum}}

\item{tau_1}{Length of first timeslice}

\item{tau_2}{Length of second timeslice}

\item{delta_ts}{Time interval between centres of the two timeslices}
}
\value{
a dataframe
}
\description{
Estimate the error on the difference between two time-slices
}
\examples{
spec.pars <- psem::GetSpecPars("Mg_Ca", T = 1e04+100, delta_t = 100,
tau_r = 100, sig.sq_a = 0.1,
seas.amp = 6, N = 10,
tau_p = 4/12,
phi_c = 0,
phi_a = pi/2, sigma.cal = 0.3)
spec.obj <- do.call(psem::ProxyErrorSpectrum, spec.pars)

ErrDiff2TimeSlices(spec.obj, 1100, 1100, 13*spec.pars$delta_t)
}
