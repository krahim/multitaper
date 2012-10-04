%     The multitaper R package
%     Multitaper and spectral analysis package for R
%     Copyright (C) 2011 Karim Rahim 
%
%     Written by Karim Rahim and Wesley Burr.
%
%     This file is part of the multitaper package for R.
%
%     The multitaper package is free software: you can redistribute it and
%     or modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 2 of the 
%     License, or any later version.
%
%     The multitaper package is distributed in the hope that it will be 
%     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
%     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
%
%     If you wish to report bugs please contact the author. 
%     karim.rahim@gmail.com
%     112 Jeffery Hall, Queen's University, Kingston Ontario
%     Canada, K7L 3N6

\name{spec.mtm}
\alias{spec.mtm}
\title{Compute and plot multitaper spectrum estimates}
\description{
   Computes and plots adaptive or nonadaptive multitaper spectrum estimates from
   contiguous time series objects.
}
\usage{
spec.mtm(timeSeries, nw=4.0, k=7, nFFT="default", taper=c("dpss"), centre=c("Slepian"),
                     dpssIN=NULL, returnZeroFreq=TRUE, Ftest=FALSE, jackknife=FALSE,
                     jkCIProb=.95, maxAdaptiveIterations=100, plot=TRUE, na.action=na.fail,
                     returnInternals=FALSE, sineAdaptive=FALSE, sineSmoothFact=0.2,
                     dtUnits=c("default"), dT=NULL, ...) 
}
\arguments{
  \item{timeSeries}{A time series of equally spaced data, this can be
created by the ts() function where deltat is specified.}
  \item{nw}{nw a positive double precision number, the time-bandwidth
    parameter.}
  \item{k}{k a positive integer, the number of tapers, often 2*nw.}
  \item{nFFT}{This function pads the data before computing the fft. nFFT
indicates the total length of the data after padding.}
  \item{taper}{Choose between dpss-based multitaper (the default,'dpss') or sine taper
  method. In the case of the sine taper, parameter nw is useless, and both Ftest
  and jackknife are forced to FALSE. The sine taper also has two specific parameters
  below.}
  \item{centre}{The time series is centred using one of three methods: expansion
    onto discrete prolate spheroidal sequences ('Slepian'), arithmetic mean ('arithMean'), 
    trimmed mean ('trimMean'), or not at all ('none'). }
  \item{dpssIN}{Allows the user to enter a dpss object which has already
    been created. This can save computation time when Slepians with
 the same bandwidth parameter and same number of tapers are used repeatedly.}
  \item{returnZeroFreq}{Boolean variable indicating if the zeroth
frequency (DC component) should be returned for all applicable arrays.}
  \item{Ftest}{Boolean variable indicating if the Ftest result should be
 computed and returned.}
  \item{jackknife}{Boolean variable indicating if jackknifed confidence
intervals should be computed and returned.}
  \item{jkCIProb}{Decimal value indicating the jackknife probability for
calculating jackknife confidence intervals. The default returns a
95\% confidence interval.}
  \item{maxAdaptiveIterations}{Maximum number of iterations in the
adaptive multitaper calculation. Generally convergence is quick, and should
require less than 100 iterations.}
  \item{plot}{Boolean variable indicating if the spectrum should be plotted.}
  \item{na.action}{Action to take if NAs exist in the data, the default
    is to fail.}
  \item{returnInternals}{Return the weighted eigencoefficients, complex mean values, 
    and so on. These are necessary for extensions to the multitaper, including 
    magnitude-squared coherence (function mtm.coh in this package).}
  \item{sineAdaptive}{In the case of using the sine taper method, choose between
  non-adaptive and adaptive taper choice.}
  \item{sineSmoothFact}{The sine taper option has an inherent smoothing parameter that
  can be set between 0.01 and 0.5. Lower values indicate smaller amounts of smoothing.}
  \item{dtUnits}{Allows indication of the units of delta-t for accurate frequency axis
  labels.}
  \item{dT}{Time step for observations. If not in seconds, dtUnits should be set to indicate
  the proper units for plot labels.}
  \item{...}{Additional parameters, such as xaxs="i" which are passed to
    the plotting function. Not all parameters are supported.} 
}
\references{
Thomson, D.J (1982)
Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
Volume \bold{70}, Number 9, pp. 1055--1096.

Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications}
Cambridge University Press. 

Riedel, K.S. and Sidorenko, A. (1995)
Minimum bias multiple taper spectral estimation. \emph{IEEE Transactions on Signal Processing}
Volume \bold{43}, Number 1, pp. 188--195.
}
\examples{
# default behaviour, dpss tapers; dT and dtUnits set to ensure axis accuracy
data(willamette)
spec.mtm(willamette, nw=4.0, k=8, dT=1.0, dtUnits="month")
spec.mtm(willamette, nw=4.0, k=8, nFFT=2048, dT=1.0, dtUnits="month")

# if you have a ts object, you can skip the dT and dtUnits parameters
will.ts <- ts(data=willamette, start=1950.75, freq=12)
spec.mtm(will.ts, nw=4.0, k=8)

# using Sine Tapers
spec.mtm(will.ts, k=10, taper="sine", sineAdaptive=FALSE)
spec.mtm(will.ts, k=10, taper="sine", sineAdaptive=TRUE, 
         maxAdaptiveIterations=100, sineSmoothFact=0.05)
}
\keyword{multitaper}
