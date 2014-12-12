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


\name{mtm.coh}
\alias{mtm.coh}
\title{Compute and plot the multitaper magnitude-squared coherence.}
\description{
   Computes and plots the adaptive multitaper spectrum estimate.
}
\usage{
mtm.coh(mtm1, mtm2, fr=NULL, tau=0, plot=TRUE, ...) 
}
\arguments{
  \item{mtm1}{An object created with spec.mtm(... ,returnInternals=TRUE).}
  \item{mtm2}{An object created with spec.mtm(... ,returnInternals=TRUE). 
    Note mtm1 and mtm2 must be created with the same frequency resolution. 
    They both must have the same values for nFFT and returnZeroFreq.} 
  \item{fr}{The frequency values for the mtm object. This can be null by
    default (which results in computation for the full frequency range) or 
    it can be a subset of frequency values.}
  \item{tau}{Phase-correction factor, if known.}
  \item{plot}{Boolean value indicating if a plot should be drawn.}
  \item{...}{Additional parameters, such as xaxs="i" which are passed through to
    the plotting function.} 
}
\references{
Thomson, DJ (1991)
Jackknifed error estimates for spectra, coherences, and transfer
functions, \emph{Advances in Spectrum Estimation} 58--113.

Thomson, D.J (1982)
Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
Volume \bold{70}, number 9, pp. 1055--1096.

Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications}
Cambridge University Press. 
}
\examples{
data(HadCRUTnh)
data(mlco2)
spec1 <- spec.mtm(HadCRUTnh, nw=5.0, k=8, plot=FALSE,
    returnInternals=TRUE, dtUnits="month", deltat=1.0)
spec2 <- spec.mtm(mlco2, nw=5.0, k=8, plot=FALSE, returnInternals=TRUE,
    dtUnits="month", deltat=1.0)
resCoh <- mtm.coh(spec1, spec2, plot=FALSE)
plot(resCoh)
plot(resCoh, cdfQuantilesTicks=1-10^(-(6:12)))
}
\keyword{multitaper}
