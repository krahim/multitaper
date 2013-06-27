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

\name{demod.dpss}
\alias{demod.dpss}
\title{Computes complex demodulates using multiple taper techniques}
\description{
  Computes complex demodulate of a given series around a given central
  frequency using multiple taper techniques. Returns amplitude, phase,
  and complex demodulate.
}
\usage{
demod.dpss(x,centreFreq,NW,blockLen,stepSize=1,wrapphase=TRUE) 
}
\arguments{
  \item{x}{Time series, required to be contiguous.}
  \item{centreFreq}{Frequency around which to demodulate.}
  \item{NW}{Parameter controlling time-bandwidth.}
  \item{blockLen}{Length of sub-block to use; demodulate is computed on 
  each block in turn.}
  \item{stepSize}{This is a proposed option that sets the index step
    size between blocks. Currently this must be set to 1 and
    changes in step size have not been implemented.}
  \item{wrapphase}{If true, routine wraps phases around +/-360 degree 
  boundaries.}
}
\references{
Thomson, D.J. (1995). The Seasons, Global Temperature, and Precession. 
\emph{Science}, Volume 268, pp. 59--68.
}
\examples{
data(CETmonthly)
nJulOff <- 1175
xd <- ts(CETmonthly[,"temp"],deltat=1/12)
demodYr <- demod.dpss(xd,centreFreq=1,NW=3,blockLen=120,stepSize=1)
phase <- demodYr$phase
offsJul <- 3*360/365 
phaseAdj <- phase
phaseAdj[1:nJulOff] <- phase[1:nJulOff] + offsJul
yr <- (time(xd)+1658)[1:length(phase)]
plot(yr, phaseAdj, type="l", lwd=2,
     ylab="Phase of the Year in Degrees",
     xlab="Gegorian calender date")
lines((1:nJulOff)/12+1659, phase[1:nJulOff], col="red")
fit <- lm( phaseAdj ~ yr)
abline(fit, lty=2, col="blue")
cat(paste("Precession Estimate: ",fit$coef[2]*60*60,digits=6," (arcseconds/yr)\n",sep=""))
}
\keyword{math}
