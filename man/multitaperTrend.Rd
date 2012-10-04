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

\name{multitaperTrend}
\alias{multitaperTrend}
\title{Estimate Linear Trend using Multitaper Techniques}
\description{
   Estimate linear trend using inverse spectrum estimation, with the 
   spectrum being computed via multitaper. This technique has improved
   spectral properties when compared to the least-squares approach. 
   Returned values from this function include the intercept, slope,
   and centered time array.
}
\usage{
multitaperTrend(xd, B, dT, t.in)
}
\arguments{
  \item{xd}{Contiguous time series to be detrended.}
  \item{B}{Bandwidth to use in estimating trend in physical units; corresponds
  to NW via equation NW=BT, where N and W are the usual Slepian definitions,
  and T is the total time elapsed, i.e. T = N*dT.}
  \item{dT}{Time step for series xd, also used in computing T.}
  \item{t.in}{Time array, used in accurately estimating the slope.}
}
\examples{
x <- 1:101
y <- 1.0 + 0.5*(x) + rnorm(n=101,mean=0,sd=2)
vars <- multitaperTrend(xd=y, B=0.05, dT=1.0, t.in=x)
plot(x,y,type="l")
lines(x,vars[[1]]+vars[[2]]*vars[[3]],type="l",col="red")
}
\keyword{math}
