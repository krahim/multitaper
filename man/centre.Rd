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


\name{centre}
\alias{centre}
\title{Centres (converts to zero-mean) the time series.}
\description{
  Centres the data using an expansion on the Slepian sequences if the bandwidth
  parameter (nw) and number of tapers (k) is specified, otherwise
  subtracts the mean or robust trimmed mean.
}
\usage{
centre(x, nw = NULL, k = NULL, deltaT = NULL, trim = 0)

}
\arguments{
  \item{x}{The data as a vector or as a time series.}
  \item{nw}{The Slepian bandwidth parameter, typically between 2.0 and 6.0.}
  \item{k}{The number of Slepian tapers used, often 2*nw.}
  \item{deltaT}{Parameter required if the data is a vector and not a time
    series, and only for the Slepian case.}
  \item{trim}{[only used if nw and k are not specified]
    The fraction (0 to 0.5) of observations to be trimmed from
    each end of `x' before the mean is computed.  Values of trim
    outside that range are taken as the nearest endpoint.}
}
\references{
  Thomson, D.J (1982)
  Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
  Volume \bold{70}, number 9, pp. 1055--1096.
  
  Slepian, D. (1978)
  Prolate spheroidal wave functions, Fourier analysis, and
  uncertainty. V--The discrete case. \emph{Bell System Technical Journal}
  Volume \bold{57}, pp. 1371--1430.
}
\examples{
data(willamette)
cent.Slepian <- centre(willamette, nw=4, k=8, deltaT=1)
cent.Trim <- centre(willamette, trim=0.2)
}
\keyword{math}
