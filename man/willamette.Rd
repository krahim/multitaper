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


% Dataset Rd file added in response to bug from Kirk Hornik Aug 12,
% 2010.

\name{willamette}
\alias{willamette}
\title{Willamette River time series}
\description{
   Willamette River time series. Each point represents the log of the
   average daily flow over a one month period from October, 1950, to
   August 1983. The sampling time is 1/12 year and the Nyquist frequency
   is 6 cycles per year. The data is from the companion code to \dQuote{Spectral
   Analysis for the Physical Applications} (1993) and was originally compiled by
   the US Geological Survey. 
}
\usage{willamette}
\format{A vector containing 395 observations}   
\references{
Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications.}
Cambridge University Press.
}
\examples{
data(willamette)
# time series object, January = year.0, December = year.917
will.ts <- ts(data=willamette, start=(1950+9/12), freq=12)
}
\keyword{datasets}
