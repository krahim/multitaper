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

\name{sineTaper}
\alias{sineTaper}
\title{Computes sine tapers}
\description{
  Computes sine tapers for use in transfer function estimation and 
  plotting. Not called from within spec.mtm.
}
\usage{
sineTaper(n, k)
}
\arguments{
  \item{n}{The data as a vector or as a time series.}
  \item{k}{The Slepian bandwidth parameter, typically between 2.0 and 6.0.}
}
\references{
  Riedel, K.S. and Sidorenko, A. (1995)
  Minimum bias multiple taper spectral estimation. \emph{IEEE Transactions
  on Signal Processing}, Volume \bold{43}, Number \bold{1}, pp. 188--195.
}
\keyword{math}
