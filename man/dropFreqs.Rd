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

\name{dropFreqs}
\alias{dropFreqs}
\alias{dropFreqs.default}
\alias{dropFreqs.spec}
\alias{dropFreqs.mtm}
\alias{dropFreqs.mtm.coh}
\title{Truncate mtm or mtm.coh Objects in Frequency}
\description{
   A utility function to truncate the frequencies in a spectral
   estimate. This utility is used before calling plot(), to increase the
   visual frequency resolution of a plot by truncating frequencies outside a
   particular band of interest. This function is not a filter, but rather
   a utilty to allow R to 'zoom' a spectrum plot to a certain frequency band.
}
\usage{
dropFreqs(spec, minFreq, maxFreq)
}
\arguments{
  \item{spec}{A spectrum object 'obj', of class spec, mtm, or mtm.coh.}
  \item{minFreq}{The lower bound for the frequency band to be retained, 
  in the same units as the obj$freq array.}
  \item{maxFreq}{The upper bound for the frequency band to be retained,
  also in the same units as the obj$freq array.}
}
\examples{
data(willamette)
mtm1 <- spec.mtm(willamette, nw=4.0, k=8, plot=FALSE, dT=1.0, dtUnits="month")
mtm2 <- dropFreqs(mtm1, 0.1, 0.4)
plot(mtm2)

# another option
plot(dropFreqs(mtm1, 0.1, 0.4))

# using sine tapers
mtm.sine <- spec.mtm(willamette, k=10, plot=FALSE, dT=1.0, dtUnits="month", 
                     taper="sine", sineAdaptive=FALSE, sineSmoothFact=0.05)
plot(dropFreqs(mtm.sine, 0.1, 0.4))                     
}
\keyword{multitaper}
