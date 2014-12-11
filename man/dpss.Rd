%     The multitaper R package
%     Multitaper and spectral analysis package for R
%     Copyright (C) 2013 Karim Rahim 
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
%     Jeffery Hall, Queen's University, Kingston Ontario
%     Canada, K7L 3N6

\name{dpss}
\alias{dpss}
\title{Compute Discrete Prolate Spheroidal Sequences}
\description{
   Compute Discrete Prolate Spheroidal (Slepian) Sequences for use as
   tapers or other applications. This function uses the tridiagonal method and exploits
   symmetry. Note the odd order tapers are normalized so that the slope at the
   centre is positive in accordance with Slepian (1978) and Thomson (1982). This
   differs from Percival and Walden (1993). This code follows section
   (8.3) of Percival and Walden (1993) using LAPACK function calls Anderson (1999).
}
\usage{
dpss(n,k,nw, returnEigenvalues=TRUE)
}
\arguments{
  \item{n}{A positive integer, typically the non-zero-padded length of the time series.}
  \item{k}{A positive integer, the number of tapers, often 2*nw for spectrum estimation purposes.}
  \item{nw}{A positive double-precision number, the time-bandwidth
    parameter.}
  \item{returnEigenvalues}{If true the appropriate eigenvalues are
    calculated and returned using the function dpssToEigenvalues. If
    FALSE, the eigenvalues returned are from the LAPACK function DSTEBZ
    using the tridiagonal. See section 8.3 of Percival and Walden (1993), or equation (13) in Slepian (1978).}
}
\value{
  \item{v}{A n by k matrix of Slepian Sequences. Each column represents
    the Slepian sequence of order k-1.}
  \item{eigen}{A length k vector of eigenvalues corresponding to
    equation (13) in Slepian (1978), or the eigenvalues of the input
    tridiagonal matrix returned from the internal call to the LAPACK function DSTEBZ.}
  }
\references{
Anderson, E. (1999).
\emph{LAPACK Users' guide (Vol. 9).}
Siam.

Percival, D.B. and Walden, A.T. (1993)
\emph{Spectral analysis for physical applications.}
Cambridge University Press.

Slepian, D. (1978)
Prolate spheroidal wave functions, Fourier analysis, and
uncertainty. V--The discrete case. \emph{Bell System Technical Journal}
Volume \bold{57}, pp. 1371--1430

Thomson, D.J (1982)
Spectrum estimation and harmonic analysis. \emph{Proceedings of the IEEE}
Volume \bold{70}, number 9, pp. 1055--1096.
}
\examples{
dpss(10,4,4.0)
dpss(100,8,5.0)
}
\keyword{math}
