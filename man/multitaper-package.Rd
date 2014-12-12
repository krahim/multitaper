\name{multitaper-package}
\alias{multitaper}
\alias{multitaper.package}
\docType{package}
\title{
    multitaper: Multitaper spectral analysis tools
}
\description{
 This package implements multitaper spectral analysis using discrete prolate spheroidal (Slepian) sequences (dpss) and the sine tapers. The Slepian sequences are efficiently and accurately obtained using the tridiagonal matrix method with LAPACK functions. The package implements Thomson's adaptive weighted multitaper spectral estimate, the jackknife standard error estimate, Thomson's Harmonic F-test, and coherence (magnitude squared and phase) calculated with Slepian sequences. Jackknife variance estimates are also provided for coherence estimates. The package uses the filtering properties of the Slepian sequences to perform complex demodulation offering a similar interface to the complex demodulation function provided with S-plus. The function calls have been designed to be consistent with R and several additional spectral analysis tools are included. This package is, and has been, used by David Thomson's research group at Queen's University since 2010.
}
\details{
\tabular{ll}{
Package: \tab multitaper\cr
Type: \tab Package\cr
Version: \tab 1.0-9\cr
Date: \tab 2014-10-05\cr
License: \tab GPL-2 \cr
}
}
\author{
Karim Rahim <karim.rahim@gmail.com> and Wesley Burr <wesley.burr@gmail.com>

Maintainer: Karim Rahim <karim.rahim@gmail.com> 
}
\references{

  Rahim, K.J. and Burr W.S. and Thomson, D.J. (2014)
  Appendix A: Multitaper R package in "Applications of Multitaper Spectral
  Analysis to Nonstationary Data" PhD diss., Queen's Univieristy, 
  pp 149--183. http://hdl.handle.net/1974/12584

  Thomson, D.J. (1982) Spectrum estimation and harmonic analysis. 
  \emph{Proceedings of the IEEE} Volume \bold{70}, number 9, pp. 1055--1096.

  Percival, D.B. and Walden, A.T. (1993)
  \emph{Spectral analysis for physical applications}
  Cambridge University Press. 

  Riedel, K.S. and Sidorenko, A. (1995)
  Minimum bias multiple taper spectral estimation. \emph{IEEE Transactions on Signal Processing}
  Volume \bold{43}, Number 1, pp. 188--195.
}
\keyword{spectrum}
\keyword{psd}
\keyword{sdf}
\keyword{multitaper}
\keyword{mtm}

