\name{multitaper-package}
\alias{multitaper}
\alias{multitaper.package}
\docType{package}
\title{
    multitaper: Multitaper spectral analysis tools
}
\description{
  Provides general functions and utilities for multitaper spectral analysis 
  using the Discrete Prolate Spheroidal (Slepian) sequences and the Sine tapers. 
  The packages implements the adaptively weighted 
  multitaper spectral estimate, with nonparametric jackknife variance estimate, as
  well as Thomson's Harmonic F-test, a multitaper version of the bivariate (magnitude 
  squared and phase) coherence with jackknife variance. The package uses the filtering
  properties of the Slepian sequences to perform complex demodulation, and it provides 
  an interface similar to the complex demodulation function provided with S-plus. 
  The function calls have been designed to be similar to R function calls for spectral 
  estimation, and the package includes several spectral analysis tools. 
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

