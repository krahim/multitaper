[![](https://www.r-pkg.org/badges/version/multitaper?color=green)](https://cran.r-project.org/package=multitaper)


test
[![](http://cranlogs.r-pkg.org/badges/grand-total/multitaper?color=green)](https://cran.r-project.org/package=multitaper)
[![](http://cranlogs.r-pkg.org/badges/last-month/multitaper?color=green)](https://cran.r-project.org/package=multitaper)
[![](http://cranlogs.r-pkg.org/badges/last-week/multitaper?color=green)](https://cran.r-project.org/package=multitaper)

multitaper
==========

R package implementing multitaper spectral estimation techniques used in time 
series analysis. In general the up-to-date version is on CRAN.

https://CRAN.R-project.org/package=multitaper 

General build instructions:

1) download and unzip to a folder called multitaper-master 
2) from the parent folder: 
  a) R CMD build multitaper-master/ 
  b) R CMD check multitaper_1.0-13.tar.gz (optional) 
  c) R CMD INSTALL multitaper_1.0-13.tar.gz 

Note: 

1) The version number 1.0-13, may change, and you will have to ensure 
the PATH variable is appropriately set. This will require gfortran. 
2) The version on CRAN is checked --as-cran using the latest development build** of R.

Please see the documentation for your Linux distribution, Xcode (or Fink) for Mac, 
or Rtools, currently maintained by Duncan Murdoch, on CRAN.

https://cran.r-project.org/bin/windows/Rtools/


** Note:
I uploaded the latest version in Oct. 2020 after checking --as-cran with the current release version of R.
Please consult CRAN policy for submitting packages.
