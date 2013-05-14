multitaper
==========

R package implementing multitaper spectral estimation techniques used in time series analysis.
This version may be slightly more update than the one on CRAN.

Build instructions for general Linux machine.

1) download and unzip to a folder called multitaper-master 2) from the parent folder: a) R CMD build multitaper-master/ b) R CMD check multitaper_1.0-6.tar.gz (optional) c) R CMD INSTALL multitaper_1.0-6.tar.gz 

Note: The version number, 1.0-6, may change, and you will have to ensure the PATH variable is appropriately set. This will require gfortran. Please see the documentation for your Linux distribution, Xcode (or Fink) for Mac, or Rtools, maintained by Duncan Murdoch, on CRAN for Windows.

Note: Mac users please see: http://r.research.att.com/tools/

