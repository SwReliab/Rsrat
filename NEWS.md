# Rsrat 1.5.0

* Refactor: Rsrat removes the metrics models.

# Rsrat 1.3.4

* Bugfix: Fix an argument 'selection' on fir.srm.nhpp.

# Rsrat 1.3.3

* Bugfix: Change a type of the detection/correction matrix.
* The argument `data` in mvfplot and dmvfplot accepts Rsrat.faultdata directly.

# Rsrat 1.3.2

* Bugfix: LLF for Exp with time data
* add get.issues for GitHub
* Matrix package is imported to represent the detection/correction matrix

# Rsrat 1.3.1

* Bugfix: LLF for Pareto with time data

# Rsrat 1.3.0

* Add penalized for poireg

# Rsrat 1.2.6

* Fix the authentication for JIRA

# Rsrat 1.2.5

* Add gfortran lib
* Use httr instead of RCurl

# Rsrat 1.2.4

* Remove gfortran lib

# Rsrat 1.2.3

* Fix bugs on truncation models
  - em for tnorm has been changed
  - ems for logis and xvmax/min use optim function at the M-step
* Fix a bug of jira
* The option printflag has been changed to `trance`.
* Fix a bug for poireg according to the changes of ems for truncation models

# Rsrat 1.2.2

* Fix a bug for winR

# Rsrat 1.2.1

* Fix a bug for faultdata, mvfplot, dmvfplot and fit.srm.nhpp
* Fix DESCRIPTION

# Rsrat 1.2.0

* Add scraping (JIRA)

# Rsrat 1.1.2

* Fix a bug for faultdata and mvfplot
* Change emsteps with Rcpp

# Rsrat 1.1.1

* Add data (dacs)
* Change interface of faultdata to use dataframe
* Modify some examples

# Rsrat 1.1.0

* Add dmvf method in NHPP
* Add mvfplot and dmvfplot

# Rsrat 1.0.1

* Fix a document of fit.srm.poireg
* Fix a bug: srms list is changed after executing fit.srm.poireg

# Rsrat 1.0.0

* Add poisson-regression-based model to handle s-metrics.

# Rsrat 0.10.0

* Add logistic-based model to handle d-metrics.
* Add fit.srm.nhpp to estimate NHPP models.

# Rsrat 0.9.3

* Add functions to NHPP (llf, inv_mvf, imtbf, cmtbf)
* Bug fix: pgumbel.min

# Rsrat 0.9.2

* Modify the document files.

# Rsrat 0.9.1

* Added a `NEWS.md` file to track changes to the package.
