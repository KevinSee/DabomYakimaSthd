
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DabomYakimaSthd

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KevinSee/DabomYakimaSthd/master?urlpath=rstudio)

This repository contains the data and code for running the **D**am
**A**dult **B**ranch **O**ccupancy **M**odel
([DABOM](https://github.com/KevinSee/DABOM))for steelhead returning to
the Yakima River.

### How to download or install

You can download the compendium as a zip from from this URL:
<https://github.com/KevinSee/DabomYakimaSthd/archive/master.zip>

Or you can install this compendium as an R package, DabomYakimaSthd,
from GitHub with:

``` r
# install.packages("devtools")
remotes::install_github("KevinSee/DabomYakimaSthd")
```

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.

# Project Notes

No hatcheries in the Yakima, but there may be some hatchery strays into
the basin

  - Run model using all fish (including hatchery strays)
  - Run model using only wild fish
  - Too few hatchery fish to run one model with different movement rates
    between hatchery and wild fish

The trap rate does change throughout the season, perhaps 2-4 times.

Trapping became consistent starting with spawn year 2011-2012.

  - Starts in early Sept
  - Ends in mid May
  - Encompasses ~ 95% of the run

Spawn year constitutes all steelhead who cross Prosser dam from July 1
to June 30

  - Same definition at Roza dam

Can compare model estimates with dam counts at Roza dam in upper Yakima

## Site information

TTN is due to go online next year There should be an array on Marza
Drain (MD+) in summer of 2020 Fish can get into MD+ without going past
TOP (but usually don’t?)

WNS: Yakima has the data for this site; not on PTAGIS. May not be
available for every year

TEAN = LMT

  - Combine older sites (e.g. NFTEAN, TEANAR, TEANMF, TEANWF)

LNR may be a single array

There are some USFW arrays between TOP and TP2 / SM1 Yakima Nation will
investigate whether there is data for those, and what the configuration
actually is
