# Getting started {#start}

## Checklist (before the course)

We will provide a temporary access to a cloud computing environment
that readily contains the available software packages. Instructions to
access the environment will be sent to the registered participants.

Setting up the system on your own computer is not required for the
course but it can be useful for later use. The required software:

* [R (version >4.1.0)](https://www.r-project.org/) 

* [RStudio](https://www.rstudio.com/products/rstudio/download/);
  choose "Rstudio Desktop" to download the latest version. Optional
  but preferred. For further details, check the [Rstudio home
  page](https://www.rstudio.com/).

* Install and load the required R packages (see Section \@ref(packages))

* After a successful installation you can start with the
  case study examples in this training material


## Support and resources

 * Online support on installation and other matters, join us at [Gitter](https://gitter.im/microbiome/miaverse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

 * Additional reading and online material are listed in Section \@ref(material).


**You can run the workflows by simply copy-pasting the examples.** For
further, advanced material, you can test and modify further examples
from the online book, and apply these techniques to your own data.

## Installing and loading the required R packages {#packages}

This section shows how to install and load all required packages into
the R session. Only uninstalled packages are installed.


```r
# List of packages that we need from cran and bioc 
cran_pkg <- c("BiocManager", "bookdown", "dplyr", "ecodist", "ggplot2", 
              "gridExtra", "kableExtra",  "knitr", "scales", "vegan", "matrixStats")
bioc_pkg <- c("yulab.utils","ggtree","ANCOMBC", "ape", "DESeq2", "DirichletMultinomial", "mia", "miaViz")

# Get those packages that are already installed
cran_pkg_already_installed <- cran_pkg[ cran_pkg %in% installed.packages() ]
bioc_pkg_already_installed <- bioc_pkg[ bioc_pkg %in% installed.packages() ]

# Get those packages that need to be installed
cran_pkg_to_be_installed <- setdiff(cran_pkg, cran_pkg_already_installed)
bioc_pkg_to_be_installed <- setdiff(bioc_pkg, bioc_pkg_already_installed)
```


```r
# If there are packages that need to be installed, installs them from CRAN
if( length(cran_pkg_to_be_installed) ) {
   install.packages(cran_pkg_to_be_installed)
}
```


```r
# If there are packages that need to be installed, installs them from Bioconductor
if( length(bioc_pkg_to_be_installed) ) {
   BiocManager::install(bioc_pkg_to_be_installed, ask = F)
}
```
 
Now all required packages are installed, so let's load them into the session.
Some function names occur in multiple packages. That is why miaverse's packages
mia and miaViz are prioritized. Packages that are loaded first have higher priority.


```r
# Reorders bioc packages, so that mia and miaViz are first
bioc_pkg <- c(bioc_pkg[ bioc_pkg %in% c("mia", "miaViz") ], 
              bioc_pkg[ !bioc_pkg %in% c("mia", "miaViz") ] ) 

# Loading all packages into session. Returns true if package was successfully loaded.
loaded <- sapply(c(bioc_pkg, cran_pkg), require, character.only = TRUE)
as.data.frame(loaded)
```

```
##                      loaded
## mia                    TRUE
## miaViz                 TRUE
## yulab.utils            TRUE
## ggtree                 TRUE
## ANCOMBC                TRUE
## ape                    TRUE
## DESeq2                 TRUE
## DirichletMultinomial   TRUE
## BiocManager            TRUE
## bookdown               TRUE
## dplyr                  TRUE
## ecodist                TRUE
## ggplot2                TRUE
## gridExtra              TRUE
## kableExtra             TRUE
## knitr                  TRUE
## scales                 TRUE
## vegan                  TRUE
## matrixStats            TRUE
```



