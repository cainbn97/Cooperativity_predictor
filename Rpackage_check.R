#!/usr/bin/Rscript

## Author: cainu5
##
## Check for R package dependencies
##
## 01/20/22

packages = c("readr", "outliers",
             "chisq.posthoc.test")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      print(paste(x, 'is not currently installed.', sep = " "), quote = FALSE)
	  # install.packages(x, dependencies = TRUE)
      # library(x, character.only = TRUE)
    }
	else
	{	
	}
  }
)