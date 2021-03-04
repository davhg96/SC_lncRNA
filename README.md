# Single cell lncRNA analysis
# Reproducibility
To reproduce this analysis first create an R project in your local machine and pull the files from this GitHub repository. Reproducibility in this project is managed by the package renv. To start working on this firs install the package:
```
if (!requireNamespace("remotes"))
  install.packages("remotes")

remotes::install_github("rstudio/renv")
```
At this point you should have all necessary files downloaded in a directory in your local machine, check that you can see a file called `rent.lock`. Set the working directory to that folder and run `renv::restore()` to install all the necessary packages and versions to reproduce this analysis.

# SRC: Analysis  files
All these scripts will access the data from a data directory located in your working directory and will create an output directory with all the output automatically.

- Functions.R: Will call all the necessary packages and load custom functions to the working environment.
- dop_ndop_fgfnoVLMC.R: Contains the differential expression analysis and plotting.
- CorrelationAndieplots.R: Will run differential expression and will plot the insertion distribution and correlation of lncRNAs and protein coding genes, as well as lncRNA length distribution.