# Overview
The following files allow one to reproduce analyses in the manuscript entitled "Species-level studies do not upscale to community-wide plant-soil feedbacks".
# Repo Contents
## Code: R code
This script contains all the data analysis and visualization code of Fig.1, Fig.2, Fig.3 and Fig.S2.
## Data  
This script contains all the data of Fig.1, Fig.2, Fig.3 and Fig.S2.

# System Requirements
## Hardware Requirements
R package requires only a standard computer with enough RAM to support the in-memory operations.
## Software Requirements
### OS Requirements
The package development version is tested on Windows and Mac systems. 

R version 4.2.2 or higher. 

# Installation Guide
## Intallation of R and Rstudio
R: https://www.r-project.org/

Rstudio: https://www.rstudio.com/

## Note: 

Install from Windows: follow the default path and click "Next".

Rstudio and R associated Settings (Rstudio->Tools->Global options)

The package should take approximately 40-60 seconds to install with vignettes on a recommended computer.

## Install software supporting R: JAVA, Rtools (windows system), Xcode(Mac system)

a) If you do not have Java installed on your computer, 
please download: https://www.java.com/zh-CN/download/ Download and install.

b) If you do not have Rtools installed on your computer, 
please download RTools from https://cran.r-project.org/bin/windows/Rtools/ to download for installation (please choose according to their own version of a R corresponding Rtools version), and check whether the installation is successful.

c) If your computer is a Mac system, 
from https://developer.apple.com/cn/xcode/, please download and install.

The package should take approximately 40-60 seconds to install with vignettes on a recommended computer.

## Package Installation
Users should install the following packages prior to library.
~~~
install.packages(c('readxl', 'openxlsx', 'ggfortify', 'ggforce', 'dplyr', 'ggrepel', 'ggpubr', 'export', 'ggspatial', 'ggplot2','cowplot', 'sf', 'multcomp', 'tidyverse',
 'car', 'lme4', 'vegan', 'FUNGuildR', 'lmerTest', 'ggbeeswarm', 'reshape2', 'DHARMa', 'glmmTMB', 'emmeans', 'broom.mixed', 'MuMIn', 'ggridges', 'party', 'caret', 'pheatmap'))
~~~

Note: ggepi & patchwork can be installed via devtools 
~~~
 devtools::install_github("lwjohnst86/ggepi")
 devtools::install_github("thomasp85/patchwork")
~~~

The package should take approximately 50-60 seconds to install with vignettes on a recommended computer.

# Demo
## Instructions to run on data: example-Fig.3 code
~~~
 
~~~
## Expected output
~~~
 
~~~
## Note: 
The dataset stored in this repository is same to the dataset in figshare( ).
Running each line of code on a recommended computer takes about 1-5 seconds, except for "Step 2: Predictability comparison among random forest models" in Fig.2's code, which takes 20-30 minutes.

# License
This project is covered under the GNU GENERAL PUBLIC LICENSE Version 2 (GPL-2).
