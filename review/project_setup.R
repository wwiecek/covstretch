# Basic packages
library(tidyverse)
library(ggpubr)
library(latex2exp)
library(shades)
library(grid)
library(assertthat)
# Optimization packages
library(nloptr)
library(odin)

# Basic plotting setting
theme_set(theme_minimal(base_size = 10))

# Parametrization and basic helper functions
source("setup/model_2vaccine.R")
source("optimisation-epi/single_run.R")
source("setup/parametrisation.R")

# Optimization functions
source("optimisation-epi/objective-functions.R")
source("optimisation-epi/nlopt_general.R")

# Visualization functions



