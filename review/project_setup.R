library(tidyverse)
library(ggpubr)
library(latex2exp)
library(shades)
library(grid)

library(nloptr)
library(odin)

theme_set(theme_minimal(base_size = 10))

source("setup/model_2vaccine.R")
source("optimisation-epi/single_run.R")
source("setup/parametrisation.R")
source("optimisation-epi/objective-functions.R")



