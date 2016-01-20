# Clear All
rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(grid)
library(reshape2)
library(parallel)

cores <- detectCores()

source("../indiff/indifference.r")
#indiff contains the indifference points of the 10 lotteries

# Grab our data
load("../data/agg-dat/Agg1Mil-S10k.Rda")
MM <- tbl_df(MM)

# Don't always want to use it all, there is tons of data
sam.prop <- .1

# Bounds for the means of data
lbound <- -1.8
ubound <- 1.4
lbound <- -1.9
ubound <- 1.55

# Certain data have rm inside the relevent range of the HL-MPL, others don't.
# Make the distiction
IN  <- MM %>%
	filter(rm > lbound, rm < ubound) %>%
	sample_frac(sam.prop) %>%
	mutate(id = 1:n())


