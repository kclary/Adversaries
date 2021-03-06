#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("test-suite.R")

configs <- read.csv("configs/all_adv_configurations_rw.csv")
test.single.config(args[[1]], configs, trials=10, all=FALSE)
