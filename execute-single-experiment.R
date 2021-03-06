#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("test-suite.R")

configs <- read.csv("configs/all_adv_configurations.csv")
test.single.config(args[[1]], configs, trials=100, all=FALSE)