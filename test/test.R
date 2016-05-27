#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(correctMultiCount)
library(rbenchmark)

data(baseCount)
data(multiCount)
df=correctCounts(baseCount, multiCount)
b = benchmark(df=correctCounts(baseCount, multiCount))
