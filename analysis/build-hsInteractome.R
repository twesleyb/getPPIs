#!/usr/bin/env Rscript

library(getPPIs)

# getPPIs does the previous three steps.
hitpredict <- getPPIs("HitPredict", taxid=9606)

