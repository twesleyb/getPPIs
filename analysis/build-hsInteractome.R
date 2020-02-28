#!/usr/bin/env Rscript

## Build a human interactome using the getPPIs package.

library(getPPIs)

hsInteractome <- getPPIs("HitPredict", taxid=9606)

save(hsInteractome,file="hsInteractome.rda")

