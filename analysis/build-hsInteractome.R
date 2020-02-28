#!/usr/bin/env Rscript

## Build a human interactome using the getPPIs package.

# Load the library.
library(getPPIs)

# Build the human interactome.
hsInteractome <- getPPIs("HitPredict", species="human")

# Save as rda object.
save(hsInteractome,file="hsInteractome.rda")

