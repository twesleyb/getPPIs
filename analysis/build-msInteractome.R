#!/usr/bin/env Rscript

## Build a human interactome using the getPPIs package.

# Load the library.
library(getPPIs)

# Build the human interactome.
msInteractome <- getPPIs("HitPredict", species="mouse")

# Save as rda object.
save(msInteractome,file="msInteractome.rda")

