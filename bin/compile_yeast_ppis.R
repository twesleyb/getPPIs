#!/usr/bin/env Rscript

# Compile yeast interactions from YeastMine.

# Load renv.
here <- getwd()
root <- dirname(here)
renv::load(root)

# Load the data.
devtools::load_all()

data(yeast_homologs) # homologs
data(yeast_interactions) # interactions

