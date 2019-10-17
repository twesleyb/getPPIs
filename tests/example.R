# Usage Example.

# Installation.
devtools::install_github("twesleyb/getPPIs")

# Load the library.
library(getPPIs)

# Load an example dataset.
data(iPSD)
head(iPSD)

# getPPIs
getPPIs(organism = "HitPredict", taxid = 10090, 

