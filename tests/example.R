# Usage Example.

# Installation.
devtools::install_github("twesleyb/getPPIs")

# Load the library.
library(getPPIs)

# Load an example dataset.
data(iPSD)
head(iPSD)

# getPPIs
hitpredict <- getPPIs(organism = "HitPredict", 
		      mygenes = iPSD, 
		      taxid = 10090, 
	              methods = "all", 
		      cutoff = 0,
	       	      file = "",
		      downloads = "../downloads"),
		      keepdata = TRUE)

