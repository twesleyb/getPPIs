#------------------------------------------------------------------------------
## Installation.
#------------------------------------------------------------------------------

devtools::install_github("twesleyb/getPPIs")

# Load the library.
library(getPPIs)

#------------------------------------------------------------------------------
## Download, compile, annotate, and subset PPIs in a single step.
#------------------------------------------------------------------------------

# Load an example dataset.
data(iPSD)

# Get PPIs amongst iPSD proteins (mouse).
hitpredict <- getPPIs(organism = "HitPredict", 
		      mygenes = iPSD, 
		      taxid = 10090, 
	              methods = "all", 
		      cutoff = 0,
		      downloads = "../downloads",
		      keepdata = TRUE,
		      saveppis = TRUE)


#------------------------------------------------------------------------------
## Building an interactome from scratch.
#------------------------------------------------------------------------------

# Defaults.
downloads <- "../downloads"
keepdata <- TRUE

# Download all PPIs in the HitPredict database.
hitpredict <- getHitPredict(organism="HitPredict", downloads, keepdata)

# Map genes to mouse homologs.
hitpredict <- getHomoloGene(hitpredict, taxid = 10090)

# Annotate hitpredict data with more descriptive method names.
# Keep all experimental evidence, no confidence score cutoff.
hitpredict <- getMethods(hitpredict, methods="all", cutoff = 0)


