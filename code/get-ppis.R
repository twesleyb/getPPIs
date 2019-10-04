#!/usr/bin/env Rscript
# Compiling Protein-Protein interctions from HitPredict.

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Imports.
suppressPackageStartupMessages({
	library(AnnotationDbi)
	library(org.Mm.eg.db)
	library(org.Hs.eg.db)
	library(org.Rn.eg.db)
	library(TBmiscr) 
})

# Directories.
here <- getwd()
if (strsplit(osVersion," ")[[1]][1]=="Windows") {
  here <- "D:/projects/iPSD-PTM/data"; setwd(here)
}
rootdir <- dirname(here)
datadir <- file.path(rootdir,"data")
downloads <- file.path(rootdir,"downloads")

#-------------------------------------------------------------------------------
## Load HitPredict interactions.
#-------------------------------------------------------------------------------

# Download HitPredict data. This will take several minutes.
url <- "http://www.hitpredict.org/download/HitPredict_interactions.txt.tgz"
gzfile <- file.path(downloads, "HitPredict_interactions.txt.tgz")
if (!file.exists(gzfile)) {
	message("Downloading PPIs from HitPredict.org...")
	download.file(url, gzfile)
} else {
	message("file already exists!")
}

# Unzip and read data.
untar(gzfile, exdir=downloads)
myfile <- file.path(datadir, "HitPredict_interactions.txt")
hitpredict <- data.table::fread(myfile, header = TRUE, skip = 5)
unlink(myfile)

# We need to insure that genes are mapped to a stable, unique identifier.
# First, replace blanks with NA.
hitpredict[hitpredict == ""] <- NA

# Subset human, mouse, and rat data.
taxids <- c(9606, 10090, 10116)
hitpredict <- subset(hitpredict, hitpredict$Taxonomy %in% taxids)

# Seperate rows with multiple Entrez values.
hitpredict <- separate_rows(hitpredict,Entrez1,sep=";")
hitpredict <- separate_rows(hitpredict,Entrez2,sep=";")

# Seperate rows with multiple Ensembl values.
hitpredict <- separate_rows(hitpredict,Ensembl1,sep=",")
hitpredict <- separate_rows(hitpredict,Ensembl2,sep=",")

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes in the HitPredict database:",
  table(c(is.na(hitpredict$Entrez1), is.na(hitpredict$Entrez2)))[2]
))

#-------------------------------------------------------------------------------
## Collect mapping data from Uniprot.
#-------------------------------------------------------------------------------

# Define a function that will download Uniprot mapping data for a given organism.
build_uniprotDB <- function(taxid){
  if (taxid == 9606) {
    message("Collecting human Uniprot mapping data!")
    db <- "HUMAN_9606"
  } else if (taxid == 10090) {
    message("Collecting mouse Uniprot mapping data!")
    db <- "MOUSE_10090"
  } else if (taxid == 10116) {
    message("Collecting rat Uniprot mapping data!")
    db <- "RAT_10116"
  } else {
    stop("Please provide valid taxid c(9606, 10090, 10116).")
  }
  # Which database should be downloaded?
  baseurl <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/"
  extension <- "_idmapping_selected.tab.gz"
  url <- paste0(baseurl,db,extension)
  gzfile <- file.path(datadir,basename(url))
  destfile <- tools::file_path_sans_ext(gzfile)
  # Download
  if (!file.exists(gzfile)) {
    message(paste0("Downloading ", basename(gzfile),"..."))
    download.file(url, gzfile, quiet = TRUE)
  } else {
    message("file already exists!")
  }
  # Unzip
  if (!file.exists(destfile)) {
    message(paste0("Unzipping ", basename(destfile),"..."))
    R.utils::gunzip(gzfile, destname = destfile)
  } else {
    message("file already exists!")
  }
  # Read data into R.
  uniprotDB <- fread(destfile, sep="\t", skip =1)
  colnames(uniprotDB) <- c("UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)",
                           "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90",
                           "UniRef50", "UniParc", "PIR", "NCBI-taxon", "MIM",
                           "UniGene", "PubMed", "EMBL", "EMBL-CDS", "Ensembl",
                           "Ensembl_TRS","Ensembl_PRO","Additional PubMed")
  # Remove unzipped file!
  unlink(destfile)
  return(uniprotDB)
}

# Get Uniprot data for all three organisms.
uniprotDB <- lapply(as.list(taxids), build_uniprotDB)
names(uniprotDB) <- as.character(taxids)

# Add Taxonomy id column and combine into a single data frame.
for (i in 1:length(uniprotDB)) {
  uniprotDB[[i]] <- tibble::add_column(uniprotDB[[i]],Taxonomy = taxids[i], .after=22)
}
uniprotDB <- do.call(rbind, uniprotDB)

# Seperate rows with multiple Ensembl IDs.
uniprotDB <- separate_rows(uniprotDB,Ensembl,sep=";")

# Create a Ensembl to Entrez map.
uniprotDB <- data.table(UniprotAC = uniprotDB$`UniProtKB-AC`,
                        UniprotID = uniprotDB$`UniProtKB-ID`,
                        Entrez = uniprotDB$`GeneID (EntrezGene)`,
                        Ensembl = uniprotDB$Ensembl)

#-------------------------------------------------------------------------------
## Map missing Entrez IDs.
#-------------------------------------------------------------------------------

# Collect all proteins.
prots <- data.table(UniprotAC = hitpredict %>% select(Uniprot1, Uniprot2) %>% stack(),
                    UniprotID = hitpredict %>% select(Name1, Name2) %>% stack(),
                    Entrez = hitpredict %>% select(Entrez1, Entrez2) %>% stack(),
                    Ensembl = hitpredict %>% select(Ensembl1, Ensembl2) %>% stack(),
                    Taxonomy = hitpredict %>% select(Taxonomy))
prots <- prots %>% select(-c(colnames(prots)[grepl("ind",colnames(prots))])) %>%
  distinct()
colnames(prots) <- c("UniprotAC","UniprotID","Entrez","Ensembl","Taxonomy")

# Collect protins with missing entrez identifiers.
missing_prots <- prots %>% filter(is.na(Entrez))

# Map Uniprot IDs to Ensembl ID with UniprotDB map.
missing_prots$Uniprot_Ensembl <- as.character(
  uniprotDB$Ensembl[match(missing_prots$UniprotID,uniprotDB$UniprotID)])

# Define a function that utilizes the AnnotationDbi mapIds() function and
# organism specific databases (e.g. org.Mm.eg.db) to map Uniprot Ids to Entrez IDS.
getEntrez <- function(missing_prots,taxid) {
	# Get gene id mapping database for a given taxid.
	if (taxid == 9606) {
		require(org.Hs.eg.db)
		database <- org.Hs.eg.db
		message("Mapping human ensembl IDs to entrez!")
	} else if (taxid == 10090) {
		require(org.Mm.eg.db)
		database <- org.Mm.eg.db
		message("Mapping mouse ensembl IDs to entrez!")
	} else if (taxid == 10116) {
		require(org.Rn.eg.db)
		database <- org.Rn.eg.db
		message("Mapping rat ensembl IDs to entrez!")
	} else {
		stop("Please provide valid taxid c(9606, 10090, 10116).")
	}
	# Get all ensembl IDs for a given taxid.
	ensembl <- missing_prots$Uniprot_Ensembl[missing_prots$Taxonomy == taxid]
	# Map Uniprot ids to entrez.
	entrez <- mapIds(database, keys = ensembl, column = "ENTREZID",
			 keytype = "ENSEMBL", multiVals = "first")
	if (class(entrez) != "list") { entrez <- as.list(entrez)}
	idx <- seq(1:length(entrez))[unlist(lapply(entrez, is.null))]
	entrez[idx] <- NA
	entrez <- unlist(entrez)
	# Map missing Entrez.
	is_missing <- is.na(missing_prots$Entrez) & missing_prots$Taxonomy == taxid
	message(paste("Missing protein identifiers:", sum(is_missing)))
	missing_prots$Entrez[is_missing] <- entrez
	check <- is.na(missing_prots$Entrez) & missing_prots$Taxonomy == taxid
	message(paste("Mapped proteins:", sum(is_missing)-sum(check)))
	return(missing_prots)
}

# Do mapping for each species.
missing_prots <- getEntrez(missing_prots,taxids[1])
missing_prots <- getEntrez(missing_prots,taxids[2])
missing_prots <- getEntrez(missing_prots,taxids[3])

# Make a map
protmap <- as.list(missing_prots$Entrez)
names(protmap) <- missing_prots$UniprotAC

is_missing1 <- is.na(hitpredict$Entrez1)
hitpredict$Entrez1[is_missing1] <- unlist(protmap[hitpredict$Uniprot1[is_missing1]])

is_missing2 <- is.na(hitpredict$Entrez2)
hitpredict$Entrez2[is_missing2] <- unlist(protmap[hitpredict$Uniprot2[is_missing2]])

# Check the missingness again.
print(paste(
  "Number of unmapped genes remaining in the HitPredict database:",
  table(c(is.na(hitpredict$Entrez1), is.na(hitpredict$Entrez2)))[2]
))

#-------------------------------------------------------------------------------
## Map Entrez gene IDs to their homologous mouse gene.
#-------------------------------------------------------------------------------

# Download and load NCBI homology gene data.
url <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
destfile <- file.path(datadir, "homologene.data")
if (!file.exists(destfile)) {
	download.file(url, destfile)
} else {
	message("file already exists!")
}
homology_data <- read.delim(destfile, header = FALSE)

# Fix column names.
# Gene ID is a genes organism specific Entrez ID.
# HID is the genes homology id.
colnames(homology_data) <- c(
  "HID", "TaxonomyID", "GeneID",
  "GeneSymbol", "ProteinGI", "ProteinAccession"
)

# Create homology map for human, mouse, and rat.
homology_data <- subset(homology_data, homology_data$TaxonomyID %in% taxids)
homology_map <- as.list(homology_data$HID)
names(homology_map) <- homology_data$GeneID

# Map Entrez gene IDs to HID.
# If you have not subsetted the data, then this operation
# will be expensive and very slow!
hitpredict$HIDA <- homology_map[as.character(hitpredict$Entrez1)]
hitpredict$HIDB <- homology_map[as.character(hitpredict$Entrez2)]

#-------------------------------------------------------------------------------
## Map HIDs to mouse Entrez using data from MGI database.
#-------------------------------------------------------------------------------

# Download mouse homology data from MGI.
url <- "http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt"
myfile <- file.path(datadir, "HGNC_homologene.rpt")
if (!file.exists(myfile)) {
  download.file(url, myfile)
} else {
  message("file already exists!")
}
mus_homology_data <- read.delim(myfile, header = TRUE, sep = "\t", row.names = NULL)

# Fix column names.
col_names <- colnames(mus_homology_data)[-1]
mus_homology_data <- mus_homology_data[, -ncol(mus_homology_data)]
colnames(mus_homology_data) <- col_names

# Map HIDs to Mus Entrez.
idx <- match(hitpredict$HIDA, mus_homology_data$HomoloGene.ID)
hitpredict$musEntrezA <- mus_homology_data$EntrezGene.ID[idx]
idx <- match(hitpredict$HIDB, mus_homology_data$HomoloGene.ID)
hitpredict$musEntrezB <- mus_homology_data$EntrezGene.ID[idx]

# Keep genes that are mapped to homologous mouse genes.
keep <- !is.na(hitpredict$musEntrezA) & !is.na(hitpredict$musEntrezB)
hitpredict <- hitpredict[keep, ]
message(paste(dim(hitpredict)[1],"mouse, human, and rat PPIs were compiled from HitPredict."))

#-------------------------------------------------------------------------------
## Add gene symbols to hitpredict data.
#-------------------------------------------------------------------------------

# Must be organism specific!
Symbol1 <- mapIds(org.Mm.eg.db, keys = as.character(hitpredict$musEntrezA), column = "SYMBOL",
                  keytype = "ENTREZID", multivals = "first")
Symbol2 <- mapIds(org.Mm.eg.db, keys = as.character(hitpredict$musEntrezB), column = "SYMBOL",
                  keytype = "ENTREZID", multivals = "first")
hitpredict <- tibble::add_column(hitpredict, Symbol1, .after = 6)
hitpredict <- tibble::add_column(hitpredict, Symbol2, .after = 7)

# Remove any remaining rows with missing entrez IDs.
check <- sum(is.na(hitpredict$musEntrezA) | is.na(hitpredict$musEntrezB)) == 0
check

# Write to file.
fwrite(hitpredict,"compiled_mus_hs_rat_2mus_hitpredict.csv")






