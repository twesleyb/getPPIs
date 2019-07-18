#' ---
#' title: Compiling protein protein interactions.
#' author: Tyler W Bradshaw
#' ---

#-------------------------------------------------------------------------------
#' ## Setup the workspace.
#-------------------------------------------------------------------------------

# Clean-up.
rm(list = ls())
dev.off()
cat("\f")
options(stringsAsFactors = FALSE)

# Set the working directory
dir <- "D:/Projects/getPPIs"
setwd(dir)

#-------------------------------------------------------------------------------
#' ## Download HitPredict interactions.
#-------------------------------------------------------------------------------

# Download HitPredict data. This will take several minutes.
url <- "http://hintdb.hgc.jp/htp/download/HitPredict_interactions.txt.tgz"
gzfile <- paste(dir, "HitPredict_interactions.txt.tgz", sep = "/")
if (!file.exists(gzfile)) {
  download.file(url, gzfile)
} else {
  print("file exists!")
}

# Unzip and read data.
untar(gzfile)
unlink(gzfile)
file <- paste(dir, "HitPredict_interactions.txt", sep = "/")
data <- read.delim(file, header = TRUE, skip = 5)
unlink(file)

# We need to insure that genes are mapped to their stable, unique
# Entrez identifier. First, replace blanks with NA.
data$Entrez1[data$Entrez1 == ""] <- NA
data$Entrez2[data$Entrez2 == ""] <- NA

# Subset human, mouse, and rat data.
taxids <- c(9606, 10090, 10116)
data <- subset(data, data$Taxonomy %in% taxids)

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# Utilize the AnnotationDbi mapIds() function and organism specific
# datases (e.g. org.Mm.eg.db) to map uniprot Ids to Entrez IDS.
uniprot <- unique(c(data$Uniprot1, data$Uniprot2))

# Mapping mouse uniprot...
entrez <- mapIds(org.Mm.eg.db,
  keys = uniprot, column = "ENTREZID",
  keytype = "UNIPROT", multiVals = "first"
)
map <- as.list(entrez)
names(map) <- uniprot

# Map missing Entrez1
is_missing <- is.na(data$Entrez1)
data$Entrez1[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Map missing Entrez2
is_missing <- is.na(data$Entrez2)
data$Entrez2[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# Mapping human uniprot ids...
library(org.Hs.eg.db)

entrez <- mapIds(org.Hs.eg.db,
  keys = uniprot, column = "ENTREZID",
  keytype = "UNIPROT", multiVals = "first"
)
map <- as.list(entrez)
names(map) <- uniprot

# Map missing Entrez1
is_missing <- is.na(data$Entrez1)
data$Entrez1[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Map missing Entrez2
is_missing <- is.na(data$Entrez2)
data$Entrez2[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Number of unmapped genes:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# Mapping rat uniprot ids...
library(org.Rn.eg.db)

entrez <- mapIds(org.Rn.eg.db,
  keys = uniprot, column = "ENTREZID",
  keytype = "UNIPROT", multiVals = "first"
)
map <- as.list(entrez)
names(map) <- uniprot

# Map missing Entrez1
is_missing <- is.na(data$Entrez1)
data$Entrez1[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Map missing Entrez2
is_missing <- is.na(data$Entrez2)
data$Entrez2[is_missing] <- map[as.character(data$Uniprot1[is_missing])]

# Number of unmapped genes remaining:
print(paste(
  "Number of unmapped genes:",
  table(c(is.na(data$Entrez1), is.na(data$Entrez2)))[2]
))

# We have ~halved the number of unmapped genes.

# Subset the interaction data, keeping just those with Uniprot ID mapped to
# human, mouse, or rat Entrez IDs.
out <- is.na(data$Entrez1) | is.na(data$Entrez2)
data <- data[!out, ]
dim(data)

#-------------------------------------------------------------------------------
#' ## Map Entrez gene IDs to their homologous mouse gene.
#-------------------------------------------------------------------------------

# Download and load NCBI homology gene data.
url <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
destfile <- paste(dir, "homologene.data", sep = "/")
download.file(url, destfile)
homology_data <- read.delim(destfile, header = FALSE)
unlink(destfile)

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
data$HIDA <- homology_map[as.character(data$Entrez1)]
data$HIDB <- homology_map[as.character(data$Entrez2)]

#-------------------------------------------------------------------------------
#' ## Map HIDs to mouse Entrez using MGI.
#-------------------------------------------------------------------------------

# Download mouse homology data from MGI.
url <- "http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt"
file <- paste(dir, "HGNC_homologene.rpt", sep = "/")
if (!file.exists(file)) {
  download.file(url, file)
} else {
  print("file exists!")
}
mus_homology_data <- read.delim(file, header = TRUE, sep = "\t", row.names = NULL)
unlink(file)

# Fix column names.
col_names <- colnames(mus_homology_data)[-1]
mus_homology_data <- mus_homology_data[, -ncol(mus_homology_data)]
colnames(mus_homology_data) <- col_names

# Map HIDs to Mus Entrez.
idx <- match(data$HIDA, mus_homology_data$HomoloGene.ID)
data$musEntrezA <- mus_homology_data$EntrezGene.ID[idx]

idx <- match(data$HIDB, mus_homology_data$HomoloGene.ID)
data$musEntrezB <- mus_homology_data$EntrezGene.ID[idx]

# Keep genes that are mapped to homologous mouse genes.
keep <- !is.na(data$musEntrezA) & !is.na(data$musEntrezB)
data <- data[keep, ]
dim(data)

#-------------------------------------------------------------------------------
# Build PPI network for proteins identified by TMT MS.
#-------------------------------------------------------------------------------

# Load proteins identified by TMT.
Rdatadir <- "D:/Documents/R/Synaptopathy-Proteomics/RData"
file <- paste(Rdatadir, "Network_and_metaModules.Rds", sep = "/")
network_data <- readRDS(file)
net <- network_data$net
meta <- network_data$meta

# Map 3 un-mapped entrez IDs by hand.
meta[is.na(meta$entrez), ]
not_mapped <- list(
  "Ndufb1|P0DN34" = 102631912,
  "F8a1|Q00558" = 14070,
  "Pc|Q05920" = 18563
)
for (i in 1:length(not_mapped)) {
  idx <- meta$protein == names(not_mapped)[i]
  meta$entrez[idx] <- not_mapped[[i]]
}
# Check.
sum(is.na(meta$entrez))

# Add gene names.
meta$gene <- mapIds(org.Mm.eg.db,
  keys = meta$entrez, column = "SYMBOL",
  keytype = "ENTREZID", multiVals = "first"
)

# Add uniprot.
meta$uniprot <- sapply(strsplit(meta$protein, "\\|"), "[", 2)

# Get PPIs
idx <- data$musEntrezA %in% meta$entrez & data$musEntrezB %in% meta$entrez
sif <- data[idx, c(18, 19, 14)]
head(sif)

# Add gene symbol annotation.
entrez <- unique(c(sif$musEntrezA, sif$musEntrezB))
symbol <- mapIds(org.Mm.eg.db,
  keys = entrez, column = "SYMBOL",
  keytype = "ENTREZID", multiVals = "first"
)
# Create map and add gene symbols.
map <- as.list(symbol)
names(map) <- entrez
sif$geneA <- map[sif$musEntrezA]
sif$geneB <- map[sif$musEntrezB]

# Remove any missing values. There should be none.
out <- is.na(sif$musEntrezA) | is.na(sif$musEntrezB)
table(out)
sif <- sif[!out, ]
dim(sif)

# Sif is a data frame of lists... fix this.
apply(sif, 2, function(x) class(x))
sif <- as.data.frame(apply(sif, 2, function(x) unlist(x)))
rownames(sif) <- NULL

# Write to file.
outdir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
file <- paste(outdir, "SIF.xlsx", sep = "/")
write.excel(sif, file)

# Create a simple node attributes file for mapping node name to gene symbol.
# foo <- data.frame("EntrezID" = sif$musEntrezA,
#                  "Symbol" = sif$geneA)

# man <- data.frame("EntrezID" = sif$musEntrezB,
#                  "Symbol" = sif$geneB)

# noa <- unique(rbind(foo,man))

# Add module id and meta module id
# map1 <- as.list(meta$module)
# names(map1) <- meta$entrez

# map2 <- as.list(meta$metaModule)
# names(map2) <- meta$entrez

# add module and meta module annotations.
# noa$Module <- map1[noa$EntrezID]
# noa$MetaModule <- map2[noa$EntrezID]
# noa$MetaModule <- paste0("MM",noa$MetaModule)

# Convert module color to hex
# noa$ModuleColor <- lapply(as.list(noa$Module),function(x) col2hex(x))

# Insure that noa is a df.
# noa <- as.data.frame(apply(noa,2,function(x) unlist(x)))

# Write to file.
# outdir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
# file <- paste(outdir,"NOA.xlsx",sep="/")
# write.excel(noa,file)

#-------------------------------------------------------------------------------
#' ## Evaluate topology of the network.
#-------------------------------------------------------------------------------

# Load the data.
dir <- "D:/Documents/R/Synaptopathy-Proteomics/Tables/Network"
file <- paste(dir, "SIF.xlsx", sep = "/")
sif <- read_excel(file, sheet = 1)

# Create a data frame with all node attributes.
nodes <- data.frame(
  Entrez = unlist(meta$entrez),
  Uniprot = unlist(meta$uniprot),
  Symbol = unlist(meta$gene)
)

# Make igraph object.
g <- graph_from_data_frame(d = sif, vertices = nodes, directed = FALSE)

# Coerce to simple graph--remove duplicate edges and self-loops.
g <- simplify(g)
is.simple(g)

# Number of nodes and edges.
length(V(g)) # All but three unmapped genes.
length(E(g))

# Number of connected components.
connected_components <- components(g)
connected_components$csize[1] # The largest connected component.

# Summarize basic properties of the graph.
mytable <- data.frame(
  Nodes = formatC(length(V(g)), format = "d", big.mark = ","),
  Edges = formatC(length(E(g)), format = "d", big.mark = ",")
)
table <- tableGrob(mytable, rows = NULL)
grid.arrange(table)

# Save
outputfigsdir <- "D:/Documents/R/Synaptopathy-Proteomics/Figures/Combined/Final_WGCNA_Analysis"
file <- paste0(outputfigsdir, "/", "PPI_Network_properties.tiff")
ggsave(file, table, width = 2, height = 2, units = "in", dpi = 300)

# Mean path length.
round(mean_distance(g, directed = FALSE), 2)

# Calculate node connectivity (degree).
connectivity <- degree(g, loops = FALSE)
head(connectivity)

# Evaluate scale free fit with WGCNA function scaleFreePlot()
plot <- ggplotScaleFreePlot(connectivity, nBreaks = 10)$ggplot
plot

# The synaptic proteome exhibits a scale free topology!

# Save as tiff.
file <- paste0(outputfigsdir, "/", "PPI_Network_ScaleFreeFit.tiff")
ggsave(file, plot, width = 3, height = 2.5, units = "in")


# Histogram of k.
plot <- qplot(connectivity,
              geom = "histogram",
              binwidth = 5,
              main = "Connectivity Histogram",
              xlab = "Connectivity (k)",
              ylab = "Frequency",
              fill = I("black"),
              col = I("black"),
              alpha = 0.2
) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "black", size = 14, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  )

plot

# Save fig.
file <- paste0(outputfigsdir, "/", outputMatName, "PPI_Histogram_Connectivity", ".tiff")
ggsave(file, plot, width = 3, height = 2.5, units = "in", dpi = 300)

