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

# Summary of interaction methods.
df <- as.data.frame(table(interactions$"Interaction Method"))
colnames(df)[1] <- "Method"
df <- arrange(df,desc(Freq))
knitr::kable(df)

# Create data.frame mapping yeast genes to their human homolog(s).
gene_map <- data.table('Yeast Gene' = homologs$"Gene Standard Name",
		       'Human Entrez' = homologs$"Entrez",
		       'Source' = homologs$Source)
gene_map <- gene_map %>% filter(`Yeast Gene` != "") %>% as.data.table()
gene_map <- gene_map %>% group_by(`Yeast Gene`) %>%
	dplyr::summarize(Entrez = paste(unique(`Human Entrez`),collapse=";"),
			 Source = paste(unique(Source),collapse=";"))

# Map Yeast genes to their human homolog.
idxA <- match(interactions$"Gene A Standard Name", gene_map$"Yeast Gene")
idxB <- match(interactions$"Gene B Standard Name", gene_map$"Yeast Gene")
interactions$"Human Homolog A" <- gene_map$"Entrez"[idxA]
interactions$"Human Homolog B" <- gene_map$"Entrez"[idxB]

# Seperate rows (yeast interactions) with multiple human homologs.
interactions <- interactions %>% 
	tidyr::separate_rows(`Human Homolog A`, sep=";") %>%
	tidyr::separate_rows(`Human Homolog B`, sep=";")

# Subset, keep interactions between human homologs.
hs_interactions <- interactions %>% 
	filter(!is.na(`Human Homolog A`) & !is.na(`Human Homolog B`))

# Map human genes to mouse homologs.
geneA <- hs_interactions$"Human Homolog A"
geneB <- hs_interactions$"Human Homolog B"
hs_interactions$"Mouse Homolog A" <- getHomologs(geneA, species="mouse")
hs_interactions$"Mouse Homolog B" <- getHomologs(geneB, species="mouse")

# Subset, keep interactions between mouse homologs.
ms_interactions <- hs_interactions %>% 
	filter(!is.na(`Mouse Homolog A`) & !is.na(`Mouse Homolog B`))

# Make column headers more clear.
colnames(ms_interactions)[1] <- "Yeast Gene A"
colnames(ms_interactions)[2] <- "Yeast Gene B"

# Save the data. 
y2h2m_interactions <- ms_interactions
datadir <- file.path(root,"data")
save(y2h2m_interactions,file=file.path(datadir,"y2h2m_interactions.rda"),version=2)
