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

# Seperate rows with multiple homologs.
interactions <- interactions %>% 
	tidyr::separate_rows(`Human Homolog A`, sep=";") %>%
	tidyr::separate_rows(`Human Homolog B`, sep=";")

# Subset, keep interactions between human homologs.
hs_interactions <- interactions %>% 
	filter(!is.na(`Human Homolog A`) & !is.na(`Human Homolog B`))

# Save the data.
subdat <- hs_interactions[sample(nrow(hs_interactions),1000),]
fwrite(subdat,"data.csv")
