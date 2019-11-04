# Parsing epilepsy genes.

# Load the genes.
library(readxl)
library(data.table)
library(dplyr)
library(getPPIs)

# Load epilepsy genes.
myfile <- list.files(pattern="Wang")
rawdat <- as.data.table(read_excel(myfile))
data <- reshape2::melt(rawdat, measure.vars=colnames(rawdat))
data <- setnames(data, old=colnames(data), new=c("gene_category","symbol"))

# Split into list of genes.
gene_list <- data %>% group_by(gene_category) %>%  group_split()
names(gene_list) <- colnames(rawdat)
gene_list <- lapply(gene_list, function(x) na.omit(x) %>% select(symbol))

# Get entrez ids with mapIDs from getPPIs.
entrez <- lapply(gene_list,function(x) 
		 mapIDs(x$symbol,from="symbol",to="entrez",species="human"))

# Loop to store entrez ids in list.
for (i in 1:length(entrez)){
	df <- gene_list[[i]]
	e <- entrez[[i]]
	names(e) <- df$symbol
	gene_list[[i]] <- e[!is.na(e)]
}	

# Map gene symbols to entrez.
epilepsy_genes <- gene_list
names(epilepsy_genes) <- c("disease_genes","neurdevelopment_associated_genes",
			   "epilepsy_related_genes","epilepsy_associated_genes",
			   "all_epilepsy_genes")

myfile <- "epilepsy_genes.rda"
save(epilepsy_genes,file=myfile)


