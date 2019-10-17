# Analysis of the mouse interactome.


# Load the data.
library(getPPIs)
data(musInteractome)

# Build a network.
mygenes <- unique(c(musInteractome$osEntrezA,musInteractome$osEntrezB))
nodes <- mygenes

  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db,
    keys = as.character(nodes),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  noa <- data.table::data.table(node = nodes, symbol = symbols)

  # Check.
  if (sum(is.na(symbols)) != 0) {
    message(paste("Warning: Unable to map",sum(is.na(symbols)), "Entrez IDs to gene symbols!"))
  }
  not_mapped <- nodes[is.na(symbols)]
# Remove from sif and noa.
sif <- filter(sif, sif$osEntrezA %in% not_mapped | sif$osEntrezB %in% not_mapped)
noa <- filter(noa, nodes %in% not_mapped)




g <- buildNetwork(musInteractome,mygenes)

# Node names.
nodes <- names(V(g))
names(nodes) <- vertex_attr(g,"symbol")

# Partition the network.
partition <- getProfile(g,partition_type="ModularityVertexPartition",nsteps=1)

#------------------------------------------------------------------------------

