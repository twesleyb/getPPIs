# Installation of organism specific mapping databases.

# Install mapping packages.
# Fixme: should check if package already exists!
installDBs <- function(){
	annotationDBs <- list(Anopheles_gambia = list(taxid=7165, database="org.Ag.eg.db", alias=NA),
			      Arabidopsis_thaliana = list(taxid=7302, database="org.At.tair.db", alias=NA),
			      Bos_taurus = list(taxid=9913, database="org.Bt.eg.db", alias="bovine"),
			    Caenorhabditis_elegans = list(taxid=6239, database="org.Ce.eg.db",alias="worm"),
			    Canis_familiaris = list(taxid=9615,database="org.Cf.eg.db",alias="dog"),
			    Drosophila_melanogaster = list(taxid=7227,database="org.Dm.eg.db",alias="fly"),
			    Danio_rerio = list(taxid=7955, database="org.Dr.eg.db",alias="zebrafish"),
			    Escherichia_coli = list(taxid=83333, database="org.EcK12.eg.db",alias=c("ecoli","k12")),
			    Escherichia_coliSakai = list(taxid=386585, database="org.EcSakai.eg.db",alias="sakai"),
			    Gallus_gallus = list(taxid=9031, database="org.Gg.eg.db",alias="chicken"),
			    Homo_sapiens = list(taxid=9606,database="org.Hs.eg.db",alias="human"),
			    Mus_musculus = list(taxid=10090,database="org.Mm.eg.db",alias="mouse"),
			    Macaca_mulatta = list(taxid=9544,database="org.Mmu.eg.db",alias="macaque"),
			    Plasmodium_falciparum = list(taxid=5833,database="org.Pf.plasmo.db",alias="malaria"),
			    Pan_troglodytes = list(taxid=9598,database="org.Pt.eg.db",alias="chimp"),
			    Rattus_norvegicus = list(taxid=10116,database="org.Rn.eg.db",alias="rat"),
			    Saccharomyces_cerevisiae = list(taxid=4932,database="org.Sc.sgd.db",alias="yeast"),
			    Sus_scrofa = list(taxid=9823,database="org.Ss.eg.db",alias="pig"),
			    Xenopus_laevis= list(taxid=8355,database="org.Xl.eg.db",alias="frog"))

	# Loop to install packages.
		packages <- unlist(sapply(annotationDBs,"[",2))
		for (package in packages){
		       if (package !%in% rownames(installed.packages())) {
			       BiocManager::install(package)
		       } else {
			       message(paste(package,"is already installed!"))
		       }
		}
}
