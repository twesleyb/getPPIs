#---------------------------------------------------------------------
#' ## Trying to answer H_A's StackOverflow question. 

getwd()
files <- list.files(pattern=".txt")

data.list <- lapply(files, function(fil) {
  scan(file=fil, what=character())
})

library(dplyr)
names(data.list) <- basename(files) %>% stringr::str_remove("\\.txt$")

str(data.list)
# List of 8
# $ GSE108363_BCGdown_D:chr [1:350] "IL1B" "IL6" "IL1A" "CCL20" ...
# $ GSE108363_BCGdown_V: chr [1:267] "IL6" "CCL20" "IL1A" "CXCL5" ...
# $ GSE108363_BCGup_D  : chr [1:250] "FABP4" "CMTM2" "FUCA1" "CD36" ...
# $ GSE108363_BCGup_V  : chr [1:429] "FCN1" "FCGR3B" "MNDA" "CPVL" ...
# $ GSE108363_MTBdown_D: chr [1:86] "CCL20" "IL1B" "IL1A" "IL6" ...
# $ GSE108363_MTBdown_V: chr [1:244] "IL1B" "IL1A" "CCL20" "IL6" ...
# $ GSE108363_MTBup_D  : chr [1:128] "FUCA1" "FGL2" "TGFBI" "CPVL" ...
# $ GSE108363_MTBup_V  : chr [1:286] "FABP4" "RNASE1" "MNDA" "CPVL" ...

intersect(data.list$GSE108363_BCGdown_D, data.list$GSE108363_BCGdown_V) %>% length
intersect(data.list$List1, data.list$List2) %>% length

sapply(data.list, length)



#-------------------------------------------------------------------------------
# Using the intersect function to see the overlaps 
#-------------------------------------------------------------------------------

files_list <- list.files(pattern=".txt")
data.file1 <- files_list[1]
data.file2 <- files_list[2]
#data.file3 <- files_list[3]

#data.file1 <- "GSE108363_BCGdown_D.txt"
#data.file2 <- "GSE108363_BCGdown_V.txt"
#data.file3 <- "GSE108363_BCGup_D.txt"
#data.file4 <- "GSE108363_BCGup_V.txt"
#data.file5 <- "GSE108363_MTBdown_D.txt"
#data.file6 <- "GSE108363_MTBdown_V.txt"
#data.file7 <- "GSE108363_MTBup_D.txt"
#data.file8 <- "GSE108363_MTBup_V.txt"

genevect1 <- scan(data.file1, what=character(), sep="\n")
genevect2 <- scan(data.file2, what=character(), sep="\n")
#genevect3 <- scan(data.file3, what=character(), sep="\n")
#genevect4 <- scan(data.file4, what=character(), sep="\n")
#genevect5 <- scan(data.file5, what=character(), sep="\n")
#genevect6 <- scan(data.file6, what=character(), sep="\n")
#genevect7 <- scan(data.file7, what=character(), sep="\n")
#genevect8 <- scan(data.file8, what=character(), sep="\n")

#filelist <- list(data.file1,data.file2,data.file3)
filelist <- list(data.file1,data.file2)

#filelist <- list(data.file1, data.file2, data.file3, data.file4, data.file5, data.file6, data.file7, data.file8)
all(sapply(filelist, file.exists))
#data.list[[genevect1]]

#-------------------------------------------------------------------------------
# read files:
#-------------------------------------------------------------------------------

gene.lists <- lapply(filelist, function(f) {
  scan(file=f, what=character())
})

#-------------------------------------------------------------------------------
# Overlaps
#-------------------------------------------------------------------------------

# Initiate an empty list and matrix for storing output of loop.
genes.overlap <- list()
nfiles <- length(gene.lists)
mx.overlap.count <- matrix(NA,nrow=nfiles)

# Generate contrasts:
contrasts <- combn(nfiles,2)

# Loop to determine intersection:
for (i in 1:dim(contrasts)[2]){
  list1 <- contrasts[1,i]
  list2 <- contrasts[2,i]
  g1 <- gene.lists[[list1]]
  g2 <- gene.lists[[list2]]
  genes.overlap[[i]] <- intersect(g1, g2)
  b <- length(genes.overlap[[i]])
  mx.overlap.count[i] <- b
}

head(genes.overlap)
