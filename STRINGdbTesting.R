## getPPIs

getwd()

# Create and set the working directory. 
root_dir <- "D:/Documents/R/getPPIs"
dir.create(root_dir)
setwd(root_dir)

# Install STRINGdb
library(STRINGdb) # Current version is 11. 
string_db <- STRINGdb$new( version="10", species=9606,score_threshold=0, input_directory="" )

data(diff_exp_example1)
head(diff_exp_example1)

example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )

hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )

string_db$get_interactions( c(tp53, atm) )
