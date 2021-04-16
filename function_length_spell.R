#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

colname_to_use_for_filtering = args[2]
df = data.table::fread(args[1], header = T)
df = as.data.frame(df)
duration = function(df){
  tmp_df = df[,c("eid",colname_to_use_for_filtering)]
  return(aggregate(df[,colname_to_use_for_filtering], by=list(eid=df$eid), function(x) {sum(x, na.rm =T)}))
}

df_out = duration(df)

write_delim( df_out, paste0("length_of_", colname_to_use_for_filtering, '_per_individual.txt' ), 
            delim = "\t", append = F, quote_escape = 'none', col_names = T)