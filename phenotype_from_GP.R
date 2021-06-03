# load library
suppressPackageStartupMessages(library(data.table))
# load table with phenotypes
df_pheno = data.table::fread('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/phenotype/gp_clinical.txt', nThread = 8)
# index the table wit hall the pheontypes
data.table::setkey(df_pheno, eid)

# function to extract the people with any of the phenotype in the list. 
people.with.pheno = function(file_csv){
  # This list return a data.table with one participant per raw and 2 columns, one columns is the participants' ID the other is a TRUE/FALSE for the phenotype
  # phenotype name is coming from the file name
  pheno_name = stringr::str_remove(basename(file_csv), ".csv")
  df_tmp = suppressWarnings(read.csv(file_csv, header = T))
  df_long = df_pheno[, .(pheno2 = df_tmp$read_codes %in% read_2), by=list(eid)]
  df_long$pheno3 = df_pheno[, .(pheno3 = df_tmp$read_codes %in% read_3), by=list(eid)]$pheno3
  df_long$pheno23 = apply(df_long[,2:3], 
                          1, 
                          function(x) sum(x) )
  df_compact = df_long[, .(pheno = sum(pheno23) > 0 ), by=list(eid)]
  colnames(df_compact) = c("eid", pheno_name)
  return(df_compact)
}
# directory to the files with phenotypes.
# one file per phenotype, contains ICD codes and read codes
dr_in = '/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/definitions/primary_care'
# outuput directory.
dr_out = '/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/returned_datasets'
# create an empty dt to append all the phenotypes
df_tmp2 = data.table::data.table(eid = unique(df_pheno$eid) )
# loop acroos all the phenotype files in the directory
for (pheno_file in list.files(dr_in, include.dirs = F, full.names = T, pattern = "*.csv"))
{
  df_tmp = people.with.pheno(pheno_file)
  data.table::setkey(df_tmp, eid)
  df_tmp2 = df_tmp2[df_tmp, on="eid"]
}
# export the phenotype dt
gp_pheno <- df_tmp2 
save(gp_pheno, 
     file = file.path(dr_out, "ukb_event_clinics.RData")
)