
library("tidyverse")
library("data.table")
df = data.table::fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/UKBb/variopath/VCF_200K_filtered_for_variopath/no_AF_filter/All_variopath_variants_in_unrelated_european_noAFfilter.tsv", nThread = 4)                       
setkey(df, VAR_ID)

length(unique(df$VAR_ID))

p1 = ggplot()+
  geom_histogram( data = as.data.table(table(df$PARTICIPANT_ID)), aes( x = N ), fill = "green", alpha = 0.3 ) +
  geom_histogram( data = as.data.table(table(df[UKB_AF <0.01, PARTICIPANT_ID])), aes( x = N ), fill = "blue", alpha = 0.3  ) + 
  geom_histogram( data = as.data.table(table(df[UKB_AF <0.001, PARTICIPANT_ID])), aes( x = N ), fill = "red", alpha = 0.3  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(1,200))
p1
ggsave("~/Desktop/Variant in UKB.png", device = "png")


pheno = readxl::read_xlsx("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/definitions/dictionary_11_March_2021.xlsx", sheet = 1)                       



dir_es = "/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/definitions/secondary_care"

list_files_scores = list.files(path = dir_es, pattern = "*.csv$", full.names = TRUE, recursive = TRUE)

# Read and create reference table with all effect sizes
data_all =  
  lapply(list_files_scores, data.table::fread, header = TRUE) %>%
  bind_rows(.id = "file")

data_all = na.omit(data_all[,c(1:5)])
data_long = data_all %>%   
  pivot_wider( id_cols = Disease, values_from = ICD10code, names_from = ICD10code, names_glue = "{ICD10code}_{.value}" )


data_long$ICD_codes = apply(data_long[2:dim(data_long)[2]], 
                            1,
                            function(x)
                            paste(x[x != "NULL"], collapse = "; ")
)


