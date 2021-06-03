library(skimr)
library(ggplot2)
library(gghighlight)
library(tidyverse)
library(reshape2)
library(plyr)
library(ggExtra)

gwas_c = data.table::fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/hi_c/GWAS/EBI_GWAS_catalogue_ALL_ASSOCIATION_20210502.tsv")

OB = filter(gwas_c, is.infinite(`OR or BETA`) == FALSE)
OB = filter(OB, is.na(`OR or BETA`) == FALSE)
OB = filter(OB, grepl(pattern = "uropean", x = OB$`INITIAL SAMPLE SIZE`))
OB$cat = ifelse(OB$SNPS %in% kvar, "p_lp","GWAS")

OB$recalc = as.numeric(
  apply(OB,
      1,
      function(x)
        if (grepl(pattern = "unit", 
                 x = x["95% CI (TEXT)"])) {
          log(as.numeric(as.character(x["OR or BETA"])))
        } else {
          x["OR or BETA"]
        }
      )
  )


# import effect sizes
dir_es = "/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/UKBb/variopath/burden_test_raremetal/output_pre_MDT/phenotypes"

list_files_scores = list.files(path = dir_es, pattern = "*.dominant.singlevar.score.txt.gz$", full.names = TRUE, recursive = TRUE)

# Read and create reference table with all effect sizes
data_all =  
  lapply(list_files_scores, data.table::fread, skip = "#CHROM", header = TRUE) %>%
  bind_rows(.id = "file")

sign_pvalue = 0.05

df_eff = data.table::dcast(data_all, 
                           `#CHROM` + POS + REF + ALT ~ file, 
                           value.var = c('ALT_EFFSIZE', "PVALUE"), function(x){max(x,na.rm = T)},
                           subset = .(PVALUE < sign_pvalue & is.na(PVALUE) == FALSE))

# return only the max squared effect size per variant
df_eff = as.data.frame(df_eff)

df_eff[sapply(df_eff, is.infinite)] <- NA

df_eff$sqr_eff = apply(df_eff[,grep(colnames(df_eff) ,pattern = "ALT_EFFSIZE*")], 
                            1, 
                            function(x){
                              max(x, na.rm = T)
                            })


# import effect sizes
dir_es = "/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/UKBb/variopath/burden_test_raremetal/output_after_MDT_clean/phenotypes"

list_files_scores = list.files(path = dir_es, pattern = "*.dominant.singlevar.score.txt.gz$", full.names = TRUE, recursive = TRUE)

# Read and create reference table with all effect sizes
data_all =  
  lapply(list_files_scores, data.table::fread, skip = "#CHROM", header = TRUE) %>%
  bind_rows(.id = "file")

sign_pvalue = 0.05

df_eff_post = data.table::dcast(data_all, 
                           `#CHROM` + POS + REF + ALT ~ file, 
                           value.var = c('ALT_EFFSIZE', "PVALUE"), function(x){max(x,na.rm = T)},
                           subset = .(PVALUE < sign_pvalue & is.na(PVALUE) == FALSE))

# return only the max squared effect size per variant
df_eff_post = as.data.frame(df_eff_post)

df_eff_post[sapply(df_eff_post, is.infinite)] <- NA

df_eff_post$sqr_eff = apply(df_eff_post[,grep(colnames(df_eff_post) ,pattern = "ALT_EFFSIZE*")], 
                       1, 
                       function(x){
                         max(x, na.rm = T)
                       })


ggplot() +
  geom_density(data = OB, mapping = aes(recalc), adjust = 3, fill = 'red',alpha=0.5) +
  geom_density(data = df_eff, mapping = aes(sqr_eff), fill = 'blue', alpha=0.5) +
  geom_density(data = df_eff_post, mapping = aes(sqr_eff), fill = 'green', alpha=0.5) +
  xlim(-10,10) +
  theme_minimal()
  