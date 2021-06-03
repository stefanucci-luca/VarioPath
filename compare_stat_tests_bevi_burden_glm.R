library(skimr)
library(ggplot2)
library(gghighlight)
library(tidyverse)
library(reshape2)
library(plyr)
library("ggpubr")
library(ggExtra)
# import effect sizes
dir_es = "/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/UKBb/variopath/burden_test_raremetal/output_pre_MDT/phenotypes"

list_files_scores = list.files(path = dir_es, pattern = "*.dominant.singlevar.score.txt.gz$", full.names = TRUE, recursive = TRUE)

# Read and create reference table with all effect sizes
data_all =  
  lapply(list_files_scores, data.table::fread, skip = "#CHROM", header = TRUE) %>%
  bind_rows(.id = "file")

skim(data_all)

sign_pvalue = 0.05

df_eff = data.table::dcast(data_all, 
                           `#CHROM` + POS + REF + ALT ~ file, 
                           value.var = c('ALT_EFFSIZE', "PVALUE"), function(x){max(x,na.rm = T)},
                           subset = .(PVALUE < sign_pvalue & is.na(PVALUE) == FALSE))

# return only the max squared effect size per variant
df_eff = as.data.frame(df_eff)

df_eff[sapply(df_eff, is.infinite)] <- NA

df_eff$sqr_eff = sqrt(apply(df_eff[,grep(colnames(df_eff) ,pattern = "ALT_EFFSIZE*")], 
                            1, 
                            function(x){
                              max(x^2, na.rm = T)
                            }))

df_eff$eff_max_raw = apply(df_eff[,grep(colnames(df_eff) ,pattern = "ALT_EFFSIZE*")], 
                            1, 
                            function(x){
                              max(x, na.rm = T)
                            })

df_eff$pvalue = apply(df_eff[,grep(colnames(df_eff) ,pattern = "PVALUE*")], 
                           1, 
                           function(x){
                             min(x, na.rm = T)
                           })

df_eff$id = paste(paste0( 'chr', df_eff$`#CHROM`), 
                  df_eff$POS, df_eff$REF, df_eff$ALT, sep="_")



eff_luan = read.csv("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luanluan/output/csvfiles/ukb_genetic_associations_dvt_pe.csv") %>% 
  filter(p < sign_pvalue)



# BeviMed
# import effect sizes
dir_es = "~/Desktop/VarioPath/bevimed/bevimed_results_20200520"

list_files_scores = list.files(path = dir_es, pattern = "*.csv$", full.names = TRUE, recursive = TRUE)

# Read and create reference table with all effect sizes
data_all_bv = rbindlist(lapply(list_files_scores, fread))
data_all_bv$id = paste(paste0( 'chr', data_all_bv$CHR), 
                       data_all_bv$POS, data_all_bv$REF, data_all_bv$ALT, sep="_")

#################################################



df_compare = data.frame("variant" = df_eff$id,
                        "eff_size_burden" = df_eff$sqr_eff,
                        "eff_burden_raw" = df_eff$eff_max_raw,
                        "pval" = df_eff$pvalue)

df_compare = merge(eff_luan, df_compare, by.x = "expvar", by.y = "variant", all.y = T)
df_compare = merge(data_all_bv, df_compare, by.x = "id", by.y = "expvar", all.x = T)



cor(df_compare$beta, df_compare$eff_burden_raw)
plot(df_compare$BeviMed, df_compare$pval)
plot(df_compare$BeviMed, df_compare$p)

View(df_compare %>% filter(pval < 0.01 & BeviMed > 0.8 & p < 0.02))

