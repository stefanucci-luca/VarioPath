
library(dplyr)
library(ggplot2)

MDTs = c('bleeding_and_coagulation',
         'hereditary_spherocytosis',
         'platelet',
         'thrombosis'
         )


#________________PATHOGENICITY CONVERSION TABLE_________________________________________
# Get additional variant filtering information
# 
googledrive::drive_download("https://docs.google.com/spreadsheets/d/17fAnlS8nNuHAxb7ZIKV4tqif7JJ-aOCVtPuTv0AQfyA/edit?usp=sharing")
pat_df = readxl::read_xlsx('Pathogenicity_JC_KM_20210126.xlsx',
                           trim_ws = T)

file.remove('Pathogenicity_JC_KM_20210126.xlsx')

pat_keep_df <- pat_df[,c("PATHOGENICITY","stay/go", "group")] %>% 
  filter(`stay/go` == 'stay')

pat_keep_df_reshaped = aggregate(PATHOGENICITY ~ group, data = pat_keep_df, c)

assign.new.pathogenicity = function(list_to_assign, conversion_df){
                            lapply(list_to_assign, 
                                   function(x)
                                    for (conv in 1:dim(conversion_df)[1] ) {
                                      if (x %in% unlist(conversion_df[conv,'PATHOGENICITY']) ) {
                                        return(conversion_df[conv,'group'])
                                      }
                                    }
                                   )
                                }

setwd('/Users/luca/Desktop/VarioPath/MDT_variants')
dir_MDT = dir('/Users/luca/Desktop/VarioPath/MDT_variants')
latest_dir = dir_MDT[startsWith(dir_MDT, 'V20')][length(dir_MDT[startsWith(dir_MDT, 'V20')])]

for (MDT in MDTs) {
  
  outdir='/Users/luca/Desktop/VarioPath/MDT_variants/plots'
  file_name_out_1 = paste("plot_genet_burden_", MDT, ".svg", sep = '')
  file_name_out_2 = paste("plot_variant_burden_", MDT, ".svg", sep = '')
  file_in = file.path( latest_dir, paste( 'MDT_for_', MDT, '.xls',
                                                  sep=''))
  
  variants = readxl::read_xls(file_in, 
                   col_names = T) %>% 
                   dplyr::select("GENE", "CHROM", "POS", "REF", "ALT", "PARTECIPANTS", 
                          "AF_ukb_calc", 'het', 'hom', 'hem', 'PATHOGENICITY','MOI_original_column')
  
  variants$counted.partecipants = rowSums(variants[,c('het', 'hom', 'hem')]) 
  variants$variatn_id = paste(variants$CHROM, variants$POS, variants$REF, variants$ALT, sep = "_")
  
  variants$new_pathogenicity = unlist(assign.new.pathogenicity(variants$PATHOGENICITY, pat_keep_df_reshaped))
  
  ggplot(variants, aes(forcats::fct_infreq(GENE), ..count.., fill = MOI_original_column)) +
    geom_bar() + 
    coord_flip() +
    guides(fill=guide_legend(ncol=1)) +
    theme(axis.text.x = element_text(size = 24, hjust = 1),
          strip.text.x = element_text(size = 24, colour = "black"),
          axis.text.y = element_text(size = 24, angle = 0), 
          legend.position = 'top',
    )
  
  ggsave(filename = file_name_out_1, 
         device = 'svg', 
         path = outdir,
         dpi = 'print', 
         height = 20, 
         width = 20,
         units = 'cm')
  
  ggplot(variants, aes(x = reorder(variatn_id, -counted.partecipants), 
                       y = counted.partecipants, 
                       fill = as.factor(new_pathogenicity) )) +
                   geom_col() +
                   facet_wrap('GENE', scales = 'free_x', ) +
                   scale_fill_manual(name = "Pathogenicity",  
                                       values = c("#1779ba",
                                                  "#767676",
                                                  "#3adb76",
                                                  "#ffae00",
                                                  "#cc4b37",
                                                  "#000000"),
                                       labels = c("2+ db agree P/LP; no conflict",
                                                  "1 db P/LP; no conflict",
                                                  "At least 1 db P/LP but conflict with VUS/DM?",
                                                  "At least 1 db P/LP but conflict with benign/likely benign/risk factor",
                                                  "No P/LP; DM? and/or VUS",
                                                  "No P/LP; DM? and likely benign/benign or likely benign/benign plus TG/EAHAD curated") 
                                       ) +
                   theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                         strip.text.x = element_text(size = 24, colour = "black"),
                         axis.text.y = element_text(size = 24),
                         legend.text = element_text(size = 20, colour = "black"),
                         legend.title = element_text(size = 24, colour = "black"),
                         legend.position = 'top',
                        )
  
  ggsave(filename = file_name_out_2, 
         device = 'svg', 
         path = outdir,
         dpi = 'print', 
         height = 60, 
         width = 90,
         units = 'cm')
  
}

df_path_gen = read.csv('/Users/luca/Desktop/VarioPath/MDT_master_spreadsheet_to_clean.tab',sep = '\t')

# Add MOI info
df_path_gen = df_path_gen %>% 
  filter(PATHOGENICITY %in% pat_keep_df$PATHOGENICITY)
#
df_path_gen$new_pathogenicity = unlist(assign.new.pathogenicity(df_path_gen$PATHOGENICITY, pat_keep_df_reshaped))


ggplot(df_path_gen, aes(forcats::fct_inseq(as.factor(new_pathogenicity)), ..count.., fill = as.factor(new_pathogenicity) )) +
  geom_bar() + 
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(name = "Pathogenicity",  
                    values = c("#1779ba",
                               "#767676",
                               "#3adb76",
                               "#ffae00",
                               "#cc4b37",
                               "#000000"),
                    labels = c("2+ db agree P/LP; no conflict",
                               "1 db P/LP; no conflict",
                               "At least 1 db P/LP but conflict with VUS/DM?",
                               "At least 1 db P/LP but conflict with benign/likely benign/risk factor",
                               "No P/LP; DM? and/or VUS",
                               "No P/LP; DM? and likely benign/benign or likely benign/benign plus TG/EAHAD curated") 
  ) +
  theme(axis.text.x = element_text(size = 24, hjust = 1),
        strip.text.x = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 24),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 24, colour = "black"),
        legend.position = 'right',
  )

ggsave(filename = 'Pathogenicity_distribution.svg', 
       device = 'svg', 
       path = outdir,
       dpi = 'print', 
       height = 20, 
       width = 60,
       units = 'cm')









variants1 = readxl::read_xls("V20210211/MDT_for_thrombosis.xls", 
                            col_names = T) %>% 
  dplyr::select("GENE", "CHROM", "POS", "REF", "ALT", "PARTECIPANTS", 
                "AF_ukb_calc", 'het', 'hom', 'hem', 'PATHOGENICITY','MOI_original_column')

variants1$counted.partecipants = rowSums(variants1[,c('het', 'hom', 'hem')]) 
variants1$variatn_id = paste(variants1$CHROM, variants1$POS, variants1$REF, variants1$ALT, sep = "_")
variants1$new_pathogenicity = unlist(assign.new.pathogenicity(variants1$PATHOGENICITY, pat_keep_df_reshaped))
variants1$mdt = "thrombosis"

variants2 = readxl::read_xls("V20210211/MDT_for_bleeding_and_coagulation.xls", 
                             col_names = T) %>% 
  dplyr::select("GENE", "CHROM", "POS", "REF", "ALT", "PARTECIPANTS", 
                "AF_ukb_calc", 'het', 'hom', 'hem', 'PATHOGENICITY','MOI_original_column')

variants2$counted.partecipants = rowSums(variants2[,c('het', 'hom', 'hem')]) 
variants2$variatn_id = paste(variants2$CHROM, variants2$POS, variants2$REF, variants2$ALT, sep = "_")
variants2$new_pathogenicity = unlist(assign.new.pathogenicity(variants2$PATHOGENICITY, pat_keep_df_reshaped))
variants2$mdt = "bleeding and coagulation"

variants3 = readxl::read_xls("V20210211/MDT_for_platelet.xls", 
                             col_names = T) %>% 
  dplyr::select("GENE", "CHROM", "POS", "REF", "ALT", "PARTECIPANTS", 
                "AF_ukb_calc", 'het', 'hom', 'hem', 'PATHOGENICITY','MOI_original_column')

variants3$counted.partecipants = rowSums(variants3[,c('het', 'hom', 'hem')]) 
variants3$variatn_id = paste(variants3$CHROM, variants3$POS, variants3$REF, variants3$ALT, sep = "_")
variants3$new_pathogenicity = unlist(assign.new.pathogenicity(variants3$PATHOGENICITY, pat_keep_df_reshaped))
variants3$mdt = "platelet"

variants = rbind(variants1, variants2, variants3)

variants = variants[which(variants$GENE %in% c("SERPINC1", "PROC", "PROS1",
                                               "MPL","GP1BB", "TUBB1",
                                               "F8","F9","VWF") ),]

ggplot(variants, aes(forcats::fct_infreq(GENE), ..count.., fill = MOI_original_column)) +
  geom_bar() + 
  coord_flip() +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(size = 24, hjust = 1),
        strip.text.x = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 24, angle = 0), 
        legend.position = 'top',
  )

ggsave(filename = file_name_out_1, 
       device = 'svg', 
       path = outdir,
       dpi = 'print', 
       height = 20, 
       width = 20,
       units = 'cm')

ggplot(variants, aes(x = reorder(variatn_id, -counted.partecipants), 
                     y = counted.partecipants, 
                     fill = as.factor(new_pathogenicity) )) +
  geom_col() +
  facet_wrap('GENE', scales = 'free', ) +
  scale_fill_manual(name = "Pathogenicity",  
                    values = c("#1779ba",
                               "#767676",
                               "#3adb76",
                               "#ffae00",
                               "#cc4b37",
                               "#000000"),
                    labels = c("2+ db agree P/LP; no conflict",
                               "1 db P/LP; no conflict",
                               "At least 1 db P/LP but conflict with VUS/DM?",
                               "At least 1 db P/LP but conflict with benign/likely benign/risk factor",
                               "No P/LP; DM? and/or VUS",
                               "No P/LP; DM? and likely benign/benign or likely benign/benign plus TG/EAHAD curated") 
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 0, angle = 90, hjust = 1),
        strip.text.x = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 24),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 24, colour = "black"),
        legend.position = 'top',
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )

ggsave(filename = file_name_out_2, 
       device = 'svg', 
       path = outdir,
       dpi = 'print', 
       height = 60, 
       width = 90,
       units = 'cm')







  