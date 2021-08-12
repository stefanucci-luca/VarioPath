
require("ggrepel")
df1 = readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/statistical analysis/Single_var_het_carriers_results_FBC_sex_sep.xlsx", sheet = 1)

df1$ID = paste(df1$CHROM, df1$POS, df1$REF, df1$ALT,
               sep = "_")

set.seed(1)
for (var in grep( glob2rx("*_effect_F"), colnames(df1), value = T ) ) {
  
  var2 = str_replace(var, pattern = "F", replacement = "M")
  
  p1 = ggplot(df1[which(
    as.numeric(unlist(df1[,grep(var,colnames(df1)) + 2])) <= 0.05 & 
    as.numeric(unlist(df1[,grep(var2,colnames(df1)) + 2])) <= 0.05),], 
       aes(
         x= as.numeric(
         unlist(df1[which(
                          as.numeric(unlist(df1[,grep(var,colnames(df1)) + 2])) <= 0.05 & 
                          as.numeric(unlist(df1[,grep(var2,colnames(df1)) + 2])) <= 0.05),
                    var])
         ), 
           y = as.numeric(
             unlist(df1[which(
               as.numeric(unlist(df1[,grep(var,colnames(df1)) + 2])) <= 0.05 & 
                 as.numeric(unlist(df1[,grep(var2,colnames(df1)) + 2])) <= 0.05),
               var2])
           ) 
         ) ) +
  geom_point() +
  geom_text_repel( data = df1[which(
      as.numeric(unlist(df1[,grep(var,colnames(df1)) + 2])) <= 0.05 & 
        as.numeric(unlist(df1[,grep(var2,colnames(df1)) + 2])) <= 0.05),], 
      aes(label = GENE), 
      size = 3.5 ) +
  geom_vline( xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0) + 
  theme_minimal() 
    ggsave(paste( "/Volumes/GoogleDrive/My Drive/Phd/VarioPath/statistical analysis/gender_plot_dragana_analysisi_20210606/",
                  var, "_forest_plot.png"
    ),
    plot = p1, device = "png", units = "cm", limitsize = FALSE, dpi = "retina", width = 10, height = 10
    )
  }

