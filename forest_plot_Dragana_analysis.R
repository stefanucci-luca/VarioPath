library(tidyverse)
df1=readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/statistical analysis/Single_var_het_carriers_results_FBC.xlsx",sheet=1)

plots = c()
outliers_feat = c()
outliers_gene = c()
outliers_val = c()
outliers_p = c()
outliers_var = c()
outliers_het = c()
outliers_hom = c()

for (var in grep( glob2rx("*_effect"), colnames(df1), value = T ) ) {

  a =  qqnorm(as.numeric(unlist(df1[,var])))
  b = lm(a$y ~ a$x)
  c = car::outlierTest(b)

  plots = c(plots, p)
  
  outliers_gene = c(outliers_gene, unlist(df1[as.numeric(names(c$rstudent)),1] ))
  outliers_val = c(outliers_val, unlist(df1[as.numeric(names(c$rstudent)),grep(var, colnames(df1))] ))
  outliers_p = c(outliers_p, unlist(df1[as.numeric(names(c$rstudent)),grep(var, colnames(df1))+2] ))
  outliers_het = c(outliers_het, unlist(df1[as.numeric(names(c$rstudent)),6] ))
  outliers_hom = c(outliers_hom, unlist(df1[as.numeric(names(c$rstudent)),7] ))
  for (var2 in 1:length(c$rstudent)) {
    outliers_feat = c(outliers_feat, var)
    outliers_var = c(outliers_var, paste(df1[as.numeric(names(c$rstudent)[var2]),2],
                                         df1[as.numeric(names(c$rstudent)[var2]),3],
                                         df1[as.numeric(names(c$rstudent)[var2]),4],
                                         df1[as.numeric(names(c$rstudent)[var2]),5],
                                         sep  = "_"
      )
    )
  }

  }


df2 = data.frame(
  "feat" = outliers_feat, 
  "gene"= outliers_gene, 
  'eff' = outliers_val, 
  'pval' = outliers_p,
  'var' = outliers_var,  
  'n_het' = outliers_het,
  'n_hom' = outliers_hom 
)

df2 = df2[sort(as.integer(df2$pval), decreasing = F),]


df1$ID = paste(df1$CHROM, df1$POS, df1$REF, df1$ALT, sep="_")

for (var in grep( glob2rx("*_effect"), colnames(df1), value = T ) ) {

df_tmp =   df1[
  which(as.numeric(unlist(df1[,grep(var, colnames(df_tmp)) + 2])) < 0.01),
  ]
  
if (dim(df_tmp)[1] > 0 ){
p1 = ggplot(data=df_tmp,
            aes(x = ID,
                y = as.numeric(unlist(df_tmp[,grep(var, colnames(df_tmp))])), 
                ymin =  as.numeric(unlist(df_tmp[,grep(var, colnames(df_tmp))])) - 
                  as.numeric(
                  unlist(df_tmp[,grep(var, colnames(df_tmp)) + 1])
                                                     ), 
                ymax = as.numeric(unlist(df_tmp[,grep(var, colnames(df_tmp))])) + 
                  as.numeric(
                  unlist(df_tmp[,grep(var, colnames(df_tmp)) + 1])
                                                     )
            )) +
  geom_pointrange(color = ggsci::pal_startrek("uniform")(7)[1]) +
  geom_hline(yintercept =0, linetype=2)+
  xlab('variant')+ ylab("beta")+
  coord_flip() +
  facet_wrap(~ GENE, nrow=50, scales = "free_y", strip.position = "left" )+
  theme_minimal()


p1 +
  theme(plot.title=element_text(size=12,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        panel.background = element_rect("transparent"), 
        strip.background = element_rect(fill = "white") )

ggsave(paste( "/Volumes/GoogleDrive/My Drive/Phd/VarioPath/statistical analysis/forest_plot_dragana_analysisi_20210606/",
              var, "_forest_plot.png"
        ),
       plot = p1, device = "png", units = "cm", limitsize = FALSE, dpi = "retina"
        )

}
}








