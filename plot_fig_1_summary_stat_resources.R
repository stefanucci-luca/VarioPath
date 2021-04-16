library(gsheet)
library(data.table)
library(ggplot2)
library(stringr)
library(scales)
library(reshape2)
library(gridExtra)
library(grid)
library(corrplot)
library(cowplot)
library(gridGraphics)
library(grid)
# Functions
# Given a vector remove all the insatnces of number in that vector
is.vec.num = 
  function(vec){
    list_res = c()
    suppressWarnings(
    for (i in 1:length(vec)) {
      if ( is.na(as.numeric(vec[i])) == TRUE ) {
        list_res = c(list_res, vec[i])
        list_collaps = paste(list_res, collapse = "|")
      }
    })
    return(list_collaps)
  }
# Take a vector as input and remove the pathogenicity taht are jsut numbers
clean.path.num =
  function(vec){
    apply( 
      as.data.frame(vec),
      1,
      function(x)
        is.vec.num(
          str_split(
            string = x ,
            pattern = "\\|", simplify = T )
        )
    )
  }
# remove duplicate from list "|" separated
rem.path.dup =
  function(vec){
    apply( 
      as.data.frame(vec),
      1,
      function(x)
        paste(
          unique(
          str_trim(
          str_split(
            string = x ,
            pattern = "\\|", simplify = T )
        )
        ), collapse = "|"
      )
    )
  }

# Import dataframes
# The VEP annotated one is useful to assign effect of a variant on the transcript
df_vep_var_variopath = fread('/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/VarioPath_VEP_annotation_only_variopath_transcripts_variants_20200615.tab')
# pathogenicity one is reportin DM/pathogenic and so on.
pathogenicity_df = fread('/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/VarioPath_variants_20200615.tsv')

# Vep annotation has some duplication (#18,345), possible sources of duplication are:
# • A variant maps 2 genes and both are in our list of selected genes.
# • Variants were called in a different way in the different databases and VEP normalise them.
# • There is a duplication that I carried from the original list.
# I'll remove the duplicated variants that are on the same gene, but keep the one that match to 2 different genes.
not_duplicated_id = which(duplicated(c(df_vep_var_variopath$Uploaded_variation, df_vep_var_variopath$Gene)) == FALSE)
df_vep_var_variopath_no_duplication = df_vep_var_variopath[c(not_duplicated_id),]

# Create a table for patogenicity per variant.
# merge df_vep_var_variopath and pathogenicity_df
## line to used #pathogenicity_df$variopath_ID = str_replace_all(pathogenicity_df$openCGA_ID, pattern = ":", replacement = "_")
# pathogenicity and df_vep_var_variopath dataframes are using two differen nomenclature for the variants. In particular:
# pathogenicity df use a structure such as : 10:100203:TA:A
# while df_vep_var_variopath would report the same variatn as 10:100203:T:-
# These lines of code convert pathogenicity structure to df_vep_var_variopath 
df = pathogenicity_df
df$mapply_dist = mapply(function(x,y) which.min(x==y),strsplit(df$REF,""),
                        strsplit(df$ALT,""))

df <- as.data.frame(dplyr::mutate_all(df, as.character))
df2 = as.data.frame(df)

for (variable in 1:dim(df)[1]) {
  if (df[variable,"mapply_dist"] == 2){
    if (str_length(df2[variable,"REF"]) > str_length(df2[variable,"ALT"]) ) {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"REF"] <- paste0(na.omit( strsplit(df[variable,"REF"],"")[[1]][df[variable,"mapply_dist"]:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"ALT"]) == 1) {
        df2[variable,"ALT"] <- "-"
      } else {
        df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][1:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      }
    } else {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][df[variable,"mapply_dist"]:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"REF"]) == 1) {
        df2[variable,"REF"] <- "-"
      } else {
        df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][1:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      }
    }
  } else if (df[variable,"mapply_dist"] > 2){
    if (str_length(df2[variable,"REF"]) > str_length(df2[variable,"ALT"]) ) {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][2:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"ALT"]) == 1) {
        df2[variable,"ALT"] <- "-"
      } else {
        df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][1:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      }
    } else {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][2:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"REF"]) == 1) {
        df2[variable,"REF"] <- "-"
      } else {
        df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][1:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      }
    }
  }
}

# Add new IDs to the pathogenicity df
pathogenicity_df$variopath_ID = paste(df2$`#CHROM`, df2$POS, df2$REF, df2$ALT, sep = "_")  
# create a new database with the relevant information from the 2 dataframes.
df_to_plot = merge(pathogenicity_df, df_vep_var_variopath_no_duplication, by.x = 'variopath_ID', by.y = 'Uploaded_variation')


ggplot(data = subset(as.data.frame(table(df_vep_var_variopath_no_duplication$IMPACT)) ), 
       mapping = aes(x = reorder(Var1, - Freq), Freq, fill= Var1)) +
  geom_col() +
  theme_minimal() +
  ggtitle("Num Impact") +
  ggsci::scale_fill_npg() + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "null", 
        axis.text = element_text(size = 12))
ggsave(filename = paste0('/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/PhD Thesis/results/figures/freq_impact_full_list_variants_', 
                         format(Sys.Date(), "%Y%m%d"),
                         ".svg"),
       device = 'svg'
       )

ggplot(data = subset(as.data.frame(table(pathogenicity_df$QUAL[pathogenicity_df$QUAL != ""])), subset = Freq > 1500 ), 
       mapping = aes(x = reorder(Var1, - Freq), Freq, fill= Var1)) +
  geom_col() +
  theme_minimal() +
  ggsci::scale_fill_npg() + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "null", 
        axis.text = element_text(size = 12))
ggsave(filename = paste0('/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/PhD Thesis/results/figures/freq_pathogenicity_full_list_variants_', 
                         format(Sys.Date(), "%Y%m%d"),
                         ".svg"),
       device = 'svg'
)


## Gene stats
# import data on the genes
df_gene = gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1mDGFXJ9zHWlSezQaKNfbyhdfQo7hlTTB50qXw2fSPN0/edit#gid=29422767")

str(df_gene)

# correct coouple of typos in the field separator
df_gene$`Chosen disease domain(s)` = str_replace_all(string = df_gene$`Chosen disease domain(s)`, 
                                                     pattern = '\\;|\\:',
                                                     replacement = "\\|")

df_gene$`Chosen disease domain(s)` = str_replace_all(string = df_gene$`Chosen disease domain(s)`, 
                                                     pattern = ' disorders',
                                                     replacement = "")
list_domain = as.vector(na.omit(unique(str_trim(unlist(str_split(string = df_gene$`Chosen disease domain(s)`,
                        pattern = "\\|"))))))

number_of_genes_per_domain = as.data.frame(table(as.vector(na.omit(str_trim(unlist(str_split(string = df_gene$`Chosen disease domain(s)`,
                                                                                                  pattern = "\\|")))))))

xlsx::write.xlsx(x = list_domain, file = '/Volumes/GoogleDrive/My Drive/Phd/VarioPath/genes/genes_and_domain_to_prep_interaction_analysis.xlsx', 
                 sheetName = 'list_of_domains', row.names = F, append = T)
xlsx::write.xlsx(x = as.data.frame(df_gene)[,c(2,14)], file = '/Volumes/GoogleDrive/My Drive/Phd/VarioPath/genes/genes_and_domain_to_prep_interaction_analysis.xlsx', 
                 sheetName = 'list_of_genes_and_domains', row.names = F, append = T)

all_domain_combinations = expand.grid(list_domain,list_domain)

matrix_dis_domain = as.data.frame(
                                str_split(string = df_gene$`Chosen disease domain(s)`,
                                pattern = "\\|", simplify = T)
                                )

all_domain_combinations$overlap = apply(all_domain_combinations,
                                                    1,
                                                    function(y)
                                                        sum(
                                                          apply( matrix_dis_domain, 
                                                                 1, 
                                                                 function(x)
                                                                   return( y[1] %in% x & y[2] %in% x ) 
                                                          )
                                                        )
                                                      )

matrix_heatmap = dcast(all_domain_combinations , Var1 ~ Var2 , value.var = "overlap")
rownames(matrix_heatmap) = matrix_heatmap[,1]
matrix_heatmap = as.matrix(matrix_heatmap[,-1])

matrix_heatmap_log = log(matrix_heatmap+1e-10)

plot1 = corrplot(matrix_heatmap_log, method = 'color', 
         type = 'upper', order = 'hclust', hclust.method = 'complete',
         is.corr = FALSE, tl.pos='l', tl.cex = 0.61, tl.col = "black",tl.srt = 45, cl.pos = "b", 
         diag = TRUE,
         col= colorRampPalette(c("blue","yellow","blue"))(50))

grid.echo()
P1 <- grid.grab()
grid.draw(P1)

order_columns = rownames(plot1)

number_of_genes_per_domain = number_of_genes_per_domain[order(match(number_of_genes_per_domain[,1],as.data.frame(rev(order_columns))[,1])),]

P2 = ggplot(data = number_of_genes_per_domain, 
       mapping = aes(x = Var1, Freq)) +
  geom_col() +
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(limits=number_of_genes_per_domain$Var1) +
  theme(axis.text = element_blank(), #element_text(angle = 0, size = 8, hjust = 0, vjust = 0), 
        legend.position = "null")

ggdraw() + cowplot::draw_plot(P1,width =0.4, height =0.75, 
                              hjust =  0,vjust =  0) + 
           cowplot::draw_plot(P2, width = 0.5, height = 0.495, 
                              hjust = -.68, vjust = -0.06) #

ggsave(filename = paste0('/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/PhD Thesis/results/figures/variopath_gene_domain_', 
                         format(Sys.Date(), "%Y%m%d"),
                         ".svg"),
       width = unit(13, 'cm'),
       height = unit(13, 'cm'),
       device = 'svg'
)


df_gene_MOI = xlsx::read.xlsx('/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/PhD Thesis/results/figures/MOI_list.xls', 
                              sheetName = 1) %>% na.omit()

ggplot(data = df_gene_MOI, 
       mapping = aes(x = reorder(Var1, - Freq), Freq, fill= Var1)) +
  geom_col() +
  theme_minimal() +
  ggsci::scale_fill_simpsons() + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "null", 
        axis.text = element_text(size = 12))
ggsave(filename = paste0('/Volumes/GoogleDrive/My Drive/Phd/Shared Luca Mattia/PhD Thesis/results/figures/variopath_gene_moi_', 
                         format(Sys.Date(), "%Y%m%d"),
                         ".svg"),
       device = 'svg'
)
