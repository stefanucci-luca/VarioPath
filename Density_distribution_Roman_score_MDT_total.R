
library(ggplot2)
data_back_ground = data.table::fread("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/protein_structure/Roman's_VarioPath_scores_20210312.tsv")

data_MDT = readxl::read_xlsx("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/protein_structure/Roman_scores_MDT_genes_20210521.xlsx",
                             sheet = 1)

plot(density(as.numeric(as.character(data_back_ground$V9)), na.rm = T))

ggplot() + 
  geom_density( as.data.frame(
    as.numeric(
      as.character(data_back_ground$V9)
      )
    ),  mapping = aes(x = as.numeric(as.character(data_back_ground$V9)))
) +
  geom_density(  data_MDT, mapping =  aes(x = as.numeric(Score)), color = "Purple" ) +
  theme_bw() +
  xlim(c(-10,10)) +
  theme( panel.background = element_rect(fill = "transparent",color = NA), # bg of the panel
         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
         panel.grid.major = element_blank(), # get rid of major grid
         panel.grid.minor = element_blank(), # get rid of minor grid
         legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
         legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
  )
