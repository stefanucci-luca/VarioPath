#########
# Add information on the probability of intolarence of Loss of function for the VarioPath genes

# Add libraries
require(dplyr)
require(readr)
require(stringr)
library(ggplot2)

# General settings
project_dir <- "~/Desktop/VarioPath/"
setwd(project_dir)

# Table of pli for all transcripts from gnomad V2 Nature paper:
# The mutational constraint spectrum quantified from variation in 141,456 humans
# Konrad J. Karczewski et al.,Nat. 2020

gnomad_pli <- read_tsv("gnomad_paper_PMID32461654_supplementary_dataset_11_full_constraint_metrics.tsv")
gnomad_pli_tx_num <- length(unique(gnomad_pli$transcript))
message("Gonmad df contains information on ", gnomad_pli_tx_num, " transcripts")

#from the file "Interpratation_descrition_supp_gnomad_paper" I extract the columns I am interested in.
columns_interest_pli <- c("gene", 
                          "transcript",
                          "pLI",
                          "pRec", # Probability that transcript falls into distribution of recessive genes (~46% o/e pLoF ratio; computed from gnomAD data
                          "pNull") # Probability that transcript falls into distribution of unconstrained genes (~100% o/e pLoF ratio; computed from gnomAD data

# Import VarioPath gene list
file_vp <- "disease_gene_list_20200915.csv"
vp_gl <- read_csv(file_vp)
vp_gl_tx_num <- length(unique(vp_gl$`chosen transcript`))
message("VarioPath df contains information on ", vp_gl_tx_num, " transcripts")

# retrieve information on the chosen transcripts for VarioPath genes.
ov_gnomad_vp <- which(gnomad_pli$transcript %in% vp_gl$`chosen transcript`)

vp_overlap_gnomad_tx <- length(unique(ov_gnomad_vp))
message(vp_overlap_gnomad_tx, " VarioPath transcripts are present in gnomad resource")

# Add pLI info to trasncript table
intersection_gnomad_variopath <- which(gnomad_pli$transcript %in% vp_gl$`chosen transcript`)
message("The VarioPath genes in gnomadV2 are: ",
        length(intersection_gnomad_variopath), 
        " out of ",
        vp_gl_tx_num, " genes")

vp_df_pli <- merge(vp_gl,
                   gnomad_pli[intersection_gnomad_variopath, columns_interest_pli],
                   by.x = "chosen transcript",
                   by.y = "transcript"
                   ) %>% 
              select("gene", "pLI", everything())

# Plot the pLI distribustion
pdf(paste0(str_replace(file_vp, ".csv", "" ) ,"_pLI_plots.pdf"), width=20/2.54, height=12/2.54)
ggplot() + 
  geom_density(data = vp_df_pli, aes( x = pLI, y = ..scaled.. , color = "red")) +
  geom_density(data = gnomad_pli, aes( x = pLI, y=..scaled.., color = "black")) +
  ylab("density")+
  xlab("probability of loss of function intolerance (pLI)") +
  scale_color_identity(name = "Transcript",
                       breaks = c("black", "red"),
                       labels = c("Gnomad V2", "VarioPath"),
                       guide = "legend") +
  theme_minimal() 
dev.off()

write_csv(vp_df_pli,
          path = paste0( str_replace(file_vp, ".csv", "" ), "_pLI_info.csv" ) 
)


