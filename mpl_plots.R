
library(data.table)
library(tidyverse)
library(skimr)

mpl_people = fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/rds-who1000-wgs10k/user/dv310/VarioPath/variopath_variants_not_AF_fitlering/All_variopath_variants_in_unrelated_european_noAFfilter_gene_name_PLT.tsv") %>% 
  filter(V11 == "MPL")
colnames(mpl_people) = c("VAR_ID", "CHROM", "POS", "REF", "ALT", "HOM", "HET", "AF", "PARTICIPANT_ID", "GT", "GENE")
setkey(mpl_people, PARTICIPANT_ID)
fbc = fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/rds-who1000-wgs10k/user/dv310/VarioPath/UKBB_FBC_tech_adj_only.tsv", nThread = 4)
setkey(fbc, subject_id)

dt_1 = mpl_people[fbc, on=c(PARTICIPANT_ID="subject_id")]

gender = fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/ukbcvo/sex_otUKBcnv.tab", nThread = 4)
gender = gender[,.(f.eid,  f.22001.0.0)]
setkey(gender, f.eid)

dt_1 = dt_1[gender, on=c(PARTICIPANT_ID="f.eid")]

skim(dt_1)

dt_1 = na.omit(dt_1, cols=c("f.22001.0.0"))

for (var in colnames(dt_1)[14:dim(dt_1)[2]] ) {
  g = ggplot() + 
    geom_density(data = dt_1[f.22001.0.0==1], aes(x = get(var), y = ..density..), fill="blue", alpha=0.2) + 
    annotate("segment",
             x = unlist(dt_1[VAR_ID == "chr1_43348956_G_A" & f.22001.0.0==1, ..var]), 
             xend = unlist(dt_1[VAR_ID == "chr1_43348956_G_A" & f.22001.0.0==1, ..var]),
             y = 0,
             yend = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8,
             colour = "blue", alpha = 0.5,
    ) +
    geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T), linetype="dashed", size=1, color = "blue")  + 
    geom_density(data = dt_1[f.22001.0.0==0], aes(x = get(var), y = -..density..),  fill= "red", alpha =0.2)  + 
    annotate("segment",
             x = unlist(dt_1[VAR_ID == "chr1_43348956_G_A" & f.22001.0.0==0, ..var]), 
             xend = unlist(dt_1[VAR_ID == "chr1_43348956_G_A" & f.22001.0.0==0, ..var]), 
             y = 0,
             yend = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8, # standard hight
             colour = "red", alpha = 0.5,
    ) +
    geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T), linetype="dashed", size=1, color = "red") +
    ggtitle(paste("MPL variants for", eval(var)))+
    coord_cartesian(xlim = c(quantile(dt_1[,..var], 0.005, na.rm = T, names = FALSE),
                             quantile(dt_1[,..var], 0.995, na.rm = T, names = FALSE))) +
    theme_minimal()
  ggsave(filename = paste0("~/Desktop/VarioPath/plot_mpl_var/plot_mpl_", var, "_chr1_43348956_G_A_20210607.pdf"), 
         device = "pdf", dpi = "retina", plot = g, units = "cm", width = 20, height = 10)
}




for (var in colnames(dt_1)[15:dim(dt_1)[2]] ) {
      g = ggplot() + 
        geom_density(data = dt_1[f.22001.0.0==1], aes(x = get(var), y = ..density..), fill="blue", alpha=0.2) + 
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43338634_G_C" & f.22001.0.0==1, ..var]), 
                 xend = unlist(dt_1[VAR_ID == "chr1_43338634_G_C" & f.22001.0.0==1, ..var]),
                 y = 0,
                 yend = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8,
                 colour = "blue", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43338634_G_C", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8)) +
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43340042_C_T" & f.22001.0.0==1, ..var]),  
                 xend = unlist(dt_1[VAR_ID == "chr1_43340042_C_T" & f.22001.0.0==1, ..var]), 
                 y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8,
                 yend = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*2, # standard hight
                 colour = "blue", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43340042_C_T", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*2)) +
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43338705_T_-" & f.22001.0.0==1, ..var]), 
                 xend = unlist(dt_1[VAR_ID == "chr1_43338705_T_-" & f.22001.0.0==1, ..var]), 
                 y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*2,
                 yend = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*3, # standard hight
                 colour = "blue", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43338705_T_-", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*3)) +
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43337929_T_A" & f.22001.0.0==1, ..var]),  
                 xend = unlist(dt_1[VAR_ID == "chr1_43337929_T_A" & f.22001.0.0==1, ..var]), 
                 y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*3,
                 yend = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*4, # standard hight
                 colour = "blue", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43337929_T_A", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*4)) +
        geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T), linetype="dashed", size=1, color = "blue")  + 
        geom_density(data = dt_1[f.22001.0.0==0], aes(x = get(var), y = -..density..),  fill= "red", alpha =0.2)  + 
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43338634_G_C" & f.22001.0.0==0, ..var]), 
                 xend = unlist(dt_1[VAR_ID == "chr1_43338634_G_C" & f.22001.0.0==0, ..var]), 
                 y = 0,
                 yend = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8, # standard hight
                 colour = "red", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43338634_G_C", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8)) +
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43340042_C_T" & f.22001.0.0==0, ..var]), 
                 xend = unlist(dt_1[VAR_ID == "chr1_43340042_C_T" & f.22001.0.0==0, ..var]), 
                 y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8,
                 yend = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*2, # standard hight
                 colour = "red", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43340042_C_T", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*2)) +
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43338705_T_-" & f.22001.0.0==0, ..var]), 
                 xend = unlist(dt_1[VAR_ID == "chr1_43338705_T_-" & f.22001.0.0==0, ..var]), 
                 y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*2,
                 yend = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*3, # standard hight
                 colour = "red", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43338705_T_-", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*3)) +
        annotate("segment",
                 x = unlist(dt_1[VAR_ID == "chr1_43337929_T_A", ..var]), 
                 xend = unlist(dt_1[VAR_ID == "chr1_43337929_T_A", ..var]),
                 y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*3,
                 yend = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*4, # standard hight
                 colour = "red", alpha = 0.5,
        ) +
        # geom_text(aes(label = "chr1_43337929_T_A", x = quantile(dt_1[,..var],na.rm = T, names = F, 0.99), 
        #               y = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8*4)) +
        geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T), linetype="dashed", size=1, color = "red") +
        ggtitle(paste("MPL variants for", eval(var)))+
        coord_cartesian(xlim = c(quantile(dt_1[,..var], 0.005, na.rm = T, names = FALSE),
                                 quantile(dt_1[,..var], 0.995, na.rm = T, names = FALSE))) +
        theme_minimal()
      ggsave(filename = paste0("~/Desktop/VarioPath/plot_mpl_var/plot_mpl_", var, "_no_var_name_20210603.pdf"), 
             device = "pdf", dpi = "retina", plot = g, units = "cm", width = 20, height = 10)
}




for (var in colnames(dt_1)[15:dim(dt_1)[2]] ) {
  g = ggplot() + 
    geom_density(data = dt_1[f.22001.0.0==1], aes(x = get(var), y = ..density..), fill="blue", alpha=0.2) + 
    annotate("segment",
             x = unlist(dt_1[VAR_ID != "No MPL variant"  & f.22001.0.0==1, ..var]), 
             xend = unlist(dt_1[VAR_ID != "No MPL variant"  & f.22001.0.0==1, ..var]),
             y = 0,
             yend = max(density(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8,
             colour = "blue", alpha = 0.5,
    ) +
    geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T), linetype="dashed", size=1, color = "blue")  + 
    geom_density(data = dt_1[f.22001.0.0==0], aes(x = get(var), y = -..density..),  fill= "red", alpha =0.2)  + 
    annotate("segment",
             x = unlist(dt_1[VAR_ID != "No MPL variant"  & f.22001.0.0==0, ..var]), 
             xend = unlist(dt_1[VAR_ID != "No MPL variant" & f.22001.0.0==0, ..var]), 
             y = 0,
             yend = -max(density(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T)$y, na.rm = T, names = FALSE)/8, # standard hight
             colour = "red", alpha = 0.5,
    ) +
    geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T), linetype="dashed", size=1, color = "red") +
    ggtitle(paste("MPL variants for", eval(var)))+
    coord_cartesian(xlim = c(quantile(dt_1[,..var], 0.005, na.rm = T, names = FALSE),
                             quantile(dt_1[,..var], 0.995, na.rm = T, names = FALSE))) +
    theme_minimal()
  ggsave(filename = paste0("~/Desktop/VarioPath/plot_mpl_var/plot_mpl_", var, "_all_variants_20210603.pdf"), 
         device = "pdf", dpi = "retina", plot = g, units = "cm", width = 20, height = 10)
}



dt_1$VAR_ID[which(is.na(dt_1$VAR_ID))] = "No MPL variant"

for (var in colnames(dt_1)[15:dim(dt_1)[2]] ) {
  g = ggplot() + 
    geom_density(data = dt_1[f.22001.0.0==1 & VAR_ID == "No MPL variant" ], aes(x = get(var), y = ..density..), fill = "blue", alpha = 0.2 ) + 
    geom_density(data = dt_1[f.22001.0.0==1 & VAR_ID != "No MPL variant" ], aes(x = get(var), y = ..density.., color = VAR_ID)) + 
    geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==1, ..var]), na.rm = T), linetype="dashed", size=1, color = "blue")  + 
    geom_density(data = dt_1[f.22001.0.0==0 & VAR_ID == "No MPL variant" ], aes(x = get(var), y = -..density..), fill = "red", alpha = 0.2 ) +
    geom_density(data = dt_1[f.22001.0.0==0 & VAR_ID != "No MPL variant" ], aes(x = get(var), y = -..density.., color = VAR_ID))  + 
    geom_vline( xintercept=median(unlist(dt_1[f.22001.0.0==0, ..var]), na.rm = T), linetype="dashed", size=1, color = "red") +
    coord_cartesian(xlim = c(quantile(dt_1[,..var], 0.005, na.rm = T, names = FALSE),
                             quantile(dt_1[,..var], 0.995, na.rm = T, names = FALSE))) +
    theme_minimal()
  ggsave(filename = paste0("~/Desktop/VarioPath/plot_mpl_var/plot_mpl_", var, "_density_20210603.pdf"), 
         device = "pdf", dpi = "retina", plot = g, units = "cm", width = 20, height = 10)
}
  

fwrite(file = "~/Desktop/VarioPath/plot_mpl_var/mpl_variant_carriers_FBC.csv", x = dt_1[VAR_ID != "No MPL variant"], sep = ",", append = F, quote = F, row.names = F, col.names = T, nThread = 4)
