
library("tidyverse")
library("data.table")
df = data.table::fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP-1/UKBb/variopath/VCF_200K_filtered_for_variopath/no_AF_filter/All_variopath_variants_in_unrelated_european_noAFfilter.tsv", nThread = 4)                       
setkey(df, VAR_ID)

length(unique(df$VAR_ID))


mpl = c('chr1_43338634_G_C',
'chr1_43338705_T_-',
'chr1_43337929_T_A',
'chr1_43340042_C_T')
setkey(df, VAR_ID)
mpl_people = df[VAR_ID %in% mpl]

fwrite(mpl_people, "~/Desktop/mdt_carriers.tab", quote = F, sep ="\t", row.names = F, col.names = T )

length(unique(df$VAR_ID))

df2=unique(df[sort(df$PARTICIPANT_ID),])
plot(hist(table(df2$PARTICIPANT_ID), breaks = c(1,2,3,4,5,6,7,8), ))
