library(forestmodel)
library(epitools)
library(epiDisplay)
library(tidyverse)
library(sjPlot)


########
# function to extract the ID information from the partecipants column
get.ids = function(partecipantID_vector){
  require(stringr)
  id = unlist(
    lapply(
      lapply(partecipantID_vector, 
             function(x)
               str_extract_all(x , pattern = "[^=| |;]+") 
      ), 
      '[[', 
      1)
  ) 
  id_clean = id[id %>% str_length() > 4]
  return(id_clean)
}
# function to count the cases at teh variant level
num.phenotype =  function(partecipantID_vector, case.list){
            apply(partecipantID_vector, 
                  1, 
                  function(x) 
                    for ( i in str_extract_all(x , pattern = "[^=| |;]+" ) ) { return(sum(i %in% case.list) )} )
            }

######## Prepare the cohort databases
# Load Luanluan phenotype extraction extraction 
load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/project/who1000-1/rds-who1000-cbrc/data/shared_luanluan_luca_UKB/Luanluan/data/ukb_dvt_pe_combined.rdata')
#dim(df)
#str(df)
########
# Import unrelated europena list to subset UKB 500K cohort
unrelated.european = read.csv('/Users/luca/Library/Group\ Containers/G69SCX94XU.duck/Library/Application\ Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/project/who1000-1/rds-who1000-cbrc/data/shared_luanluan_luca_UKB/unrelated_ethnic_groups/unrelated_UKB_200_ethnicity.txt.somalier-ancestry_not_1KGP_EUR.tsv',
                              sep = '\t', header = F) %>% 
  select(V1)
# dim(unrelated.european)
# str(unrelated.european)
########
# constraint gourp to unrelated european and remove all the others eid
df_eur_un = df[which(df$eid %in% unrelated.european$V1),]
# dim(df_eur_un) == dim(unrelated.european)
# str(df_eur_un)
# Cut unnecessary column (for the moment) to speed up the analysis.
df.eur.un.small = df_eur_un[,c("eid", "sex", "bmi", "ages", "smallbin", "vte")]
# Be sure that VTE is just the definition of ICD codes (i.e columns ps_pe and ps_dvt )
df.eur.un.small$vte = apply(df_eur_un[,c('ps_pe','ps_dvt')], 
                            1, 
                            function(x)
                            if ( sum(x) > 0) {
                              return(1)
                            } else {
                              return(0)
                            }
                            )
# Import the MDT thrombosis list to extract the case IDs 
cases = readxl::read_xls('~/Desktop/VarioPath/MDT_variants/V20210211/MDT_for_thrombosis.xls') 
cases_id = get.ids(cases$PARTECIPANTS)
# Keep only eur unrelated cases
cases_id = cases_id[cases_id %in% unrelated.european$V1]

######################### GENE BURDEN TEST ##############################################
# add variatn label
df.eur.un.small$has_variant = as.factor(
  unlist(
    lapply(df.eur.un.small$eid, 
           function(x)
             if (x %in% cases_id){
               return("1")
             } else {
               return("0")
             }
    )
  )
)
# Add gene information (i.e. in what gene is the variant)
df.eur.un.small$gene = 'None'

df.eur.un.small$gene = unlist( # unlist remove the lapply list structure
  paste( # some people can have more than one variant. 
    # This means that they return 2 or more genes. paste() collapse them.
    lapply(df.eur.un.small$eid, 
           function(x)
             if (x %in% cases_id) {
               cases$GENE[ grep(pattern = x, x = as.vector(cases$PARTECIPANTS)) ]
             } else {
               'None'
             }
    ),
    sep = " ")
)

for (gene in unique(df.eur.un.small$gene) ) {
  if ( gene != 'None' ) {
    print(gene)
    glm1 = glm(vte~has_variant, 
               family=binomial(link = "logit"), 
               data=df.eur.un.small[df.eur.un.small$gene == gene | df.eur.un.small$gene == 'None' ,])
    pdf(file = paste('Desktop/VarioPath/phenotype analysis/thrombosis_OR/',
                     gene, 'or.pdf', sep = ""), title = gene )
        print(forest_model(glm1))
    dev.off()
  }
}

######################### VARIANT TEST ##################################################
# select the cases ID (i.e. df.eur.un.small$vte == 1)
VTE.cases = df.eur.un.small[df.eur.un.small$vte == 1,'eid']
# Use num of phenotype to extract the number of cases that have the variant and have also the phenotype
cases$phenotype  = num.phenotype(cases[,'PARTECIPANTS'], VTE.cases)
# the number of cases that have variants but no phenotype. calculated as: n_cases_total_with_variant - n_cases_with_variant_and_variatns 
cases$no_phenotype = cases$male + cases$female - cases$phenotype

# The number of cases that don't have phenoytpe and don't have variant. (only unrelated european cohort)
no.pheno_no.variant = dim (df.eur.un.small[df.eur.un.small$gene == 'None' & 
                                        df.eur.un.small$has_variant == '0' & 
                                        df.eur.un.small$vte == '0',])[1]
# The number of cases that have the phenoytpe but don't have variants. (only unrelated european cohort)
pheno_no.variant = dim (df.eur.un.small[df.eur.un.small$gene == 'None' & 
                                             df.eur.un.small$has_variant == '0' & 
                                             df.eur.un.small$vte == '1',])[1]
# Calculate OR per variant
# OR ratio are calculated using the fisher exact test. the formula that it use is (a/b) / (c/d)
# a = variant and phenotype
# b = variant and no phenotype
# c = no variant and phenotype
# d = no variant and no phenoytpe
cases$or=0
cases$p.value=1
cases$conf_low=0
cases$conf_high=0
for (variant in 1:dim(cases)[1]) {
  pheno_variant = cases$phenotype[variant]
  no.pheno_variant = cases$no_phenotype[variant]
  # construct the 2x2 contingency table
  matrix2.2 <- matrix(c(pheno_variant, pheno_no.variant, no.pheno_variant, no.pheno_no.variant ), 
                      nrow = 2, 
                      ncol = 2, 
                      byrow = TRUE)
  rownames(matrix2.2) <- c("phenotype", "no_phenotype")
  colnames(matrix2.2) <- c("variant", "no_variant")
  # stat test. i.e. Fisher exact test
  test_stat = fisher.test(matrix2.2, alternative = 'greater', conf.int = T, or = TRUE, simulate.p.value = TRUE, B = 1000)
  # test_stat = chisq.test(matrix2.2, correct = TRUE, rescale.p = TRUE, simulate.p.value = TRUE)
  # test_stat2 = oddsratio.wald(matrix2.2, correction = T)
  cases$or[variant] = test_stat$estimate # calculated in fisher test
  cases$conf_low[variant] = test_stat$conf.int[1] 
  cases$conf_high[variant] = test_stat$conf.int[2]
  cases$p.value[variant] = test_stat$p.value
  }

# conf interval plot
# extract only the interesting ccolumns
cas_to_plot = cases[cases$or > 1 & cases$or != 'Inf' ,c("GENE","CHROM","POS","REF","ALT",
         "PARTECIPANTS","AF_ukb_calc","het",
         "hom","hem","male","female","phenotype",
         "no_phenotype","or","p.value","conf_low","conf_high")] %>% 
  arrange(or, desc = TRUE)
  
# create IDs for the y olumn of the plot
cas_to_plot$ID = paste(cas_to_plot$GENE,cas_to_plot$CHROM,
                       cas_to_plot$POS,cas_to_plot$REF,cas_to_plot$ALT, sep = "_")
# plot and save
ggplot(cas_to_plot, mapping = aes(y=ID, x = or, xmin=conf_low, xmax=conf_high))+
  geom_point() +
  geom_errorbarh(height=.05) +
  scale_x_continuous(limits=c(-2.5, 40) )


sinfo = sessionInfo()
save.image("~/Desktop/VarioPath/phenotype analysis/OR_thrombosis_preliminary_results_20210213.RData")

########################### END OF THE SCRIPT #################################################


