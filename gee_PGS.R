#####################################################################################################################
# This script is run to perform GEE association analysis between Mullen scores and polygenic scores
# using the BATSS sample for the paper: 
# The latent structure of emerging cognitive abilities: an infant twin study.
# Bussu G., Taylor M., Tammimies K., Ronald A., Falck-Ytter T.
#
# coded by Giorgia Bussu, June 2022
# 
####################################
# The script need as input developmental and polygenic scores from the selected dataset in the long format.
#
# The script gives as output GEE parameters from the different models tested, with corrected p-values.
#
#####################################################################################################################

rm(list=ls())

library(drgee)
library(gt)
library(tidyverse)

# set a working directory
setwd('your working directory path')
list.files() 

# import and check data

pgs_data<-read.csv('PGSdata_matched_final.csv',header=T,sep=',')
names(pgs_data); dim(pgs_data)

pgs_data<-pgs_data[-1]

data <- read.csv(file='phenotypic_data_resRAW_final_ELCwGM.csv',header=T,sep=',')
names(data); dim(data)

data<-data[-1]

# match and merge datasets
data<-data[match(pgs_data$id,data$id),]

data_pgstwin<-data.frame(cbind(data,pgs_data))

#####################################################################################################################
#### Data tidying and scaling
#####################################################################################################################

# remove redundant columns
indx_to_remove<-c(2:16, 23:26, 45:48)
data_pgstwin<-data_pgstwin[-indx_to_remove]

# scale PGS data
vars <- colnames(data_pgstwin)[c(8:25)]
data_pgstwin[vars]<-lapply(data_pgstwin[vars],as.numeric)
data_pgstwin[vars]<-lapply(data_pgstwin[vars],scale)

#####################################################################################################################
#### GEE models
#####################################################################################################################


#### Neurodevelopmental and psychiatric conditions

asd_sum<-summary( gee( formula = ELCstd ~ prs_asd + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
              data = data_pgstwin, clusterid = "Twinpair", cond = F) )

adhd_sum<-summary( gee( formula = ELCstd ~ prs_adhd + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
              data = data_pgstwin, clusterid = "Twinpair", cond = F) )

bip_sum<-summary( gee( formula = ELCstd ~ prs_bip + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
              data = data_pgstwin, clusterid = "Twinpair", cond = F) )

mdd_sum<-summary( gee( formula = ELCstd ~ prs_mdd + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
              data = data_pgstwin, clusterid = "Twinpair", cond = F) )

scz_sum<-summary( gee( formula = ELCstd ~ prs_scz + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
              data = data_pgstwin, clusterid = "Twinpair", cond = F) )

results <- matrix(0,7,3)
results<-data.frame(results)

results[1,]<-asd_sum$coefficients[2,c(1,2,4)]
results[2,]<-adhd_sum$coefficients[2,c(1,2,4)]
results[3,]<-bip_sum$coefficients[2,c(1,2,4)]
results[4,]<-mdd_sum$coefficients[2,c(1,2,4)]
results[5,]<-scz_sum$coefficients[2,c(1,2,4)]

#### IQ and educational attainment

ea_sum<-summary( gee( formula = ELCstd ~ prs_ea + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
                       data = data_pgstwin, clusterid = "Twinpair", cond = F) )

iq_sum<-summary( gee( formula = ELCstd ~ prs_iq + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11,
                       data = data_pgstwin, clusterid = "Twinpair", cond = F) )

results[6,]<-ea_sum$coefficients[2,c(1,2,4)]
results[7,]<-iq_sum$coefficients[2,c(1,2,4)]

p_psy<-p.adjust(results[1:5,3],method="fdr")
p_iq<-p.adjust(results[6:7,3],method="fdr")

results$p_fdr<-c(p_psy,p_iq)

#####################################################################################################################
#### Save and print results
#####################################################################################################################

results<-cbind(c('Autism','ADHD','Bipolar','MDD','SCZ','EduAttain','IQ'),results)
names(results)<-c('PGS test','Beta','Std','p','p_FDR')

# save results as csv and print table
write.csv(results,'GEE_PGS_MSEL_results.csv',row.names=F)
gt(results) %>% gtsave('results_GEE_PGS_ELC_BATSS.docx')

