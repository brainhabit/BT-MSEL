#####################################################################################################################
# This script is run to prepare data for twin modelling (wide format)
# using the BATSS sample for the paper: 
# The latent structure of emerging cognitive abilities: an infant twin study.
# Bussu G., Taylor M., Tammimies K., Ronald A., Falck-Ytter T.
#
# coded by Giorgia Bussu, June 2022
# 
####################################
# The script needs as input developmental and demographic scores from the entire BATSS sample in the long format.
#
# The script gives as output the selected dataset based on inclusion criteria, in wide (twin_data.csv) and long format (phenotypic_data.csv).
#
#####################################################################################################################

rm(list=ls())

### set a working directory (otherwise create a new R project):

setwd('your folder')
list.files() 

### import and check data:

data_msel <- readxl::read_xlsx('Mullen.xlsx')
names(data_msel);dim(data_msel)

data_demo <- readxl::read_xlsx('Demographics_5m.xlsx')
names(data_demo);dim(data_demo)

data_exclusion<-readxl::read_xlsx('BabyTwins background summary and exclusion.xlsx')
names(data_exclusion); dim(data_exclusion)

### checking match across files

length(which(data_msel$`Participant EASE Code`!=data_demo$`Participant EASE Code`| data_msel$`Participant EASE Code`!=data_exclusion$Code | data_demo$`Participant EASE Code`!=data_exclusion$Code))
mm_indx<-which(data_msel$`Participant EASE Code`!=data_demo$`Participant EASE Code`| data_msel$`Participant EASE Code`!=data_exclusion$Code | data_demo$`Participant EASE Code`!=data_exclusion$Code)

# MSEL file had twins swapped within 5 pairs -> now matching to demographic and exclusion files
msel_matchingindx<-match(data_msel$`Participant EASE Code`,data_exclusion$Code)
data_msel<-data_msel[msel_matchingindx,]

#####################################################################################################################
### exclusion based on general criteria
#####################################################################################################################

excl_indx<-which(data_exclusion$`Exclusion version A`==1)

### additional exclusion criteria: exposure to Swedish

# select those who have Swedish as 2nd or 3rd language
lang2_indx<-which(data_demo$`Heard language 2`=='Swedish')
lang3_indx<-which(data_demo$`Heard language 3`=='Swedish')

# translate % exposure to language to numeric variable
data_demo$language2_home<-as.numeric(gsub("%","",data_demo$`Heard language 2- Proportion Heard by Child at home`))
data_demo$language2_else<-as.numeric(gsub("%","",data_demo$`Heard language 2- Proportion heard by child elsewhere`))

data_demo$language3_home<-as.numeric(gsub("%","",data_demo$`Heard language 3- Proportion Heard by Child at home`))
data_demo$language3_else<-as.numeric(gsub("%","",data_demo$`Heard language 3- Proportion heard by child elsewhere`))

# select those with exposure<10% NB here the indx goes over those who have Swedish as 2nd/3rd language, not whole sample
excl_indx_lang2<-which((data_demo$language2_else[lang2_indx]+data_demo$language2_home[lang2_indx])<0.1)
excl_indx_lang3<-which((data_demo$language3_else[lang3_indx]+data_demo$language3_home[lang3_indx])<0.1)

# exclusion indx reported to the whole sample
excl_indx_lang<-c(lang2_indx[excl_indx_lang2],lang3_indx[excl_indx_lang3])

# entire list exclusion criteria
to_exclude<-c(excl_indx,excl_indx_lang)

# clean data files from excluded participants
data_demo_clean<-data_demo[-to_exclude,]
data_msel_clean<-data_msel[-to_exclude,]
data_exclusion_clean<-data_exclusion[-to_exclude,]

# add dataset background info
data_background<- readxl::read_xlsx('Background_5m.xlsx')

backmatch<-match(data_background$Kod,data_msel_clean$`Participant EASE Code`)
data_background_matched<-data_background[which(!is.na(backmatch)),]

# double check background match
head(cbind(data_background_matched$Kod,data_msel_clean$`Participant EASE Code`))
which(data_background_matched$Kod!=data_msel_clean$`Participant EASE Code`)

# include only valid msel
indx_valid<-which(data_msel_clean$`Data validity code`=='Valid')

table(data_msel_clean$`Data validity code`)

indx_invalid<-which(data_msel_clean$`Data validity code`=='Invalid')
indx_unknown<-which(data_msel_clean$`Data validity code`=='Not-known')

data_msel_clean<-data_msel_clean[indx_valid,]
data_demo_clean<-data_demo_clean[indx_valid,]
data_background_matched<-data_background_matched[indx_valid,]
data_exclusion_clean<-data_exclusion[indx_valid,]

# add data on pregnancy
data_pregnancy<- readxl::read_xlsx('pregnancy_time_data.xlsx')

indx_preg<-match(data_msel_clean$`Participant EASE Code`,data_pregnancy$ID)

preg_week<-data_pregnancy$Gweeks[indx_preg]
preg_day<-data_pregnancy$Gdays[indx_preg]
preg_day[which(is.na(preg_day))]<-0

preg_term<-preg_week*7+preg_day

#####################################################################################################################
### dataset included in our analyses
#####################################################################################################################

data<-data.frame(cbind(data_msel_clean$`Participant EASE Code`,data_exclusion_clean$Zygosity,data_msel_clean$Gender,data_msel_clean$`Age at Date of Assessment (days)`,data_msel_clean$GM_raw_score,data_msel_clean$FM_raw_score,data_msel_clean$VR_raw_score,data_msel_clean$RL_raw_score,data_msel_clean$EL_raw_score,data_msel_clean$`Early Learning Composite Score`,preg_term,data_demo_clean$`Bio Mum Age`,data_demo_clean$`Bio Dad Age`,data_demo_clean$`A. Highest level of education`,data_demo_clean$`B. Highest level of education`,data_background_matched$`F16. Ungef?r hur h?g ?r familjens* gemensamma m?nadsinkomst innan skatt (l?n och andra ers?ttningar/bidrag)?`,data_demo_clean$TWAB,data_demo_clean$`Twin pair no.`))
names(data)<-c('id','zygosity','sex','age','GM','FM','VR','RL','EL','ELC','term_age','Mum_age','Dad_age','A_edu','B_edu','family_income','TWAB','Twinpair')

# parental education level to MAX within family
data$A_edu_num<-as.numeric(as.factor(data$A_edu))
data$A_edu_num[which(data$A_edu_num==4)]<-6
data$A_edu_num[which(data$A_edu_num==5)]<-4
data$A_edu_num[which(data$A_edu_num==6)]<-5

data$B_edu_num<-as.numeric(as.factor(data$B_edu))
data$B_edu_num[which(data$B_edu_num==4)]<-6
data$B_edu_num[which(data$B_edu_num==5)]<-4
data$B_edu_num[which(data$B_edu_num==6)]<-5

data$edu_mean<-rowMeans(cbind(data$A_edu_num,data$B_edu_num),na.rm = T)


# parental age to mean across mum and dad (NA.RM=TRUE for 3 missing Dad age)
data$Mum_age<-as.numeric(data$Mum_age)
data$Dad_age<-as.numeric(data$Dad_age)
data$parental_age<-rowMeans(cbind(data$Mum_age,data$Dad_age),na.rm = T)

#####################################################################################################################
##### basic data tidying
#####################################################################################################################

names(data);dim(data)

data<-data[,-c(12:15,19,20)]

table(data$zygosity)

# check how many pairs are incomplete by making a variable that counts 
#   the frequency of their pair ID:

data <- merge(data,data.frame(table(Twinpair=data$Twinpair)),by='Twinpair')
table(data$Freq) # check how many pair IDs appear only once 

# carry out the exclusion + make sure you keep singletons as NA:

data$Twinpair[which(data$Freq==1)]
data_test<-rbind(data[1:5,],data[5:76,],data[76:91,],data[91:114,],data[114:115,],data[115:168,],data[168:177,],data[177:236,],data[236:361,],data[361:474,],data[474:521,],data[521:567,])
data_test <- merge(data_test,data.frame(table(id=data_test$id)),by='id')
naindx<-matrix(which(data_test$Freq.y==2),nrow=2)[1,]

data_test$GM[naindx]<-NA
data_test$FM[naindx]<-NA
data_test$VR[naindx]<-NA
data_test$RL[naindx]<-NA
data_test$EL[naindx]<-NA
data_test$ELC[naindx]<-NA
data_test$TWAB[naindx[which(data_test$TWAB[naindx]==1)]]<-2
data_test$TWAB[naindx[-which(data_test$TWAB[naindx]==1)]]<-1

data<-data_test

dim(data) 

# remove the frequency variable:

data <- data[-c(17,18)]

#####################################################################################################################
##### transform data for twin analysis
#####################################################################################################################

# binary sex: 0=Females; 1=Males.
data$sex[which(data$sex=='Male')]<-1
data$sex[which(data$sex=='Female')]<-0
data$sex<-as.numeric(data$sex)

# ordinal discrete income
data$income<-as.numeric(as.factor(data$family_income))
data$income[which(data$income==2)]<-12
data$income[which(data$income==1)]<-2
data$income[which(data$income==11)]<-1
data$income[which(data$income==12)]<-11

data<-data[,-c(13)]

# numeric variables
vars <- colnames(data)[c(2,5:15)]
data[vars]<-lapply(data[vars],as.numeric)

### check the distribution of the Mullen scales:

vars <- colnames(data)[c(6:11)]
lapply(data[vars],psych::describe)

### control for the effects of covariates: sex, age, parental age, income, education.

data$age<-scale(data$age)
data$term_age<-scale(data$term_age)
data$parental_age<-scale(data$parental_age)
data$income<-scale(data$income)
data$edu_mean<-scale(data$edu_mean)


resGM <- resid(lm(scale(GM)~sex+age+term_age,data=data,na.action=na.exclude))
resFM <- resid(lm(scale(FM)~sex+age+term_age,data=data,na.action=na.exclude))
resVR <- resid(lm(scale(VR)~sex+age+term_age,data=data,na.action=na.exclude))
resRL <- resid(lm(scale(RL)~sex+age+term_age,data=data,na.action=na.exclude))
resEL <- resid(lm(scale(EL)~sex+age+term_age,data=data,na.action=na.exclude))

#resELC<- resid(lm(ELC~sex+age+term_age,data=data,na.action=na.exclude))
rawComposite<-data$GM+data$FM+data$VR+data$RL+data$EL
rawComposite<-15*scale(rawComposite)+100

resRAWC <- resid(lm(scale(rawComposite)~sex+age+term_age,data=data,na.action=na.exclude))

# standardize the residuals and add them to the data:

data$GMstd <- resGM
data$FMstd <- resFM
data$VRstd <- resVR
data$RLstd <- resRL
data$ELstd <- resEL
data$ELCstd <- resRAWC

# check that the standardization worked:

stand <- colnames(data)[c(17:22)]
lapply(data[stand],psych::describe)

#####################################################################################################################
######## Randomize twin order
#####################################################################################################################

### give each person a random number
set.seed(2022)
npairs <- nrow(data)/2
rand <- c(rep(4,npairs),rep(5,npairs))
data$rand <- sample(rand) # randomly choose 1 number for each person

### divide the twins into two subsets

twin1 <- subset(data,TWAB==1)
twin2 <- subset(data,TWAB==2)

### flip twin 2's random number

randVar <- data.frame(cbind(twin1$Twinpair,twin1$rand))
colnames(randVar) <- c('Twinpair','rand')

twin2 <- twin2[-c(23)]
twin2 <- merge(twin2,randVar,by='Twinpair')

# check
table(twin1$rand);table(twin2$rand)

# for twin 1, convert the random numbers to 0s and 1s:

twin1$rand[twin1$rand==4] <- 0
twin1$rand[twin1$rand==5] <- 1

# and then recode twin 2, but in the opposite direction to twin 1:

twin2$rand[twin2$rand==4] <- 1
twin2$rand[twin2$rand==5] <- 0

# check:

table(twin1$rand);table(twin2$rand) 

# now join twin 1 and twin 2 back together:

data <- data.frame(rbind(twin1,twin2))
names(data);dim(data)

# sort the data by pairnr and check the first few rows:

data <- data[order(data[,2]),]
head(data)

# split the data by random number:

twin1 <- subset(data,rand==1)
twin2 <- subset(data,rand==0)

# check that the dimensions are the same:

dim(twin1);dim(twin2)

# add '1' to the end of the variables for twin 1 and '2' for twin 2:

twin1labs <- paste(colnames(twin1),1,sep='')
twin2labs <- paste(colnames(twin2),2,sep='')

names(twin1) <- twin1labs
names(twin2) <- twin2labs

# combine data so that 1 pair each line
dataD <- data.frame(cbind(twin1,twin2))

# remove unused variables:

dataD <- dataD[-c(6:16,23:39,46)]

# relabel a few variables:

names(dataD)[names(dataD)=='Twinpair1'] <- 'Twinpair'
names(dataD)[names(dataD)=='id1'] <- 'id'
names(dataD)[names(dataD)=='zygosity1'] <- 'zygosity'
names(dataD)[names(dataD)=='sex1'] <- 'sex'
names(dataD)[names(dataD)=='age1'] <- 'age'

#####################################################################################################################
#### Save data
#####################################################################################################################

write.csv(data,'phenotypic_data.csv')
write.csv(dataD,'twin_data.csv')
