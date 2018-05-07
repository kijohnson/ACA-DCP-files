########################################################################################
########################################################################################
##### R Code to Analyze the SEER Data to Assess Effect of ACA-DCP on AYA Cancer Outcomes (Barnes)
########################################################################################
########################################################################################

#####################################################################
#### Prepare SEER DATA
#####################################################################

## Read in data
data <- read.table("~/export4SEER/export4.txt",sep=";",header=T,quote = "")

## make R objects for a couple objects (for defining a couple more exclusions not specified in the SEER*STAT application)
seer.diagnoseyear <- data$X.Year.of.diagnosis.
seer.cancertype <- data$X.AYA.site.recode.WHO.2008.
seer.stage <- data$X.Derived.AJCC.Stage.Group..6th.ed..2004... 
seer.age <- data$X.Age.at.diagnosis. # equivalent to the data$X.Age.recode.with.single.ages.and.85.. variable

## Define additional exclusions
removes.year <- which(seer.diagnoseyear == 2010)
removes.stage <- which((seer.stage=="0" | seer.stage=="0a" | seer.stage=="0is") & seer.cancertype !="8.5.2 Carcinoma of bladder")
removes.age <- which(seer.age<19 | seer.age==26 | seer.age>29)  ####### THIS IS NEW

## remove the observations that met the exclusion criteria
data <- data[-unique(c(removes.year,removes.stage,removes.age)),]

## define even more exlusions (to be removed later--need to make a larger set of data for insurance analyses)
seer.age <- data$X.Age.at.diagnosis. # equivalent to the data$X.Age.recode.with.single.ages.and.85.. variable
seer.cancertype <- data$X.AYA.site.recode.WHO.2008.
seer.stage <- data$X.Derived.AJCC.Stage.Group..6th.ed..2004... #data$X.Derived.AJCC.Stage.Group..7th.ed..2010...
removes.ages <- NULL
removes.cancertypes <- NULL # taking out CNS tumors, Leukemia
for(i in 1:nrow(data)){
  chopped <- strsplit(as.character(seer.cancertype[i]),split=NULL)[[1]]
  if((chopped[1]=="1" | chopped[1]=="3") & chopped[2]==".") removes.cancertypes <- c(removes.cancertypes,i)
}
removes.stage <- c(which(is.na(seer.stage)))
removes.seer <- removes <- unique(c(removes.ages, removes.cancertypes,removes.stage))


### make R objects for the covariates of interest using the variable names from the dataset
seer.age <- data$X.Age.at.diagnosis. 
seer.diagnoseyear <- data$X.Year.of.diagnosis.
seer.diagnosemonth <- data$X.Month.of.diagnosis.
seer.cancertype <- data$X.AYA.site.recode.WHO.2008.
seer.stage <- data$X.Derived.AJCC.Stage.Group..6th.ed..2004...
seer.stage2 <- data$X.Derived.AJCC.Stage.Group..7th.ed..2010...
seer.alive <- data$X.Vital.status.recode..study.cutoff.used..
seer.alive2 <- data$X.SEER.cause.specific.death.classification
seer.survmonths <- data$X.Survival.months.
seer.insurance <- data$X.Insurance.Recode..2007...
seer.marital <- data$X.Marital.status.at.diagnosis.
seer.race <- data$X.Race.recode..W..B..AI..API..
seer.hispanic <- data$X.NHIA.Derived.Hisp.Origin.
seer.sex <- data$X.Sex.
seer.state <- data$X.State.
seer.residence <- data$X.Rural.Urban.Continuum.Code.2013.  
seer.residence2 <- data$X.Rural.Urban.Continuum.Code.2003.
seer.countyeducation <- data$X.....High.school.education.ACS.2008.12. 
seer.countyeducation.old <- data$X.....High.school.education.ACS.2007.11.  ## for sensitivity analyses
seer.countyeducation.new <- data$X.....High.school.education.ACS.2011.15.  ## for sensitivity analyses
seer.countyincome <- data$X.Median.household.income..in.tens..ACS.2008.12. 
seer.countyincome.old <- data$X.Median.household.income..in.tens..ACS.2007.11. ## for sensitivity analyses
seer.countyincome.new <- data$X.Median.household.income..in.tens..ACS.2011.15. ## for sensitivity analyses
seer.site <- data$X.CS.Schema...AJCC.6th.Edition.

### Modify the classes within the variables and apply desired formatting
marital <- (ifelse(seer.marital=="Married (including common law)",1,0)+ifelse(seer.marital=="Unknown",2,0)) # married=1, unmarried=0, unk=2
marital[which(marital==2)] <- NA
race <- (ifelse(seer.race=="Black",1,0) + ifelse(seer.race=="Asian or Pacific Islander",2,0) + ifelse(seer.race=="American Indian/Alaska Native",3,0) + ifelse(seer.race=="Unknown",4,0)) # white=0, black=1, Asian =2, American Indian=3, unk  =4
race[which(race==4)] <- NA
hispanic <- ifelse(seer.hispanic=="Non-Spanish-Hispanic-Latino",0,1) # hispanic=1
sex <- ifelse(seer.sex=="Female",1,0) # female=1
age <- seer.age
diagnoseyear <- seer.diagnoseyear
diagnosemonth <- seer.diagnosemonth
diagnoseyear2 <- diagnoseyear + ifelse(diagnosemonth %in% c("April","May","June"),0.25,0) + ifelse(diagnosemonth %in% c("July","August","September"),0.5,0) + ifelse(diagnosemonth %in% c("October","November","December"),0.75,0)
state <- seer.state
residence  <- (ifelse(seer.residence=="Counties in metropolitan areas ge 1 million pop" | seer.residence=="Counties in metropolitan areas of 250,000 to 1 million pop" | seer.residence=="Counties in metropolitan areas of lt 250 thousand pop",1,0) + 
                 ifelse(seer.residence=="Unknown/missing/no match" | seer.residence=="Unknown/missing/no match (Alaska or Hawaii - Entire State)",2,0))
residence[which(residence==2)] <- NA
raceethnicity.han <- 1+ifelse(race==0 & hispanic==0,-1,0) + ifelse(race==1 & hispanic==0,1,0) + ifelse(hispanic==1,2,0) # 0=non-hispanic white, 1=non-hispanic other,  2=non-hispanic black. 3=hispanic 
countyeducation <- as.numeric(as.character(seer.countyeducation))
countyeducation.han <- (ifelse(countyeducation>=1090,1,0) + ifelse(countyeducation>=1410,1,0) + ifelse(countyeducation>=2070,1,0))  ## 
intervention.group <- ifelse(age<26,1,0)
post.aca <- ifelse(diagnoseyear>2010,1,0)
countyincome <- as.numeric(as.character(seer.countyincome))
countyincome.even <- ifelse(countyincome>=quantile(countyincome,0.25,na.rm=T),1,0) + ifelse(countyincome>=quantile(countyincome,0.5,na.rm=T),1,0) + ifelse(countyincome>=quantile(countyincome,0.75,na.rm=T),1,0) # 5125 5754 6935

### A few more variable derivations for sensitivity analyses
residence2  <- (ifelse(seer.residence2=="Counties in metropolitan areas ge 1 million pop" | seer.residence2=="Counties in metropolitan areas of 250,000 to 1 million pop" | seer.residence2=="Counties in metropolitan areas of lt 250 thousand pop",1,0) + 
                  ifelse(seer.residence2=="Unknown/missing/no match" | seer.residence2=="Unknown/missing/no match (Alaska or Hawaii - Entire State)",2,0))
residence2[which(residence2==2)] <- NA
residence2[which(diagnoseyear>2010)] <- residence[which(diagnoseyear>2010)]
raceethnicity <- ifelse(race==0 & hispanic==0,1,0) + ifelse(race==1 & hispanic==0,2,0) + ifelse(hispanic==0 & race==4,3,0) + ifelse(hispanic==1,4,0)  # 0=non-hispanic other, 1=non-hispanic white, 2=non-hispanic black. 3=non-hispanic unk, 4=hispanic 
countyeducation.old <- as.numeric(as.character(seer.countyeducation.old))
countyeducation.new <- as.numeric(as.character(seer.countyeducation.new))
countyeducation2 <- countyeducation.old
countyeducation2[which(diagnoseyear>2010)] <- countyeducation.new[which(diagnoseyear>2010)]
countyeducation2.han <- (ifelse(countyeducation2>=1090,1,0) + ifelse(countyeducation2>=1410,1,0) + ifelse(countyeducation2>=2070,1,0))  ## 
countyincome.old <- as.numeric(as.character(seer.countyincome.old))
countyincome.new <- as.numeric(as.character(seer.countyincome.new))
countyincome2 <- countyincome.old*1.0522
countyincome2[which(diagnoseyear>2010)] <- countyincome.new[which(diagnoseyear>2010)]
countyincome.kim <- ifelse(countyincome>=3800,1,0) + ifelse(countyincome>=4800,1,0) + ifelse(countyincome>=6300,1,0)
countyincome2.kim <- ifelse(countyincome2>=3800,1,0) + ifelse(countyincome2>=4800,1,0) + ifelse(countyincome2>=6300,1,0)
countyincome2.even <- ifelse(countyincome2>=quantile(countyincome2,0.25,na.rm=T),1,0) + ifelse(countyincome2>=quantile(countyincome2,0.5,na.rm=T),1,0) + ifelse(countyincome2>=quantile(countyincome2,0.75,na.rm=T),1,0) 

#### Group cancer types (consistent with 2016 Han et al.)
cervix <- (seer.site == "Cervix")
cancertype <- as.character(seer.cancertype)
leukemias <- braintumors <- NULL
for(i in 1:nrow(data)){
  chopped <- strsplit(as.character(seer.cancertype[i]),split=NULL)[[1]]
  if((chopped[1]=="3") & chopped[2]==".") braintumors <- c(braintumors,i)
  if((chopped[1]=="1") & chopped[2]==".") leukemias <- c(leukemias,i)
}
cancertype[braintumors] <- "Brain"
cancertype[leukemias] <- "Leukemia"
cancertype[which(cancertype == "8.5.4 Carcinoma of cervix and uterus" & cervix)] <- "8.5.4 Carcinoma of cervix"
cancertype[which(cancertype %in% c("5.1 Fibromatous neoplasms","5.2 Rhabdomyosarcoma","5.3.1.1 Specified (excluding Kaposi sarcoma)","5.3.2 Unspecified soft tissue sarcoma"))] <- "5.0 Soft Tissue Sarcomas (excluding Kaposi sarcoma)" 
cancertype[which(cancertype %in% c("8.2.1 Nasopharyngeal carcinoma","8.2.2 Other sites in lip, oral cavity and pharynx","8.2.3 Nasal cav,mid ear,sinus,larynx,ill-def head/neck"))] <- "8.2 Carcinoma of head and neck"
cancertype[which(cancertype %in% c("4.1 Osteosarcoma","4.2 Chondrosarcoma","4.3 Ewing tumor","4.4 Other specified and unspecified bone tumors"))] <- "4.0 Osseous and chondromatous neoplasms"
cancertype[which(cancertype != "8.1 Thyroid carcinoma" & cancertype != "6.1 Germ cell and trophoblastic neoplasms of gonads" &
                   cancertype != "7.1 Melanoma" & cancertype != "8.4 Carcinoma of breast" & cancertype != "2.1 Non-Hodgkin lymphoma" &
                   cancertype != "2.2 Hodgkin lymphoma" & cancertype != "8.5.4 Carcinoma of cervix" & cancertype != "8.6.1 Carcinoma of colon and rectum" &
                   cancertype != "5.0 Soft Tissue Sarcomas (excluding Kaposi sarcoma)" & cancertype != "8.5.1 Carcinoma of kidney" &
                   cancertype != "8.2 Carcinoma of head and neck" & cancertype !="4.0 Osseous and chondromatous neoplasms" & 
                   cancertype != "8.5.3 Carcinoma of gonads" & cancertype != "8.6.2 Carcinoma of stomach" & cancertype != "8.3 Carcinoma of trachea,bronchus, and lung" & cancertype != "Brain" &  cancertype != "Leukemia")] <- "99.9.9.9 Other"
cancertype[which(cancertype=="8.4 Carcinoma of breast" & sex==0)] <- "99.9.9.9 Other"
cancertype <- as.factor(cancertype)
cancertype2 <- cancertype
cancertype2[which(cancertype=="Brain" | cancertype=="Leukemia")] <- "99.9.9.9 Other" ## have to do this to make imputation work.
cancertype2 <- droplevels(cancertype2)

### Format the response variables (insurance, stage, survival)
alive <- ifelse(seer.alive=="Alive",1,0) # 1=alive ## death from all causes
alive2 <- ifelse(seer.alive2=="Alive or dead of other cause",1,0) ## death from cancer
alive2[which(seer.alive2=="N/A not first tumor")] <- NA
insurance <- (ifelse(seer.insurance=="Any Medicaid",1,0) + ifelse(seer.insurance=="Insured" | seer.insurance=="Insured/No specifics",2,0) + ifelse(seer.insurance=="Insurance status unknown",9,0)) # 0=no insurance, 1 = medicaid, 2 = other insurance, 9=unk
insurance[which(insurance==9)] <- NA
stage <- numeric(length(seer.stage))
for(i in 1:nrow(data)){
  chopped <- strsplit(as.character(seer.stage[i]),split=NULL)[[1]]
  val <- sum(chopped=="I")
  if("V" %in% chopped) val <- 4
  if(is.na(seer.stage[i]) | seer.stage[i]=="UNK Stage") val <- NA
  else if(seer.stage[i]=="OCCULT") val <- 4
  stage[i] <- val
}
stage[which(stage==0)] <- 1 ## Combining stage 0 with 1 (since the only stage 0 are for bladder, and there are only a couple hundred cases)
survmonths <- as.numeric(as.character(seer.survmonths))

latestage <- ifelse(stage==4,1,0)
earlystage <- ifelse(stage<=1,1,0)

## A different stage derivation for a sensitivity analysis (uses AJCC 7 for post-2010 cases)
stage2 <- stage 
for(i in which(diagnoseyear>=2010)){
  chopped <- strsplit(as.character(seer.stage2[i]),split=NULL)[[1]]
  val <- sum(chopped=="I")
  if("V" %in% chopped) val <- 4
  if(is.na(seer.stage2[i]) | seer.stage2[i]=="UNK Stage") val <- NA
  else if(seer.stage2[i]=="OCCULT") val <- 4
  stage2[i] <- val
}
stage2[which(stage2==0)] <- 1

## Combine derived variables into a dataset
data.insurance <- data.insuranceno <- data.frame(earlystage=earlystage,insurance=as.factor(insurance),binarysurvival=rep(1,length(insurance)),latestage=latestage, intervention.group, post.aca, 
                                                 marital=as.factor(marital), race=as.factor(raceethnicity.han), countyeducation=as.factor(countyeducation.han), 
                                                 countyincome=as.factor(countyincome.even), residence=as.factor(residence), age, sex, diagnoseyear, diagnoseyear2, cancertype, 
                                                 alive, alive2, survmonths, stage=as.factor(stage),race.old=as.factor(race),hispanic.old=hispanic,stage.cov=as.factor(stage),insurance.cov=as.factor(insurance),
                                                 stage2=as.factor(stage2),stage2.cov=as.factor(stage2),countyincome2=as.factor(countyincome2.even),countyeducation2=as.factor(countyeducation2.han),
                                                 countyincome.even=as.factor(countyincome.kim),countyincome2.even=as.factor(countyincome2.kim),residence2=as.factor(residence2),statecounty=as.factor(data$X.State.county.),state=as.factor(state))

## Generating additional survival information for multiple imputation. see White & Royston 2009. Imputing missing covariate values for the Cox model. Statistics in Medicine. DOI: 10.1002/sim.3618
library(survival)
my.surv <- Surv(data.insurance$survmonths, 1-data.insurance$alive2)
my.fit  <- summary(survfit(my.surv ~ 1))
h.sort.of <- my.fit$n.event / my.fit$n.risk
H.tilde   <- cumsum(h.sort.of)
data.insurance$H.T <- rep(NA,nrow(data.insurance))
for(i in 1:length(my.fit$time)){
  data.insurance$H.T[which(data.insurance$survmonths>=my.fit$time[i])] <- H.tilde[i]
}
data.insuranceno$H.T <- data.insurance$H.T

## Mark observations removed from stage at diagnosis analyses (i.e. differentiate missing due to not applicable staging [mark with "9"] and simply missing [mark with "NA"])
data.insurance$stage <- data.insurance$stage.cov <-   data.insuranceno$stage <- data.insuranceno$stage.cov <- as.character(data.insurance$stage)
data.insurance$stage[removes] <- 9
data.insurance$stage <- data.insurance$stage.cov <-   data.insuranceno$stage <- data.insuranceno$stage.cov <- as.factor(data.insurance$stage)
data.insurance$stage2 <- data.insurance$stage2.cov <-   data.insuranceno$stage2 <- data.insuranceno$stage2.cov <- as.character(data.insurance$stage2)
data.insurance$stage2[removes] <- 9
data.insurance$stage2 <- data.insurance$stage2.cov <-   data.insuranceno$stage2 <- data.insuranceno$stage2.cov <- as.factor(data.insurance$stage2)

data.han <- data.noimpute <- droplevels(data.insurance[-removes,])

data.insurance.seer <- data.insuranceno.seer <- data.insurance
data.han.seer <- data.han

#####################################################################
############## PREPARE NCDB DATA
#####################################################################

data <- read.csv("~/Downloads/NCDB1198ALLv2.csv") # 2015 submission data

### GENERATE CANCER TYPES

# importing the AYA Recode definitions as downloaded from the SEER website: https://seer.cancer.gov/ayarecode/aya-who2008.html
recodes <- read.table("~/Downloads/ayarecodewho2008.txt",sep = ";",header=T)
recodes <- recodes[-which(is.na(recodes[,5])),]
recodes[which(recodes[,1]==""),1] <- recodes[which(recodes[,1]=="")-1,1]
rownames(recodes) <- seq(1,nrow(recodes))
recodes[21:27,1] <- recodes[20,1]
recodes$Primary.Site2 <- gsub("C", "", recodes[,3])
for(i in 1:ncol(recodes)){
  recodes[,i] <- as.character(recodes[,i])
}

# function to expand the ranges and stuff; prep format for later
expander2 <- function(x){
  x <- as.character(x)
  for(i in 1:length(x)){
    newx <- NULL # I'll be adding all the values (everything that was separated by commas, including the values denoted in the ranges)
    newVec <- strsplit(x[i], ",")[[1]] # this will make a character vector, separating elements from where the commas were
    expanders <- grep("-",newVec) # find all the elements that represent a range of values
    for(j in 1:length(newVec)){
      if(j %in% expanders){
        subVec <- strsplit(newVec[j], "-")[[1]] # the newVec[j] element represents a range, and this just splits the range into 2 elements from where the "-" was
        additions <- seq(as.numeric(subVec[1]),as.numeric(subVec[2])) # create a sequence between (and including) those numbers
        newx <- c(newx,additions) 
      }
      else newx <- c(newx,newVec[j])
    }
    x[i] <- paste(newx,collapse=",")
  }
  return(x)
}

# expanding the histology and primary site information based on the function above
histology <- expander2(recodes[,4])
primarysite <- expander2(recodes[,6])
# don't need to expand the behavior since it's already in right format

# the rd (recode definition) dataframe has all the pertinent information.
rd <- data.frame(CancerType=as.character(recodes[,1]),Behavior=as.character(recodes[,2]),PrimarySite=as.character(primarysite),Histology=as.character(histology))
for(i in 1:ncol(rd)){
  rd[,i] <- as.character(rd[,i])
}
## I'm importing updated cancer type names for compatibility with code I've already written (the names don't quite match-funny spaces and stuff; shouldn't affect analyses though)
rd.names <- read.csv("cancer-type-definition.csv")[,-1]
rd[,1] <- as.character(rd.names[,1])

data$PRIMARY_SITE <- as.numeric(gsub("C", "", data$PRIMARY_SITE)) # get rid of the "C" in primary site and convert to numeric; makes it way easier in loop below

cancertype <- as.character(rep(NA,nrow(data)))
for(i in 1:(nrow(rd)-1)){
  bs <- as.numeric(strsplit(rd$Behavior[i], ",")[[1]])
  hs <- as.numeric(strsplit(rd$Histology[i], ",")[[1]])
  ps <- as.numeric(strsplit(rd$PrimarySite[i], ",")[[1]])
  cancertype[which(data$BEHAVIOR %in% bs & data$HISTOLOGY %in% hs & data$PRIMARY_SITE %in% ps)] <- rd$CancerType[i]
}
cancertype[which(is.na(cancertype))] <- "Unclassified and Non-Malignant"

data$Cancer.Type.AYA.Recode.WHO.2008 <- cancertype # the "AYA Recode" version of cancer type

## formatting for the analysis (with cancer type definitions from the Han et al. 2016 paper)

seer.cancertype <- cancertype
## need to define carcinoma of the cervix (separated from uterus in the Han analysis but not in the normal AYA recode)
# following definitions taken from https://staging.seer.cancer.gov/cs/schema/02.05.50/cervix/?breadcrumbs=(~schema_list~),(~view_schema~,~cervix~),(~view_input~,~cervix~,~ssf25~)
cervix <- (data$BEHAVIOR==3 & data$HISTOLOGY %in% c(8000:9136,9141:9582,9700:9701) & data$PRIMARY_SITE %in% c(530,531,538,539))

leukemias <- braintumors <- NULL
for(i in 1:nrow(data)){
  chopped <- strsplit(as.character(seer.cancertype[i]),split=NULL)[[1]]
  if((chopped[1]=="3") & chopped[2]==".") braintumors <- c(braintumors,i)
  if((chopped[1]=="1") & chopped[2]==".") leukemias <- c(leukemias,i)
}

cancertype[braintumors] <- "Brain"
cancertype[leukemias] <- "Leukemia"
cancertype[which(cancertype == "8.5.4 Carcinoma of cervix and uterus" & cervix)] <- "8.5.4 Carcinoma of cervix"
cancertype[which(cancertype %in% c("5.1 Fibromatous neoplasms","5.2 Rhabdomyosarcoma","5.3.1.1 Specified (excluding Kaposi sarcoma)","5.3.2 Unspecified soft tissue sarcoma"))] <- "5.0 Soft Tissue Sarcomas (excluding Kaposi sarcoma)" 
cancertype[which(cancertype %in% c("8.2.1 Nasopharyngeal carcinoma","8.2.2 Other sites in lip, oral cavity and pharynx","8.2.3 Nasal cav,mid ear,sinus,larynx,ill-def head/neck"))] <- "8.2 Carcinoma of head and neck"
cancertype[which(cancertype %in% c("4.1 Osteosarcoma","4.2 Chondrosarcoma","4.3 Ewing tumor","4.4 Other specified and unspecified bone tumors"))] <- "4.0 Osseous and chondromatous neoplasms" 
cancertype[which(cancertype != "8.1 Thyroid carcinoma" & cancertype != "6.1 Germ cell and trophoblastic neoplasms of gonads" &
                   cancertype != "7.1 Melanoma" & cancertype != "8.4 Carcinoma of breast" & cancertype != "2.1 Non-Hodgkin lymphoma" &
                   cancertype != "2.2 Hodgkin lymphoma" & cancertype != "8.5.4 Carcinoma of cervix" & cancertype != "8.6.1 Carcinoma of colon and rectum" &
                   cancertype != "5.0 Soft Tissue Sarcomas (excluding Kaposi sarcoma)" & cancertype != "8.5.1 Carcinoma of kidney" &
                   cancertype != "8.2 Carcinoma of head and neck" & cancertype !="4.0 Osseous and chondromatous neoplasms" & 
                   cancertype != "8.5.3 Carcinoma of gonads" & cancertype != "8.6.2 Carcinoma of stomach" & cancertype != "8.3 Carcinoma of trachea,bronchus, and lung" & cancertype != "Brain" &  cancertype != "Leukemia")] <- "99.9.9.9 Other"
cancertype[which(cancertype=="8.4 Carcinoma of breast" & data$SEX==1)] <- "99.9.9.9 Other"

data$Cancer.Type.Justin.Recode <- cancertype ## consistent with the Han et al. 2016 paper

###### Prepare the rest of the variables

## make R objects for a couple objects (for defining a couple more exclusions not specified in the SEER*STAT application)
seer.diagnoseyear <- data$YEAR_OF_DIAGNOSIS #data$year_of_diagnosis
seer.cancertype <- data$Cancer.Type.AYA.Recode.WHO.2008 
seer.stage <- data$ANALYTIC_STAGE_GROUP
seer.stage[which(seer.stage=="" | seer.stage=="5")] <- 9
seer.age <- data$AGE

## Define exclusions
removes.year <- which(seer.diagnoseyear == 2010 | seer.diagnoseyear<2007 | seer.diagnoseyear>2014)
removes.stage <- which((seer.stage==0) & seer.cancertype !="8.5.2 Carcinoma of bladder")
removes.age <- which(seer.age<19 | seer.age==26 | seer.age>29)
tab <- table(data$PUF_FACILITY_ID,data$YEAR)
bads <- which(rowSums(tab[,-c(1,2,3,7,12)]>0)<7)
excludes <- rownames(tab)[bads]
removes.facilities <- which(data$PUF_FACILITY_ID %in% excludes)
removes.notprimaries <- which(data$SEQUENCE_NUMBER>1)

## remove the observations that met the exclusion criteria
data <- data[-unique(c(removes.year,removes.stage,removes.age,removes.facilities,removes.notprimaries)),]

## define more exlusions (to be removed later--need to make a larger set of data for insurance analyses)
seer.cancertype <- data$Cancer.Type.Justin.Recode
seer.stage <- data$ANALYTIC_STAGE_GROUP
seer.stage[which(seer.stage=="5")] <- 9
seer.age <- data$AGE
seer.diagnoseyear <- data$YEAR_OF_DIAGNOSIS
removes.ages <- NULL
removes.cancertypes <- which(seer.cancertype=="Brain" | seer.cancertype=="Leukemia") # taking out CNS tumors, Leukemia
removes.stage <- c(which(seer.stage==8))
removes.ncdb <- removes <- unique(c(removes.ages, removes.cancertypes,removes.stage))

### Format the variables to match (as much as possible) the derivations from the SEER analysis
seer.stage[which(seer.stage==0)] <- 1
if(sum(seer.stage==5 | seer.stage==9)>0) seer.stage[which(seer.stage==5 | seer.stage==9)] <- NA
seer.stage[removes] <- NA
seer.race <- data$RACE
seer.race[which(seer.race==99)] <- NA
seer.race[which(seer.race>2)] <- 3
seer.race <- seer.race-1
seer.hispanic <- data$SPANISH_HISPANIC_ORIGIN
seer.hispanic[which(seer.hispanic==9)] <- NA
seer.hispanic[which(seer.hispanic>0)] <- 1
seer.cancertype <- data$Cancer.Type.Justin.Recode
seer.countyincome <- data$MED_INC_QUAR_12-1
seer.countyeducation <- data$NO_HSD_QUAR_12-1
seer.residence <- data$UR_CD_13
seer.residence <- ifelse(seer.residence<4,1,0)
seer.insurance <- data$INSURANCE_STATUS
seer.insurance[which(seer.insurance %in% c(3,4))] <- 1
seer.insurance[which(seer.insurance==9)] <- NA
seer.insurance <- 3-seer.insurance 
seer.insurance[which(seer.insurance==3)] <- 0 #0=no insurance, 1 = medicaid, 2 = other insurance, 9=unk
seer.sex <- data$SEX-1 # female=1
seer.survmonths <- data$DX_LASTCONTACT_DEATH_MONTHS
seer.alive <- data$PUF_VITAL_STATUS
intervention.group <- ifelse(seer.age<26,1,0)
post.aca <- ifelse(seer.diagnoseyear>=2010,1,0)
raceethnicity.han <- 1+ifelse(seer.race==0 & seer.hispanic==0,-1,0) + ifelse(seer.race==1 & seer.hispanic==0,1,0) + ifelse(seer.hispanic==1,2,0) # 0=non-hispanic white, 1=non-hispanic other,  2=non-hispanic black. 3=hispanic 

## Create a dataframe with all the derived variables
data.insurance <- data.insuranceno <- data.frame(earlystage=ifelse(seer.stage==1,1,0),insurance=as.factor(seer.insurance),latestage=ifelse(seer.stage==4,1,0), intervention.group, post.aca, 
                                                 marital=rep(1,length(post.aca)), race=as.factor(raceethnicity.han), countyeducation=as.factor(seer.countyeducation), 
                                                 countyincome=as.factor(seer.countyincome), residence=as.factor(seer.residence), age=seer.age, sex=seer.sex, diagnoseyear=seer.diagnoseyear, cancertype=as.factor(seer.cancertype), 
                                                 alive=seer.alive, alive2=seer.alive, survmonths=seer.survmonths, stage=as.factor(seer.stage),race.old=as.factor(seer.race),hispanic.old=seer.hispanic,stage.cov=as.factor(seer.stage),insurance.cov=as.factor(seer.insurance),
                                                 stage2=as.factor(seer.stage),stage2.cov=as.factor(seer.stage),countyincome2=as.factor(seer.countyincome),countyeducation2=as.factor(seer.countyeducation),
                                                 countyincome.even=as.factor(seer.countyincome),countyincome=as.factor(seer.countyincome),residence2=as.factor(seer.residence),binarysurvival=rep(1,length(post.aca)))

## Mark observations removed from stage at diagnosis analyses (i.e. differentiate missing due to not applicable staging [mark with "9"] and simply missing [mark with "NA"])
data.insurance$stage <- data.insurance$stage.cov <-   data.insuranceno$stage <- data.insuranceno$stage.cov <- as.character(data.insurance$stage)
data.insurance$stage[removes] <- 9
data.insurance$stage <- data.insurance$stage.cov <-   data.insuranceno$stage <- data.insuranceno$stage.cov <- as.factor(data.insurance$stage)
data.insurance$stage2 <- data.insurance$stage2.cov <-   data.insuranceno$stage2 <- data.insuranceno$stage2.cov <- as.character(data.insurance$stage2)
data.insurance$stage2[removes] <- 9
data.insurance$stage2 <- data.insurance$stage2.cov <-   data.insuranceno$stage2 <- data.insuranceno$stage2.cov <- as.factor(data.insurance$stage2)

## Generating additional survival information for multiple imputation. see White & Royston 2009. Imputing missing covariate values for the Cox model. Statistics in Medicine. DOI: 10.1002/sim.3618
data.survival <- data.survivalno <- data.insurance
library(survival)
my.surv <- Surv(data.survival$survmonths, 1-data.survival$alive2)
my.fit  <- summary(survfit(my.surv ~ 1))
h.sort.of <- my.fit$n.event / my.fit$n.risk
H.tilde   <- cumsum(h.sort.of)
#H.tilde   <- c(H.tilde, tail(H.tilde, 1))
data.survival$H.T <- rep(NA,nrow(data.survival))
for(i in 1:length(my.fit$time)){
  data.survival$H.T[which(data.survival$survmonths>=my.fit$time[i])] <- H.tilde[i]
}
data.insurance$H.T <- data.survivalno$H.T <- data.survival$H.T

data.insurance$stage <- data.insurance$stage.cov
data.insurance$earlystage <- ifelse(as.numeric(data.insurance$stage)<=1,1,0)
data.insurance$latestage <- ifelse(as.numeric(data.insurance$stage)==4,1,0)
data.insurance$insurance <- data.insurance$insurance.cov

data.han <- data.noimpute <- droplevels(data.insurance[-removes,])

data.insurance.ncdb <- data.insuranceno.ncdb <- data.insurance
data.han.ncdb <- data.han

######## ANALYSES (without imputation)

## include basic DID calculations
## DID plus robust SE
## then stratify by cancer type

###############
## DID estimates, no multiple imputation (note that the estimates given in paper utilize multiple imputation; see below for imputation procedure)
###############

###### INSURANCE (Figure 2a, Supplemental Table 1)

# More only private insurance or uninsured
data.insurance.ncdb$insurance.resp1 <- ifelse(data.insurance.ncdb$insurance==0,0,1)
data.insurance.seer$insurance.resp1 <- ifelse(data.insurance.seer$insurance==0,0,1)
keep.data.ncdb <- which(data.insurance.ncdb$insurance== 0 | data.insurance.ncdb$insurance== 2)
keep.data.seer <- which(data.insurance.seer$insurance== 0 | data.insurance.seer$insurance== 2)
# Medicaid or uninsured compared to private insurance
data.insurance.seer$insurance.resp2 <- ifelse(data.insurance.seer$insurance==2,1,0)
data.insurance.ncdb$insurance.resp2 <- ifelse(data.insurance.ncdb$insurance==2,1,0)
## Unadjusted estimates
ins1.ncdb.u <- lm(insurance.resp1~post.aca*intervention.group,data=data.insurance.ncdb[keep.data.ncdb,])
ins1.seer.u <- lm(insurance.resp1~post.aca*intervention.group,data=data.insurance.seer[keep.data.seer,])
ins2.ncdb.u <- lm(insurance.resp2~post.aca*intervention.group,data=data.insurance.ncdb)
ins2.seer.u <- lm(insurance.resp2~post.aca*intervention.group,data=data.insurance.seer)
## Adjusted estimates
ins1.ncdb.a <- lm(insurance.resp1~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome,data=data.insurance.ncdb[keep.data.ncdb,])
ins1.seer.a <- lm(insurance.resp1~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome,data=data.insurance.seer[keep.data.seer,])
ins2.ncdb.a <- lm(insurance.resp2~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome,data=data.insurance.ncdb)
ins2.seer.a <- lm(insurance.resp2~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome,data=data.insurance.seer)
## Robust SEs
library(lmtest)
library(sandwich)
results.ins1.ncdb.u <- as.matrix(coeftest(ins1.ncdb.u, vcov.=vcovHC))
results.ins1.seer.u <- as.matrix(coeftest(ins1.seer.u, vcov.=vcovHC))
results.ins2.ncdb.u <- as.matrix(coeftest(ins2.ncdb.u, vcov.=vcovHC))
results.ins2.seer.u <- as.matrix(coeftest(ins2.seer.u, vcov.=vcovHC))
results.ins1.ncdb.a <- as.matrix(coeftest(ins1.ncdb.a, vcov.=vcovHC))
results.ins1.seer.a <- as.matrix(coeftest(ins1.seer.a, vcov.=vcovHC))
results.ins2.ncdb.a <- as.matrix(coeftest(ins2.ncdb.a, vcov.=vcovHC))
results.ins2.seer.a <- as.matrix(coeftest(ins2.seer.a, vcov.=vcovHC))
## the "post.aca*intervention.group" coefficient is the DID (always the last one in the matrix)

## Confidence interval estimation
ci.compute <- function(x,df=NULL){
  if(is.null(df)){
    df <- Inf
  }
  x[1] + c(-1,1)*qt(0.975,df)*x[2]
}
c(results.ins1.ncdb.u[nrow(results.ins1.ncdb.u),c(1,4)],ci.compute(results.ins1.ncdb.u[nrow(results.ins1.ncdb.u),c(1,2)],ins1.ncdb.u$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins1.seer.u[nrow(results.ins1.seer.u),c(1,4)],ci.compute(results.ins1.seer.u[nrow(results.ins1.seer.u),c(1,2)],ins1.seer.u$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins2.ncdb.u[nrow(results.ins2.ncdb.u),c(1,4)],ci.compute(results.ins2.ncdb.u[nrow(results.ins2.ncdb.u),c(1,2)],ins2.ncdb.u$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins2.seer.u[nrow(results.ins2.seer.u),c(1,4)],ci.compute(results.ins2.seer.u[nrow(results.ins2.seer.u),c(1,2)],ins2.seer.u$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins1.ncdb.a[nrow(results.ins1.ncdb.a),c(1,4)],ci.compute(results.ins1.ncdb.a[nrow(results.ins1.ncdb.a),c(1,2)],ins1.ncdb.a$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins1.seer.a[nrow(results.ins1.seer.a),c(1,4)],ci.compute(results.ins1.seer.a[nrow(results.ins1.seer.a),c(1,2)],ins1.seer.a$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins2.ncdb.a[nrow(results.ins2.ncdb.a),c(1,4)],ci.compute(results.ins2.ncdb.a[nrow(results.ins2.ncdb.a),c(1,2)],ins2.ncdb.a$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound
c(results.ins2.seer.a[nrow(results.ins2.seer.a),c(1,4)],ci.compute(results.ins2.seer.a[nrow(results.ins2.seer.a),c(1,2)],ins2.seer.a$df.residual)) # DID estimate, p-value, 95% CI lower bound & upper bound

## the parallel trends assumption can be tested by substituting the full data (in code above; i.e. data.insurance.ncdb or data.insurance.seer) with the following data:
data.insurance.ncdb.assumption <- data.insurance.ncdb[which(data.insurance.ncdb$diagnoseyear==2007 | data.insurance.ncdb$diagnoseyear==2009),]
data.insurance.ncdb.assumption$post.aca <- ifelse(data.insurance.ncdb.assumption$diagnoseyear==2009,1,0) ## a "placebo" ACA-DCP effect
data.insurance.seer.assumption <- data.insurance.seer[which(data.insurance.seer$diagnoseyear==2007 | data.insurance.seer$diagnoseyear==2009),]
data.insurance.seer.assumption$post.aca <- ifelse(data.insurance.seer.assumption$diagnoseyear==2009,1,0) ## a "placebo" ACA-DCP effect

## Estimates stratified by cancer type can be generated by subsetting the data in the lm() lines above
# example for colorectal cancer:
# instead of "data.insurance.seer", use "data.insurance.seer[which(data.insurance.seer$cancertype=="8.6.1 Carcinoma of colon and rectum"),]"
# a for() loop can also easily be implemented

###### STAGE AT DIAGNOSIS (Figure 2b-c, Supplemental Table 2)

## Unadjusted estimates
estage.ncdb.u <- lm(earlystage~post.aca*intervention.group,data=data.han.ncdb)
estage.seer.u <- lm(earlystage~post.aca*intervention.group,data=data.han.seer)
lstage.ncdb.u <- lm(latestage~post.aca*intervention.group,data=data.han.ncdb)
lstage.seer.u <- lm(latestage~post.aca*intervention.group,data=data.han.seer)
## Adjusted estimates
estage.ncdb.a <- lm(earlystage~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome+cancertype,data=data.han.ncdb)
estage.seer.a <- lm(earlystage~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome+cancertype,data=data.han.seer)
lstage.ncdb.a <- lm(latestage~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome+cancertype,data=data.han.ncdb)
lstage.seer.a <- lm(latestage~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome+cancertype,data=data.han.seer)

## Perform robust SE estimation, CI computation, parallel trends assumption testing, and stratified cancer type estimates

####### SURVIVAL (Supplemental Table 3)
data.insurance.ncdb$cens <- 1-data.insurance.ncdb$alive 
data.insurance.seer$cens <- 1-data.insurance.seer$alive ## or use data.insurance.seer$alive2 for cancer specific survival (instead of overall survival)
## Unadjusted estimates
survival.ncdb.u <- coxph(Surv(survmonths,cens)~post.aca*intervention.group,data=data.insurance.ncdb,x=T,robust=T)
survival.seer.u <- coxph(Surv(survmonths,cens)~post.aca*intervention.group,data=data.insurance.seer,x=T,robust=T)
ret.ncdb.u <- summary(survival.ncdb.u)$coefficients
ret.seer.u <- summary(survival.seer.u)$coefficients
c(exp(ret.ncdb.u[nrow(ret.ncdb.u),1]),exp(ci.compute(ret.ncdb.u[nrow(ret.ncdb.u),c(1,(ncol(ret.ncdb.u)-2))])),ret.ncdb.u[nrow(ret.ncdb.u),ncol(ret.ncdb.u)]) # hazard ratio estimate, 95% CI lower & upper bound, p-value
c(exp(ret.seer.u[nrow(ret.seer.u),1]),exp(ci.compute(ret.seer.u[nrow(ret.seer.u),c(1,(ncol(ret.seer.u)-2))])),ret.seer.u[nrow(ret.seer.u),ncol(ret.seer.u)]) # hazard ratio estimate, 95% CI lower & upper bound, p-value
## Adjusted estimates
survival.ncdb.a <- coxph(Surv(survmonths,cens)~post.aca*intervention.group+residence+countyeducation+countyincome+sex+race+strata(cancertype),data=data.insurance.ncdb,x=T,robust=T)
survival.seer.a <- coxph(Surv(survmonths,cens)~post.aca*intervention.group+residence+countyeducation+countyincome+sex+marital+race+strata(cancertype),data=data.insurance.seer,x=T,robust=T)
ret.ncdb.a <- summary(survival.ncdb.a)$coefficients
ret.seer.a <- summary(survival.seer.a)$coefficients
c(exp(ret.ncdb.a[nrow(ret.ncdb.a),1]),exp(ci.compute(ret.ncdb.a[nrow(ret.ncdb.a),c(1,(ncol(ret.ncdb.a)-2))])),ret.ncdb.a[nrow(ret.ncdb.a),ncol(ret.ncdb.a)]) # hazard ratio estimate, 95% CI lower & upper bound, p-value
c(exp(ret.seer.a[nrow(ret.seer.a),1]),exp(ci.compute(ret.seer.a[nrow(ret.seer.a),c(1,(ncol(ret.seer.a)-2))])),ret.seer.a[nrow(ret.seer.a),ncol(ret.seer.a)]) # hazard ratio estimate, 95% CI lower & upper bound, p-value


###########################################################
#### Multiple Imputation
## The values presented in the paper (figures and tables) are based on MI rather than complete cases (like shown above)
###########################################################

## Obtaining imputations (18 imputations-the approximate percent of observations with 1 or more missing values)
library(Hmisc)
## For the stage at diagnosis analyses
set.seed(1) ## these seeds were used in the analyses and can be used to exactly reproduce the results
impute_arg.seer <- aregImpute(~ marital + countyeducation + countyincome + residence + age + sex + diagnoseyear+  cancertype + race.old + hispanic.old + stage.cov + insurance.cov + alive2 + survmonths, data = data.han.seer, n.impute = 18, burnin=3,boot.method="approximate bayesian") 
set.seed(1)
impute_arg.ncdb <- aregImpute(~ (countyeducation + countyincome + residence + sex + race.old + hispanic.old) * post.aca*intervention.group + cancertype + stage.cov + insurance.cov + alive + H.T + survmonths, data = data.han.ncdb, n.impute = 18, burnin=3,boot.method="approximate bayesian") 
## For the insurance analyses
set.seed(1)
impute_arg2.seer <- aregImpute(~ marital + countyeducation + countyincome + residence + age + sex + diagnoseyear+  cancertype + race.old + hispanic.old + stage.cov + insurance.cov + alive2 + survmonths, data = data.insurance.seer, n.impute = 18, burnin=3,boot.method="approximate bayesian") # 
set.seed(1)
impute_arg2.ncdb <- aregImpute(~ (countyeducation + countyincome + residence + sex + race.old + hispanic.old) * post.aca*intervention.group + cancertype + stage.cov + insurance.cov + alive + H.T + survmonths, data = data.insurance.ncdb, n.impute = 18, burnin=3,boot.method="approximate bayesian") 

##### Analyses

## For SEER

results.seer <- list() # an object to store the DID results from the regression models
for(i in 1:18){
  #### Setting up data
  ## re-setting the data
  data.insurance.seer <- data.insuranceno.seer
  
  data.insurance.seer$stage <- data.insurance.seer$stage.cov <-   data.insuranceno.seer$stage <- data.insuranceno.seer$stage.cov <- as.character(data.insurance.seer$stage)
  data.insurance.seer$stage[removes.seer] <- 9
  data.insurance.seer$stage <- data.insurance.seer$stage.cov <-   data.insuranceno.seer$stage <- data.insuranceno.seer$stage.cov <- as.factor(data.insurance.seer$stage)
  
  data.han.seer <- data.noimpute <- droplevels(data.insurance.seer[-removes.seer,])
  
  ## putting in the values for the stage at diagnoses analyses
  data.han.seer$marital[which(is.na(data.noimpute$marital))] <- impute_arg$imputed$marital[,i]-1
  data.han.seer$countyeducation[which(is.na(data.noimpute$countyeducation))]  <- impute_arg$imputed$countyeducation[,i]-1
  data.han.seer$countyincome[which(is.na(data.noimpute$countyincome))]  <- impute_arg$imputed$countyincome[,i]-1
  data.han.seer$residence[which(is.na(data.noimpute$residence))]  <- impute_arg$imputed$residence[,i]-1
  data.han.seer$race.old[which(is.na(data.noimpute$race.old))]  <- impute_arg$imputed$race.old[,i]-1
  data.han.seer$stage.cov[which(is.na(data.noimpute$stage.cov))] <- impute_arg$imputed$stage.cov[,i]
  data.han.seer$insurance.cov[which(is.na(data.noimpute$insurance.cov))] <- impute_arg$imputed$insurance.cov[,i]-1
  data.han.seer$survmonths[which(is.na(data.noimpute$survmonths))] <- impute_arg$imputed$survmonths[,i]
  data.han.seer$alive2[which(is.na(data.noimpute$alive2))] <- impute_arg$imputed$alive2[,i]
  raceethnicity.han <- 1+ifelse(data.han.seer$race.old==0 & data.han.seer$hispanic.old==0,-1,0) + ifelse(data.han.seer$race.old==1 & data.han.seer$hispanic.old==0,1,0) + ifelse(data.han.seer$hispanic.old==1,2,0) # 0=non-hispanic white, 1=non-hispanic other, 2=non-hispanic black. 3=hispanic 
  data.han.seer$race <- as.factor(raceethnicity.han)
  data.han.seer$stage <- data.han.seer$stage.cov
  data.han.seer$earlystage <- ifelse(as.numeric(data.han.seer$stage)<=1,1,0)
  data.han.seer$latestage <- ifelse(as.numeric(data.han.seer$stage)==4,1,0)
  data.han.seer$insurance <- data.han.seer$insurance.cov
  
  ## put in values for insurance analyses
  data.insurance.seer$marital[which(is.na(data.insuranceno.seer$marital))] <- impute_arg2$imputed$marital[,i]-1
  data.insurance.seer$countyeducation[which(is.na(data.insuranceno.seer$countyeducation))]  <- impute_arg2$imputed$countyeducation[,i]-1
  data.insurance.seer$countyincome[which(is.na(data.insuranceno.seer$countyincome))]  <- impute_arg2$imputed$countyincome[,i]-1
  data.insurance.seer$residence[which(is.na(data.insuranceno.seer$residence))]  <- impute_arg2$imputed$residence[,i]-1
  data.insurance.seer$race.old[which(is.na(data.insuranceno.seer$race.old))]  <- impute_arg2$imputed$race.old[,i]-1
  data.insurance.seer$stage.cov[which(is.na(data.insuranceno.seer$stage.cov))] <- impute_arg2$imputed$stage.cov[,i]
  data.insurance.seer$insurance.cov[which(is.na(data.insuranceno.seer$insurance.cov))] <- impute_arg2$imputed$insurance.cov[,i]-1
  data.insurance.seer$survmonths[which(is.na(data.insuranceno.seer$survmonths))] <- impute_arg2$imputed$survmonths[,i]
  data.insurance.seer$alive2[which(is.na(data.insuranceno.seer$alive2))] <- impute_arg2$imputed$alive2[,i]
  raceethnicity.han <- 1+ifelse(data.insurance.seer$race.old==0 & data.insurance.seer$hispanic.old==0,-1,0) + ifelse(data.insurance.seer$race.old==1 & data.insurance.seer$hispanic.old==0,1,0) + ifelse(data.insurance.seer$hispanic.old==1,2,0) # 0=non-hispanic white, 1=non-hispanic other, 2=non-hispanic black. 3=hispanic 
  data.insurance.seer$race <- as.factor(raceethnicity.han)
  data.insurance.seer$stage <- data.insurance.seer$stage.cov
  data.insurance.seer$earlystage <- ifelse(as.numeric(data.insurance.seer$stage)<=1,1,0)
  data.insurance.seer$latestage <- ifelse(as.numeric(data.insurance.seer$stage)==4,1,0)
  data.insurance.seer$insurance <- data.insurance.seer$insurance.cov
  
  ### Analyses (only adjusted insurance, combining Medicaid and uninsured, is shown for simplicity; but any/all of the above analyses can be included here)
  mat <- matrix(NA,nrow=6,ncol=4)
  colnames(mat) <- c("estimate","SE","df","p-value")
  rownames(mat) <- c("Insurance1","Insurance2","Early Stage","Late Stage","All Cause Survival","Cancer Survival")
  ## Insurance
  data.insurance.seer$insurance.resp1 <- ifelse(data.insurance.seer$insurance==0,0,1)
  keep.data.seer <- which(data.insurance.seer$insurance== 0 | data.insurance.seer$insurance== 2)
  ins1.seer.a <- lm(insurance.resp1~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome,data=data.insurance.seer) # regression model estimation
  results.ins1.seer.a <- as.matrix(coeftest(ins1.seer.a, vcov.=vcovHC)) ## robust SE estimation
  mat[1,c(1,2,4)] <- results.ins1.seer.a[nrow(results.ins1.seer.a),c(1,2,4)]
  mat[1,3] <- ins1.seer.a$df.residual
  
  data.insurance.seer$insurance.resp2 <- ifelse(data.insurance.seer$insurance==2,1,0)
  ins2.seer.a <- lm(insurance.resp2~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome,data=data.insurance.seer) # regression model estimation
  results.ins2.seer.a <- as.matrix(coeftest(ins2.seer.a, vcov.=vcovHC)) ## robust SE estimation
  mat[2,c(1,2,4)] <- results.ins2.seer.a[nrow(results.ins2.seer.a),c(1,2,4)]
  mat[2,3] <- ins2.seer.a$df.residual
  
  ## Stage at diagnosis
  estage.seer.a <- lm(earlystage~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome+cancertype,data=data.han.seer) # regression model estimation
  results.estage.seer.a <- as.matrix(coeftest(estage.seer.a, vcov.=vcovHC)) ## robust SE estimation
  mat[3,c(1,2,4)] <- results.estage.seer.a[nrow(results.estage.seer.a),c(1,2,4)]
  mat[3,3] <- estage.seer.a$df.residual
  
  lstage.seer.a <- lm(latestage~post.aca*intervention.group+race+sex+residence+marital+countyeducation+countyincome+cancertype,data=data.han.seer) # regression model estimation
  results.lstage.seer.a <- as.matrix(coeftest(lstage.seer.a, vcov.=vcovHC)) ## robust SE estimation
  mat[4,c(1,2,4)] <- results.lstage.seer.a[nrow(results.lstage.seer.a),c(1,2,4)]
  mat[4,3] <- lstage.seer.a$df.residual
  
  # all cause survival
  data.insurance.seer$cens <- 1-data.insurance.seer$alive
  survival.seer.a <- coxph(Surv(survmonths,cens)~post.aca*intervention.group+residence+countyeducation+countyincome+sex+marital+race+strata(cancertype),data=data.insurance.seer,x=T,robust=T)
  ret.seer.a <- summary(survival.seer.a)$coefficients
  mat[5,] <- c(ret.seer.a[nrow(ret.seer.a),1],ret.seer.a[nrow(ret.seer.a),(ncol(ret.seer.a)-2)],Inf,ret.seer.a[nrow(ret.seer.a),ncol(ret.seer.a)]) ## Note that it's not exponentiated yet (we exponentiate after combining the imputed results)
  
  # cancer survival
  data.insurance.seer$cens <- 1-data.insurance.seer$alive2
  survival.seer.a <- coxph(Surv(survmonths,cens)~post.aca*intervention.group+residence+countyeducation+countyincome+sex+marital+race+strata(cancertype),data=data.insurance.seer,x=T,robust=T)
  ret.seer.a <- summary(survival.seer.a)$coefficients
  mat[6,] <- c(ret.seer.a[nrow(ret.seer.a),1],ret.seer.a[nrow(ret.seer.a),(ncol(ret.seer.a)-2)],Inf,ret.seer.a[nrow(ret.seer.a),ncol(ret.seer.a)]) ## Note that it's not exponentiated yet (we exponentiate after combining the imputed results)
  
  #### Additional analyses/models, like for parallel trends assumption testing, or analyses stratified by cancer type, could also be added here
  
  ## store model results in the list object
  results.seer[[i]] <- mat
}



## For NCDB

results.ncdb <- list()

for(i in 1:18){
  ## re-setting the data (since some of it was misspecified with the sensitivity analyses)
  data.insurance.ncdb <- data.insuranceno.ncdb
  
  data.insurance.ncdb$stage <- data.insurance.ncdb$stage.cov <-   data.insuranceno.ncdb$stage <- data.insuranceno.ncdb$stage.cov <- as.character(data.insurance.ncdb$stage)
  data.insurance.ncdb$stage[removes.ncdb] <- 9
  data.insurance.ncdb$stage <- data.insurance.ncdb$stage.cov <-   data.insuranceno.ncdb$stage <- data.insuranceno.ncdb$stage.cov <- as.factor(data.insurance.ncdb$stage)
  data.insurance.ncdb$stage2 <- data.insurance.ncdb$stage2.cov <-   data.insuranceno.ncdb$stage2 <- data.insuranceno.ncdb$stage2.cov <- as.character(data.insurance.ncdb$stage2)
  data.insurance.ncdb$stage2[removes.ncdb] <- 9
  data.insurance.ncdb$stage2 <- data.insurance.ncdb$stage2.cov <-   data.insuranceno.ncdb$stage2 <- data.insuranceno.ncdb$stage2.cov <- as.factor(data.insurance.ncdb$stage2)
  
  data.han.ncdb <- data.noimpute <- droplevels(data.insurance.ncdb[-removes.ncdb,])
  
  ## putting in the values for the stage at diagnosis analyses
  data.han.ncdb$countyeducation[which(is.na(data.noimpute$countyeducation))]  <- impute_arg$imputed$countyeducation[,i]
  data.han.ncdb$countyincome[which(is.na(data.noimpute$countyincome))]  <- impute_arg$imputed$countyincome[,i]
  data.han.ncdb$residence[which(is.na(data.noimpute$residence))]  <- impute_arg$imputed$residence[,i]-1
  data.han.ncdb$race.old[which(is.na(data.noimpute$race.old))]  <- impute_arg$imputed$race.old[,i]-1
  data.han.ncdb$hispanic.old[which(is.na(data.noimpute$hispanic.old))]  <- impute_arg$imputed$hispanic.old[,i]-1
  data.han.ncdb$stage.cov[which(is.na(data.noimpute$stage.cov))] <- impute_arg$imputed$stage.cov[,i]
  data.han.ncdb$insurance.cov[which(is.na(data.noimpute$insurance.cov))] <- impute_arg$imputed$insurance.cov[,i]-1
  data.han.ncdb$survmonths[which(is.na(data.noimpute$survmonths))] <- impute_arg$imputed$survmonths[,i]
  data.han.ncdb$alive2[which(is.na(data.noimpute$alive2))] <- impute_arg$imputed$alive2[,i]
  raceethnicity.han <- 1+ifelse(data.han.ncdb$race.old==0 & data.han.ncdb$hispanic.old==0,-1,0) + ifelse(data.han.ncdb$race.old==1 & data.han.ncdb$hispanic.old==0,1,0) + ifelse(data.han.ncdb$hispanic.old==1,2,0) # 0=non-hispanic white, 1=non-hispanic other, 2=non-hispanic black. 3=hispanic 
  data.han.ncdb$race <- as.factor(raceethnicity.han)
  data.han.ncdb$stage <- data.han.ncdb$stage.cov
  data.han.ncdb$earlystage <- ifelse(as.numeric(data.han.ncdb$stage)<=1,1,0)
  data.han.ncdb$latestage <- ifelse(as.numeric(data.han.ncdb$stage)==4,1,0)
  data.han.ncdb$insurance <- data.han.ncdb$insurance.cov
  
  ## put in values for the insurance and survival analyses
  data.insurance.ncdb$countyeducation[which(is.na(data.insuranceno.ncdb$countyeducation))]  <- impute_arg2$imputed$countyeducation[,i]
  data.insurance.ncdb$countyincome[which(is.na(data.insuranceno.ncdb$countyincome))]  <- impute_arg2$imputed$countyincome[,i]
  data.insurance.ncdb$residence[which(is.na(data.insuranceno.ncdb$residence))]  <- impute_arg2$imputed$residence[,i]-1
  data.insurance.ncdb$race.old[which(is.na(data.insuranceno.ncdb$race.old))]  <- impute_arg2$imputed$race.old[,i]-1
  data.insurance.ncdb$hispanic.old[which(is.na(data.insuranceno.ncdb$hispanic.old))]  <- impute_arg2$imputed$hispanic.old[,i]-1
  data.insurance.ncdb$insurance.cov[which(is.na(data.insuranceno.ncdb$insurance.cov))] <- impute_arg2$imputed$insurance.cov[,i]-1
  data.insurance.ncdb$survmonths[which(is.na(data.insuranceno.ncdb$survmonths))] <- impute_arg2$imputed$survmonths[,i]
  data.insurance.ncdb$alive2[which(is.na(data.insuranceno.ncdb$alive2))] <- impute_arg2$imputed$alive2[,i]
  raceethnicity.han <- 1+ifelse(data.insurance.ncdb$race.old==0 & data.insurance.ncdb$hispanic.old==0,-1,0) + ifelse(data.insurance.ncdb$race.old==1 & data.insurance.ncdb$hispanic.old==0,1,0) + ifelse(data.insurance.ncdb$hispanic.old==1,2,0) # 0=non-hispanic white, 1=non-hispanic other, 2=non-hispanic black. 3=hispanic 
  data.insurance.ncdb$race <- as.factor(raceethnicity.han)
  data.insurance.ncdb$stage <- data.insurance.ncdb$stage.cov
  data.insurance.ncdb$earlystage <- ifelse(as.numeric(data.insurance.ncdb$stage)<=1,1,0)
  data.insurance.ncdb$latestage <- ifelse(as.numeric(data.insurance.ncdb$stage)==4,1,0)
  data.insurance.ncdb$insurance <- data.insurance.ncdb$insurance.cov
  
  ### Analyses (only adjusted insurance, combining Medicaid and uninsured, is shown for simplicity; but any/all of the above analyses can be included here)
  mat <- matrix(NA,nrow=6,ncol=4)
  colnames(mat) <- c("estimate","SE","df","p-value")
  rownames(mat) <- c("Insurance1","Insurance2","Early Stage","Late Stage","All Cause Survival","Cancer Survival")
  ## Insurance
  data.insurance.ncdb$insurance.resp1 <- ifelse(data.insurance.ncdb$insurance==0,0,1)
  keep.data.ncdb <- which(data.insurance.ncdb$insurance== 0 | data.insurance.ncdb$insurance== 2)
  ins1.ncdb.a <- lm(insurance.resp1~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome,data=data.insurance.ncdb) # regression model estimation
  results.ins1.ncdb.a <- as.matrix(coeftest(ins1.ncdb.a, vcov.=vcovHC)) ## robust SE estimation
  mat[1,c(1,2,4)] <- results.ins1.ncdb.a[nrow(results.ins1.ncdb.a),c(1,2,4)]
  mat[1,3] <- ins1.ncdb.a$df.residual
  
  data.insurance.ncdb$insurance.resp2 <- ifelse(data.insurance.ncdb$insurance==2,1,0)
  ins2.ncdb.a <- lm(insurance.resp2~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome,data=data.insurance.ncdb) # regression model estimation
  results.ins2.ncdb.a <- as.matrix(coeftest(ins2.ncdb.a, vcov.=vcovHC)) ## robust SE estimation
  mat[2,c(1,2,4)] <- results.ins2.ncdb.a[nrow(results.ins2.ncdb.a),c(1,2,4)]
  mat[2,3] <- ins2.ncdb.a$df.residual
  
  ## Stage at diagnosis
  estage.ncdb.a <- lm(earlystage~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome+cancertype,data=data.han.ncdb) # regression model estimation
  results.estage.ncdb.a <- as.matrix(coeftest(estage.ncdb.a, vcov.=vcovHC)) ## robust SE estimation
  mat[3,c(1,2,4)] <- results.estage.ncdb.a[nrow(results.estage.ncdb.a),c(1,2,4)]
  mat[3,3] <- estage.ncdb.a$df.residual
  
  lstage.ncdb.a <- lm(latestage~post.aca*intervention.group+race+sex+residence+countyeducation+countyincome+cancertype,data=data.han.ncdb) # regression model estimation
  results.lstage.ncdb.a <- as.matrix(coeftest(lstage.ncdb.a, vcov.=vcovHC)) ## robust SE estimation
  mat[4,c(1,2,4)] <- results.lstage.ncdb.a[nrow(results.lstage.ncdb.a),c(1,2,4)]
  mat[4,3] <- lstage.ncdb.a$df.residual
  
  # all cause survival
  data.insurance.ncdb$cens <- 1-data.insurance.ncdb$alive
  survival.ncdb.a <- coxph(Surv(survmonths,cens)~post.aca*intervention.group+residence+countyeducation+countyincome+sex+race+strata(cancertype),data=data.insurance.ncdb,x=T,robust=T)
  ret.ncdb.a <- summary(survival.ncdb.a)$coefficients
  mat[5,] <- c(ret.ncdb.a[nrow(ret.ncdb.a),1],ret.ncdb.a[nrow(ret.ncdb.a),(ncol(ret.ncdb.a)-2)],Inf,ret.ncdb.a[nrow(ret.ncdb.a),ncol(ret.ncdb.a)]) ## Note that it's not exponentiated yet (we exponentiate after combining the imputed results)
  mat[6,] <- c(ret.ncdb.a[nrow(ret.ncdb.a),1],ret.ncdb.a[nrow(ret.ncdb.a),(ncol(ret.ncdb.a)-2)],Inf,ret.ncdb.a[nrow(ret.ncdb.a),ncol(ret.ncdb.a)]) ## Just a placeholder to match the SEER output (no cancer survival in NCDB)
  
  #### Additional analyses/models, like for parallel trends assumption testing, or analyses stratified by cancer type, could also be added here
  
  ## store model results in the list object
  results.ncdb[[i]] <- mat
}

#### Combining the imputed results
## A function to recreate the confidence interval and p-value from the combined (across imputations) estimate and SE
summary3 <- function(mat,row=seq(1,nrow(mat))){
    j <- 0
    matter <- matrix(NA,ncol=4,nrow=length(row))
    for(i in row){
      j <- j+1
      vv <- mat[i,2]
      coef <- mat[i,1]
      df <- mat[i,3]
      z <- coef/vv
      p <- pt(z,df)
      if(is.na(p)) p <- p
      else if(p>0.5) p <- 1-p
      p <- 2*p
      ci <- coef + c(-1,1)*qt(0.975,df)*vv
      matter[j,] <- c(coef,ci,p)
    }
    rownames(matter) <- rownames(mat)
    return(matter)
}
## A function to combine the estimates and SEs across imputations
adjust.multiple.imputation.simple <- function(results.list){
  ref <- results.list[[1]]
  coefmat <- semat <- dfsmat <- NULL
  for(k in 1:length(results.list)){
    ses <- results.list[[k]][,c(1,2,3)]
    coefmat <- rbind(coefmat,as.numeric(ses[,1]))
    semat <- rbind(semat,as.numeric(ses[,2]))
    dfsmat <- rbind(dfsmat,as.numeric(ses[,3]))
  }
  library(Amelia)
  newvals <- mi.meld(coefmat,semat)
  results <- cbind(t(newvals$q.mi),t(newvals$se.mi),dfsmat[1,])
  rownames(results) <- rownames(ref)
  newests <- summary3(results)
  return(newests)
}

mi.combined.ncdb <- adjust.multiple.imputation.simple(results.ncdb)
mi.combined.seer <- adjust.multiple.imputation.simple(results.seer)
mi.combined.ncdb[c(5,6),1:3] <- exp(mi.combined.ncdb[c(5,6),1:3]) ## exponentiate survival results to generate hazard ratio and 95% CI
mi.combined.seer[c(5,6),1:3] <- exp(mi.combined.seer[c(5,6),1:3]) ## exponentiate survival results to generate hazard ratio and 95% CI
mi.combined.ncdb
mi.combined.seer


#################################################
##### Figure code - generates PDFs of the plots
#################################################

## Note that we did NOT exclude observations from 2010 when making the insurance and stage figures (but the data prepared above does not include any 2010 observations)
## the survival figures still exclude 2010

#### Functions to facilitate figure generation
meany <- function(x,y){
  means <- numeric(length(table(y)))
  for(i in 1:length(means)){
    means[i] <- mean(x[which(y==names(table(y))[i])],na.rm=T) # calculate the mean (between 0 and 1) for each category (e.g. 19-25y pre-ACA, 27-29 post-ACA)
  }
  names(means) <- names(table(y))
  return(means)
}

ny <- function(x,y){
  means <- numeric(length(table(y)))
  for(i in 1:length(means)){
    means[i] <- sum(y==names(table(y))[i])
    means[i] <- means[i] - sum(is.na(x[which(y==names(table(y))[i])]))
  }
  names(means) <- names(table(y))
  return(means)
}

## 'internal' function for estimating the lines/trends and CIs and putting them on a plot
binary.plotter <- function(x,y,lt=1,cols="black",pc=20,add=F,min=0,max=1,namer="",noaxis=F,errwidth=0.02,offset=0,cols2=NULL,xlims=NULL,lw=1,cx=1){
  library(Hmisc)
  grpmeans <- meany(x,y)
  ses <- sqrt(grpmeans*(1-grpmeans)/ny(x,y)) ## standard error of a proportion
  sd <- ses*qnorm(0.975) ## the one-sided length of a 95% confidence interval (to be added to and subtracted from estimate)
  grpmeans <- grpmeans*100 # convert proportion to percent
  sd <- sd*100 # change length of confidence interval to match
  d <- data.frame(x=as.numeric(names(table(y)))+offset,y=grpmeans,sd=sd) # a dataframe for easy calling
  if(is.null(xlims)) xlims <- c(min(d$x)-offset,max(d$x)+offset) # function can automatically determine graph limits on x axis
  if(add==F) plot(d$x,d$y,lty=lt,col=cols,xlab="",ylab=paste(namer," (%)",sep=""),type="n",pch=pc,ylim=c(min,max)*100,xaxt='n',xlim=xlims,cex.lab=cx,cex.axis=cx)
  else points(d$x,d$y,lty=lt,col=cols,type="n",pch=pc)
  with(data=d,expr = errbar(x,y, y+sd, y-sd, add=T, pch=pc, cap=errwidth,errbar.col=cols2,col=cols,lty=lt))
  lines(d$x,d$y,lty=lt,col=cols,type="l",pch=pc,lwd=lw)
}

### Function to generate the plots (the way to specify all the options)
trend.plot <- function(data,covariate,overlay=F,lims=NULL,pc=20,location="bottomleft",ylab=NULL,offset=0,offset.extra=0,cols2=NULL,xlims=NULL,lw=1,cx=1,diagnoseyear2=F,legendit=T){
  if(diagnoseyear2) data$diagnoseyear <- data$diagnoseyear2
  colnames(data)[which(colnames(data)==covariate)] <- "pgy"
  a <- 1
  if(is.null(cols2)){
    a <- 0
    cols2 <- c("blue4","green")
  }
  else if(length(cols2)==1){
    a <- 1
    cols2 <- c(cols2,cols2)
  }
  if(legendit==F) a <- 1
  if(covariate=="insurance") data$pgy <- (as.numeric(data$pgy)-1)/2
  # a little bit of redundancy here (these are calculated in another function too), but used in this 'wrapper' function to help determine automatic ranges for y axis
  vals1 <- meany(data$pgy[which(data$intervention.group==0)],data$diagnoseyear[which(data$intervention.group==0)])
  vals2 <- meany(data$pgy[which(data$intervention.group==1)],data$diagnoseyear[which(data$intervention.group==1)])
  sdvals1 <- sqrt(vals1*(1-vals1)/ny(data$pgy[which(data$intervention.group==0)],data$diagnoseyear[which(data$intervention.group==0)]))*1.96
  sdvals2 <- sqrt(vals2*(1-vals2)/ny(data$pgy[which(data$intervention.group==1)],data$diagnoseyear[which(data$intervention.group==1)]))*1.96
  vals <- c(vals1+sdvals1,vals1-sdvals1,vals2+sdvals2,vals2-sdvals2) 
  if(is.null(lims)==F) vals=lims # if lims (y limits for range of y axis) is specified it will override the automatic selection
  binary.plotter(data$pgy[which(data$intervention.group==0)],data$diagnoseyear[which(data$intervention.group==0)],cols="green",min=min(vals),max=max(vals),namer=ylab,add=overlay,pc=pc,offset=offset.extra+offset,cols2=cols2[2],xlims=xlims,lw=lw,cx=cx)
  binary.plotter(data$pgy[which(data$intervention.group==1)],data$diagnoseyear[which(data$intervention.group==1)],add=T,cols="blue4",pc=pc,offset=offset.extra-offset,cols2=cols2[1],xlims=xlims,lw=lw,cx=cx)
  axis(1,at=c(2007,2008,2009,2010,2011,2012,2013,2014),labels=c(2007,2008,2009,2010,2011,2012,2013,2014),cex.axis=cx)
  if(a==0) legend(location,c("19-25","27-29"),col=c("blue4","green"),pch=20,lwd=lw,lty=0,cex=cx)
}

## Similar to above, but assuming overlaying results from a separate dataset (so adding two additional trends--intervention & control)
trend.plot2 <- function(data,covariate,overlay=T,lims=NULL,pc=20,location="bottomleft",ylab=NULL,offset=0,offset.extra=0,cols2=NULL,xlims=NULL,lw=1,cx=1,diagnoseyear2=F){
  if(diagnoseyear2) data$diagnoseyear <- data$diagnoseyear2
  colnames(data)[which(colnames(data)==covariate)] <- "pgy"
  overlay <- T
  a <- 1
  if(is.null(cols2)){
    a <- 0
    cols2 <- c("blue4","green")
  }
  else if(length(cols2)==1){
    a <- 1
    cols2 <- c(cols2,cols2)
  }
  if(covariate=="insurance") data$pgy <- (as.numeric(data$pgy)-1)/2
  # a little bit of redundancy here (these are calculated in another function too), but used in this 'wrapper' function to help determine automatic ranges for y axis
  vals1 <- meany(data$pgy[which(data$intervention.group==0)],data$diagnoseyear[which(data$intervention.group==0)])
  vals2 <- meany(data$pgy[which(data$intervention.group==1)],data$diagnoseyear[which(data$intervention.group==1)])
  sdvals1 <- sqrt(vals1*(1-vals1)/ny(data$pgy[which(data$intervention.group==0)],data$diagnoseyear[which(data$intervention.group==0)]))*1.96
  sdvals2 <- sqrt(vals2*(1-vals2)/ny(data$pgy[which(data$intervention.group==1)],data$diagnoseyear[which(data$intervention.group==1)]))*1.96
  vals <- c(vals1+sdvals1,vals1-sdvals1,vals2+sdvals2,vals2-sdvals2) 
  if(is.null(lims)==F) vals=lims # if lims (y limits for range of y axis) is specified it will override the automatic selection
  binary.plotter(data$pgy[which(data$intervention.group==0)],data$diagnoseyear[which(data$intervention.group==0)],lt=2,cols="green",min=min(vals),max=max(vals),namer=ylab,add=overlay,pc=pc,offset=offset.extra+offset,cols2=cols2[2],xlims=xlims,lw=lw,cx=cx)
  binary.plotter(data$pgy[which(data$intervention.group==1)],data$diagnoseyear[which(data$intervention.group==1)],lt=2,add=T,cols="blue4",pc=pc,offset=offset.extra-offset,cols2=cols2[1],xlims=xlims,lw=lw,cx=cx)
  axis(1,at=c(2007,2008,2009,2010,2011,2012,2013,2014),labels=c(2007,2008,2009,2010,2011,2012,2013,2014),cex.axis=cx)
  if(a==0) legend(location,c("19-25","27-29","SEER","NCDB"),col=c("blue4","green","black","black"),pch=c(20,20,20,1),lwd=lw,lty=c(0,0,1,2),cex=cx)
}

# a function that automatically names and makes pdfs using the previous functions
trend.plot.pdf2 <- function(data,data.ncdb,covariate,overlay=F,lims=NULL,pc=16,location="bottomleft",ylab=NULL,offset=0.0,offset.extra=0,cols2=NULL,xlims=NULL,lw=1,cx=1,diagnoseyear2=F,cancertype=NULL){
  namer <- paste(covariate,"-",cancertype,"-",paste(lims,collapse="-"),"-trendplot.pdf",sep="")
  pdf(namer)
  par(mar=c(5,6,4,1))
  if(is.null(cancertype)) trend.plot(data,covariate,overlay,lims,pc,location,ylab,offset,offset.extra,cols2,xlims,lw,cx,diagnoseyear2,legendit=F)
  else trend.plot(data[which(data$cancertype==cancertype),],covariate,overlay,lims,pc,location,ylab,offset,offset.extra,cols2,xlims,lw,cx,diagnoseyear2,legendit=F)
  if(is.null(cancertype)) trend.plot2(data.ncdb,covariate,overlay,lims,1,location,ylab,offset,offset.extra,cols2,xlims,lw,cx,diagnoseyear2)
  else trend.plot2(data.ncdb[which(data.ncdb$cancertype==cancertype),],covariate,overlay,lims,1,location,ylab,offset,offset.extra,cols2,xlims,lw,cx,diagnoseyear2)
  dev.off()
}

## Create Kaplan-Meier curves for the intervention and control groups, before and after the DCP
survival.plot <- function(data.kim2,allcause=F,ylims=c(0.8,1),xlims=c(0,100),lw=3,cx=1.8,legendit=T){
  library(survival)
  if(allcause==F) data.kim2$cens <- 1-data.kim2$alive2 # cause-specific survival (i.e. died from cancer, did not die from cancer)
  else data.kim2$cens <- 1-data.kim2$alive # all-cause survival (i.e. died, not dead)
  pre.control <- which(data.kim2$post.aca==0 & data.kim2$intervention.group==0)
  post.control <- which(data.kim2$post.aca==1 & data.kim2$intervention.group==0)
  pre.intervention <- which(data.kim2$post.aca==0 & data.kim2$intervention.group==1)
  post.intervention <- which(data.kim2$post.aca==1 & data.kim2$intervention.group==1)
  
  ### Creating different survival objects/curves for each pre-/post-ACA and age group
  g1.KM <- survfit(Surv(survmonths,cens)~1,type="kaplan-meier",data=data.kim2[pre.intervention,])
  g2.KM <- survfit(Surv(survmonths,cens)~1,type="kaplan-meier",data=data.kim2[pre.control,])
  g3.KM <- survfit(Surv(survmonths,cens)~1,type="kaplan-meier",data=data.kim2[post.intervention,])
  g4.KM <- survfit(Surv(survmonths,cens)~1,type="kaplan-meier",data=data.kim2[post.control,])
  
  ymax <- 1
  ymin <- 0
  g1.KM.time <- c(0,g1.KM$time)
  g2.KM.time <- c(0,g2.KM$time)
  g3.KM.time <- c(0,g3.KM$time)
  g4.KM.time <- c(0,g4.KM$time)
  g1.KM.surv <- c(1,g1.KM$surv)
  g2.KM.surv <- c(1,g2.KM$surv)
  g3.KM.surv <- c(1,g3.KM$surv)
  g4.KM.surv <- c(1,g4.KM$surv)
  
  plot(g1.KM.time,g1.KM.surv,ylim=ylims,xlim=xlims,type="s",lty=6, xlab="Time (months)", ylab="Survival",main="",col="blue4",lwd=lw,cex.lab=cx,cex.axis=cx)
  #plot(g1.KM.time,g1.KM.surv,ylim=c(.94,1),xlim=c(0,100),type="s",lty=1, xlab="Time (months)", ylab="Survival",main="",col="blue4",lwd=3,cex.lab=2,cex.axis=2)
  lines(g2.KM.time,g2.KM.surv,lty=6,type="s",col="green",lwd=lw)
  lines(g3.KM.time,g3.KM.surv,lty=1,type="s",col="blue4",lwd=lw)
  lines(g4.KM.time,g4.KM.surv,lty=1,type="s",col="green",lwd=lw)
  if(legendit) legend("topright",c("19-25","27-29","Pre-ACA","Post-ACA"),lty=c(1,1,6,1),col=c("blue4","green","black","black"),bty="o",cex=1.5,lwd=lw)
  
}

## generate PDF using previous function
survival.plot.pdf <- function(data.kim2,allcause=F,ylims=c(0.8,1),xlims=c(0,100),lw=3,cx=1.8,cancertype=NULL){
  if(allcause) namer <- paste("survival-allcause","-",cancertype,"-",paste(ylims,collapse="-"),"-plot.pdf",sep="")
  else namer <- paste("survival-specific","-",cancertype,"-",paste(ylims,collapse="-"),"-plot.pdf",sep="")
  pdf(namer)
  par(mar=c(5,6,4,1))
  if(is.null(cancertype)==F) survival.plot(data.kim2[which(data.kim2$cancertype==cancertype),],allcause,ylims,xlims,lw,cx)
  else survival.plot(data.kim2,allcause,ylims,xlims,lw,cx)
  dev.off()
}

## note that the filenames for the PDF files generated for similar analyses for NCDB and SEER may be the same (so a file may be overwritten)

## Plotting insurance (Figures 1a-b)
data.insurance.seer$insurance.2 <- ifelse(data.insurance.seer$insurance==2,1,0)
data.insurance.ncdb$insurance.2 <- ifelse(data.insurance.ncdb$insurance==2,1,0)
data.insurance.seer$insurance <- ifelse(data.insurance.seer$insurance==2,1,0)
data.insurance.ncdb$insurance <- ifelse(data.insurance.ncdb$insurance==2,1,0)
trend.plot.pdf2(data.insurance.seer,data.insurance.ncdb,"insurance.2",ylab="Private Insurance",location="topright",cx=1.8,lw=3,lims=c(0.5,1))
trend.plot.pdf2(data.insurance.seer[which(data.insurance.seer$insurance==0 | data.insurance.seer$insurance==2),],data.insurance.ncdb[which(data.insurance.ncdb$insurance==0 | data.insurance.ncdb$insurance==2),],"insurance",ylab="Private Insurance (Excluding Medicaid)",location="bottomright",cx=1.8,lw=3,lims=c(0.5,1))

## Plotting stage at diagnosis, including some stratified by cancer type (Figures 1c-f)
trend.plot.pdf2(data.han.seer,data.han.ncdb,"earlystage",ylab="Early Stage",location="topright",cx=1.8,lw=3,lims=c(0.5,1))
trend.plot.pdf2(data.han.seer,data.han.ncdb,"latestage",location="topright",ylab="Late Stage",cx=1.8,lw=3,lims=c(0,0.5))
trend.plot.pdf2(data.han.seer,data.han.ncdb,"latestage",location="topright",ylab="Late Stage",cx=1.8,lw=3,lims=c(0,1),cancertype="8.6.1 Carcinoma of colon and rectum")
trend.plot.pdf2(data.han.seer,data.han.ncdb,"earlystage",location="topright",ylab="Early Stage",cx=1.8,lw=3,lims=c(0,1),cancertype="8.6.1 Carcinoma of colon and rectum")

## Plotting survival, including some stratified by cancer type (Figure 3)
survival.plot.pdf(data.insurance.ncdb,ylims=c(0.75,1),xlims = c(0,120),allcause=T) ## overall survival
survival.plot.pdf(data.insurance.seer,ylims=c(0.75,1),xlims = c(0,100)) ## cancer survival
survival.plot.pdf(data.insurance.seer,ylims=c(0.75,1),xlims = c(0,100),allcause=T) ## overall survival
survival.plot.pdf(data.insurance.seer,allcause=T,cancertype="8.6.1 Carcinoma of colon and rectum",ylims=c(0.5,1))
survival.plot.pdf(data.insurance.ncdb,allcause=T,cancertype="8.6.1 Carcinoma of colon and rectum",ylims=c(0.5,1),xlims=c(0,120))

#################################################
##### Code for generating the Table 1 (summarize characteristics of study populations)
#################################################

## Some functions to help
total.maker <- function(mat){
  mat <- cbind(rowSums(mat),mat)
  colnames(mat)[1] <- "Total"
  return(mat)
}

percent.adder <- function(mat){
  coly <- colSums(mat)
  cols <- matrix(rep(coly,nrow(mat)),byrow=T,ncol=ncol(mat))
  mat2 <- mat/cols
  newmat <- paste(as.character(mat)," (",as.character(round(mat2,3)*100),")",sep="")
  newmat <- matrix(newmat,ncol=ncol(mat),nrow=nrow(mat))
  colnames(newmat) <- colnames(mat)
  rownames(newmat) <- rownames(mat)
  return(newmat)
}

doit <- function(mat){
  mat <- total.maker(mat)
  mat <- percent.adder(mat)
  return(mat)
}

data.summarize <- function(data.kim,out.csv=T){
  groups <- ifelse(data.kim$post.aca==0 & data.kim$intervention.group==0,2,0) + 
    ifelse(data.kim$post.aca==1 & data.kim$intervention.group==1,1,0) + 
    ifelse(data.kim$post.aca==1 & data.kim$intervention.group==0,3,0) # 0 = pre, younger. 2=pre, older. 1=post, younger. 3=post, older
  grps <- as.matrix(table(groups))
  est <- table(data.kim$earlystage,groups,useNA="ifany")
  est <- matrix(est,ncol=ncol(est))
  ins <- table(data.kim$insurance,groups,useNA="ifany")
  ins <- matrix(ins,ncol=ncol(ins))
  bins <- table(data.kim$binarysurvival,groups,useNA="ifany")
  bins <- matrix(bins,ncol=ncol(bins))
  lst <- table(data.kim$latestage,groups,useNA="ifany")
  lst <- matrix(lst,ncol=ncol(lst))
  mar <- table(data.kim$marital,groups,useNA="ifany")
  mar <- matrix(mar,ncol=ncol(mar))
  re <- table(data.kim$race,groups,useNA="ifany")
  re <- matrix(re,ncol=ncol(re))
  ce <- table(data.kim$countyeducation,groups,useNA="ifany")
  ce <- matrix(ce,ncol=ncol(ce))
  coi <- table(data.kim$countyincome,groups,useNA="ifany")
  coi <- matrix(coi,ncol=ncol(coi))
  res <- table(data.kim$residence,groups,useNA="ifany")  ## what to do with the unknowns??
  res <- matrix(res,ncol=ncol(res))
  ag <- table(data.kim$age,groups,useNA="ifany")
  ag <- matrix(ag,ncol=ncol(ag))
  se <- table(data.kim$sex,groups,useNA="ifany")
  se <- matrix(se,ncol=ncol(se))
  diy <- table(data.kim$diagnoseyear,groups,useNA="ifany")
  diy <- matrix(diy,ncol=ncol(diy))
  can <- table(data.kim$cancertype,groups,useNA="ifany")
  can <- matrix(can,ncol=ncol(can))
  al <- table(data.kim$alive2,groups,useNA="ifany")
  al <- matrix(al,ncol=ncol(al))
  #table(data.kim$survmonths,groups)
  st <- table(data.kim$stage,groups,useNA="ifany")
  st <- matrix(st,ncol=ncol(st))
  colnames(est) <- colnames(ins) <- colnames(bins) <- colnames(lst) <- colnames(mar) <- colnames(re) <- 
    colnames(ce) <- colnames(coi) <- colnames(res) <- colnames(ag) <- colnames(se) <- colnames(diy) <- 
    colnames(can) <- colnames(al) <- colnames(st) <- rownames(grps) <- c("Pre-ACA, 21-25","Pre-ACA, 27-29", "Post-ACA, 21-25","Post-ACA, 27-29")
  rownames(est) <- c("Stage>1","Stage<=1","NA")
  rownames(ins) <- c("None","Medicaid","Private","NA")
  #rownames(bins) <- c("Dead by 12 mo","Alive at 12 mo","NA")
  rownames(lst) <- c("stage<4","stage=4","NA")
  rownames(mar) <- c("Unmarried","Married","NA")
  rownames(re) <- c("Non-Hispanic White","Non-Hispanic Other","Non-Hispanic Black","Hispanic","NA")
  rownames(ce) <- c("<10.9",">=10.9 & <14.1", ">=14.1 & <20.7", ">20.7","NA")
  rownames(coi) <- c("<38",">=38 & <43", ">=43 & <63", ">63","NA")
  rownames(res) <- c("Non-metro","Metro","NA")
  rownames(ag) <- c("19","20","21","22","23","24","25","27","28","29")
  rownames(se) <- c("Male","Female")
  rownames(diy) <- c("2007","2008","2009","2011","2012","2013","2014")
  rownames(can) <- sort(as.character(unique(data.kim$cancertype)))
  cat(dim(al))
  #rownames(al) <- c("Death from Cancer","Alive, or death from other cause","NA")
  rownames(al) <- c("Death from Cancer","Alive, or death from other cause")
  if(nrow(st)==6) rownames(st) <- c("I","II","III","IV","Excluded","NA")
  else rownames(st) <- c("I","II","III","IV","NA")
  est <- doit(est)
  ins <- doit(ins)
  bins <- doit(bins)
  lst <- doit(lst)
  mar <- doit(mar)
  se <- doit(se)
  re <- doit(re)
  ce <- doit(ce)
  coi <- doit(coi)
  res <- doit(res)
  ag <- doit(ag)
  diy <- doit(diy)
  can <- doit(can)
  al <- doit(al)
  st <- doit(st)
  
  if(out.csv){
    grps <- doit(t(grps))
    grps1 <- est1 <- ins1 <- bins1 <- lst1 <- mar1 <- re1 <- ce1 <- coi1 <- res1 <- se1 <- ag1 <- diy1 <- can1 <- al1 <- st1 <- data.frame(matrix(NA,ncol=5,nrow=1))
    colnames(grps1) <- colnames(est1) <- colnames(ins1) <- colnames(bins1) <- colnames(lst1) <- colnames(mar1) <- colnames(re1) <- colnames(ce1) <- colnames(coi1) <- colnames(res1) <- colnames(ag1) <- colnames(diy1) <- colnames(can1) <- colnames(al1) <- colnames(st1) <- colnames(se1) <-  c("Total","Pre-ACA, 21-25","Pre-ACA, 27-29", "Post-ACA, 21-25","Post-ACA, 27-29")
    rownames(grps1) <- c( "Groups: Data by Intervention and Period")
    rownames(est1) <- c("Early Stage by Groups")
    rownames(ins1) <- c("Insurance by Groups")
    rownames(bins1) <- c("Survival by Groups")
    rownames(lst1) <- c("Late Stage by Groups")
    rownames(mar1) <- c("Marital Status by Groups")
    rownames(se1) <- c("Sex by Groups")
    rownames(re1) <- c("Race/Ethnicity by Groups")
    rownames(ce1) <- c("County Education (% w/o high school) by Groups")
    rownames(coi1) <- c("County income (thousands) by Groups")
    rownames(res1) <- c("Residence by Groups")
    rownames(ag1) <- c("Age by Groups")
    rownames(diy1) <- c("Year Diagnosed by Groups")
    rownames(can1) <- c("Cancer Type by Groups")
    rownames(al1) <- c("Cancer Survival by Groups")
    rownames(st1) <- c("Stage by Groups")
    #matr <- data.frame(rbind(grps1,grps,ins1,ins,est1,est,lst1,lst,st1,st,al1,al,bins1,bins,se1,se,mar1,mar,re1,re,ce1,ce,coi1,coi,res1,res,can1,can,ag1,ag,diy1,diy))
    matr <- data.frame(rbind(grps1,grps,ins1,ins,st1,st,al1,al,se1,se,mar1,mar,re1,re,ce1,ce,coi1,coi,res1,res,can1,can))
    write.csv(matr,"data_summary.csv")
  }
}
data.insurance.ncdb$marital <- rep(c(0,1,NA),length=nrow(data.insurance.ncdb)) # for easy compatibility with code; there isn't marital information in NCDB
data.summarize(data.insurance.ncdb)
data.summarize(data.insurance.seer)


### DID Portions of the tables

ci.compute <- function(x,df=NULL){
  if(is.null(df)){
    df <- Inf
  }
  x[1] + c(-1,1)*qt(0.975,df)*x[2]
  #else x[1] + c(-1,1)*qnorm(0.975)*x[2]
}

ci <- function(x,df=Inf){
  return(c(x[2],df))
}

did <- function(outcome,post.aca,intervention.group,ret="summary",boot=F){
  x1 <- outcome[which(post.aca==1 & intervention.group==1)]
  x0 <- outcome[which(post.aca==0 & intervention.group==1)]
  y1 <- outcome[which(post.aca==1 & intervention.group==0)]
  y0 <- outcome[which(post.aca==0 & intervention.group==0)]
  d1 <-  (mean(x1,na.rm=T)-mean(x0,na.rm=T))*100
  d2 <- (mean(y1,na.rm=T)-mean(y0,na.rm=T))*100
  d1c <-  c(mean(x1,na.rm=T),mean(x0,na.rm=T))*100
  d2c <- c(mean(y1,na.rm=T),mean(y0,na.rm=T))*100
  dd <- d1-d2
  mod <- lm(outcome~post.aca*intervention.group)
  results <- as.matrix(coeftest(mod, vcov.=vcovHC))
  estp <- results[nrow(results),c(1,4)]
  cis <- ci.compute(results[nrow(results),c(1,2)],mod$df.residual)
  estp <- estp*100
  cis <- cis*100
  if(ret=="estimates"){
    return(c(d1,d2,dd))
  }
  else return(c(paste(round(d1,1)," (",round(d1c[2],1),", ",round(d1c[1],1),")",sep=""),paste(round(d2,1)," (",round(d2c[2],1),", ",round(d2c[1],1),")",sep=""),round(dd,1),paste(round(estp[1],2)," (",round(cis[1],2),", ",round(cis[2],2),")",sep="")))
}

did.summarize <- function(data,outcome){
  mat <- matrix(NA,nrow=length(unique(data$cancertype))+1,ncol=4)
  mat[1,] <- did(outcome,data$post.aca,data$intervention.group)
  l <- 1
  for(i in sort(unique(data$cancertype))){
    l <- l+1
    mat[l,] <- did(outcome[which(data$cancertype==i)],data$post.aca[which(data$cancertype==i)],data$intervention.group[which(data$cancertype==i)])
  }
  rownames(mat) <- c("Main",as.character(sort(unique(data$cancertype))))
  colnames(mat) <- c("Intervention","Control","DID","ModDID")
  return(mat)
}

did.summarize.all <- function(data.insurance,data.stage){
  di2 <- did.summarize(data.insurance,ifelse(data.insurance$insurance==2,1,0))
  di1 <- did.summarize(data.insurance[-which(data.insurance$insurance==1),],ifelse(data.insurance$insurance[-which(data.insurance$insurance==1)]==2,1,0))
  de <- did.summarize(data.stage,data.stage$earlystage)
  dl <- did.summarize(data.stage,data.stage$latestage)
  de <- rbind(de,NA,NA)
  dl <- rbind(dl,NA,NA)
  rownames(de)[seq(nrow(de)-1,nrow(de))] <- rownames(di1)[seq(nrow(de)-1,nrow(de))]
  rownames(dl)[seq(nrow(de)-1,nrow(de))] <- rownames(di1)[seq(nrow(de)-1,nrow(de))]
  return(cbind(de,dl,di1,di2))
}

library(lmtest)
library(sandwich)
matc.ncdb <- did.summarize.all(data.insurance.ncdb,data.han.ncdb)
matc.seer <- did.summarize.all(data.insurance.seer,data.han.seer)
write.csv(matc.ncdb,"DCP-cancertype-DID-NCDB.csv")
write.csv(matc.seer,"DCP-cancertype-DID-SEER.csv")


##########################
## Code for Figure 2
##########################

library(ggplot2)
library(readxl)  #for reading in excel file
test_ins<- read_excel("~/testb.xlsx", sheet=2) #insurance data (rows ordered)
test_estage<-read_excel("~/testb.xlsx", sheet=1) #early stage data (rows ordered)
test_lstage<-read_excel("~/testb.xlsx", sheet=3) #late stage data (rows ordered)
names(test_ins)<-c("x", "y","ylo", "yhi", "Population") 
names(test_estage)<-c("x", "y","ylo", "yhi", "Population") 
names(test_lstage)<-c("x", "y","ylo", "yhi", "Population") 

vals <- test_ins[ test_ins$Population == 'SEER', ]
#xvals <- vals[with(vals, order(y)), ]$x #### this reorders by value of SEER y
xvals <- vals$x ##### this keeps the order from the excel sheets
test_ins$xseery <- factor(test_ins$x, levels = xvals)

# uses pointrangeh instead to fix legend, which works to make sideways
# changes order in legend with guide
# uses new var created above for ordering by SEER subgroup value
# colors from color brewer will show up different in printed/copied grayscale
library(ggstance)
library(stringr)
p1<-ggplot(test_ins, aes(x=y, y=xseery, color=Population))+ 
  geom_pointrangeh(aes(xmin=ylo, xmax=yhi),
                   position = position_dodgev(height=0.5), alpha=.9, fatten=3) + 
  scale_color_manual(values=c("#fdae61","#2b83ba")) + geom_vline(xintercept = 0, linetype=2)+
  xlab("") +
  ylab("") + 
    scale_x_continuous(limits=c(-40,+42)) +
  theme_minimal() +
  theme(legend.position="none",axis.text.y=element_text(face=ifelse(levels(test_ins$xseery)=="Overall","bold","plain"),
                                 color = ifelse(levels(test_ins$xseery)=="Overall","#2b83ba","black")))  +
  guides(fill=guide_legend(reverse=TRUE), 
            colour=guide_legend(reverse=TRUE))
p1

# reorder x by value of y in SEER subset
vals <- test_estage[ test_estage$Population == 'SEER', ]
#xvals <- vals[with(vals, order(y)), ]$x #### this reorders by value of SEER y
xvals <- vals$x ##### this keeps the order from the excel sheets
test_estage$xseery <- factor(test_estage$x, levels = xvals)


# uses pointrangeh instead to fix legend, which works to make sideways
# changes order in legend with guide
# uses new var created above for ordering by SEER subgroup value
# colors from color brewer will show up different in printed/copied grayscale
library(ggstance)
library(stringr)
p2<-ggplot(test_estage, aes(x=y, y=xseery, color=Population))+ 
  geom_pointrangeh(aes(xmin=ylo, xmax=yhi),
                   position = position_dodgev(height=0.5), alpha=.9, fatten=3) + 
  scale_color_manual(values=c("#fdae61","#2b83ba")) +  
   geom_vline(xintercept = 0, linetype=2)+
  xlab("") +
  ylab("") + 
  scale_x_continuous(limits=c(-40,+42)) +
  theme_minimal() +
  theme(legend.position="none", axis.text.y=element_text(face=ifelse(levels(test_estage$xseery)=="Overall","bold","plain"),
                                 color = ifelse(levels(test_estage$xseery)=="Overall","#2b83ba","black"))) +
  guides(fill=guide_legend(reverse=TRUE), 
            colour=guide_legend(reverse=TRUE))
p2

# reorder x by value of y in SEER subset
vals <- test_lstage[ test_lstage$Population == 'SEER', ]
#xvals <- vals[with(vals, order(y)), ]$x #### this reorders by value of SEER y
xvals <- vals$x ##### this keeps the order from the excel sheets
test_lstage$xseery <- factor(test_lstage$x, levels = xvals)


# uses pointrangeh instead to fix legend, which works to make sideways
# changes order in legend with guide
# uses new var created above for ordering by SEER subgroup value
# colors from color brewer will show up different in printed/copied grayscale
library(ggstance)
library(stringr)
p3<-ggplot(test_lstage, aes(x=y, y=xseery, color=Population))+ 
  geom_pointrangeh(aes(xmin=ylo, xmax=yhi),
                   position = position_dodgev(height=0.5), alpha=.9, fatten=3) + 
  scale_color_manual(values=c("#fdae61","#2b83ba")) +  
   geom_vline(xintercept = 0, linetype=2)+
  xlab("Difference in Differences") +
  ylab("") + 
  scale_x_continuous(limits=c(-40,+42)) +
  theme_minimal() +
  theme(legend.position="bottom", axis.text.y=element_text(face=ifelse(levels(test_lstage$xseery)=="Overall","bold","plain"),
                                 color = ifelse(levels(test_lstage$xseery)=="Overall","#2b83ba","black"))) +
  guides(fill=guide_legend(reverse=TRUE), 
            colour=guide_legend(reverse=TRUE))
p3





