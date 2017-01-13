library (ggplot2)
library (ROCR)
library(gridExtra)
library(grid)
library(caret)
library(rpart)
require(plyr)

setwd ("/Users/Vlad/projects/project-VSP/")


ClinVar <-read.csv ("Data/ClinVar_built.csv",  sep = ",")
ClinVar <- ClinVar[ which(ClinVar$GeneSymbol=='MAPT'),]
ClinVar <- rename (ClinVar, c("Mutation_3" = "MUTATION"))

MAPT_MutPred <-read.csv ("Data/MAPT_MutPredScores.csv", sep = ",")


#df <- df[,colSums(is.na(df))<nrow(df)]

MAPT_dbNSFP <- read.csv (file="Data/MAPT_dbNSFPa_output.csv", stringsAsFactors=FALSE, na.strings=c(NA,"NA"," NA"), sep =",")

MAPT_dbNSFP <- rename(MAPT_dbNSFP, c("MutationAssessor_variant"="MUTATION"))


Merged_DF<-merge(ClinVar,MAPT_MutPred, by="MUTATION", all =FALSE)


#read data files
#PROV<-read.csv(PROV_DATA, sep =',')
#PROV_transformed<-reshape(PROV, idvar="position", v.names = "PROV_Scores", timevar ="Index", varying = c("A", "C", "D", "E","F", "G", "H", "I", "K", "L", "M", "N","P", "Q", "R", "S", "T", "V", "W", "Y","Del"), direction="long")
#colnames (PROV_transformed)[2] <-"AA_POSITION"
#colnames (PROV_transformed)[4] <-"PROVEAN_SCORE"

Merged_DF1<-merge(df,PKD1_dbNSFP, by=c("genename","MUTATION"), all =FALSE)
Merged_DF2<-merge(df,PKD2_dbNSFP, by=c("genename","MUTATION"), all =FALSE)
Merged_DF <-rbind(Merged_DF1, Merged_DF2)
Merged_DF <-Merged_DF[!duplicated(Merged_DF[c("MUTATION")]),]
#Merged_DF$GVGD<-Merged_DF$GD-(Merged_DF$GV)*1.5


#Merged_DF$FATHMM_score<-1/Merged_DF$FATHMM_score  
#Merged_DF<-cbind(as.character(Merged_DF$MUTATION),apply(Merged_DF[,2:dim(Merged_DF)[2]],2,as.numeric))
#colnames(Merged_DF)[1]<-"MUTATION"

#Merged_DF[is.na(Merged_DF)]<-0
write.csv(Merged_DF, file = "Merged_DF.csv")