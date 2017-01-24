#BalchLab 2016
#Creates a geostatistical map of genes using ExAC data

library (fields)
library (sp)
library (gstat)
library (ggplot2)
require (reshape2)
library (survival)
library (ggthemes)
library (scales)
library ('dplyr')
library (mice)
library (randomForest)
library(magrittr)

#Set variables 


setwd ("/Users/Vlad/projects/project-VSP/")
myvars <- c("GENE_ID", "AA_POS", "MUTATION", "ALLELE.COUNT", "Pathogenic")
PROV_DATA <- "Data/MAPT_PROVEANScores.csv"
EXAC_DATA <- "Data/MAPT_ExACScores.csv"
MUTPRED_DATA <- "Data/MAPT_MutPredScores.csv"
AmylSurf_DATA <- "Data/MAPT_Amyl_consurf.csv"
SavePlot <-"/VSP-images/PKD2_MutPred-EXAC.png"
AC_Limit <- 1000
GraphTitle <- c("MAPT MutPpred - ExAC, Log transformed")

AmylSurf <- read.csv (AmylSurf_DATA,  sep = ",")
                
ClinVar <-read.csv ("Data/ClinVar_built.csv",  sep = ",")
ClinVar <- ClinVar[ which(ClinVar$GeneSymbol=='MAPT'),]
colnames (ClinVar)[10] <-"MUTATION"
colnames (ClinVar)[5] <-"GENE_ID"
colnames (ClinVar)[11] <-"AA_POS"
ClinVar["ALLELE.COUNT"] <- 0.10
ClinVar["Pathogenic"] <- "P"

#read data files
EXAC<-read.csv(EXAC_DATA)
EXAC["Pathogenic"] <- ""
EXAC <- EXAC[myvars]
ClinVar <-ClinVar[myvars]

#PROV<-read.csv(PROV_DATA, sep =',')
#PROV_transformed<-reshape(PROV, idvar="position", v.names = "PROV_Scores", timevar ="Index", varying = c("A", "C", "D", "E","F", "G", "H", "I", "K", "L", "M", "N","P", "Q", "R", "S", "T", "V", "W", "Y","Del"), direction="long")
MutPred<-read.csv(MUTPRED_DATA)
#colnames (PROV_transformed)[2] <-"AA_POSITION"
#colnames (PROV_transformed)[4] <-"PROVEAN_SCORE"
Merged_DF <- rbind(EXAC, ClinVar)
Merged_DF<-merge(Merged_DF,MutPred, by=c("MUTATION"), all =FALSE)



Merged_DF[is.na(Merged_DF)]<-0

#modify amino acid notation to drop .p (ExAC format)
EXAC$Ptein.Consequence<-gsub("p.","", as.character(EXAC$amino_acid_change))

#adjust prediction scores - so that the range is similar to AA coordinates
Merged_DF$PROVEAN_SCORE <- Merged_DF$PROVEAN_SCORE*-100
Merged_DF$MUTPRED.Score<-round(Merged_DF$MUTPRED.Score*100, 0)
Merged_DF$MUTPRED.Score<-rescale(Merged_DF$MUTPRED.Score, to =c(1,800)) #TODO:needs to be auto set to length of protein
Merged_DF$ALLELE.COUNT.log<-log(Merged_DF$ALLELE.COUNT)


Merged_DF<-merge(Merged_DF,AmylSurf, by=c("AA_POS"), all =FALSE)


Merged_DF$a4v<-rescale(Merged_DF$a4v, to =c(1,700))

write.csv(Merged_DF, file = "Data/TAU_Merged_Data.csv")


#set max allele count to a value best suited for analysis
#TODO:Figure out how to do this algorithmically
#Merged_DF$ALLELE.COUNT[Merged_DF$ALLELE.COUNT>AC_Limit]<-AC_Limit

#Spacial Interpolation using amino acid index and prediction score as x and y axis
OP      <- par( mar=c(2,2,2,2))
coordinates(Merged_DF) <- c("AA_POS","a4v")
plot(Merged_DF, pch=16, cex=((Merged_DF$ALLELE.COUNT.log-1)/200))
text(Merged_DF, as.character(Merged_DF$ALLELE.COUNT.log), pos=3, col="grey", cex=0.8)

# Create an empty grid where n is the total number of cells
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

grd              <- as.data.frame(spsample(Merged_DF, "regular", n=60000))
names(grd)       <- c("AA_POS", "a4v")
coordinates(grd) <- c("AA_POS", "a4v")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# Interpolate the surface using a power value of 2 (idp=2.0)
Merged_DF.idw <- idw(ALLELE.COUNT.log~1,Merged_DF,newdata=grd,idp=3.0)
# Plot the raster and the sampled points
#OP      <- par( mar=c(6,2,2,2))
png(filename='MAPT_Aggregation.png')

par(mar = c(4, 3, 2,4 ), mgp = c(1, 0.5, 0))

image(Merged_DF.idw,"var1.pred",col=terrain.colors(50), xlab="AA POSITION", ylab="Aggregation", mgp = c(1, 4, 0))
x.tick.number <- 10
y.tick.number <- 10

at <- seq(1, max(Merged_DF$AA_POS), length.out=x.tick.number)
axis(1, at = Merged_DF$AA_POSITION, labels = Merged_DF$AA_POSITION, las =1, pos=0)
                             #round(EXAC$MUT_PRED, 0)
at <- seq(1, max(Merged_DF$ALLELE.COUNT.log), length.out=y.tick.number)
axis(2, at = at, labels = NULL , las = 1, pos=0)

contour(Merged_DF.idw,"var1.pred", add=TRUE, nlevels=10, col="#656565")
box()
plot(Merged_DF, add=TRUE, pch=16, cex=0.5)
#text(coordinates(Merged_DF), as.character(round(Merged_DF$Pathogenic,1)), pos=4, cex=0.8, col="blue")
text(coordinates(Merged_DF), as.character(Merged_DF$Pathogenic,1), pos=4, cex=0.8, col="red")
parameters<-par(OP)
str(OP)
title(main=GraphTitle, font.main =4)
dev.off()



# #Sals function to strip amino acid notation for variants... Maybe put in a separate file? 
# split_cols<-function(foo){
#   out<-vector()
#   for(i in 1:length(foo)){
#     
#     vec<-unlist(strsplit(foo[i],split =""))
#     
#     first<-vec[1]
#     last<-vec[length(vec)]
#     mid<-paste(vec[2:(length(vec)-1)],collapse="")
#     
#     int_row<-c(first,mid,last)
#     out=rbind(out,int_row)
#   }
#   return(out)
#   
# }

#Attempt at Kriging, work in progress
# EXAC %>% as.data.frame %>% 
#   ggplot(aes(AA_POSITION, PROVEAN_SCORE)) + geom_point(aes(size=ALLELE.COUNT), color="blue", alpha=3/4) + 
#   ggtitle("Zinc Concentration (ppm)") + coord_equal() + theme_bw()
# 
# class(EXAC)
# coordinates(EXAC) <- ~AA_POSITION + PROVEAN_SCORE
# class(EXAC)
# bbox(EXAC)
# coordinates(EXAC)
# proj4string(EXAC)
# identical(bbox(EXAC), EXAC@bbox)
# identical(coordinates(EXAC), EXAC@coords)
# EXAC@data%>%glimpse
# EXAC%>%as.data.frame%>%glimpse
# 
# lzn.vgm <- variogram(log(ALLELE.COUNT)~1, EXAC) # calculates sample variogram values 
# 
# lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph",800,-200)) # fit model
# 
# plot(lzn.vgm, lzn.fit)
# 
# 
# 
# # to compare, recall the bubble plot above; those points were what there were values for. this is much more sparse
# plot1 <- EXAC%>% as.data.frame %>%
#   ggplot(aes(AA_POSITION, PROVEAN_SCORE)) + geom_point(size=1) + coord_equal() + 
#   ggtitle("Points with measurements")
# plot1
# 
# 
# x<-c(1:968)
# y<-c(-163:804)
# x_name<-"AA_POSITION"
# y_name<-"PROVEAN_SCORE"
# require(reshape2)
# EXAC.grid<-melt(data.frame(x,y))
# 
# colnames(EXAC.grid)<-c(x_name,y_name)
# 
# # this is clearly gridded over the region of interest
# plot2 <- EXAC.grid %>% as.data.frame %>%
#   ggplot(aes(AA_POSITION, PROVEAN_SCORE)) + geom_point(size=1) + coord_equal() + 
#   ggtitle("Points at which to estimate")
# 
# coordinates(EXAC) <- ~ AA_POSITION + PROVEAN_SCORE # step 3 above
# lzn.kriged <- krige(log(ALLELE.COUNT) ~ 1, EXAC, EXAC.grid, model=lzn.fit)
# 
# library(gridExtra)
# grid.arrange(plot1, plot2, ncol = 2)
# 
# 
# p + geom_text(angle=45) 

# 
# 
# 
# #Stripping protein variant notation into REF - ALT - AA_POSITION
# PROVEAN<-split_cols(EXAC$Ptein.Consequence)
# write.csv(PROVEAN,"PKD2_provean.csv")
# 
# PKD2_muts<-as.character(exac_PKD2$Conseqence)
# 
# PKD2_PROVEAN<-split_cols(PKD2_muts) 
# 
# write.csv(PKD2_PROVEAN,"PKD2_provean.csv")
# 
# mut_pred<-read.csv("gainullin_mutpred.csv")
# 
# mut_pred$mutation<-sub("^","p.",mut_pred$mutation)
# 
# 
# PKD2_mut_pred <-subset(mut_pred, ID=="PKD2_HUMAN")
# 
# 
# #Not sure what this is anymore
# PKD2_tran<-read.csv("PKD2_var_tran.csv")
# NameList<-c("amino_acid_change","effect")
# idx<-match(NameList, names(exac2))
# idx<-sort(c(idx-1,idx))
# exac_sub<-exac2[,idx]
# 
# 
# drops<-c("Alternate", "Annotation")
# exac_sub<-exac_sub[,!(names(exac_sub) %in% drops)]
# exac_sub<-unique(exac_sub)
# PKD2_exac_forMutpred<-exac_sub[!names(exac_sub) %in% 'effect']
# 
# write.csv(PKD2_exac_forMutpred,'PKD2_HUMAN_all.csv')
# 
# drops<-c("transcript_name", "gene", "effect")
# PKD2_tran<-PKD2_tran[,!(names(PKD2_tran) %in% drops)]
# 
# merged_exac_PKD2 <-merge(x=PKD2_tran, y=exac_sub, by = "amino_acid_change", all = TRUE)
# merged_exac_PKD2$effect[is.na(merged_exac_PKD2$effect)]<-0
# 
# 
