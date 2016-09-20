#BalchLab 2016
#Creates a geostatistical map of genes using ExAC data

library(fields)
library(sp)
library(gstat)

library (ggplot2)
require (reshape2)
library (survival)
library ('ggthemes')
library ('scales')
library ('dplyr')
library ('mice')
library ('randomForest')
library(scales) # for "comma"
library(magrittr)

#Set variables 


setwd ("D:/Balclab/projects/project-VSP/")

#read data files
exac_ENSP<-read.csv("Data/Exac_PKD2_ENSP.csv")
PROV<-read.csv("Data/PKD2PROVEANScores.csv", sep ='')
PROV_transformed<-reshape(PROV, idvar="position", v.names = "PROV_Scores", timevar ="Index", varying = c("A", "C", "D", "E","F", "G", "H", "I", "K", "L", "M", "N","P", "Q", "R", "S", "T", "V", "W", "Y","Del"), direction="long")
colnames (PROV_transformed)[2] <-"AA_POSITION"
colnames (PROV_transformed)[4] <-"PROVEAN_SCORE"
Merged_DF_PKD2<-merge(exac_ENSP,PROV_transformed, by=c("AA_POSITION", "PROVEAN_SCORE"), all =TRUE)
Merged_DF_PKD2[is.na(Merged_DF_PKD2)]<-0

#modify amino acid notation to drop .p (ExAC format)
exac_ENSP$Ptein.Consequence<-gsub("p.","", as.character(exac_ENSP$amino_acid_change))

#adjust prediction scores - so that the range is similar to AA coordinates
Merged_DF_PKD2$PROVEAN_SCORE <- Merged_DF_PKD2$PROVEAN_SCORE*-100
exac_ENSP$MUT_PRED<-round(exac_ENSP$MUT_PRED*100, 0)
exac_ENSP$MUT_PRED<-rescale(exac_ENSP$MUT_PRED, to =c(1,1000))

#set max allele count to a value best suited for analysis 
Merged_DF_PKD2$Allele.Count[Merged_DF_PKD2$Allele.Count>5]<-5


exac_CFTR<-read.csv("CFTR_Exac.csv")
exac_CFTR$Conseqence <- as.character(exac_CFTR$Conseqence)
exac_CFTR$PROVEAN_SCORE <- exac_CFTR$PROVEAN_SCORE*120
exac_CFTR$Allele.Count[exac_CFTR$Allele.Count>10]<-10


#Spacial Interpolation for  CFTR
coordinates(exac_CFTR) <- c("AA_POSITION","PROVEAN_SCORE")
plot(exac_CFTR, pch=16, ,cex=( (exac_CFTR$Allele.Count-1)/250))
text(exac_CFTR, as.character(exac_CFTR$Allele.Count), pos=3, col="grey", cex=0.8)

# Create an empty grid where n is the total number of cells
grd_cf              <- as.data.frame(spsample(exac_CFTR, "regular", n=10000))
names(grd_cf)       <- c("AA_POSITION", "PROVEAN_SCORE")
coordinates(grd_cf) <- c("AA_POSITION", "PROVEAN_SCORE")
gridded(grd_cf)     <- TRUE  # Create SpatialPixel object
fullgrid(grd_cf)    <- TRUE  # Create SpatialGrid object

# Interpolate the surface using a power value of 2 (idp=2.0)
exac_CFTR.idw <- idw(Allele.Count~1,exac_CFTR,newdata=grd_cf,idp=2.0)

# Plot the raster and the sampled points
OP      <- par( mar=c(2,2,2,2))
image(exac_CFTR.idw,"var1.pred",col=terrain.colors(20))
contour(exac_CFTR.idw,"var1.pred", add=TRUE, nlevels=10, col="#656565")
plot(exac_CFTR, add=TRUE, pch=16, cex=0.5)
text(coordinates(exac_CFTR), as.character(round(exac_CFTR$AA_POSITION,1)), pos=4, cex=0.8, col="blue")
par(OP)

#Spacial Interpolation using amino acid index and prediction score as x and y axis
#OP      <- par( mar=c(2,2,2,2))
coordinates(Merged_DF_PKD2) <- c("AA_POSITION","PROVEAN_SCORE")
plot(Merged_DF_PKD2, pch=16, ,cex=( (Merged_DF_PKD2$Allele.Count-1)/200))
text(Merged_DF_PKD2, as.character(Merged_DF_PKD2$Allele.Count), pos=3, col="grey", cex=0.8)

# Create an empty grid where n is the total number of cells
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

grd              <- as.data.frame(spsample(Merged_DF_PKD2, "regular", n=60000))
names(grd)       <- c("AA_POSITION", "PROVEAN_SCORE")
coordinates(grd) <- c("AA_POSITION", "PROVEAN_SCORE")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# Interpolate the surface using a power value of 2 (idp=2.0)
Merged_DF_PKD2.idw <- idw(Allele.Count~1,Merged_DF_PKD2,newdata=grd,idp=3.0)
# Plot the raster and the sampled points
#OP      <- par( mar=c(6,2,2,2))
par(mar = c(4, 3, 2,4 ), mgp = c(1, 0.5, 0))
image(Merged_DF_PKD2.idw,"var1.pred",col=terrain.colors(50), xlab="AA POSITION", ylab="PROVEAN", mgp = c(1, 4, 0))

axis(1, at = exac_ENSP$AA_POSITION, labels = exac_ENSP$AA_POSITION, las =1, pos=0)
                                #round(exac_ENSP$MUT_PRED, 0)
x.tick.number <- 10
at <- seq(1, max(Merged_DF_PKD2.idw$PROVEAN_SCORE), length.out=x.tick.number)
axis(2, at = at, labels = NULL , las = 1, pos=0)

contour(Merged_DF_PKD2.idw,"var1.pred", add=TRUE, nlevels=10, col="#656565")
#box()
plot(Merged_DF_PKD2, add=TRUE, pch=16, cex=0.5)
#text(coordinates(exac_ENSP), as.character(round(exac_ENSP$AA_POSITION,1)), pos=4, cex=0.8, col="blue")
parameters<-par(OP)
str(OP)
title(main="PKD2 - PROVEAN Exac", font.main =4)



# Interpolate the surface using a power value of 2 (idp=2.0)


# Plot the raster and the sampled points
OP      <- par( mar=c(0,0,0,0))
image(exac_ENSP.idw,"var1.pred",col=terrain.colors(20), axes = TRUE)
axis(1, at = seq(100, 800, by =100))
axis(2, at = seq(100, 800, by =100))

contour(exac_ENSP.idw,"var1.pred", add=TRUE, nlevels=10, col="#656565", )
plot(dat, add=TRUE, pch=16, cex=0.5)
text(coordinates(dat), as.character(round(dat$Z,1)), pos=4, cex=0.8, col="blue")
par(OP)



#Attempt at Kriging, work in progress
exac_ENSP %>% as.data.frame %>% 
  ggplot(aes(AA_POSITION, PROVEAN_SCORE)) + geom_point(aes(size=Allele.Count), color="blue", alpha=3/4) + 
  ggtitle("Zinc Concentration (ppm)") + coord_equal() + theme_bw()

class(exac_ENSP)
coordinates(exac_ENSP) <- ~AA_POSITION + PROVEAN_SCORE
class(exac_ENSP)
bbox(exac_ENSP)
coordinates(exac_ENSP)
proj4string(exac_ENSP)
identical(bbox(exac_ENSP), exac_ENSP@bbox)
identical(coordinates(exac_ENSP), exac_ENSP@coords)
exac_ENSP@data%>%glimpse
exac_ENSP%>%as.data.frame%>%glimpse

lzn.vgm <- variogram(log(Allele.Count)~1, exac_ENSP) # calculates sample variogram values 

lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph",800,-200)) # fit model

plot(lzn.vgm, lzn.fit)



# to compare, recall the bubble plot above; those points were what there were values for. this is much more sparse
plot1 <- exac_ENSP%>% as.data.frame %>%
  ggplot(aes(AA_POSITION, PROVEAN_SCORE)) + geom_point(size=1) + coord_equal() + 
  ggtitle("Points with measurements")
plot1


x<-c(1:968)
y<-c(-163:804)
x_name<-"AA_POSITION"
y_name<-"PROVEAN_SCORE"
require(reshape2)
exac_ENSP.grid<-melt(data.frame(x,y))

colnames(exac_ENSP.grid)<-c(x_name,y_name)

# this is clearly gridded over the region of interest
plot2 <- exac_ENSP.grid %>% as.data.frame %>%
  ggplot(aes(AA_POSITION, PROVEAN_SCORE)) + geom_point(size=1) + coord_equal() + 
  ggtitle("Points at which to estimate")

coordinates(exac_ENSP) <- ~ AA_POSITION + PROVEAN_SCORE # step 3 above
lzn.kriged <- krige(log(Allele.Count) ~ 1, exac_ENSP, exac_ENSP.grid, model=lzn.fit)

library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)


p + geom_text(angle=45) 

#Sals function to strip amino acid notation for variants... Maybe put in a separate file? 
split_cols<-function(foo){
  out<-vector()
  for(i in 1:length(foo)){
    
    vec<-unlist(strsplit(foo[i],split =""))
    
    first<-vec[1]
    last<-vec[length(vec)]
    mid<-paste(vec[2:(length(vec)-1)],collapse="")
    
    int_row<-c(first,mid,last)
    out=rbind(out,int_row)
  }
  return(out)
  
}


#Stripping protein variant notation into REF - ALT - AA_POSITION
PROVEAN<-split_cols(exac_ENSP$Ptein.Consequence)
write.csv(PROVEAN,"PKD2_provean.csv")

CFTR_muts<-as.character(exac_CFTR$Conseqence)

CFTR_PROVEAN<-split_cols(CFTR_muts) 

write.csv(CFTR_PROVEAN,"CFTR_provean.csv")

mut_pred<-read.csv("gainullin_mutpred.csv")

mut_pred$mutation<-sub("^","p.",mut_pred$mutation)


PKD2_mut_pred <-subset(mut_pred, ID=="PKD2_HUMAN")


#Not sure what this is anymore
PKD2_tran<-read.csv("PKD2_var_tran.csv")
NameList<-c("amino_acid_change","effect")
idx<-match(NameList, names(exac2))
idx<-sort(c(idx-1,idx))
exac_sub<-exac2[,idx]


drops<-c("Alternate", "Annotation")
exac_sub<-exac_sub[,!(names(exac_sub) %in% drops)]
exac_sub<-unique(exac_sub)
PKD2_exac_forMutpred<-exac_sub[!names(exac_sub) %in% 'effect']

write.csv(PKD2_exac_forMutpred,'PKD2_HUMAN_all.csv')

drops<-c("transcript_name", "gene", "effect")
PKD2_tran<-PKD2_tran[,!(names(PKD2_tran) %in% drops)]

merged_exac_PKD2 <-merge(x=PKD2_tran, y=exac_sub, by = "amino_acid_change", all = TRUE)
merged_exac_PKD2$effect[is.na(merged_exac_PKD2$effect)]<-0


