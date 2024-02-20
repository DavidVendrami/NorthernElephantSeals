##### msats #####

library(adegenet)
library(dplyr)

# The following function to add panels letters was found in some blog, but I do not rememebr who wrote it and I can't find it anymore! 
# Please get in touch if this looks like the function you wrote and I'll happily acknowledge that's yours. 
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.02,0.9),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user") + 0.1
  text(labels=label[1], x=this.x, y=this.y, xpd=T, adj=0,...)
}


data<-read.structure('msats_genotypes.stru', n.ind=219, n.loc=24, onerowperind = T, col.lab = 1, col.pop=0,col.others=0,row.marknames=0)
# How many genotypes are there? 219
#
# How many markers are there? 24
#
# Which column contains labels for genotypes ('0' if absent)? 1
#
# Which column contains the population factor ('0' if absent)? 0
#
# Which other optional columns should be read (press 'return' when done)? 1: 
#
# Which row contains the marker names ('0' if absent)? 0
#
# Are genotypes coded by a single row (y/n)? y
#
# Converting data from a STRUCTURE .stru file to a genind object... 

labs<-read.table('msats_dataset_description.txt',h=T)
labs$Animal_ID<-gsub('ES-','ES',labs$Animal_ID)
cols2<-c("#E41A1CFF","#4DAF4AFF","#984EA3FF","#FF7F00FF","#a65628","gold")

x.data <- tab(data, freq=TRUE, NA.method="mean")
pca.data<-dudi.pca(x.data, center = TRUE, scale = FALSE, scannf = FALSE, nf = 10)

pops<-data.frame(Animal_ID=row.names(data@tab))
pops<-left_join(pops,labs,by='Animal_ID')

### SNPs ###
snps<-read.PLINK("NES_notImp.raw",n.cores=1)
bla<-read.table("74_Samples_sMLH.txt",h=T)
ch<-read.table("NES_notImp.raw",h=T)
pca1<-glPca(snps,parallel=FALSE, nf = 10)

info<-ch[,1:6]
info$FID<-bla$Class
info$Multiclass<-bla$Primary_classification
cols<-c()
cols[1:74]<-"#777777"
info$cols<-cols
info$cols[info$FID=="Case"]<-"#6495ED"

# Multiclass classification
pops$Primary_classification[pops$Primary_classification=="Primary_bacterial_infection" | pops$Primary_classification=="Secondary_bacterial_infection"] <- "Bacteria"
#data@pop <- as.factor(pops$Primary_classification)
#s.class(pca.data$li, fac=pop(data),col=transp(cols2,.6),axesel=FALSE, cstar=0, cpoint=3, grid=F,label=c())

colsM<-c()
colsM[1:219]<-"gold" # Cong_defect = gold
pops$colsM<-colsM
pops$colsM[pops$Primary_classification == "Bacteria"] <- "#4DAF4AFF" # blue = case
pops$colsM[pops$Primary_classification == "Malnutrition"] <- "#984EA3FF"
pops$colsM[pops$Primary_classification == "Otostrongylis"] <- "#E41A1CFF"
pops$colsM[pops$Primary_classification == "Protozoa"] <- "#a65628"
pops$colsM[pops$Primary_classification == "Trauma"] <- "#FF7F00FF"

bla$Primary_classification[bla$Primary_classification=="Primary_bacterial_infection" | bla$Primary_classification=="Secondary_bacterial_infection"] <- "Bacteria"
colsM<-c()
colsM[1:74]<-"gold" # Cong_defect = gold
info$colsM<-colsM
info$colsM[bla$Primary_classification == "Bacteria"] <- "#4DAF4AFF" # blue = case
info$colsM[bla$Primary_classification == "Malnutrition"] <- "#984EA3FF"
info$colsM[bla$Primary_classification == "Otostrongylis"] <- "#E41A1CFF"
info$colsM[bla$Primary_classification == "Protozoa"] <- "#a65628"
info$colsM[bla$Primary_classification == "Trauma"] <- "#FF7F00FF"

pdf("PCA_plots_Multiclass.pdf",width=9, height=9)
par(mfrow=c(2,2))
plot(pca.data$li$Axis1, pca.data$li$Axis2, xlab="PC1", ylab="PC2", main="", pch = 20, cex = 2.5, col=transp(pops$colsM,.6),
 	 axes=F)
axis(1,c(seq(-1.5,1.5,by=0.5)))
axis(2,c(seq(-1.5,1.5,by=0.5)))
box(bty="l")
put.fig.letter(label="(a) Microsatellites PC1 and PC2", location="topleft",cex=1)

plot(pca.data$li$Axis3, pca.data$li$Axis4, xlab="PC3", ylab="PC4", main="", pch = 20, cex = 2.5, col=transp(pops$colsM,.6),
 	 axes=F)
axis(1,c(seq(-1.5,1.5,by=0.5)))
axis(2,c(seq(-1.5,1.5,by=0.5)))
box(bty="l")
put.fig.letter(label="(b) Microsatellites PC3 and PC4", location="topleft",cex=1)

plot(pca1$scores[,1], pca1$scores[,2], xlab="PC1", ylab="PC2", main="", pch = 20, cex = 2.5, col=transp(info$colsM,.6),
	 axes=F)
axis(1,c(seq(-10,25,by=5)))
axis(2,c(seq(-15,15,by=5)))
box(bty="l")
put.fig.letter(label="(c) SNPs PC1 and PC2", location="topleft",cex=1)

plot(pca1$scores[,3], pca1$scores[,4], xlab="PC3", ylab="PC4", main="", pch = 20, cex = 2.5, col=transp(info$colsM,.6),
	 axes=F)
axis(1,c(seq(-10,25,by=5)))
axis(2,c(seq(-15,15,by=5)))
box(bty="l")
put.fig.letter(label="(d) SNPs PC1 and PC2", location="topleft",cex=1)
dev.off()
