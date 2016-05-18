library(tidyr)
data <- read.csv("LesionEccentricity.HEM.csv")
#remove threshhold rows
data2 <- data[c(-1, -2, -3, -4),]
#make two columns for Chrom and Pos
data3 <- separate(data2, X.1, into = c("Chrom", "Pos"))
head(data3)

#order Chrom values, then Pos values
#must make Pos values numeric object
data3$Pos <- as.numeric(data3$Pos)
data4 <- data3[order(data3$Chrom, data3$Pos), ]

#write csv without thresholds
write.csv(data4, "LesionEccentricity.HEM.PlottingFormat.csv")

##################### DETERMINING DATA STRUCTURE FOR ECCENTRICITY
library(vioplot)
data <- read.csv("GWASphenotypesALL.csv")
eccentricity  <- as.data.frame(c(data[,72:73],data[,75:76]))
yellow  <- as.data.frame(c(data[,47:48],data[,50:51]))
green  <- as.data.frame(c(data[,32:33],data[,35:36]))
size  <- as.data.frame(c(data[,87:88],data[,90:91]))
 str(eccentricity)
boxplot(data$Lesion.0.m.eccentricity_Apple517)
boxplot(data$Lesion.0.m.eccentricity_B05.10)
boxplot(data$Lesion.0.m.eccentricity_Supersteak)
boxplot(data$Lesion.0.m.eccentricity_UKRazz)

vioplot(data$Lesion.0.m.eccentricity_Apple517)
vioplot(data$Lesion.0.m.eccentricity_B05.10)
vioplot(data$Lesion.0.m.eccentricity_Supersteak)
vioplot(data$Lesion.0.m.eccentricity_UKRazz)

#check distribution
library(fitdistrplus)
plotdist(data$Lesion.0.m.eccentricity_Apple517, demp = TRUE, hist = TRUE) 
         #main = "Lesion.0.m.eccentricity_Apple517" )
plotdist(data$Lesion.0.m.eccentricity_B05.10)
plotdist(data$Lesion.0.m.eccentricity_Supersteak)
plotdist(data$Lesion.0.m.eccentricity_UKRazz)
plotdist(data$Lesion.Size_Apple517, demp = TRUE)

plotdist(data$Lesion.0.m.eccentricity_Apple517, "pois", para=list(lambda = mean(data$Lesion.0.m.eccentricity_Apple517)))

###########################   WRITING POSITIONS BEYOND 95 THRESHOLD   ########################################
HEM.thresh <- read.csv("YlwLesionProportion.HEM.csv")#RF-this is just a reorganized file of LesionSize.HEM.csv.
#just need HEM.csv because in conversion to AGI, we'll need format like "II.51768" (not plotting format)

sig.HEM <- data.frame()
for(i in 5:dim(HEM.thresh)[1]) {     #when first four rows are threshold values
    if(i %% 1000 == 0) {print(paste(i, "--", Sys.time()))}
    if(any(abs(as.numeric(HEM.thresh[i,3:6]))-abs(as.numeric(HEM.thresh[1,3:6])) > 0)) { # 95% threshold is HEM.thresh[1,3:6], 99% threshold is HEM.thresh[3,3:6]
        sig.HEM <- (rbind(sig.HEM,HEM.thresh[i,]))
    }
}
colnames(sig.HEM) <- colnames(HEM.thresh)
write.csv(rbind(HEM.thresh[1,],sig.HEM),"LesionEccentricity.HEM.95Sig.csv")

#99th percentile
sig.HEM2  <- data.frame()
for(i in 5:dim(HEM.thresh)[1]) {
    if(i %% 1000 ==0) {print(paste(i, "--", Sys.time()))}
    if(any((abs(as.numeric(HEM.thresh[i,3:6])) - abs(as.numeric(HEM.thresh[3,3:6])))>0))
       { sig.HEM2 <- (rbind(sig.HEM2,HEM.thresh[i,]))}
}
colnames(sig.HEM2) <- colnames(HEM.thresh)
write.csv(rbind(HEM.thresh[3,],sig.HEM2), "LesionEccentricity.HEM.99Sig.csv")

#special loop for lesion eccentricity, looking only at Apple517 column
sig.HEM3  <- data.frame()
for(i in 5:dim(HEM.thresh)[1]){
    if(i %% 1000 == 0) {print(paste(i, "--", Sys.time()))}
    if( (abs(as.numeric(HEM.thresh[i,3])) - as.numeric(HEM.thresh[3,3]) ) >0){
        sig.HEM3 <- rbind(sig.HEM3, HEM.thresh[i,])
    }
}
colnames(sig.HEM3) <- colnames(HEM.thresh)
write.csv(rbind(HEM.thresh[3,],sig.HEM3), "LesionEccentricity.HEM.3.99Sig.csv")


####################################    MATCHING SNP TO AGI VIA JASON'S MAP    ##########################################
library(plyr)
map <- read.csv("SNPAGImap.csv")   #load SNP-AGI conversion df
sigsnps <- read.csv("YlwLesionProportion.HEM.csv")
colnames(map)
colnames(sigsnps)
sigsnps <- sigsnps[, -1] 
colnames(sigsnps)[1] <- "SNP" 
mergedsigsnps <- merge(map, sigsnps, by = "SNP")

write.csv(mergedsigsnps, "YlwLesionProportion.HEM.AGI.csv")

#separate chrom from bp
dat <- read.csv("YlwLesionProportion.HEM.AGI.csv")
library(tidyr)
dats <- separate(dat, SNP, into = c("chr", "pos"))
write.csv(dats, "YlwLesionProportion.HEM.AGI.sep.csv")


sigsnps  <- read.csv("LesionEccentricity.HEM.99Sig.csv") #load list of significant SNPs

sigsnps <- sigsnps[, -1] #remove an unnecessary column from the 95sig snps df

#check column names from both df so we can merge by SNP position
colnames(map)
colnames(sigsnps)
#give appropriate colname to SNP column so that merge function will work
colnames(sigsnps)[1] <- "SNP"  
#merge the two dfs by SNP id
mergedsigsnps <- merge(map, sigsnps, by = "SNP")
#check dimensions - merged should have one less because it removed threshold row
dim(mergedsigsnps)
dim(sigsnps)
#take out rows where there is no gene
justgenes <- mergedsigsnps[!is.na(mergedsigsnps$AGI),]

#create df with two columns that counts how many time each gene appears in the list
newbs <- ddply(justgenes, .(AGI), nrow)
#rename counts column
colnames(newbs)[2] <- "SNPsInGene"
#include the info from the new df above using merge by AGI
genecount <- merge(justgenes, newbs, by = "AGI")
#subset the df picking out only genes where more than one SNP showed up as significant
morethanone <- subset(genecount, SNPsInGene > 1)
#create list of unique gene names (morethanone df has duplicate telling which SNP was significant)
genes <- unique(morethanone$AGI)
write.csv(morethanone, "LesionEccentricity.HEM.99Sig.genelist.csv")
write.csv(genes, "LesionEccentricity.HEM.99Sig.justgenes.csv")

################ 99 Sig SNPS #########################
map <- read.csv("SNPAGImap.csv")   #load SNP-AGI conversion df
sigsnps  <- read.csv("LesionEccentricity.HEM.95Sig.csv") #load list of significant SNPs

sigsnps <- sigsnps[, -1] #remove an unnecessary column from the 95sig snps df
sigsnps <- sigsnps[, -1] #remove another unnecessary col

#check column names from both df so we can merge by SNP position
colnames(map)
colnames(sigsnps)
#give appropriate colname to SNP column so that merge function will work
#make sure first column of "sigsnps" is the chrom with snp position eg. "I.48475"
colnames(sigsnps)[1] <- "SNP"  
#merge the two dfs by SNP id
mergedsigsnps <- merge(map, sigsnps, by = "SNP")
#check dimensions - merged should have one less because it removed threshold row
dim(mergedsigsnps)
dim(sigsnps)
#take out rows where there is no gene
justgenes <- mergedsigsnps[!is.na(mergedsigsnps$AGI),]
library(plyr)
#create df with two columns that counts how many time each gene appears in the list
newbs <- ddply(justgenes, .(AGI), nrow)
#rename counts column
colnames(newbs)[2] <- "SNPsInGene"
#include the info from the new df above using merge by AGI
genecount <- merge(justgenes, newbs, by = "AGI")
#subset the df picking out only genes where more than one SNP showed up as significant
morethanone <- subset(genecount, SNPsInGene > 1)
#create list of unique gene names (morethanone df has duplicate rows for each significant SNP's separate position within the gene)
genes <- unique(morethanone$AGI)
write.csv(morethanone, "LesionEccentricity.HEM.95Sig.genelist.csv")
write.csv(genes, "LesionEccentricity.HEM.95Sig.justgenes.csv")

#################################################################################################################
                                ## NOW TO COMPARE GENE LISTS ##
#################################################################################################################

yellow <- read.csv("YlwLesionProportion.HEM.99Sig.justgenes.csv")
dim(yellow) #686
green <- read.csv("LesionGrn.mm.2.HEM.99Sig.justgenes.csv")
dim(green) #1352
area <- read.csv("LesionArea.mm.2.HEM.99Sig.justgenes.csv")
dim(area) #864


yellow95  <- read.csv("YlwLesionProportion.HEM.95Sig.justgenes.csv")
dim(yellow95) #3890 (31567 before removing SNPs in more than one gene)
green95 <- read.csv("LesionGrn.mm.2.HEM.95Sig.justgenes.csv")
dim(green95) #3504 (28909 before removing SNPs in more than one gene)
area95 <- read.csv("LesionArea.mm.2.HEM.95Sig.justgenes.csv")
dim(area95) #4464 (37780 before removing SNPs in more than one gene)

#extract vectors of each gene list
str(green)
areavec <- as.vector(area$x)#99th
grnvec <- as.vector(green$x)#99th
ylwvec <- as.vector(yellow$x)#99th
area95vec <- as.vector(area95$x)
grn95vec <- as.vector(green95$x)
ylw95vec <- as.vector(yellow95$x)

# #allthree
# star <- c(areavec, grnvec)
# overlap <- star[duplicated(star)]
# length(overlap)
# over2 <- c(overlap, ylwvec)
# over3 <- over2[duplicated(over2)]
# length(over3)

#take two 99 vectors and find common genes - doesn't exclude genes shared by all 3 - important to maintain these orders eg. areavec[areavec%in%ylwvec] - see below
#NOTE: %in% function returns binary vector the length of the vector preceding %in%, so this vector preceding this symbol must be the same vector which is subsetted to obtain an accurate result
AY <- areavec[areavec%in%ylwvec]
length(AY) #96 genes shared between AY
AG <- areavec[areavec%in%grnvec]
length(AG) #180 genes shared between AG
GY  <- grnvec[grnvec%in%ylwvec]
length(GY) #145 genes shared between GY

#write genes 99Sig for two phenotypes (all isolates) into a csv
write.csv(AY, "SizeYlw99SigSharedGenes.csv")
write.csv(AG, "SizeGrn99SigSharedGenes.csv")
write.csv(GY, "GrnYlw99SigSharedGenes.csv")

#########NOW make a list of genes shared by two phenos (99th) but also excluding genes found above 95th for third pheno (using AY, AG, GY from above) **(also removes any genes that will show up in allthree list)
ayx <- c(AY, grn95vec) #make a list of AY (genes shared by area99 and yellow99), and grn95
ayx2 <- ayx[duplicated(ayx)] #see which genes show up twice
length(ayx2) #725
ayx3 <- AY[! AY %in% ayx2] #remove those genes that showed up in both AY and grn95 (remaining genes will be above 99 for area, yellow, but below 95 for green)
write.csv(ayx3,"AY99belowG95.csv")

agx <- c(AG, ylw95vec) #make a list of AY (genes shared by area99 and yellow99), and grn95
agx2 <- agx[duplicated(agx)] #see which genes show up twice
length(agx2) #121
agx3 <- AG[! AG %in% agx2] #remove those genes that showed up in both AY and grn95 (remaining genes will be above 99 for area, yellow, but below 95 for green)
write.csv(agx3,"AG99belowY95.csv")

gyx <- c(GY, area95vec) #make a list of AY (genes shared by area99 and yellow99), and grn95
gyx2 <- gyx[duplicated(gyx)] #see which genes show up twice
length(gyx2) #85 duplicates
gyx3 <- GY[! GY %in% gyx2] #remove those genes that showed up in both GY and grn95 from GY (remaining genes will be above 99 for area, yellow, but below 95 for green)
write.csv(gyx3,"GY99belowA95.csv")


#what's shared by all three?
allthree <- AY[AY%in%AG]
length(allthree)
#33

#Now, make exclusive lists for each pheno (genes that show up in one pheno's 99th list, but NOT in the other two phenos 95th list) - removed 95th SNPs from other phenos one at a time
#yellow exclusive
#remove gen95 from ylw99
y <- c(ylwvec, grn95vec) #ylwvec starts at 686 length, y at 4190
y2 <- y[duplicated(y)]
length(y2) #duplicates are 322 in length
y3 <- ylwvec[! ylwvec %in% y2] #remove from ylw99 those genes that appeared in both ylw99 and y2 (duplicates)
#now remove area95 from above
y4 <- c(y3, area95vec)
y5 <- y4[duplicated(y4)]
length(y5) #167
y6 <- y3[! y3 %in% y5] #leaves us 197 genes from 686 that we started with
write.csv(y6,"yellowExclusive.csv")

#green exclusive
y <- c(grnvec, ylw95vec) #grnvec starts at 1352 length
y2 <- y[duplicated(y)]
y3 <- grnvec[! grnvec %in% y2]
#now remove area95 from above
y4 <- c(y3, area95vec)
y5 <- y4[duplicated(y4)]
length(y5) #299
y6 <- y3[! y3 %in% y5] #leaves us 429 genes from 1352 that we started with
write.csv(y6,"greenExclusive.csv")

#area exclusive
y <- c(areavec, ylw95vec) #areavec starts at 864 length
y2 <- y[duplicated(y)]
y3 <- areavec[! areavec %in% y2] #y3 is 441 in length
#now remove grn95 from above
y4 <- c(y3, grn95vec)
y5 <- y4[duplicated(y4)]
length(y5) #141
y6 <- y3[! y3 %in% y5] #leaves us 300 genes from 864 that we started with
write.csv(y6,"areaExclusive.csv")



#BELOW is the old method - where I took genes 99sig for one pheno and removed genes that showed up above the 99th for the other two phenos
# #read lists of genes according to phenotype and rm those that don't strictly apply to one and only one phenotype
# size <- as.vector(read.csv("size.csv")[,1])
# green <- as.vector(read.csv("green.csv")[,1])
# yellow <- as.vector(read.csv("yellow.csv")[,1])
# GY <- as.vector(read.csv("GrnYlw99SigSharedGenes.csv")[,2])
# GS <- as.vector(read.csv("SizeGrn99SigSharedGenes.csv")[,2])
# YS <- as.vector(read.csv("SizeYlw99SigSharedGenes.csv")[,2])

# #remove duplicates and write csv with genes exclusive to one pheno
# y <- c(yellow, YS, GY)
# y2 <- y[duplicated(y)]
# length(y2)
# y3 <- y[! y %in% y2]
# write.csv(y3,"yellowExclusive.csv")
# 
# g <- c(green, GS, GY)
# g2 <- g[duplicated(g)]
# g3 <- g[! g %in% g2]
# write.csv(g3, "greenExclusive.csv")
# 
# s <- c(size, GS, YS)
# s2 <- s[duplicated(s)]
# s3 <- s[! s %in% s2]
# write.csv(s3, "sizeExclusive.csv")

############### SEARCHING GENE LISTS FOR MATCHES WITH MUTANT LISTS ####################

allthree <- as.vector(shared$x[1:33])
ag <- as.vector(read.csv("AG99belowY95.csv")$x)
ay <- as.vector(read.csv("AY99belowG95.csv")$x)
gy <- as.vector(read.csv("GY99belowA95.csv")$x)
a <- as.vector(read.csv("areaExclusive.csv")$x)
y <- as.vector(read.csv("yellowExclusive.csv")$x)
g <- as.vector(read.csv("greenExclusive.csv")$x)

mich <- as.vector(read.csv("MutantsMichelle.csv")[,1])
baohua <- as.vector(read.csv("TDNA_Lines_Baohua_20160122.csv")[,1])
denby  <- as.vector(read.csv("Denby_Homozygous_Lines.csv")$AGI)

#gene jason validated
jason <- c("AT1G01060","AT1G14690","AT1G22070","AT1G31260","AT2G25540","AT2G32150","AT2G39940","AT3G15500","AT3G26830","AT3G52430","AT4G01860","AT4G01880","AT4G01883","AT4G14368","AT4G14400","AT4G14420","AT4G17010","AT4G39350")


#insert lists/vector you'd like to compare
all2 <- c(jason, genes)

#find duplicates
dupes <- all2[duplicated(all2)]
length(dupes)
dupes




