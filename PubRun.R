setwd("~/OneDrive - University College London/PhD/TropicH")
options(scipen = 999) # stop scientific notation globally (restore to default with options(scipen=0))
library(data.table)
# Load matrix of plot SADs in each continent with plot level attributes
africasadp <-read.csv('AfSADp.csv',header=T)  # Africa80p.csv OG, AfSADp.csv NEW POST-AFPD
amazonsadp <-read.csv('AmSADp.csv',header=T) # Amazon80p.csv OG, AmSADp.csv NEW POST-TNRS AND GERARDO
asiasadp <-read.csv('Asia80fulp1.csv',header=T) # now without BBT plot, but <0.9ha anyway so doesn't matter
#####DATA PREPARATION
####NOTE FOR THESE FUNCTIONS DON'T REMOVE INDET FROM DATASET BUT DO IN CORRECTION FUNCTIONS - CHECK WHICH!!!
amazon1 <- amazonsadp[amazonsadp$PlotSize>=0.9,] # cut out small plots
amazon2 <- amazon1[,5:ncol(amazon1)] # remove plot level attributes
amazon2 <- amazon2[,colSums(amazon2!=0)>0] # remove empty columns
africa1 <- africasadp[africasadp$PlotArea>=0.9,]
africa2 <- africa1[,5:ncol(africa1)]
africa2 <- africa2[,colSums(africa2!=0)>0] 
asia1 <- asiasadp[asiasadp$Size>=0.9,]
asia2 <- (asia1[,6:ncol(asia1)])
asia2 <- asia2[,colSums(asia2!=0)>0] 

#### RAREFACTION (JUST HYPERDOM 50%) 
# NOTE IDEALLY COMBINE RAREFACTION FUNCTION WITH RAREFACTIONPERC FUNCTION SO DON'T HAVE TO RUN EFFECTIVELY TWICE
# ONLY THING RAREFACTIONPERC IS MISSING OF NOTE IS FISHERS ALPHA SO JUST ADD THAT IN AND ALL GOOD
# ALTHOUGH SUBSEQUENT PLOTTING AND OTHER CODE WILL BE AFFECTED
asiar <- completeRare(asia2,repeats=200)
amazonr <- completeRare(amazon2,repeats=200)
africar <- completeRare(africa2,repeats=200)
write.csv(asiar,'asiar.csv',row.names = F)
write.csv(amazonr,'amazonr.csv',row.names = F)
write.csv(africar,'africar.csv',row.names = F)

## 50% RAREFACTION RESULTS
resasrar <- asiar[nrow(asiar),]
resamrar <- amazonr[which.min(abs(amazonr$Stems - sum(asia2))),] # sum(asia2) = 77587
resafrar <- africar[which.min(abs(africar$Stems - sum(asia2))),]
resrar <- rbind.data.frame('Africa'=resafrar,'Amazon'=resamrar,'Asia'=resasrar)
resrar$`H#Min` <- resrar$`H#`-1.96*resrar$`H#SD`
resrar$`H#Max` <- resrar$`H#`+1.96*resrar$`H#SD`
resrar$`H%Min` <- resrar$`H%`-1.96*resrar$`H%SD`
resrar$`H%Max` <- resrar$`H%`+1.96*resrar$`H%SD`
resrar$TSMin <- resrar$TS-1.96*resrar$TSSD
resrar$TSMax <- resrar$TS+1.96*resrar$TSSD
resrar$AlphaMin <- resrar$AlphaMLE-1.96*resrar$AlphaSE
resrar$AlphaMax <- resrar$AlphaMLE+1.96*resrar$AlphaSE
resrar$GammaMin <- resrar$GammaMLE-1.96*resrar$GammaSE
resrar$GammaMax <- resrar$GammaMLE+1.96*resrar$GammaSE
resrar0 <- resrar[c('H#','H#Min','H#Max','TS','TSMin','TSMax','H%','H%Min','H%Max','AlphaMLE','AlphaMin','AlphaMax','GammaMLE','GammaMin','GammaMax')]
resrarf <- rbind.data.frame(rep('Rarefaction Results',15),resrar0)

# RAREFACTION (10%-90%)
asiarperc <- completeRareperc(asia2,repeats=200)
amazonrperc <- completeRareperc(amazon2,repeats=200)
africarperc <- completeRareperc(africa2,repeats=200)
write.csv(asiarperc,'asiarperc.csv',row.names = F)
write.csv(amazonrperc,'amazonrperc.csv',row.names = F)
write.csv(africarperc,'africarperc.csv',row.names = F)
asiarperc <- read.csv('asiarperc.csv')
amazonrperc <- read.csv('amazonrperc.csv')
africarperc <- read.csv('africarperc.csv')

## 10-90% RAREFACTION RESULTS TABLE - PUB ED TABLES 2&3 PREP
resasrarp1 <- asiarperc[nrow(asiarperc),]
resasrarp1$h9SD <- 178.827725631529
resasrarp1$perc9SD <- 8.25788810223842
resamrarp1 <- amazonrperc[which.min(abs(amazonrperc$stems - sum(asia2))),] # sum(asia2) = 77587
resafrarp1 <- africarperc[which.min(abs(africarperc$stems - sum(asia2))),]
resrarp1 <- rbind.data.frame('Africa'=resafrarp1,'Amazon'=resamrarp1,'Asia'=resasrarp1)
write.csv(resrarp1,'resrarp1.csv',row.names = T)
# RESULTS COMPILED IN RESRARP1.XLSX WITH FORMATTING 

# EMPIRICAL RESULTS
africarperc[nrow(africarperc),]
amazonrperc[nrow(amazonrperc),]
asiarperc[nrow(asiarperc),]

#### HYPERDOM IDENTITY ANALYSIS - PUB SUP FIG 1,2,3
# custom divisions over which to draw samples by combining below in excel
# from IdentityPrep.xlsx
library(clipr) # allows use of 'read_clip_tbl' function which pastes data from clipboard copied from excel into DF
x1 <- read_clip_tbl(header=F)
x2 <- read_clip_tbl(header=F)
x3 <- read_clip_tbl(header=F)
# Africa hyperdominance subsampling
hypersub1 <- hyperSubProps(africa2,sequence = x1$V1) # hypersub11 before
write.csv(hypersub1,'afhypsub.csv')
# Amazon hyperdominance subsampling
hypersub2 <- hyperSubProps(amazon2,sequence = x2$V1) # hypersub21
write.csv(hypersub2,'amhypsub.csv')
# Asia hyperdominance subsampling
hypersub3 <- hyperSubProps(asia2,sequence = x3$V1) #hypersub31
write.csv(hypersub3,'ashypsub.csv')
# COMPILED INTO SUP123.XLSX


## EMPIRICAL AND UNCORRECTED EXTRAPOLATED 50% HYPERDOM RESULTS
# set row defining sampling level to which results are referring
samplevel <- c(rep('Empirical Results',9),rep('LS-Fitted Empirical Results',4),rep('LS-Extrapolated Results',11),rep('Pseudo-Corrected LS-Extrapolated Results 50%',3),rep('UnCorrected LS-Extrapolated Results 10-90%',60))
# run empirical and extrapolated master functions on each of the regional sets
afexh1 <- extrapH(africa2,91516984378) # 138206082133 total trees including non-sampled countries
amexh1 <- extrapH(amazon2,330935460291) # 405515779149
asexh1 <- extrapH(asia2,217394931190) # 254242453402
# make dataframe out of regional results and sampling identifier column
exh <- rbind.data.frame('Level'=samplevel,'Africa'=afexh1,'Amazon'=amexh1,'Asia'=asexh1)
names(exh) <- c('H#','TS','H%','Alpha','min','max','Stems','Plots','%ID','H#','TS','H%','Stems','H#','min','max','TS','min','max','H%','min','max','Total Stems (LS)','Total Stems (Crowther)','Fit Error','H#','H%','Species','h1','h2','h3','h4','h5','h6','h7','h8','h9','perc1','perc2','perc3','perc4','perc5','perc6','perc7','perc8','perc9','stems','SpeciesMin','h1Min','h2Min','h3Min','h4Min','h5Min','h6Min','h7Min','h8Min','h9Min','perc1Min','perc2Min','perc3Min','perc4Min','perc5Min','perc6Min','perc7Min','perc8Min','perc9Min','stemsMin','SpeciesMax','h1Max','h2Max','h3Max','h4Max','h5Max','h6Max','h7Max','h8Max','h9Max','perc1Max','perc2Max','perc3Max','perc4Max','perc5Max','perc6Max','perc7Max','perc8Max','perc9Max','stemsMax')

# below version compatible with extrapolation correction
exh1 <- rbind.data.frame('Africa'=afexh1,'Amazon'=amexh1,'Asia'=asexh1)
names(exh1) <- c('H#','TS','H%','Alpha','min','max','Stems','Plots','%ID','H#','TS','H%','Stems','H.mean','H.min','H.max','S.mean','S.min','S.max','Prop.mean%','Prop.min','Prop.max','Total Stems (LS)','Tot.t','Fit Error','H#','H%','Species','h1','h2','h3','h4','h5','h6','h7','h8','h9','perc1','perc2','perc3','perc4','perc5','perc6','perc7','perc8','perc9','stems','SpeciesMin','h1Min','h2Min','h3Min','h4Min','h5Min','h6Min','h7Min','h8Min','h9Min','perc1Min','perc2Min','perc3Min','perc4Min','perc5Min','perc6Min','perc7Min','perc8Min','perc9Min','stemsMin','SpeciesMax','h1Max','h2Max','h3Max','h4Max','h5Max','h6Max','h7Max','h8Max','h9Max','perc1Max','perc2Max','perc3Max','perc4Max','perc5Max','perc6Max','perc7Max','perc8Max','perc9Max','stemsMax')
exh1$Tot.a <- c(184930623,739643061,542157318) # add column delineating forest region areas (just area of WWF ecoregion with at least one plot in data)
exh1$Region <- c('Africa', 'Amazonia', 'Southeast Asia')
exhp <- exh1[c('Region','H.mean','H.min','H.max','S.mean','S.min','S.max','Prop.mean%','Prop.min','Prop.max','Tot.t','Tot.a','h1','h2','h3','h4','h5','h6','h7','h8','h9','perc1','perc2','perc3','perc4','perc5','perc6','perc7','perc8','perc9','h1Min','h2Min','h3Min','h4Min','h5Min','h6Min','h7Min','h8Min','h9Min','perc1Min','perc2Min','perc3Min','perc4Min','perc5Min','perc6Min','perc7Min','perc8Min','perc9Min','h1Max','h2Max','h3Max','h4Max','h5Max','h6Max','h7Max','h8Max','h9Max','perc1Max','perc2Max','perc3Max','perc4Max','perc5Max','perc6Max','perc7Max','perc8Max','perc9Max')]


####NOTE FOR THIS FUNCTION HAVE TO REMOVE INDET FROM DATASET AND PREDATASET, DO NOT REMOVE THEM FROM NON-CORRECTION FUNCTIONS
# delete Indet species from fittings to get 'predataset's
africai <- subset( africa1, select = -Indet )
amazoni <- subset( amazon1, select = -Indet )
asiai <- subset( asia1, select = -c(Indet,Dataset) )# delete unnecessary FP v Slik dataset column so same as other sets
# # Make dataframes in correct format for Prado scripts
africa <- data.frame('species' = names(africai[,5:length(africai)]),'N.ind' = colSums(africai[,5:length(africai)]),'N.plots' = colSums(africai[,5:length(africai)]!=0))
amazon <- data.frame('species' = names(amazoni[,5:length(amazoni)]),'N.ind' = colSums(amazoni[,5:length(amazoni)]),'N.plots' = colSums(amazoni[,5:length(amazoni)]!=0))
asia <- data.frame('species' = names(asiai[,5:length(asiai)]),'N.ind' = colSums(asiai[,5:length(asiai)]),'N.plots' = colSums(asiai[,5:length(asiai)]!=0))

afextrap <- extrapPradoH(region = 'Africa',predataset = africai, dataset = africa, uncex = exhp, minSpecies = 2e3,maxSpecies = 1e4,repetitions = 250)
amextrap <- extrapPradoH(region = 'Amazonia',predataset = amazoni,dataset = amazon, uncex = exhp, minSpecies = 1e4,maxSpecies = 2.5e4,repetitions = 250)
asextrap <- extrapPradoH(region = 'Southeast Asia',predataset = asiai, dataset = asia, uncex = exhp,minSpecies = 1e4,maxSpecies = 2.5e4,repetitions = 250)
extrap <- rbind.data.frame('Africa'=afextrap,'Amazonia'=amextrap,'Southeast Asia'=asextrap)
names(extrap) <- c('H#','H#Min','H#Max','H1','H1Min','H1Max','H2','H2Min','H2Max','H3','H3Min','H3Max','H4','H4Min','H4Max','H5','H5Min','H5Max','H6','H6Min','H6Max','H7','H7Min','H7Max','H8','H8Min','H8Max','H9','H9Min','H9Max','TS','TSMin','TSMax','H%','H%Min','H%Max','H%1','H%1Min','H%1Max','H%2','H%2Min','H%2Max','H%3','H%3Min','H%3Max','H%4','H%4Min','H%4Max','H%5','H%5Min','H%5Max','H%6','H%6Min','H%6Max','H%7','H%7Min','H%7Max','H%8','H%8Min','H%8Max','H%9','H%9Min','H%9Max')
extrap[subset(names(extrap),!grepl('Min|Max',names(extrap)))]

## 50% HYPERDOM RAREFIED, EMPIRICAL AND CORRECTED EXTRAPOLATED RESULTS TABLE
# DATA FOR FIG1 PIE CHARTS RAREFY AND EXTRAP - TO PYTHON 
# PUB TABLE 1 RAREFIED 50%HYPERDOM
# Combine rarefied and empirical, extrapolated results together
reesults <- cbind.data.frame(resrarf,rbind.data.frame('Level' = rep('Corrected Extrapolated'), extrap))
write.csv(reesults,'reesults.csv')
reesults <- read.csv('reesults.csv')
extrap <- reesults[2:4,17:79]
names(extrap) <- c('H#','H#Min','H#Max','H1','H1Min','H1Max','H2','H2Min','H2Max','H3','H3Min','H3Max','H4','H4Min','H4Max','H5','H5Min','H5Max','H6','H6Min','H6Max','H7','H7Min','H7Max','H8','H8Min','H8Max','H9','H9Min','H9Max','TS','TSMin','TSMax','H%','H%Min','H%Max','H%1','H%1Min','H%1Max','H%2','H%2Min','H%2Max','H%3','H%3Min','H%3Max','H%4','H%4Min','H%4Max','H%5','H%5Min','H%5Max','H%6','H%6Min','H%6Max','H%7','H%7Min','H%7Max','H%8','H%8Min','H%8Max','H%9','H%9Min','H%9Max')
afextrap <- as.numeric(as.character(unname(unlist(extrap[1,]))))
amextrap <- as.numeric(as.character(unname(unlist(extrap[2,]))))
asextrap <- as.numeric(as.character(unname(unlist(extrap[3,]))))

##### PUB FIG 3: plot H% against H% threshold for 10-90% definition of H% at the rarefied level of Asia dataset size and extrapolated to regional level
png('Fig31.png',res = 150, width = 1180, height = 1080)
plot(seq(9.5,89.5,by=10),resasrarp1[12:20], col = '#0000FF', pch=19,ylim= c(0,90),xlab='Dominance Threshold (Percentage of Stems Accounted For)', ylab ='Dominant Percentage (Dominant Proportion of Total Species)')
points(seq(10.5,90.5,by=10),resamrarp1[12:20],col = '#00FFFF',pch=19)
points(seq(10,90,by=10),resafrarp1[12:20],col = '#FF00FF',pch=19) # EEDD62 #EFE255
arrows(seq(9.5,89.5,by=10), as.numeric(resasrarp1[12:20]-1.96*resasrarp1[32:40]), seq(9.5,89.5,by=10), as.numeric(resasrarp1[12:20]+1.96*resasrarp1[32:40]), length=0.05, angle=90, code=3,col='#0000FF')
arrows(seq(10.5,90.5,by=10), as.numeric(resamrarp1[12:20]-1.96*resamrarp1[32:40]), seq(10.5,90.5,by=10), as.numeric(resamrarp1[12:20]+1.96*resamrarp1[32:40]), length=0.05, angle=90, code=3,col = '#00FFFF')
arrows(seq(10,90,by=10), as.numeric(resafrarp1[12:20]-1.96*resafrarp1[32:40]), seq(10,90,by=10), as.numeric(resafrarp1[12:20]+1.96*resafrarp1[32:40]), length=0.05, angle=90, code=3,col='#FF00FF')
# extrapolated
points(seq(9.5,89.5,by=10),asextrap[seq(37,62,3)],pch=18, col = '#0000FF')
points(seq(10.5,90.5,by=10),amextrap[seq(37,62,3)],pch=18,col='#00FFFF')
points(seq(10,90,by=10),afextrap[seq(37,62,3)],pch=18,col = '#FF00FF')
# arrows(seq(9.5,89.5,by=10), as.numeric(asextrap[seq(38,63,3)]), seq(9.5,89.5,by=10), as.numeric(asextrap[seq(39,64,3)]), length=0.05, angle=90, code=3,col='#0000FF')
# arrows(seq(10.5,90.5,by=10), as.numeric(amextrap[seq(38,63,3)]), seq(10.5,90.5,by=10), as.numeric(amextrap[seq(39,64,3)]), length=0.05, angle=90, code=3,col = '#00FFFF')
# arrows(seq(10,90,by=10), as.numeric(afextrap[seq(38,63,3)]), seq(10,90,by=10), as.numeric(afextrap[seq(39,64,3)]), length=0.05, angle=90, code=3,col='#FF00FF')
axis(1, at = as.numeric(seq(10, 90, by = 10)))
# https://r-charts.com/base-r/grid/
axis(2, at = as.numeric(seq(0, 90, by = 10)),tck=1,lty = 2, col = "gray")
legend("topleft",c("Amazonia Rarefied", "Africa Rarefied","Southeast Asia Rarefied","Amazonia Extrapolated", "Africa Extrapolated","Southeast Asia Extrapolated"), col=c("#00FFFF","#FF00FF","#0000FF",'#00FFFF','#FF00FF',"#0000FF"),pch = c(19,19,19,18,18,18),bty = 'n')
dev.off()

# # BELOW NOT RUN AS ASIA DATASET NOT CHANGED SO UNNECESSARY
# # PUB ED FIG 1: Total Species
# png('PubTSAsia.png',res = 200, width = 1180, height = 1080)
# plot(asiar$Stems,asiar$TS,col = 'red', xlim=c(0,100000),xlab='Number of Stems',ylab='Total Number of Species')
# points(asiar2$Stems,asiar2$TS,col = 'purple')
# legend("topleft",c("Asia All",'Asia >0.9ha'),lty=1, col=c("red",'purple'))
# dev.off()
# # PUB ED FIG 1: No. Hyperdominants
# png('PubHAsia.png',res = 200, width = 1180, height = 1080)
# plot(asiar$Stems,asiar$`H#`,col = 'red', xlim=c(0,100000),xlab='Number of Stems',ylab='Number of Hyperdominants')
# points(asiar2$Stems,asiar2$`H#`,col = 'purple')
# legend("topleft",c("Asia All",'Asia >0.9ha'),lty=1, col=c("red",'purple'))
# dev.off()

# PUB ED FIG 4: RADs and Preston Plot LS fits
# Africa output
png('PrestonAf.png',res = 200, width = 1480, height = 1080)
fitPlot(africa2,color='#FF00FF')
dev.off()
png('RADAf.png',res = 200, width = 1480, height = 1080)
fitPlot1(africa2,color='#FF00FF')
dev.off()
# Amazon output
png('PrestonAm.png',res = 200, width = 1480, height = 1080)
fitPlot(amazon2,color='#00FFFF')
dev.off()
png('RADAm.png',res = 200, width = 1480, height = 1080)
fitPlot1(amazon2,color='#00FFFF')
dev.off()
# Asia output
png('PrestonAs.png',res = 200, width = 1480, height = 1080)
fitPlot(asia2,color='#0000FF')
dev.off()
png('RADAs.png',res = 200, width = 1480, height = 1080)
fitPlot1(asia2,color='#0000FF')
dev.off()

aff <- read.csv('AffiliCompilation.csv')
as.character(aff$AffiliationNos)[1]
