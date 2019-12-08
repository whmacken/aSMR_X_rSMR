####This script contains the working steps for estimating SI values. This includes modelling aSMR from rSMR,
####building polynomial models of standardised SI over edatopic grids, and predicting max SI based on temperature data
###Kiri Daust, July 2018

###Model aSMR <-> rSMR
.libPaths("E:/R packages351")
rm(list=ls())
require(tcltk)
require(foreach)
require(ggplot2)
require (plyr)
require(dplyr) ## this package may cause issues with count function of plyr
require(magrittr)
require(data.table)
require (reshape)
require(reshape2)
require(gridExtra)
require(mgcv)
require(rgl)
require(deldir)
require(stringr)
require(janitor)
require(modeest)
require(expss)

####Function to create 3d barplots from Stack Overflow####################
binplot.3d <- function(x, y, z, alpha=1, topcol="#ff0000", sidecol="#aaaaaa", linecol="#000000")
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  x1 <- c(rep(c(x[1], x[2], x[2], x[1]), 3), rep(x[1], 4), rep(x[2], 4))
  z1 <- c(rep(0, 4), rep(c(0, 0, z, z), 4))
  y1 <- c(y[1], y[1], y[2], y[2], rep(y[1], 4), rep(y[2], 4), rep(c(y[1], y[2], y[2], y[1]), 2))
  x2 <- c(rep(c(x[1], x[1], x[2], x[2]), 2), rep(c(x[1], x[2], rep(x[1], 3), rep(x[2], 3)), 2))
  z2 <- c(rep(c(0, z), 4), rep(0, 8), rep(z, 8))
  y2 <- c(rep(y[1], 4), rep(y[2], 4), rep(c(rep(y[1], 3), rep(y[2], 3), y[1], y[2]), 2))
  rgl.quads(x1, z1, y1, col=rep(sidecol, each=4), alpha=alpha)
  rgl.quads(c(x[1], x[2], x[2], x[1]), rep(z, 4), c(y[1], y[1], y[2], y[2]), col=rep(topcol, each=4), alpha=1) 
  rgl.lines(x2, z2, y2, col=linecol)
}

barplot3d <- function(z, alpha=1, col="#ff0000", scale=1)
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  z <- as.matrix(z)
  xy <- dim(z)
  x <- seq(xy[1])
  y <- seq(xy[2])
  z <- z / max(z, na.rm=TRUE) * max(x, y) * scale
  for (i in x) 
  {
    for (j in y) 
    {
      binplot.3d(c(i, i+1), c(j, j+1), z[i,j], alpha=alpha, topcol=col)
    }
  }
}
wireframe(zfit)

layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9)),widths=c(1,1,1), heights =c(1,1,1), respect=TRUE)
par(mai = c(0.5, 0.5, 0.2, 0.2)) #speSEfies the margin size in inches
######### change names function
changeNames  <- function(x, old, new){
  result <- vector("numeric", length(x))
  for (i in 1:length(x)){
    code <- x[i]
    index <- match(code, old)
    result[i] <- new[index]
  }
  return(result)
}


###Create rSMR -> aSMR crosswalk######################
wd <- tk_choose.dir(); setwd(wd)

####This part includes removing monthly CMD based on surplus (doesn't make much difference)#####
allDat <- fread("BECv11_100Pt_Normal_1961_1990MSY.csv", data.table = FALSE) ###Climate data 1961-90for representative points per BGC 
allDat <- fread("AlbertaSNR_500Pt_Normal_1961_1990MSY.csv", data.table = FALSE) ###Climate data 1961-90for representative points per BGC 
allDat <- fread("USA_July12_TrainingPts_2_Normal_1961_1990MSY.csv", data.table = FALSE) ###Climate data 1961-90for representative points per BGC 

temp <- allDat[,grep("CMD",colnames(allDat))]
###work on checking for excess summer moister leading to Moist zonal sites
#temp2 <- allDat[,grep("PPT0", colnames(allDat))]
temp2 <- allDat[c("ID1", "ID2", "PPT05","PPT06","PPT07","PPT08","PPT09" )]
temp2$Moist <- rowSums(temp2[,3:7] >90)
temp2$Moist <- ifelse(temp2$Moist > 4, "6",
                      ifelse(temp2$Moist >3, "5.5",""))
temp2 <- temp2[c("ID1","ID2","Moist")]
temp2 <- temp2[temp2$Moist != "",] # only those points calculated to be moist
#####################################
allDat <- allDat[,c("ID2","PPT_at","PPT_wt","PAS")]
allDat <- cbind(allDat,temp)
#allDat <- cbind(allDat,temp2[c("Moist")])
#allDat <- cbind(allDat,temp,temp2) ### use this one if Moist indicator is completed above
allDat$PPT.dorm <- allDat$PPT_at + allDat$PPT_wt
CMD <- aggregate( . ~ ID2, allDat, mean) ##
Moist <- aggregate(ID1 ~ ID2 + Moist, temp2, length) # count of number of sites in each BGC meeting moist criteria
BGC_Special <- fread("BGC_Special.csv", data.table = FALSE)
CMD <- merge(CMD,BGC_Special, by.x = "ID2")

####To adjust in zones with prolonged snowpack
CMD$CMDKiri <- ifelse(CMD$Special == "snow", CMD$CMD07+CMD$CMD08+CMD$CMD09, 
                 CMD$CMD02+CMD$CMD03+CMD$CMD04+CMD$CMD05+CMD$CMD06+CMD$CMD07+CMD$CMD08+CMD$CMD09)
CMD$CMD <- CMD$CMDKiri

CMD <- CMD[,c("ID2","CMD","PPT.dorm")]
CMD$Def <- 300 - CMD$PPT.dorm ###adds deficit from incomplete recharge in dormant season. 
#----------------------
Adjust <- fread("CMD_adjustments.csv", data.table=getOption("datatable.fread.datatable", TRUE))## temporary adjustment of CMD based on vegetation analysis
CMD <- merge(CMD,Adjust, by.x = "ID2")
CMD$CMD <- CMD$CMD + CMD$CMD_Adjust
#----------------------
CMD$Def[CMD$Def < 0] <- 0
CMD$CMD <- CMD$CMD + CMD$Def
CMD <- CMD[,c("ID2","CMD")]

###_____________________________________________####
###Option with just with ppt#######
allDat <- fread("AlbertaSNR_500Pt_Normal_1961_1990MSY.csv", data.table = FALSE)
allDat <- allDat[,c("ID2","PPT_at","PPT_wt","CMD")]
allDat$PPT.dorm <- allDat$PPT_at + allDat$PPT_wt
CMD <- aggregate(cbind(PPT.dorm, CMD) ~ ID2, allDat, mean)###Mean by BGC
CMD$Def <- 300 - CMD$PPT.dorm ###adds deficit from incomplete recharge in dormant season. 
CMD$Def[CMD$Def < 0] <- 0
CMD$CMD <- CMD$CMD + CMD$Def
CMD <- CMD[,c("ID2","CMD")]
###_____________________________________________####

###for each wetter rSMR, previous CMD is divided by 2
for (i in 1:3){
  CMD[,2+i] <- CMD[,1+i]/2
}
colnames(CMD) <- c("BGC","rSMR4","rSMR5","rSMR6","rSMR7")
CMD <- CMD[,c(1,3:5,2)]

###for each drier rSMR, previous CMD + 100
for (i in 1:4){
  CMD[,length(CMD)+1] <- CMD[,length(CMD)] + 100
}
colnames(CMD)[6:9] <- c("rSMR3","rSMR2","rSMR1","rSMR0")

CMD <- CMD[,order(colnames(CMD))]## creates full grid of CMD values by BGC by rSMR class

#####________Creates boxplot of CMD ranges by expert rSMR4 for initial model and ID of outliers______________#
#expGrid <- read.csv("ExpertGrid.csv")
expGrid <- read.csv("ExpertGrid_adjusted.csv")
expGrid$BGC <- gsub("[[:space:]]","",expGrid$BGC)
colnames(expGrid)[-1] <- paste(colnames(expGrid)[-1],"_E", sep = "")
CMDrange <- merge(CMD, expGrid, by = "BGC")
CMDrange$rSMR4_E <- as.factor(CMDrange$rSMR4_E)
p <- ggplot(CMDrange, aes(rSMR4_E, rSMR4))+
  geom_boxplot(stat = "boxplot",  varwidth=TRUE) +
  geom_point(shape = 21, fill = "red", size = 1)+
  #geom_jitter(width = 0) +
  xlab ("aSMR Class")+
  ylab ("CMD")+
  geom_text (aes(label= BGC), size = 1, hjust = -1, position = position_dodge(width=1) )#
ggsave("CMD_ranges.pdf", plot = p, dpi = "print", device = "pdf", width = 15, height = 15, units = "cm")
ggsave("CMD_ranges.png", plot = p, dpi = "print", device = "png", width = 15, height = 15, units = "cm")
print(p)
#______________________________________________________________________________##

#############################################Now calculate aSMR with ruleset
rules <- read.csv("aSMR_Rules_FullStep.csv") ### import rules on range of CMD which equates to aSMR
rules <- read.csv("aSMR_Rules_HalfStep.csv") ### import rules on range of CMD which equates to aSMR
aSMRClass <- function(x){
  for(i in 1:length(ruleSelect$CMD)){
    if(x < ruleSelect$CMD[i+1]){
      x <- ruleSelect$aSMR[i]
      break
    }
  }
  return(x)
}

###Calculate values based on rules###
test <- foreach(SMR = colnames(CMD)[-1], .combine = cbind) %do% {
  temp <- CMD[,SMR]
  if(SMR == "rSMR7"){
    ruleSelect <- rules[rules$SMRLevel == 7,-1]
  }else if(SMR == "rSMR6"){
    ruleSelect <- rules[rules$SMRLevel == 6,-1]
  }else if(SMR == "rSMR5"){
    ruleSelect <- rules[rules$SMRLevel == 5,-1]
  }else{
    ruleSelect <- rules[rules$SMRLevel == 0,-1]
  }
  out <- sapply(temp,FUN = aSMRClass)
  out
}


test <- as.data.frame(test)
test <- cbind(CMD$BGC, test)
colnames(test) <- colnames(CMD)
test$BGC <- gsub("[[:space:]]","",test$BGC)
SMRCross <- melt(test) ###aSMR lookup
colnames(SMRCross) <- c("BGC", "rSMR", "aSMRC")
write.csv(test, file= "modelled_USAv11_rSMR_aSMR_grid_HalfStep.csv")
save(SMRCross,file = "rSMR_aSMR_CalcList.RData")
load ("rSMR_aSMR_CalcList.RData")

###Code to test quality of aSMR crosswalk table as compared to Expert Grid
### Jump to next section for comparison to vegetation
expGrid <- read.csv("ExpertGrid_adjusted.csv") ## or original
expGrid$BGC <- gsub("[[:space:]]","",expGrid$BGC)
rSMRListexp <- melt(expGrid, id.vars= c("BGC"), variable_name = "rSMR", value.name = "aSMRE")
colnames(rSMRListexp) <- c("BGC", "rSMR", "aSMRE")
rSMR_exp_calc <- merge(SMRCross,rSMRListexp, by = c("BGC", "rSMR")  )
aSMRdiff <- rSMR_exp_calc
aSMRdiff$diff <- aSMRdiff$aSMRE - aSMRdiff$aSMRC
aSMRDiffCount2 <- aggregate(BGC ~ diff, aSMRdiff, length)
aSMRdiffmat <- cast(aSMRdiff, Species ~ aSMRC) #ignore message here casts to matrix
d <- ggplot(aSMRdiff, aes(diff))+
  geom_histogram(binwidth = .25, fill = "dodgerblue3") +
  xlab ("Difference Expert-Modelled aSMR Class")+
  ylab ("Count")
plot(d)
ggsave("aSMR difference between expert and climate model2.jpeg", plot = d, dpi = "print", device = "jpeg", width = 6, height = 3, units = "in")
dev.off() 
save(SMRCross,file = "rSMR_aSMR_Both.RData")
colnames(expGrid)[-1] <- paste(colnames(expGrid)[-1],"_E", sep = "")
comp <- merge(expGrid, test, by = "BGC") 

for(i in 1:8){
  comp[,length(comp) + 1] <- comp[,i+1] - comp[,i+9]
}

comp <- comp[,-(2:17)] ###outputs matrix showing difference between expert and modelled aSMR
colnames(comp) <- colnames(test)
hist(comp$SMR4, col = "purple")

for(i in 1:8){
  hist(comp[,i+1], col = "purple", main = colnames(current)[i+1])
}


comp$Zone <- gsub("[[:lower:]]|[[:digit:]]","",comp$BGC) # adds Zone to comparison
zoneQual <- aggregate(rSMR4 ~ Zone, comp, mean)
barplot(zoneQual$rSMR4, names.arg = zoneQual$Zone, width=4) # bar plot of difference between expert and modelled by zone

compLong <- melt(comp[,-10])
colnames(compLong) <- c("BGC","rSMR","Diff")
write.csv(compLong, file= "rSMR_aSMR_differencelist_halfstep.csv") # output for overall stats
compLong <- compLong[abs(compLong$Diff) >= 1,] # set difference level to report below
temp <- melt(expGrid)
colnames(temp) <- c("BGC","rSMR","Expert")
temp$rSMR <- gsub("r","",temp$rSMR)
compLong <- merge(compLong, temp, by = c("BGC","rSMR"), all.x = TRUE) ### list of 1 or greater differences
comp2 <- comp
write.csv(comp2, file= "compared_rSMR_aSMR_grid_halfstep300rules2b100.csv")
comp2[comp2 == 0.5 | comp2 == -0.5] <- 0 # set differences of .5 to zero
write.csv(comp2, file= "compared_rSMR_aSMR_grid_fullstep.csv")
### outputs PDFs barcharts by zone comparing difference
for(j in unique(comp$Zone)){
  pdf(paste(j,"fullstep.pdf", sep = ""))
  current <- comp2[comp2$Zone == j,]
  layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9)),widths=c(1,1,1), heights =c(1,1,1), respect=TRUE)
  par(mai = c(0.5, 0.5, 0.2, 0.2)) #speSEfies the margin size in inches
  for(i in 1:8){
    hist(current[,i+1], breaks = 7, col = "purple", main = colnames(current)[i+1])
  }
  dev.off()
}

###Output for AllBGCs combined - set comp2 above for either half step or full step
#for(j in unique(comp$Zone)){
  pdf("AllBGCHalfstep.pdf")
  current <- comp
  layout = layout  #  layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9)),widths=c(1,1,1), heights =c(1,1,1), respect=TRUE)
  par = par #(mai = c(0.5, 0.5, 0.2, 0.2)) #speSEfies the margin size in inches
   for(i in 1:8){
    hist(current[,i+1], breaks = 7, col = "purple", main = colnames(current)[i+1])
  }
  dev.off() 

  ##################Calculating the estimated aSMR value of species in the BECMaster
  ####Create Plots BGC & edatopic position data frame
  load("rSMR_aSMR_Both.RData")
    load("VegDat_Clean.RData")
 PlotCount <- count (vegData, c("PlotNumber"))
 #PlotCount2 <- aggregate(Species ~ PlotNumber, vegData, length)
 PlotCount3 <- as.data.frame (unique(vegData$PlotNumber)) 
 #vegData <- read.table("BecMaster15VegData.txt", header = TRUE)
  
  codeCross <- read.csv("TreeCodeCrosswalk.csv", stringsAsFactors = FALSE)
  
  #vegData <- separate(vegData, Species, c("Species","Type"), "-", remove = TRUE)
  #---------Option for Trees only
  #  vegData <- vegData[vegData$Type %in% c(1,2),]
  #--------------------------
  BGCLookup <- read.csv("BGCLookup.csv", stringsAsFactors = FALSE)
  BGCLookup <- BGCLookup[,c(3,12)]
  colnames(BGCLookup)[1] <- "PlotNumber"
  
  ##import edatopic data
  plotEnv <- read.csv("KiriEnvDat.csv", stringsAsFactors = FALSE)# generated from Vpro BECMaster plot|BGC|SMR|SNR
  plotEnv <- plotEnv[plotEnv$NutrientRegime %in% c("A","B","C","D","E"),]
  plotEnv <- plotEnv[plotEnv$MoistureRegime %in% c(0,1,2,3,4,5,6,7,8),]
  plotEnv <- plotEnv[,-2]
  plotEnv <- merge(plotEnv,BGCLookup, by = "PlotNumber", all.x = TRUE)
  plotEnv <- plotEnv[plotEnv$BGC_LABEL != "",]
  plotEnv <- plotEnv[plotEnv$BGC_LABEL != " ",]
  plotEnv$BGC_noSpace <- gsub("[[:space:]]","",plotEnv$BGC_LABEL)
  plotEnv$BGC_Zone <- substring(plotEnv$BGC_LABEL, 1, 4)
  colnames (plotEnv) [2] <- "rSMR"
###load the aSMR to rSMR crosswalk
load( "rSMR_aSMR_Both.RData")
colnames(rSMR_exp_calc)[1] <- "BGC_noSpace"
rSMR_exp_calc$rSMR <- gsub("rSMR", "", rSMR_exp_calc$rSMR)
#rSMR_exp_calc[1:4] <- lapply(rSMR_exp_calc, factor)
#    SMRCross$rSMR <- gsub("rSMR", "", SMRCross$rSMR)
# colnames(SMRCross) [1:3] <- c("BGC_noSpace", "rSMR", "aSMR")
Env2 <- merge(plotEnv, rSMR_exp_calc, by = c("BGC_noSpace", "rSMR"),  all.x = TRUE)
plotEnv2 <- Env2[, c(3,1,5,6,2,7,8)] # colnames(plotEnv)[4] <- "Unit"
 #   modBGC <- read.csv("ModelledBGC.csv", stringsAsFactors = FALSE)
##Count the number of plots in each aSMR category
aSMRCount <- count (plotEnv2, c("aSMRC")) #count of species occurrences by modelled aSMR 1/2 step
    save(plotEnv2,file = "PlotEnv.RData")
    load( "PlotEnv.RData")
  ############Add BGC and aSMR info to vegetation list data
Veg_Env <- merge(vegData, plotEnv2, by = "PlotNumber", all.x = TRUE)
Veg_Env <- Veg_Env[!is.na(Veg_Env$aSMRC),]
Veg_Env <- Veg_Env[Veg_Env$BGC_LABEL != "",]
#Veg_Env$aSMR <- as.factor(Veg_Env$aSMRC)
#Veg_Env2 <- group_by(Veg_Env, Species, aSMR)
Veg_Env3 <- plyr::count (Veg_Env, c("Species", "aSMRC")) #count of species occurrences by modelled aSMR 1/2 step
Veg_Env4 <- aggregate(Cover ~ Species + aSMRC, Veg_Env, sum) #sum of covers of species occurrences by modelled aSMR 1/2 step
#Veg_Env3 <- table (Veg_Env, c("Species", "aSMR"), show_na = FALSE)
############summarize species occurrences by aSMR category for determining indicator values
Spp_aSMR <- cast(Veg_Env3, Species ~ aSMRC) #ignore message here casts to matrix
Spp_aSMRcov <- cast(Veg_Env4, Species ~ aSMRC) #ignore message here casts to matrix

write.csv(Spp_aSMR, file= "Species_by_aSMR_grid_halfstep.csv")
write.csv(Spp_aSMRcov, file= "Speciescover_by_aSMR_grid_halfstep.csv")

#### build some summaries to look at species distribution by calculated aSMR
###Start Kiri##############################
####calculate percent
Spp_CountPerc <- Spp_aSMR[,colnames(Spp_aSMR) != "NA"] ##remove NA
Spp_CountPerc$Total <- rowSums(Spp_CountPerc[,-1])
Spp_CountPerc[,!colnames(Spp_CountPerc) %in% c("Species","Total")] <- (Spp_CountPerc[,!colnames(Spp_CountPerc) %in% c("Species","Total")]/Spp_CountPerc$Total)*100

###Now for %cover
Spp_aSMRcov[is.na(Spp_aSMRcov)] <- 0 ####fill NA with 0
Spp_CovPerc <- Spp_aSMRcov
Spp_CovPerc$Total <- rowSums(Spp_CovPerc[,-1], na.rm = TRUE)
Spp_CovPerc[,!colnames(Spp_CountPerc) %in% c("Species","Total")] <- (Spp_CovPerc[,!colnames(Spp_CountPerc) %in% c("Species","Total")]/Spp_CovPerc$Total)*100

####Max, Min and Median values
#install.packages("modeest") ##package to calculate mode
#library(modeest)

temp <- melt(Spp_aSMRcov, id.vars = "Species")###to long format
tempCount <- Spp_aSMR

#tempCount <- melt(tempCount, id.vars = "Species")
tempCount <- tempCount[!is.na(tempCount$Species),]##remove NA species
tempCount[is.na(tempCount)] <- 0
tempCount [-1][tempCount[-1] < 6] <- 0###remove aSMR counts with few occurrences
tempCount <-cbind(tempCount,rowSums(tempCount[,-c(1)]))
colnames(tempCount) [18] <- c("Sum")
tempCount <- tempCount[tempCount$Sum >9,]###remove species with few occurrences
SppList <- as.character(unique(tempCount$Species))

####Create stats for each species
stats <- foreach(Spp = SppList, .combine = rbind) %do% {
  SppSub <- temp[temp$Species == Spp,]###subset for species
  mode <- SppSub$aSMRC[SppSub$value == max(SppSub$value)] ###Mode is aSMR at max value
  SppSub$Integral <- cumsum(SppSub$value) ###cumulative sum
  SppSub$StandValue <- (SppSub$value/sum(SppSub$value))*1000 ###standardise value
  q1 <- 0.1*sum(SppSub$value)
  q2 <- 0.9*sum(SppSub$value)
  q10 <- SppSub$aSMRC[SppSub$Integral >= q1][1] ###aSMR value at 10% percentile
  q90 <- SppSub$aSMRC[SppSub$Integral >= q2][1]###and 90%
  
  ###Caclulate standard deviation by repeating aSMR value by standardised cover value
  s2 <- sd(c(rep(0.5,SppSub$StandValue[1]),rep(1,SppSub$StandValue[2]),rep(1.5,SppSub$StandValue[3]),rep(2,SppSub$StandValue[4]),
           rep(2.5,SppSub$StandValue[5]),rep(3,SppSub$StandValue[6]),rep(3.5,SppSub$StandValue[7]),rep(4,SppSub$StandValue[8]),
           rep(4.5,SppSub$StandValue[9]),rep(5,SppSub$StandValue[10]),rep(5.5,SppSub$StandValue[11]),rep(6,SppSub$StandValue[12]),
           rep(6.5,SppSub$StandValue[13]),rep(7,SppSub$StandValue[14]),rep(7.5,SppSub$StandValue[15]),rep(8,SppSub$StandValue[16])))

  CountSub <- tempCount[tempCount$Species == Spp,]
  numPlots <- sum(CountSub$value)
  
  out <- data.frame(Species = Spp, Q10 = q10, Mode = mode, Q90 = q90, SD = s2, NumPlots = CountSub$Sum)
  out
}

write.csv(stats, "Species_aSMR_Statistics.csv")

###Build Graphics
choices <- as.character(tempCount$Species[rowSums(tempCount[,-1]) > 200]) ###only allow selection of species with certain abundance cutoff
Spp <- select.list(choices, graphics = TRUE, multiple = TRUE) ###Select which species to plot (can choose multiple)
#Spp <- c("ABIEAMA", "ABIEGRA","ABIELAS","LARIOCC", "PINUCON", "PINUPON","PSEUMEN","PICEENE","THUJPLI","TSUGHET")
Spp <- c("PETAFRI", "MAIADIL","OPLOHOR")#Moist5
Spp <- c("LYSIAME", "MENYTRI",  "SPHAFUS")
temp <- as.data.frame(Spp_aSMRcov[Spp_aSMRcov$Species %in% Spp,])##subset for selected species
rownames(temp) <- temp$Species
temp <- temp[,-1]
temp <- as.data.frame(t(temp))
temp$aSMR <- rownames(temp)

multSpp <- melt(temp, id.vars = "aSMR")###Create data set for doing multiple curves at once
colnames(multSpp)[2:3] <- c("Species","Cover")

temp$Total <- rowSums(temp[,-length(temp)]) ###Sum together if multiple species selected


####Create bar plot of cover by aSMR
ggplot(temp, aes(x = aSMR, y = Total))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(y = "Cover", title = paste(Spp, collapse = ", "))

####fit gaussian spp curves and plot on same graph######################
pdf("TreeSpeciesaSMR_curves.pdf")
jpeg("TreeSpeciesaSMR_curves.jpeg", width = 700, height = 350)
plot(0,0, type = "n", xlim = c(1,8), ylim = c(0,0.3), xlab = "aSMR", ylab = "Standardised Cover")
cols <- rainbow(length(Spp))
for(i in 1:length(Spp)){
  CurrSpp <- Spp[i]
  datSub <- multSpp[multSpp$Species == CurrSpp,]
  
  dat <- data.frame(x=as.numeric(as.character(datSub$aSMR)), r=datSub[,3])
  dat$r <- dat$r/sum(dat$r)###have to standardise
  mu <- dat$x[dat$r == max(dat$r)] ###start optimisation at maximum so doesn't get stuck
  
  g.fit <- try(nls(r ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=mu,sigma=1,k=0.5) , data = dat, nls.control(maxiter = 1000, minFactor = 0.000001)))###Non-linear least-squares regression using gaussian formula
  if(isTRUE(class(g.fit)=="try-error")){ ####if it doesn't work, try again now setting the mean to the maximum value (usually 8)
    points(r ~ x, data = dat, pch = 4, col = cols[i])
    dat$r <- jitter(dat$r)
    g.fit <- nls(r ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(sigma=1,k=0.5) , data = dat, nls.control(maxiter = 1000, minFactor = 0.000001))###Non-linear least-squares regression using gaussian formula
    v <- summary(g.fit)$parameters[,"Estimate"] ###v contains the three parameters for the best fit curve
    plot(function(x) v[2]*exp(-1/2*(x-mu)^2/v[1]^2),col=cols[i],add=T,xlim=range(dat$x))
  }else{
    v <- summary(g.fit)$parameters[,"Estimate"] ###v contains the three parameters for the best fit curve
    points(r ~ x, data = dat, pch = 4, col = cols[i])
    plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),col=cols[i],add=T,xlim=range(dat$x))
  }
}
legend("topright", Spp, col = cols, pch = 15, cex = 0.8)
dev.off() 



####End Kiri#########

######## Convert smr distribution into indicator values
##########These calculations could be augmented by other metrics such as SD  and min/max limits
Spp_aSMR[is.na(Spp_aSMR)] <- 0 ## ignore error messages here
Spp_aSMR <-cbind(Spp_aSMR,rowSums(Spp_aSMR[,-c(1)]))
colnames(Spp_aSMR) [19] <- c("Sum")
Spp_aSMR <- Spp_aSMR[Spp_aSMR$Sum >9,] ## remove species with few occurrences
Spp_aSMR2 <- adorn_percentages(Spp_aSMR[,-c(19)])
#Spp_aSMR2 <- adorn_pct_formatting(Spp_aSMR, digits = 1) # function does not work properly
Spp_aSMR2$TD0 <- Spp_aSMR2$'0'
Spp_aSMR2$TD1 <- Spp_aSMR2$TD0 + Spp_aSMR2$'0.5'+ Spp_aSMR2$'1'
Spp_aSMR2$TD2 <- Spp_aSMR2$TD1 + Spp_aSMR2$'1.5'+ Spp_aSMR2$'2'
Spp_aSMR2$TD3 <- Spp_aSMR2$TD2 + Spp_aSMR2$'2.5'+ Spp_aSMR2$'3'
Spp_aSMR2$TD4 <- Spp_aSMR2$TD3 + Spp_aSMR2$'3.5'+ Spp_aSMR2$'4'
Spp_aSMR2$Fresh <- Spp_aSMR2$'4.5'+ Spp_aSMR2$'5' + Spp_aSMR2$'5.5'
Spp_aSMR2$Moist <- Spp_aSMR2$'6'+ Spp_aSMR2$'6.5' + Spp_aSMR2$'7' 
Spp_aSMR2$Wet <- Spp_aSMR2$'7.5' + Spp_aSMR2$'8' 
###converts the aSMRs into indicator groups
Spp_aSMR2$SMR_Ind <-   ifelse(Spp_aSMR2$TD0>=.25, "TD0", 
                            ifelse (Spp_aSMR2$TD1 > 0.25, "TD1",  
                                ifelse (Spp_aSMR2$TD2 > 0.25, "TD2",
                                  ifelse (Spp_aSMR2$TD3 > 0.25, "TD3",
                                    ifelse (Spp_aSMR2$TD4 > 0.25, "TD4",
                        ifelse(Spp_aSMR2$Wet>=.25, "W", 
                                        ifelse ((Spp_aSMR2$Moist +Spp_aSMR2$Wet) > 0.40, "M",
                                              ifelse (Spp_aSMR2$Fresh > 0.25, "F", ""))))))))
Spp_aSMR2$SMR_Ind <- ifelse ((Spp_aSMR2$Moist + Spp_aSMR2$Wet)>0.10 & Spp_aSMR2$TD4 >0.10, "X",Spp_aSMR2$SMR_Ind)                                                  
write.csv(Spp_aSMR2, file= "Spp_SMR_CalculationMatrix.csv")

Spp_Ind <- Spp_aSMR2[,c(1,27)] # create simple list of calculated indicator value

###update indicator list with exceptions file Spp_SMRInd_Exceptions.csv
SMRexcept <- read.csv("Spp_SMRInd_Exceptions.csv")
SMRexcept <- SMRexcept[SMRexcept$Exceptions!="",]
SMRexcept <- SMRexcept[,c(2,4)]
SMRexcept$Exceptions <- as.character (SMRexcept$Exceptions)
Spp_Ind <- merge(Spp_Ind, SMRexcept, by = "Species", all=TRUE)
Spp_Ind$SMR_Ind <-   ifelse(!is.na (Spp_Ind$Exceptions), Spp_Ind$Exceptions,Spp_Ind$SMR_Ind)
Spp_Ind <- Spp_Ind[,c(1,2)]
write.csv(Spp_Ind, file= "Spp_SMRInd.csv")
save(Spp_Ind,file = "Spp_aSMR_Ind.RData")

#_____WILL _____________________________###
###########Create smr indicator summaries by plot
load("Spp_aSMR_Ind.RData")


#Spp_Ind2 <- unlist(lapply(Spp_Ind$SMR_Ind, tolower))
vegDataSMR <- merge(vegData, Spp_Ind, by.x = "Species", all.x = TRUE) ##add indicator value to species
##reduce to sum of indicator classes by plot
PlotInd <- dcast(vegDataSMR, PlotNumber ~  SMR_Ind, value.var = "Cover",fun.aggregate = sum, na.rm = TRUE)#roll up by indicator value
vegDataSMR2 <- merge(PlotInd, plotEnv2, by.x = "PlotNumber", all.x = TRUE) ##add ENV info
save(vegDataSMR2,file = "PlotIndicatorSum.RData")
write.csv(vegDataSMR2, file= "PlotIndicatorSum.csv")
############add estimated aSMR from vegetation calculation

vegDataSMR2$aSMR_Veg <- 
  ifelse (vegDataSMR2$VW > 20, "9VW",
      ifelse((vegDataSMR2$W + vegDataSMR2$VW)>20 & (vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3 + vegDataSMR2$TD4 + vegDataSMR2$F)/(vegDataSMR2$W + vegDataSMR2$VW) <1,"8W",
          ifelse (((vegDataSMR2$M) >30 & (vegDataSMR2$W + vegDataSMR2$VW) > 10), "7VM",
                ifelse (((vegDataSMR2$M + vegDataSMR2$W + vegDataSMR2$VW) > 10), "6M", 
ifelse(((vegDataSMR2$TD0>=20) & (vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3 + vegDataSMR2$TD4 + vegDataSMR2$F) < 20), "0XD", 
    ifelse ((((vegDataSMR2$TD0 + vegDataSMR2$TD1)>=20) & (vegDataSMR2$TD2 + vegDataSMR2$TD3 + vegDataSMR2$TD4 + vegDataSMR2$F) < 20), "1ED",
            ifelse (((((vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2) > 10) & (vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2)/(vegDataSMR2$TD3 + vegDataSMR2$TD4 + vegDataSMR2$F) > 1)), "2VD",
#            ifelse ((((vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2) > 20) & (vegDataSMR2$TD3 + vegDataSMR2$TD4 + vegDataSMR2$F) < 20), "2VD",
#            ifelse ((((vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3) > 20) & (vegDataSMR2$TD4 + vegDataSMR2$F) < 20), "3MD",      
            ifelse (((((vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3) > 10) & (vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3)/ (vegDataSMR2$TD4 + vegDataSMR2$F) >.5)), "3MD",      
                  ifelse (((vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3 + vegDataSMR2$TD4 > 10) & (vegDataSMR2$M + vegDataSMR2$W) < 2), "4SD",
                          ifelse (((vegDataSMR2$TD0 + vegDataSMR2$TD1 + vegDataSMR2$TD2 + vegDataSMR2$TD3 + vegDataSMR2$TD4 > 3) & (vegDataSMR2$F + vegDataSMR2$M + vegDataSMR2$W) < 10), "4SD","5F"))))))))))
                      #  
                                        #ifelse (((vegDataSMR2$M + vegDataSMR2$W) > 5) & vegDataSMR2$TD4 >30, "SD", "F")))))))) 
Plot_aSMR <- vegDataSMR2[,c("PlotNumber", "aSMR_Veg")]                                      

#####Summarize plot aSMR by hierarchical unit
SUhier <- read.csv("HierarchyTestExportSBSmc2.csv", stringsAsFactors = FALSE)
colnames(SUhier)[1:12]=c("PlotNumber", "1Region", "2Class", "3Order", "4SubOrder", "5Alliance", "6SubAlliance", "7Association", "8SubAssociation", "9Facies", "10Working", "11SiteSeries")
SUhierALL <-SUhier
level <- "11SiteSeries"
SUhier <- SUhier[,c("PlotNumber", level)]
Plot_aSMR <- merge(Plot_aSMR, SUhier, by = "PlotNumber", all.x = FALSE)
VegDataSMR3 <-merge(vegDataSMR2, SUhier, by = "PlotNumber", all.x = FALSE)
write.csv(VegDataSMR3, file = "SBSmc2_Plot_aSMR count.csv")
colnames(Plot_aSMR) [3] <- "SiteSeries"
SS_aSMR <- dcast (Plot_aSMR, SiteSeries ~  aSMR_Veg, fun.aggregate = length)
write.csv(SS_aSMR, file = "SBSmc2_SiteSeries_aSMR count.csv")
#add function to add the expert and modelled aSMR to the veg assessment
#SS_aSMR_compare <- merge(SS_aSMR, , by = "BGC", all.x = FALSE)  


######################Old code


###Fit gaussian curve to total data (or one species) (thanks to https://stats.stackexchange.com/questions/83022/how-to-fit-data-that-looks-like-a-gaussian)
#dat <- data.frame(x=as.numeric(as.character(temp$aSMR)), r=temp$Total)
#dat$r <- dat$r/sum(dat$r)###have to standardise

#g.fit <- nls(r ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=3,sigma=1,k=0.5) , data = dat)###Non-linear least-squares regression using gaussian formula

#v <- summary(g.fit)$parameters[,"Estimate"] ###v contains the three parameters for the best fit curve
#plot(r ~ x, data = dat, xlab = "aSMR", ylab = "Standardised Cover", main = paste(Spp, collapse = ", "))
#plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),col=2,add=T,xlim=range(dat$x)) ###plot distribution using optimised parameters

####Not used##################
#getStats <- function(data){
#  out <- as.data.frame(data$Species)
#  out$Max <- apply(data[,2:17], 1, FUN = max)
#  out$Min <- apply(data[,2:17], 1, FUN = min)
#  out$Mean <- apply(data[,2:17], 1, FUN = mean)
#  out$Median <- apply(data[,2:17], 1, FUN = median)
#  out$Mode <- apply(data[,2:17], 1, FUN = venter, k = 4)
#  colnames(out)[1] <- "Species"
#  return(out)
#}
##############################################3