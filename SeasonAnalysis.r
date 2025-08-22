#libraries
library(readxl)
library(MCMCglmm)
library(multcomp)

###load useful functions 
std.error  <- function(x) {sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))} #calculate standard error

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
} #rescales a vector to fall between 0 and 1, used for plotting

bufferX <- function(x,p) { 
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	} #finds the range of a vector and adds a proportion to either side, used for plotting

###read in dataframe for duckweed image data, make useful variables
mappeddat <- read.csv("~/Seasonal Duckweed Mbio/AlexDuckNewData.csv",header=T)
mappeddat$colsq <- (mappeddat$column)^2
mappeddat$numrow <- as.numeric(as.factor(mappeddat$row))
mappeddat$rowsq <- (mappeddat$numrow)^2
mappeddat$pergreen <- mappeddat$g/(mappeddat$r + mappeddat$g + mappeddat$b)
mappeddat$colint <- (1-(mappeddat$mean/255))*100 ##reverse the ImageJ "mean" metric so that it is color intensity out of 100 rather then whiteness out of 255
mappeddat$aggregation <- mappeddat$area/mappeddat$perim
alextrt <- read.csv("~/Seasonal Duckweed Mbio/AlexDatTRT.csv",header=T)
old <- sort(unique(alextrt$Micr))
newm <- c(1,2,3,4,6,7,8,NA,1,2,3,4,5,6,7,8)
seasonV <- c("Spr","Sum","Sum","Sum","Fall","Fall","Winter", NA, "Spr","Sum","Sum","Sum","Fall","Fall","Fall","Winter")
alextrt$nummicr <- sapply(1:nrow(alextrt),function(z) newm[which(old==alextrt$Micr[z])])
alextrt$summicr <- ifelse(alextrt$nummicr %in% c(2,3,4),"summer","notsummer")
alextrt$season <- sapply(1:nrow(alextrt),function(z) seasonV[which(old==alextrt$Micr[z])])
alextrt$springmicr <- ifelse(alextrt$season=="Spr","spring","notspring")
### old names are in order of AT's table 1 in capstone thesis. specifically:
### "LRD1" "LRD2" "LRD3" "LRD4" "LRD5" "LRD6" "LRD7" 
###corresponds to
### 12 may 2022, 10 jun 2022, 12 july 2022, 9 aug 2022, 16 oct 2022, 11 nov 2022, 18 dec 2022
####"n" is no microbes inoculated
###"TD1"  "TD2"  "TD3"  "TD4"  "TD5"  "TD6"  "TD7"  "TD8" 
###corresponds to
### 7 may 2022, 4 june 2022, 8 july 2022, 5 aug 2022, 10 sep 2022, 13 oct 2022, 5 november ,15 dec 2022
 
 
#read in raw od data from final timepoint, change column names to match contents
raw_od <- as.data.frame(read_excel("~/Seasonal Duckweed Mbio/Undergrad Plate Data- OBrien_combined.xlsx",sheet=1))
colnames(raw_od)[1:2] <- c("Plate","Row")

#convert final timepoint OD data into a long format data frame
startod <- which(colnames(raw_od)=="1")
endod <- which(colnames(raw_od)=="12")
raw_od$wavelength <- sapply(1:nrow(raw_od), function(z) strsplit(raw_od$Read[z],":")[[1]][[2]])
oddat <- data.frame( od = unlist(raw_od[,startod:endod]), 
					Plate = rep(raw_od$Plate,times=12),
					NumRow = as.numeric(as.factor(rep(raw_od$Row,times=12))),
					Col = rep(c(1:12),each = nrow(raw_od)), 
					Wavelength = rep(raw_od$wavelength,times=12) )
oddat$wellid <- paste(oddat$Plate,oddat$NumRow,oddat$Col,sep="_")
#getting rid of  outliers with optical density >1 for either wavelength
oddat$od[oddat$od>1]<-NA
alextrt$wellid <- paste(alextrt$Plate,alextrt$Row,alextrt$Col,sep="_")
oddatshort <- oddat[oddat$Wavelength == "600" & oddat$wellid%in%alextrt$wellid,] 
colnames(oddatshort)[colnames(oddatshort)=="od"] <- "od600"
oddatshort$Wavelength <- NULL
oddatshort$od450 <- oddat$od[oddat$Wavelength == "450" & oddat$wellid%in%alextrt$wellid]
oddatsort <- t(sapply(1:nrow(alextrt), function(z) oddatshort[which(oddatshort$wellid==alextrt$wellid[z]) ,]))
#sum(oddatsort[,"wellid"]==alextrt$wellid)
 
###combine data and treatments into one dataframe
#first keep uninoculated treatment
datandtrt_wuninoc <- cbind(alextrt,mappeddat)
datandtrt_wuninoc$od600 <- unlist(oddatsort[,"od600"])
datandtrt_wuninoc$od450 <- unlist(oddatsort[,"od450"])
datandtrt_wuninoc$inoc <- ifelse(datandtrt_wuninoc$Micr=="n","NoMicr","PlusMicr")
datandtrt_wuninoc$urep <- paste(datandtrt_wuninoc$Plant, datandtrt_wuninoc$rep, sep=".")
write.csv(datandtrt_wuninoc, "~/Seasonal Duckweed Mbio/MappedGrowth_trt_data.csv",row.names=F)

#check normality
shapiro.test(datandtrt_wuninoc$area)
shapiro.test(datandtrt_wuninoc$pergreen)
shapiro.test(datandtrt_wuninoc$colint)
shapiro.test(datandtrt_wuninoc$aggregation)
#allows us to continue, treating data as normal, W values very close to 1
shapiro.test(datandtrt_wuninoc$od600)
shapiro.test(log(datandtrt_wuninoc$od600))
shapiro.test(datandtrt_wuninoc$od450)
shapiro.test(log(datandtrt_wuninoc$od450))

datandtrt_wuninoc$lnod450 <- log(datandtrt_wuninoc$od450)
datandtrt_wuninoc$lnod600 <- log(datandtrt_wuninoc$od600)

#here, remove uninoculated treatment, as not relevant to most tests
datandtrt <- datandtrt_wuninoc[datandtrt_wuninoc$Micr!="n",]



#####Models and plots for effects of inoculation on plant measurements
#random effect of plant, as started at different sized frond units based on plant genotype; and are different genotypes
#random effect of urep, which refers to the location on the plate
InocMod <- MCMCglmm(sqmm~inoc, rand = ~ Plant + urep ,verbose=F,data=datandtrt_wuninoc,nitt=100000,thin=10,burnin=1000,pr=T)
inocmean <- tapply(datandtrt_wuninoc$sqmm, datandtrt_wuninoc$inoc,mean)
inocSE <- tapply(datandtrt_wuninoc$sqmm, datandtrt_wuninoc$inoc,std.error)
png("~/Seasonal Duckweed Mbio/Area_inoc.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(inocmean~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(inocmean-inocSE,inocmean+inocSE)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=inocmean+inocSE,y1=inocmean-inocSE,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("No microbes","Microbes"),las=2)
mtext(expression("Duckweed frond area, mm"^2), side=2, 2.5)
text(1.5,7.55,"n.s.",cex=1.25)
dev.off()
#
InocModCol <- MCMCglmm(colint~inoc, rand = ~ Plant + urep ,verbose=F,data=datandtrt_wuninoc,nitt=100000,thin=10,burnin=1000,pr=T)
inocmeanCol <- tapply(datandtrt_wuninoc$colint, datandtrt_wuninoc$inoc,mean,na.rm=T)
inocSECol <- tapply(datandtrt_wuninoc$colint, datandtrt_wuninoc$inoc,std.error)
png("~/Seasonal Duckweed Mbio/Col_inoc.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(inocmeanCol~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(inocmeanCol-inocSECol,inocmeanCol+inocSECol)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=inocmeanCol+inocSECol,y1=inocmeanCol-inocSECol,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("No microbes","Microbes"),las=2)
mtext("Frond color intensity %", side=2, 2.5)
text(1.5,70.4,"*",cex=2)
dev.off()
#
InocModGrn <- MCMCglmm(pergreen~inoc, rand = ~ Plant + urep ,verbose=F,data=datandtrt_wuninoc,nitt=100000,thin=10,burnin=1000,pr=T)
inocmeanGrn <- 100*tapply(datandtrt_wuninoc$pergreen, datandtrt_wuninoc$inoc,mean,na.rm=T)
inocSEGrn <- 100*tapply(datandtrt_wuninoc$pergreen, datandtrt_wuninoc$inoc,std.error)
png("~/Seasonal Duckweed Mbio/Grn_inoc.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(inocmeanGrn~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(inocmeanGrn-inocSEGrn,inocmeanGrn+inocSEGrn)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=inocmeanGrn+inocSEGrn,y1=inocmeanGrn-inocSEGrn,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("No microbes","Microbes"),las=2)
mtext("Frond greenness %", side=2, 2.5)
text(1.5,44.15,"*",cex=2)
dev.off()
#
InocModAgg <- MCMCglmm(aggregation~inoc, rand = ~ Plant + urep ,verbose=F,data=datandtrt_wuninoc,nitt=100000,thin=10,burnin=1000,pr=T)
inocmeanAgg <- tapply(datandtrt_wuninoc$aggregation, datandtrt_wuninoc$inoc,mean,na.rm=T)
inocSEAgg <- tapply(datandtrt_wuninoc$aggregation, datandtrt_wuninoc$inoc,std.error)
png("~/Seasonal Duckweed Mbio/Agg_inoc.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(inocmeanAgg~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(inocmeanAgg-inocSEAgg,inocmeanAgg+inocSEAgg)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=inocmeanAgg+inocSEAgg,y1=inocmeanAgg-inocSEAgg,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("No microbes","Microbes"),las=2)
mtext("Frond aggregation", side=2, 2.5)
text(1.5,17.5,"n.s.",cex=1.25)
dev.off()
#
InocModOD450 <- MCMCglmm(lnod450~inoc, rand = ~ Plant + urep ,verbose=F,data=datandtrt_wuninoc,nitt=100000,thin=10,burnin=1000,pr=T)
inocmeanod450 <- tapply(datandtrt_wuninoc$lnod450, datandtrt_wuninoc$inoc,mean,na.rm=T)
inocSEod450 <- tapply(datandtrt_wuninoc$lnod450, datandtrt_wuninoc$inoc,std.error)
png("~/Seasonal Duckweed Mbio/od450_inoc.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(inocmeanod450~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(inocmeanod450-inocSEod450,inocmeanod450+inocSEod450)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=inocmeanod450+inocSEod450,y1=inocmeanod450-inocSEod450,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("No microbes","Microbes"),las=2)
mtext("Log of optical density, 450 nm", side=2, 2.5)
text(1.5,-2.85,"*",cex=2)
dev.off()
#
InocModOD600 <- MCMCglmm(lnod600~inoc, rand = ~ Plant + urep ,verbose=F,data=datandtrt_wuninoc,nitt=100000,thin=10,burnin=1000,pr=T)
inocmeanod600 <- tapply(datandtrt_wuninoc$lnod600, datandtrt_wuninoc$inoc,mean,na.rm=T)
inocSEod600 <- tapply(datandtrt_wuninoc$lnod600, datandtrt_wuninoc$inoc,std.error)
png("~/Seasonal Duckweed Mbio/od600_inoc.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(inocmeanod600~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(inocmeanod600-inocSEod600,inocmeanod600+inocSEod600)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=inocmeanod600+inocSEod600,y1=inocmeanod600-inocSEod600,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("No microbes","Microbes"),las=2)
mtext("Log of optical density, 600 nm", side=2, 2.5)
text(1.5,-2.992,"*",cex=2)
dev.off()

#percent changes
inocmeanslist <- list(inocmean,inocmeanCol,inocmeanGrn,inocmeanAgg,inocmeanod450,inocmeanod600)
PercentDiffInoc <- round(100*(unlist(lapply(1:length(inocmeanslist), function(z)
			(inocmeanslist[[z]][2]-inocmeanslist[[z]][1])/inocmeanslist[[z]][1]))),digits=1)
names(PercentDiffInoc) <- c("area","color intensity","greenness","aggregation","od450","od600")
PercentDiffInoc

#####Models and plots for effects of summer vs not summer microbes on plant measurements
#random effect of plant, as started at different sized frond units based on plant genotype; and are different genotypes
#random effect of urep, which refers to the location on the plate
SummMicr <- MCMCglmm(sqmm~summicr, rand = ~ Plant + urep ,verbose=F,data=datandtrt,nitt=100000,thin=10,burnin=1000,pr=T)
SummMicrG <- MCMCglmm(pergreen~summicr, rand = ~ Plant + urep ,verbose=F,data=datandtrt,nitt=100000,thin=10,burnin=1000,pr=T)
SummMicrCI <- MCMCglmm(colint~summicr, rand = ~ Plant +urep ,verbose=F,data=datandtrt,nitt=100000,thin=10,burnin=1000,pr=T)
SummMicrAgg <- MCMCglmm(aggregation~summicr, rand = ~ Plant +urep ,verbose=F,data=datandtrt,nitt=100000,thin=10,burnin=1000,pr=T)
SummMicrod450 <- MCMCglmm(lnod450~summicr, rand = ~ Plant +urep ,verbose=F,data=datandtrt,nitt=100000,thin=10,burnin=1000,pr=T)
SummMicrod600 <- MCMCglmm(lnod600~summicr, rand = ~ Plant +urep ,verbose=F,data=datandtrt,nitt=100000,thin=10,burnin=1000,pr=T)
summean <- tapply(datandtrt$sqmm, datandtrt$summicr,mean)
sumSE <- tapply(datandtrt$sqmm, datandtrt$summicr,std.error)
png("~/Seasonal Duckweed Mbio/Area_sum.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(summean~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(summean-sumSE,summean+sumSE)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=summean+sumSE,y1=summean-sumSE,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("Other","Summer"),las=2)
mtext(expression("Duckweed frond area, mm"^2), side=2, 2.5)
mtext("Microbe Season",side=1, line =5)
text(1.5,7.85,"n.s.",cex=1.25)
dev.off()
#
summeanGreen <- 100*tapply(datandtrt$pergreen, datandtrt$summicr,mean,na.rm=T)
sumSEGreen <- 100*tapply(datandtrt$pergreen, datandtrt$summicr,std.error)
png("~/Seasonal Duckweed Mbio/Green_sum.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(summeanGreen~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(summeanGreen-sumSEGreen,summeanGreen+sumSEGreen)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=summeanGreen+sumSEGreen,y1=summeanGreen-sumSEGreen,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("Other","Summer"),las=2)
mtext("Frond greenness %", side=2, 2.5)
mtext("Microbe Season",side=1, line =5)
text(1.5,44.31,"p<0.1",cex=1.25)
dev.off()
#
summeanColor <- tapply(datandtrt$colint, datandtrt$summicr,mean,na.rm=T)
sumSEColor <- tapply(datandtrt$colint, datandtrt$summicr,std.error)
png("~/Seasonal Duckweed Mbio/ColInt_sum.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(summeanColor~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(summeanColor-sumSEColor,summeanColor+sumSEColor)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=summeanColor+sumSEColor,y1=summeanColor-sumSEColor,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("Other","Summer"),las=2)
mtext("Frond color intensity %", side=2, 2.5)
mtext("Microbe Season",side=1, line =5)
text(1.5,71.2,"*",cex=2)
dev.off()
#
summeanagg <- tapply(datandtrt$aggregation, datandtrt$summicr,mean,na.rm=T)
sumSEagg <- tapply(datandtrt$aggregation, datandtrt$summicr,std.error)
png("~/Seasonal Duckweed Mbio/Aggregation_sum.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(summeanagg~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(summeanagg-sumSEagg,summeanagg+sumSEagg)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=summeanagg+sumSEagg,y1=summeanagg-sumSEagg,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("Other","Summer"),las=2)
mtext("Frond aggregation", side=2, 2.5)
mtext("Microbe Season",side=1, line =5)
text(1.5,18.25,"*",cex=2)
dev.off()
#
summeanod450 <- tapply(datandtrt$lnod450, datandtrt$summicr,mean,na.rm=T)
sumSEod450 <- tapply(datandtrt$lnod450, datandtrt$summicr,std.error)
png("~/Seasonal Duckweed Mbio/od450_sum.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(summeanod450~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(summeanod450-sumSEod450,summeanod450+sumSEod450)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=summeanod450+sumSEod450,y1=summeanod450-sumSEod450,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("Other","Summer"),las=2)
mtext("Log of optical density, 450 nm", side=2, 2.5)
mtext("Microbe Season",side=1, line =5)
text(1.5,-2.787,"*",cex=2)
dev.off()
#
summeanod600 <- tapply(datandtrt$lnod600, datandtrt$summicr,mean,na.rm=T)
sumSEod600 <- tapply(datandtrt$lnod600, datandtrt$summicr,std.error)
png("~/Seasonal Duckweed Mbio/od600_sum.png",height=4,width=3, units = "in", res = 1200)
par(mar=c(6,4,1,1))
plot(summeanod600~c(1,2),pch=21,cex=3,xlim=c(0,3),xaxt="n", ylab="",xlab="",
	ylim=bufferX(range(c(summeanod600-sumSEod600,summeanod600+sumSEod600)),0.1), bg=rgb(0.8,0.8,0.8))
arrows(c(1,2), y0=summeanod600+sumSEod600,y1=summeanod600-sumSEod600,length=0,lwd=3)
axis(side=1,at=c(1,2),labels=c("Other","Summer"),las=2)
mtext("Log of optical density, 600 nm", side=2, 2.5)
mtext("Microbe Season",side=1, line =5)
text(1.5,-2.947,"*",cex=2)
dev.off()

#percent changes
summmeanslist <- list(summean,summeanColor,summeanGreen,summeanagg,summeanod450,summeanod600)
PercentDiffSumm <- round(100*(unlist(lapply(1:length(summmeanslist), function(z)
			(summmeanslist[[z]][2]-summmeanslist[[z]][1])/summmeanslist[[z]][1]))),digits=1)
names(PercentDiffSumm) <- c("area","color intensity","greenness","aggregation","od450","od600")
PercentDiffSumm

######correlation analysis/exploration, figure S2
resp_vars <- colnames(datandtrt)%in%c("sqmm", "pergreen" , "colint" , "aggregation" ,"lnod450" , "lnod600")
respdat <- datandtrt_wuninoc[,resp_vars]
respdat_micr <- datandtrt[,resp_vars]
bwr <- colorRampPalette(c(rgb(1,0,0),rgb(1,1,1),rgb(0,0,1)))
resp_cors <- cor(respdat[,c(1,4,2,3,5,6)],use="complete.obs")#note order not the same as resp_vars
RespNames <- c("Area", "Aggregation" , "%Greenness" , "%Color Intensity" ,"ln(OD450)" , "ln(OD600)")#note order not the same as resp_vars
vals <- seq(from = 0, to = 1, length.out=ncol(resp_cors))
xvals <- matrix(rep(vals,times=6),nrow=6,byrow=T)
yvals <- matrix(rep(vals,times=6),nrow=6,byrow=F)
png("~/Seasonal Duckweed Mbio/Correlations.png",height=5,width=5, units = "in", res = 1200)
par(mar=c(7,7,1,1))
image(resp_cors,zlim=c(-1,1),col=bwr(100),xaxt="n",yaxt="n")
	text(xvals,yvals,round(resp_cors,digits=1))
	axis(at=c(vals),side=1,labels=RespNames,las=2)
	axis(at=c(vals),side=2,labels=RespNames,las=2)
dev.off()


#sums of squares
ssbyvar <- function(response,category.vec){ #sums of squares function
                means <- tapply(response,category.vec,mean,na.rm=T) #take the means by category
                ssresid <- sum(sapply(sort(unique(category.vec)), function(z) sum( (response[category.vec==z] - means[names(means)==z])^2,na.rm=T ))) #square of difference of each datapoint from its associated treatment mean (residual variation)
                sstot <- sum((response-mean(response,na.rm=T))^2,na.rm=T) #square of difference of each datapoint from the grand mean (total variation)
                sst <- (sstot-ssresid) # total variation - residual variation = treatment variation
                return(sst/sstot) # treatment variance as a fraction of total variation
                }

ssPlant <- 100*sapply(c(1,4,2,3,5,6), function(z) ssbyvar(respdat_micr[,z],datandtrt$Plant) )
round(ssPlant)
ssMicr <- 100*sapply(c(1,4,2,3,5,6), function(z) ssbyvar(respdat_micr[,z],datandtrt$Micr) ) - ssPlant
# ssPlant <- 100*sapply(c(1,4,2,3,5,6), function(z) ssbyvar(respdat[,z],datandtrt_wuninoc$Plant) )
# ssMicr <- 100*sapply(c(1,4,2,3,5,6), function(z) ssbyvar(respdat[,z],paste(datandtrt_wuninoc$Micr,datandtrt_wuninoc$Plant)) ) - ssPlant
png("~/Seasonal Duckweed Mbio/VarExplbyTrts_Inoc.png",height=4,width=4, units = "in", res = 1200)
par(mar=c(7,4,1,1))
bg <- barplot(rbind(ssPlant,ssMicr)[,c(1,4,2,3,5,6)],ylim=c(0,40),col=c(rgb(0,0,0),rgb(0.75,0.75,0.75)))
	axis(at=bg, labels=RespNames,side=1,las=2)
	legend(bg[1],40,c("Site","Microbes (within site)"),fill=c(rgb(0,0,0),rgb(0.75,0.75,0.75)),bty="n")
	mtext("%Variation explained", side=2,line=2.5)
dev.off()



# #####figure for all the treatment means, each response variable. Figure S3
# 
allmn <- tapply(datandtrt_wuninoc$sqmm, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),mean)
allSE <- tapply(datandtrt_wuninoc$sqmm, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),std.error)

allmnaggr <- tapply(datandtrt_wuninoc$aggregation, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),mean,na.rm=T)
allSEaggr <- tapply(datandtrt_wuninoc$aggregation, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),std.error)
allmngrn <- tapply(datandtrt_wuninoc$pergreen, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),mean,na.rm=T)
allSEgrn <- tapply(datandtrt_wuninoc$pergreen, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),std.error)
allmncolint <- tapply(datandtrt_wuninoc$colint, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),mean,na.rm=T)
allSEcolint <- tapply(datandtrt_wuninoc$colint, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),std.error)

allmnod450 <- tapply(datandtrt_wuninoc$lnod450, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),mean,na.rm=T)
allSEod450 <- tapply(datandtrt_wuninoc$lnod450, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),std.error)
allmnod600 <- tapply(datandtrt_wuninoc$lnod600, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),mean,na.rm=T)
allSEod600 <- tapply(datandtrt_wuninoc$lnod600, paste(datandtrt_wuninoc$Plant,datandtrt_wuninoc$Micr),std.error)

alldates <- c("12 May", "10 Jun", "12 Jul", "9 Aug", "16 Oct", "11 Nov", "18 Dec", "None", "None",
			   "7 May", "4 June", "8 Jul", "5 Aug", "10 Sep", "13 Oct", "5 Nov", "15 Dec")

png("~/Seasonal Duckweed Mbio/allmeans.png",height=10,width=5, units = "in", res = 1200)
par(mar=c(0,4.5,0,1))
par(oma=c(5,0,1,0))
par(mfrow=c(6,1))
ylims <- bufferX(range(c(allmn-allSE,allmn+allSE)),0.1)
plot(allmn~c(1:17),pch=NA,xlim=c(0,18),xaxt="n", ylab="",xlab="",ylim=ylims)
polygon(c(7.5,9.5,9.5,7.5), c(-500,-500,200,200),col=rgb(0,0,0,alpha=0.15), border =NA)
points(allmn~c(1:17),pch=21,cex=3, bg=rgb(0.8,0.8,0.8))
abline(v=8.5,lty=3)
arrows(c(1:17), y0=allmn+allSE,y1=allmn-allSE,length=0,lwd=3)
text(0,11.25,"LaRoche Pond",adj=c(0,0.5))
text(18,11.25,"Thompson Farm",adj=c(1,0.5))
mtext(expression("Area, mm"^2), side=2, 2.5)
#
ylims <- bufferX(range(c(allmnaggr-allSEaggr,allmnaggr+allSEaggr)),0.1)
plot(allmnaggr~c(1:17),pch=NA,xlim=c(0,18),xaxt="n", ylab="",xlab="",ylim=ylims)
polygon(c(7.5,9.5,9.5,7.5), c(-500,-500,200,200),col=rgb(0,0,0,alpha=0.15), border =NA)
points(allmnaggr~c(1:17),pch=21,cex=3, bg=rgb(0.8,0.8,0.8))
abline(v=8.5,lty=3)
arrows(c(1:17), y0=allmnaggr+allSEaggr,y1=allmnaggr-allSEaggr,length=0,lwd=3)
mtext("Aggregation", side=2, 2.5)
#
ylims <- bufferX(range(c(allmngrn-allSEgrn,allmngrn+allSEgrn)),0.15)
plot(allmngrn~c(1:17),pch=NA,xlim=c(0,18),xaxt="n", ylab="",xlab="",ylim=ylims)
polygon(c(7.5,9.5,9.5,7.5), c(-500,-500,200,200),col=rgb(0,0,0,alpha=0.15), border =NA)
points(allmngrn~c(1:17),pch=21,cex=3, bg=rgb(0.8,0.8,0.8))
abline(v=8.5,lty=3)
arrows(c(1:17), y0=allmngrn+allSEgrn,y1=allmngrn-allSEgrn,length=0,lwd=3)
mtext("%Greenness", side=2, 2.5)
#
ylims <- bufferX(range(c(allmncolint-allSEcolint,allmncolint+allSEcolint)),0.15)
plot(allmncolint~c(1:17),pch=NA,xlim=c(0,18),xaxt="n", ylab="",xlab="",ylim=ylims)
polygon(c(7.5,9.5,9.5,7.5), c(-500,-500,200,200),col=rgb(0,0,0,alpha=0.15), border =NA)
points(allmncolint~c(1:17),pch=21,cex=3, bg=rgb(0.8,0.8,0.8))
abline(v=8.5,lty=3)
arrows(c(1:17), y0=allmncolint+allSEcolint,y1=allmncolint-allSEcolint,length=0,lwd=3)
mtext("%Color Intensity", side=2, 2.5)
#
ylims <- bufferX(range(c(allmnod450-allSEod450,allmnod450+allSEod450)),0.1)
plot(allmnod450~c(1:17),pch=NA,xlim=c(0,18),xaxt="n", ylab="",xlab="",ylim=ylims)
polygon(c(7.5,9.5,9.5,7.5), c(-500,-500,200,200),col=rgb(0,0,0,alpha=0.15), border =NA)
points(allmnod450~c(1:17),pch=21,cex=3, bg=rgb(0.8,0.8,0.8))
abline(v=8.5,lty=3)
arrows(c(1:17), y0=allmnod450+allSEod450,y1=allmnod450-allSEod450,length=0,lwd=3)
mtext("ln(OD450)", side=2, 2.5)
#
ylims <- bufferX(range(c(allmnod600-allSEod600,allmnod600+allSEod600)),0.1)
plot(allmnod600~c(1:17),pch=NA,xlim=c(0,18),xaxt="n", ylab="",xlab="",ylim=ylims)
polygon(c(7.5,9.5,9.5,7.5), c(-500,-500,200,200),col=rgb(0,0,0,alpha=0.15), border =NA)
points(allmnod600~c(1:17),pch=21,cex=3, bg=rgb(0.8,0.8,0.8))
abline(v=8.5,lty=3)
arrows(c(1:17), y0=allmnod600+allSEod600,y1=allmnod600-allSEod600,length=0,lwd=3)
mtext("ln(OD600)", side=2, 2.5)
axis(side=1,las=2,labels=alldates, at=c(1:17))
dev.off()

### "LRD1" "LRD2" "LRD3" "LRD4" "LRD5" "LRD6" "LRD7" 
###corresponds to
### 12 may 2022, 10 jun 2022, 12 july 2022, 9 aug 2022, 16 oct 2022, 11 nov 2022, 18 dec 2022
####"n" is no microbes inoculated
###"TD1"  "TD2"  "TD3"  "TD4"  "TD5"  "TD6"  "TD7"  "TD8" 
###corresponds to
### 7 may 2022, 4 june 2022, 8 july 2022, 5 aug 2022, 10 sep 2022, 13 oct 2022, 5 november ,15 dec 2022


