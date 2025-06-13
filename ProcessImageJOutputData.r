

MapToWells <- function(dat,map,firstcol, sumcols, meancols){ #columns labeled "roi", and "image"
	n <- nrow(map)
	outdata <- matrix(NA,nrow=n, ncol=(length(sumcols)+length(meancols)+1))
	for(i in 1:n){
		rois <- map[i,firstcol:ncol(map)] 
		p <- as.character(map$image[i])
		welldat <- dat[as.character(dat$image)==p & dat$roi%in%rois,]# 
		welldat.sums <- colSums(welldat[,sumcols])
		welldat.means <- colMeans(welldat[,meancols],na.rm=T)
		wellstats <- c(nrow(welldat),welldat.sums,welldat.means)
		outdata[i,] <- wellstats
		}
	mappeddata <- cbind(map[,1:(firstcol-1)],outdata)
	colnames(mappeddata) <- c(colnames(map)[1:(firstcol-1)],"particles",colnames(dat)[sumcols],colnames(dat)[meancols])
	return(mappeddata)
}

concatarea <- as.data.frame(read_excel("~/Seasonal Duckweed Mbio/NewDuckweedPlateDat.xlsx" ,sheet=1))
concatmap <- as.data.frame(read_excel("~/Seasonal Duckweed Mbio/NewDuckweedPlateMap.xlsx" ,sheet=1))

mappeddat <- MapToWells(concatarea,concatmap,5, c(3,7), c(4,8:14))
#convert pixel area into sq mm.
#the standard width of a 96-well plate is 85.4 millimeters
imagestats <- as.data.frame(read_excel("~/Alex experiment/ImageSettings.xlsx" ,sheet=1))
plateratio <- imagestats$plate_width/85.4
mappeddat$sqmm <- sapply(1:nrow(mappeddat), function(z) mappeddat$area[z]/(plateratio[imagestats$image==mappeddat$image[z]]^2))


write.csv(mappeddat,"~/Seasonal Duckweed Mbio/DuckweedData.csv")

