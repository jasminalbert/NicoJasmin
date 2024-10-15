data_location <- "./data/"

cleandat <- function(path){
	dat <- read.csv(path, skip=1)
	datid <- dat[,1:6]
	datOpen <- data.frame(apply(dat[,7:11],2,as.numeric))
	datClosed <- data.frame(apply(dat[,12:16],2,as.numeric))
	datList <- list(open=datOpen, closed=datClosed)
	for (trt in names(datList)){
		datList[[trt]] <- cbind(datid, trt,datList[[trt]])
		naRows <- apply(is.na(datList[[trt]][,8:12]),1,all)
		datList[[trt]] <- datList[[trt]][!naRows,]
		nnames <- substring(names(datList[[trt]])[8:12], 2)
		names(datList[[trt]])[8:12] <- nnames
		
		goodS <- datList[[trt]]$seedsGood
		badS <- datList[[trt]]$seedsBad
		totS <- datList[[trt]]$seedsTotal
		noSum <- is.na(totS)
		totS[noSum] <- goodS[noSum]+badS[noSum]
		summand <- !is.na(goodS)&!is.na(badS)
		sumcheck <- totS[summand]==goodS[summand]+badS[summand]
		if(!all(sumcheck)){stop("check sums")}
		datList[[trt]]$seedsTotal <- totS
		}
	return(datList)
}


getStats <- function(datList){
	datList <- c(datList, list(both=rbind(datList$open, datList$closed)))
	res <- datList
	seedcounttype <- names(datList[[1]])[9:11]
	for (trt in names(datList)){
		data <- datList[[trt]]
		useRows <- !is.na(data$mass)&!is.na(data$seedsTotal)
		data <- data[useRows,]
		stats <- data.frame()
		for (seeds in seedcounttype){
			#seedres <- seeds
			trend <- lm(data[,seeds]~data$mass)
			m <- trend$coefficients["data$mass"]
			b <- trend$coefficients["(Intercept)"]
			xran <- range(data$mass)
			x <- seq(xran[1], xran[2], 0.001)
			y <- (m * x) + b
			N <- nrow(data)

			p_value <- summary(trend)$coefficients[2,4]
			rsq <- summary(trend)$r.squared
			title <- paste0(unique(data$S), collapse="")
			title <- unique(paste0(title, data$Sp))
			title <- paste(title, trt, seeds, sep="-")
			plot(data$mass, data[,seeds], xlab = "mass", ylab = "seeds total", ylim=range(data[,seeds]), col="darkblue", pch=21, bg="lightgrey")
			title(main=title, line=0.12, cex.main=1.19)
			lines(x, y, col="darkgrey")
			text(0,max(data[,seeds]) * .95, labels = paste("y =", round(m, digits = 2), "x +", round(b, digits = 2), "; p =", signif(p_value,3), "; \nR^2 =", signif(rsq,3), "; N=",N), adj = 0)
	
			seedstat <- c(m=m, b=b, p=p_value, R=rsq, N=N)
			stats <- rbind(stats,seedstat)	
		}
		names(stats) <- c("m","b","p","R","N")
		rownames(stats) <- seedcounttype
		res[[trt]] <- stats
	}
	return(res)
}

predSeeds <- function(dat,stats,seeds="seedsTotal", trt="both", pred=T){
	#fill in total seeds
	#for rows that dont have it
	noTs <- is.na(dat[,seeds])
	if (pred){
		check <- all(!is.na(dat$mass[noTs]))
		if(!check){stop("missing mass values")}
		stats <- c(m=stats[[trt]][seeds,"m"], b=stats[[trt]][seeds,"b"])
		x <- dat$mass[noTs]
		dat$seedsTotal[noTs] <- stats["m"]*x+stats["b"]
	}
	dat$est <- ifelse(noTs,1,0) 
	dat$spb <- dat[,seeds]/dat$blooms #some NA blooms to check
	dat$id <- apply(dat[,1:6],1,paste, collapse='')
	return(dat)
}

#datList <- cleandat(path)
#getStats(datList)
#Lauren meeting 101524: log scale, same axis, use dif trt model 
#one figure 

fig_loc <- "./figures/"
fig101124 <- paste0(fig_loc,"fig101124.pdf")
fig101524 <- paste0(fig_loc,"fig101524.pdf")
file_names <- list.files(data_location)
file_names <- file_names[-4]

res <- as.list(file_names)
names(res) <- file_names
datList <- res
pdf(fig101124)
par(mfrow=c(3,3), mar=c(1,1.5,1.5,0), oma=c(1,1,.5,.5), mgp=c(2,.25,0), tcl=-.2, xpd=F)
for (f in file_names){
	path <- paste0(data_location,f)
	#mass<-ifelse(f==file_names[4],F,T)
	datList[[f]] <- cleandat(path)
	if(!f==file_names[4]){
		res[[f]] <- getStats(datList[[f]])
	} else{res[[f]] <- NA}
	#res[[f]] <- ifelse(f==file_names[4],NA,getStats(datList[[f]]))
};
names(res) <- substring(file_names,13,nchar(file_names)-4)
names(datList) <- names(res)
datList$CLAAMO$open <- rbind(datList$GLTCCLAAMO$open, datList$WWCLAAMO$open)
datList$CLAAMO$closed <- rbind(datList$GLTCCLAAMO$closed, datList$WWCLAAMO$closed)
res$CLAAMO <- getStats(datList$CLAAMO)
dev.off()
#both seeds good is best
#probably should combine GLTC and WW

str(datList)
#combine open and closed
datList <- datList[-which(names(datList)%in%c("WWCLAAMO","GLTCCLAAMO"))]
datList <- lapply(datList, function(X){rbind(X$open, X$closed)})

datListp <- datList
spb <- datList
for(sp in names(datList)){
	p <- ifelse(sp=="WWCLAPUR",F,T)
	datListp[[sp]] <- predSeeds(datList[[sp]], res[[sp]], pred=p)
	spb[[sp]] <- reshape(datListp[[sp]][,c("id","trt","spb")], direction="wide", idvar="id", timevar="trt")
	naVals <- apply(is.na(spb[[sp]][,2:3]),1,any)
	spb[[sp]] <- spb[[sp]][!naVals,]
	barplot(as.matrix(spb[[sp]][,2:3]))
}

#separate by mix type




#if (dat$S[1]=="GLTC" & dat$Sp[1]=="CLAPUR"){
	#	dat <- dat[-63,]}
#datOpen <- cbind(datid, datOpen)
	#colnames <- names(datOpen)
	#colnames[7:11] <- substring(colnames[7:11], 2)
	#datClosed <- cbind(datid, datClosed)
	#names(datClosed) <- colnames; names(datOpen) <- colnames
#if(!is.na(datList[[d]]$seedsGood) & !is.na(datList[[d]]$seedsBad)){
			#sumcheck<-datList[[d]]$seedsTotal==datList[[d]]$seedsGood+datList[[d]]$seedsBad 
		#}
	#datOpen$seedsTotal <- datOpen$seedsGood + datOpen$seedsBad
	#datClosed$seedsTotal <- datClosed$seedsGood + datClosed$seedsBad
	#datOpen <- datOpen[!is.na(datOpen$mass)&!is.na(datOpen$seedsTotal),]
	#datClosed <- datClosed[!is.na(datClosed$mass)&!is.na(datClosed$seedsTotal),]
	
	#return(list(open=datOpen, closed=datClosed))
#	useRows <- !is.na(datList[[d]]$mass)&!is.na(datList[[d]]$seedsTotal)
#		if (!mass){useRows<-!is.na(datList[[d]]$seedsTotal)}
#		datList[[d]] <- datList[[d]][useRows,]
