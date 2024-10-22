setwd("/Users/jasminalbert/Documents/HALLETTLAB/NicoJasmin")
data_location <- "./data/"
filename <- "bigData_workingCopy - CCDATA.csv"
data <- read.csv(paste0(data_location,filename))
str(data)
unique(data$Sp); unique(data$M)
specs <- c("CLAAMO","CLAPUR")
data <- data[data$Sp%in%specs & data$M!=0,]
idend <- which(names(data)=="Fc")
censtart <- which(names(data)=="CLAAMO")
cenend <- which(names(data)=="ERILAN")
end <- which(names(data)=="rop")
select <- c(1:idend,censtart:end)
cenmat <- data[,censtart:cenend]
data <- data[,select]
iddat <- data[,1:idend]
head(cenmat)
cenmat <- apply(cenmat,2, as.numeric)
dens <- rowSums(cenmat, na.rm=T)
rich <- rowSums(cenmat>0, na.rm=T)
cendat <- cbind(iddat,cenmat,dens,rich,data$rop)
head(cendat)
str(cendat)
all((dens==0)==(rich==0))
cendat <- cendat[!dens==0,]

dmax <- max(dens);rmax <- max(rich)
dbrks <- seq(0,round(dmax/10)*10,5);rbrks <- seq(0,rmax,1)
par(mfrow=c(1,3),mar=c(1,1,1,0.5), oma=c(2,2,0.5,0),mgp=c(1,0.1,0), tcl=-0.1)
hist(cendat$dens,breaks=dbrks)
hist(cendat$rich,breaks=rbrks)
plot(cendat$dens, cendat$rich, xlim=c(0,dmax),ylim=c(0,rmax), col=cols[specs], pch=19)
par(mfrow=c(2,3))
for(sp in specs){
	dat.sp <- cendat[cendat$Sp==sp,]
	hist(dat.sp$dens, dbrks,main=paste(sp, "nbhd density"))
	hist(dat.sp$rich, rbrks,main=paste(sp, "nbhd richness"))
	plot(dat.sp$dens, dat.sp$rich, xlim=c(0,dmax),ylim=c(0,rmax))
}
mixes <- unique(data$M)
gs <- unique(data$genSpec)
data$gs <- as.factor(data$genSpec)
plot(data$gs, data$M, cex.axis=0.5)
specs <- as.factor(specs)

alph <- 0.5
cacol <- c(col2rgb("violetred")/255);cacol <- rgb(cacol[1],cacol[2],cacol[3],alph)
cpcol <- c(col2rgb("mediumpurple")/255);cpcol <- rgb(cpcol[1],cpcol[2],cpcol[3],alph)
cols <- c(cacol,cpcol)
cendat$Sp <- as.factor(cendat$Sp)
divar <- c("dens","rich")
brks <- list(dens=dbrks, rich=rbrks)
pdf("nbhdfigs.pdf", height=6, width=8)
par(mfrow=c(2,3),mar=c(.5,.5,2,0.5), oma=c(2,1,0.5,0),mgp=c(1,0.1,0), tcl=-0.1)
for(m in gs){
	dat.m <- cendat[cendat$genSpec==m,]
	for (v in divar){
		hist(dat.m[,v], brks[[v]],main="",xlab="")	
		title(main=paste(m, "nbhd", v), line=0)
		for(sp in specs){
			hist(dat.m[,v][dat.m$Sp==sp],brks[[v]],add=T, col=cols[specs==sp],border=cols[specs=sp])
		}	
	}
	plot(dat.m$dens, dat.m$rich, type="n", xlim=c(1,dmax),ylim=c(0,rmax), xlab="", ylab="")
	title(xlab="nbhd density", ylab="nbhd richness", line=0.95, xpd=NA)
	abline(h=1:6, col="lightgrey")
	points(dat.m$dens, jitter(dat.m$rich,1.2),col=cols[dat.m$Sp], pch=19)
};dev.off()

#merge
data <- readRDS("bloomseeds2.RDS")
iddat <- cendat[,1:idend]
cendat$id <- apply(iddat, 1, paste, collapse='')
adat <- merge(cendat, data)
names(adat)[names(adat)%in%c("data$rop")]<- c("rop")
adat$tbo <- adat$rop+adat$blooms.open
adat$tbc <- adat$rop+adat$blooms.closed
adat$tb <- adat$rop+adat$blooms.closed+adat$blooms.open
adat$tso <- adat$open*adat$tb
adat$tsc <-adat$closed*adat$tb

tests <- expand.grid(trt=c("tso","tsc"), v=divar, stringsAsFactors=F)
testsL <- as.list(apply(tests,1,paste,collapse=''))
names(testsL) <- unlist(testsL)
rownames(tests) <- names(testsL)
pdf("seedsNbhdReg.pdf", height=6,width=7)
par(mfrow=c(2,2),mar=c(.5,.5,2,1), oma=c(2,2,0.5,0),mgp=c(1,0.1,0), tcl=-0.1)
for (t in names(testsL)){
	testin <- tests[t,]
	nbhd <- as.character(testin[,"v"])
	ts <- as.character(testin[,"trt"])
	testsL[[t]]$model <- summary(lm(adat[,ts]~adat[,nbhd]))
	x<- seq(0,max(adat[,nbhd], na.rm=T),0.1)
	m <- testsL[[t]]$model$coefficients[2,1]
	b <- testsL[[t]]$model$coefficients[1,1]
	p <- testsL[[t]]$model$coefficients[2,4]
	y <- m*x+b
	testsL[[t]]$y <- y
	testsL[[t]]$x <- x
	plot(adat[,nbhd],adat[,ts],col=cols[adat$Sp], pch=19, cex=0.8, bty='l')
	title(xlab=nbhd, ylab="", line=0.9, xpd=NA)
	lines(x,y, col="black", lty=2, lwd=2.5)
	text(x=max(adat[,nbhd], na.rm=T), y=max(adat[,ts],na.rm=T)*.97, labels=signif(p,3), adj=c(1,0))
}; 
title(main=c("open"), outer=T, adj=c(0.25), line=-1.4, font.main=1, cex.main=1.6)
title(main=c("closed"), outer=T, adj=c(0.75), line=-1.4, font.main=1, cex.main=1.6)
title(ylab="total seeds", outer=T, xpd=NA, cex.lab=1.3,line=0.8)
dev.off()
#odens <- lm(adat$tso~adat$dens)



plot(adat$dens,adat$tso, col=cols[adat$Sp], pch=19, cex=0.5)
plot(y=adat$tsc,x=adat$dens, col=cols[adat$Sp], pch=19, cex=0.5)
plot(y=adat$tso, x=adat$rich, col=cols[adat$Sp], pch=19, cex=0.5)
plot(y=adat$tsc,x=adat$rich, col=cols[adat$Sp], pch=19, cex=0.5)
#adat[adat$tsc>5000,]

lmpred <- function(x,y){
	model <-  summary(lm(y~x))
	m <- model$coefficients[2,1]
	b <- model$coefficients[1,1]
	p <- model$coefficients[2,4]
	x<- seq(0,max(x, na.rm=T),0.1)
	y<-m*x+b
	return(list(m,b,p,x,y))
}
pred <- function(max,m,b){
	x<- seq(0,max,0.1)
	y<-m*x+b
	if(length(y)<1){y<-rep(NA,length(x))}
	return(list(m,b,x,y))
}

cores <- data.frame()
for (s in specs){
	for (m in mixes){
		select <- adat$M==m & adat$Sp==s
		dat.tmp <- adat[select,]
		if (nrow(dat.tmp)>1){
			for (t in rownames(tests)){
				nbhd <- dat.tmp[,tests[t,"v"]]
				seeds <- dat.tmp[,tests[t,"trt"]]
				tests[t,]$rho <- cor(nbhd,seeds,"na.or.complete")
				model <- summary(lm(seeds~nbhd))
				tests[t,]$m <- model$coefficients[2,1]
				tests[t,]$b <- model$coefficients[1,1]
				tests[t,]$p <- model$coefficients[2,4]
			}
			cores <- rbind(cores,data.frame(s,m,tests))
		}
	}
}

pdf("seedsNbhdSpM.pdf", height=6,width=7)
par(mfrow=c(2,2),mar=c(.5,.5,2,1), oma=c(2,2,0.5,0),mgp=c(1,0.1,0), tcl=-0.1)
for (m in mixes){
	dat.tmp <- adat[adat$M==m,]
	gs <- unique(dat.tmp$genSpec)
	res.tmp <- cores[cores$m==m & cores$v=="dens",]
	maxy <- max(c(dat.tmp$tso,dat.tmp$tsc),na.rm=T)
	maxx <- max(dat.tmp$dens,na.rm=T)
	plot(dat.tmp$dens, dat.tmp$tso,col=cols[dat.tmp$Sp], pch=19, cex=0.8, bty='l', ylim=c(0,maxy),main=paste("mix",m,gs,"open"))
	res <-pred(maxx, res.tmp[res.tmp$trt=="tso"&res.tmp$s=="CLAAMO",]$m.1, res.tmp[res.tmp$trt=="tso"&res.tmp$s=="CLAAMO",]$b)
	y <- res[[4]];x<-res[[3]]
	p <- res.tmp[res.tmp$trt=="tso"&res.tmp$s=="CLAAMO",]$p
	lines(x,y,lwd=2,lty=2,col=cacol)
	text(maxx,maxy,labels=paste('p=',signif(p,3)),adj=c(1,1),col=cacol)
	res <-pred(maxx, res.tmp[res.tmp$trt=="tso"&res.tmp$s=="CLAPUR",]$m.1, res.tmp[res.tmp$trt=="tso"&res.tmp$s=="CLAPUR",]$b)
	y <- res[[4]];x<-res[[3]]
	p <- res.tmp[res.tmp$trt=="tso"&res.tmp$s=="CLAPUR",]$p
	lines(x,y,lwd=2,lty=2,col=cpcol)
	text(maxx,maxy*.9,labels=paste('p=',signif(p,3)),adj=c(1,1),col=cpcol)
	plot(dat.tmp$dens, dat.tmp$tsc,col=cols[dat.tmp$Sp], pch=19, cex=0.8, bty='l',ylim=c(0,maxy),main=paste("mix",m,gs,"closed"))
	res <-pred(maxx, res.tmp[res.tmp$trt=="tsc"&res.tmp$s=="CLAAMO",]$m.1, res.tmp[res.tmp$trt=="tsc"&res.tmp$s=="CLAAMO",]$b)
	y <- res[[4]];x<-res[[3]]
	p <- res.tmp[res.tmp$trt=="tsc"&res.tmp$s=="CLAAMO",]$p
	lines(x,y,lwd=2,lty=2,col=cacol)
	text(maxx,maxy,labels=paste('p=',signif(p,3)),adj=c(1,1),col=cacol)
	res <-pred(maxx, res.tmp[res.tmp$trt=="tsc"&res.tmp$s=="CLAPUR",]$m.1, res.tmp[res.tmp$trt=="tsc"&res.tmp$s=="CLAPUR",]$b)
	y <- res[[4]];x<-res[[3]]
	p <- res.tmp[res.tmp$trt=="tsc"&res.tmp$s=="CLAPUR",]$p
	lines(x,y,lwd=2,lty=2,col=cpcol)
	text(maxx,maxy*.9,labels=paste('p=',signif(p,3)),adj=c(1,1),col=cpcol)
}
dev.off()












