rm(list=c())
fig_loc <- "./figures/"
figseedspbloom <- paste0(fig_loc,"seedspbloom101624.pdf")
data <- readRDS("bloomseeds.RDS")
str(data)
names(data)[names(data)%in%c("spb.open","spb.closed")]<- c("open","closed")
data[data$M==0,]$genSpec <- "SS"
data <- data[-which(data$Sp=="CLAAMO" & data$M==6),]

data$S <- factor(data$S)
idvar <-c("S","St","Sp","Fc")
tblue <- rgb(0,0,0.5,0.5);tgreen <- rgb(0,0.5,0,0.5)
tcols <- c(tblue,tgreen)#data.frame(tcols[data$S],data$S)
sp <- unique(data$Sp)
pdf(figseedspbloom)
par(mar=c(1,1,1,0.5), oma=c(2,2,0.5,0),mgp=c(1,0.1,0), tcl=-0.1)
for (s in sp){
	dat.tmp <- data[data$Sp==s,]
	mat.sp <- dat.tmp[,c("open","closed")]
	mixes <- unique(dat.tmp$M)
	maxdat <- max(mat.sp, na.rm=T)
	par(mfrow=c(2,2))
	for (m in mixes){
		#dat.sp.m <- dat.tmp[dat.tmp$M==m,]
		#mat <- dat.sp.m[,c("spb.open","spb.closed")]
		mat <- mat.sp[dat.tmp$M==m,]
		p <- t.test(mat[,1],mat[,2], paired=T)$p.value
		gs <- unique(dat.tmp[dat.tmp$M==m,]$genSpec)
		site <- dat.tmp[dat.tmp$M==m,]$S
		sp.m.gs <- paste(s,m,gs, sep="-")
		dat.x <- jitter(rep(1:2,each=nrow(mat)),0.4)
		boxplot(mat, notch=T, xlab="seeds per bloom", ylim=c(0,maxdat))
		points(x=dat.x,y=unlist(mat), col=tcols[site], pch=19, cex=0.5)
		title(main=sp.m.gs, cex.main=1, line=-0.9)
		text(x=1.5,y=maxdat*.93,labels=paste0("p=",signif(p,3)))
	};title(ylab="seeds per bloom", xlab="pollination treatment", outer=T, line=0.3)
	dat.tmp$mgs <- apply(cbind(dat.tmp$M, dat.tmp$genSpec),1,paste,collapse="-")
	m.dat <- reshape(dat.tmp, direction="wide", idvar=idvar, timevar="mgs", drop=c("genSpec","id","open","closed","M"))
	m.mat <- m.dat[,!names(m.dat)%in%idvar]
	dat.x.m <- jitter(rep(1:ncol(m.mat), each=nrow(m.mat)),0.4)
	par(mfrow=c(1,1), xpd=NA)
	boxplot(m.mat)
	points(dat.x.m,unlist(m.mat), col=tcols[m.dat$S],pch=19, cex=0.5)
	title(main=s, cex.main=1, line=-0.9)
	title(xlab="mix id", ylab="open minus closed seeds per bloom", line=1)
};dev.off()

#anova(dif~genSpec, data=dat.tmp[!is.na(dat.tmp$dif),])
	#summary(aov(dat.tmp$dif~dat.tmp$genSpec))
	
	
#regression for site, neighbourhood density curve, neighbourhood richness	

