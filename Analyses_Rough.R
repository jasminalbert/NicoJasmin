# rough draft of analyses---I'll clean it up later
library(ggplot2)
library(dplyr)
library(lme4)
#install.packages("lme4")
install.packages("lmerTest")
library(lmerTest)

#READ: All this stuff is just setting up data frames, can ignore unless one desires to see how i am cleaning & joining data


GLTC_Clapur <- read.csv("~/Documents/LabThings/NicoJasmin/Data/SeedsData - GLTCCLAPUR.csv", skip = 1)
GLTC_Clapur$OseedsTotal == (GLTC_Clapur$OseedsGood + GLTC_Clapur$OseedsBad)
GLTC_Clapur$CseedsTotal == (GLTC_Clapur$CseedsGood + GLTC_Clapur$CseedsBad)
GLTC_Clapur$OseedsTotal <- (GLTC_Clapur$OseedsGood + GLTC_Clapur$OseedsBad)
GLTC_Clapur$CseedsTotal <- (GLTC_Clapur$CseedsGood + GLTC_Clapur$CseedsBad)
GLTC_Clapur$Omass <- as.numeric(GLTC_Clapur$Omass)
GLTC_Clapur$Cmass <- as.numeric(GLTC_Clapur$Cmass)
Omass = (GLTC_Clapur$Omass)
Cmass = (GLTC_Clapur$Cmass)
Omass_present = GLTC_Clapur[!is.na(Omass)&!is.na(GLTC_Clapur$OseedsTotal),]
Cmass_present = GLTC_Clapur[!is.na(Cmass)&!is.na(GLTC_Clapur$CseedsTotal),]

Cmass_present
#open mass -> seed conversion
trend = lm(Omass_present$OseedsTotal ~ Omass_present$Omass)
summary(trend)
str(trend)
trend$coefficients
class(trend$coefficients)
names(trend$coefficients)
m <- trend$coefficients["Omass_present$Omass"]
b <- trend$coefficients["(Intercept)"]
RangeOmass <- range(Omass_present$Omass)
#x <- seq(0, 0.16, 0.01)
#y <- (m * x) + b

for (i in 1:nrow(GLTC_Clapur)) {
  
  if (is.na(GLTC_Clapur$OseedsTotal[i]) & !is.na(GLTC_Clapur$Omass[i])) {
    GLTC_Clapur$OseedsTotal[i] <- (m * (GLTC_Clapur$Omass[i])) + b
  }
}


#closed mass -> seed conversion
trend_C = lm(Cmass_present$CseedsTotal ~ Cmass_present$Cmass)
summary(trend_C)
str(trend_C)
trend_C$coefficients
class(trend_C$coefficients)
names(trend_C$coefficients)
m_C <- trend_C$coefficients["Cmass_present$Cmass"]
b_C <- trend_C$coefficients["(Intercept)"]
RangeCmass <- range(Cmass_present$Cmass)
#x_C <- seq(0, 0.16, 0.01)
#y_C <- (m * x) + b

for (i in 1:nrow(GLTC_Clapur)) {
  
  if (is.na(GLTC_Clapur$CseedsTotal[i]) & !is.na(GLTC_Clapur$Cmass[i])) {
    GLTC_Clapur$CseedsTotal[i] <- (m_C * (GLTC_Clapur$Cmass[i])) + b_C
  }
}
GLTC_Clapur

WW_Clapur <- read.csv("~/Documents/LabThings/NicoJasmin/Data/SeedsData - WWCLAPUR.csv", skip = 1)

WW_Clapur

GLTC_Colligra <- read.csv("~/Documents/LabThings/NicoJasmin/Data/SeedsData - GLTCCOLLIGRA.csv", skip = 1)
WW_Colligra <- read.csv("~/Documents/LabThings/NicoJasmin/Data/SeedsData - WWCOLLIGRA.csv", skip = 1)

for (i in 1:nrow(WW_Colligra)) {
  if (!is.na(WW_Colligra$CseedsGood[i]) & is.na(WW_Colligra$CseedsBad[i])) {
    WW_Colligra$CseedsBad[i] <- 0
  }
}

for (i in 1:nrow(WW_Colligra)) {
  if (!is.na(WW_Colligra$OseedsGood[i]) & is.na(WW_Colligra$OseedsBad[i])) {
    WW_Colligra$OseedsBad[i] <- 0
  }
}

#TODO: do empty inflorescence to seeds conversion for the open samples of WW & GLTC Colligra here




GLTC_Colligra$CseedsTotal <- (GLTC_Colligra$CseedsGood + GLTC_Colligra$CseedsBad)
WW_Colligra$CseedsTotal <- (WW_Colligra$CseedsGood + WW_Colligra$CseedsBad)
GLTC_Colligra$OseedsTotal <- (GLTC_Colligra$OseedsGood + GLTC_Colligra$OseedsBad)
WW_Colligra$OseedsTotal <- (WW_Colligra$OseedsGood + WW_Colligra$OseedsBad)

Cleandata <- function(spec) {
  # Remove rows with NA values in Oseeds or Cseeds
  spec <- spec[!is.na(spec$OseedsTotal) & !is.na(spec$CseedsTotal), ]
  
  spec$polldiff <- ((spec$OseedsTotal / spec$Oblooms) - (spec$CseedsTotal / spec$Cblooms))
  spec <- spec[!is.na(spec$polldiff), ]
  return(spec)
}

GLTC_Clapur <- Cleandata(GLTC_Clapur)
WW_Clapur <- Cleandata(WW_Clapur)
GLTC_Colligra <- Cleandata(GLTC_Colligra)
WW_Colligra <- Cleandata(WW_Colligra)

#dropping a column for WW clapur so i can join df's
WW_Clapur <- subset(WW_Clapur, select = -notesClosed)
GLTC_Clapur <- subset(GLTC_Clapur, select = -X)
GLTC_Clapur <- subset(GLTC_Clapur, select = -X.1)

WW_Clapur
GLTC_Clapur

####
Clapur <- rbind(GLTC_Clapur, WW_Clapur)
Colligra <- rbind(GLTC_Colligra, WW_Colligra)
####

cdata <- read.csv("~/Documents/LabThings/NicoJasmin/Data/bigData_workingCopy - CCDATA.csv")
specs <- c("COLLIGRA","CLAPUR")
cdata <- cdata[cdata$Sp%in%specs,]

idend <- which(names(cdata)=="Fc")
realstart <- which(names(cdata) == "S")
censtart <- which(names(cdata)=="CLAAMO")
cenend <- which(names(cdata)=="ERILAN")
cenmat <- cdata[,censtart:cenend]

cdata_no_bs <- cdata[,realstart:cenend]

select <- c(1:idend,censtart:cenend)
cdata <- cdata[,select]
iddat <- cdata[,1:idend]

cenmat <- apply(cenmat,2, as.numeric)
dens <- rowSums(cenmat, na.rm=T)
rich <- rowSums(cenmat>0, na.rm=T)

cendat <- cbind(iddat,cenmat,dens,rich)
cendat[is.na(cendat)] <- 0

Match <- function(spec, censusdat) {
  spec$dens <- NA
  spec$rich <- NA
  spec$CLAAMO <- NA
  spec$CLAPUR <- NA
  spec$COLLIGRA <- NA
  spec$COLLOGRA <- NA
  spec$GILCAP <- NA
  spec$MADGRA <- NA
  spec$MICGRA <- NA
  spec$EPIBRA <- NA
  spec$COLHET <- NA
  spec$PHAHET <- NA
  spec$clarkia <- NA
  spec$EPIANG <- NA
  spec$ERILAN <- NA
  for (i in 1:nrow(spec)) {
    for (j in 1:nrow(cendat)) {
      if (spec$St[i] == censusdat$St[j] & spec$M[i] == censusdat$M[j] & 
          spec$Sp[i] == censusdat$Sp[j] & spec$Fc[i] == censusdat$Fc[j]) {
        spec$dens[i] <- censusdat$dens[j]
        spec$rich[i] <- censusdat$rich[j]
        spec$CLAAMO[i] <- censusdat$CLAAMO[j]
        spec$CLAPUR[i] <- censusdat$CLAPUR[j]
        spec$COLLIGRA[i] <- censusdat$COLLIGRA[j]
        spec$COLLOGRA[i] <- censusdat$COLLOGRA[j]
        spec$GILCAP[i] <- censusdat$GILCAP[j]
        spec$MADGRA[i] <- censusdat$MADGRA[j]
        spec$MICGRA[i] <- censusdat$MICGRA[j]
        spec$EPIBRA[i] <- censusdat$EPIBRA[j]
        spec$COLHET[i] <- censusdat$COLHET[j]
        spec$PHAHET[i] <- censusdat$PHAHET[j]
        spec$clarkia[i] <- censusdat$clarkia[j]
        spec$EPIANG[i] <- censusdat$EPIANG[j]
        spec$ERILAN[i] <- censusdat$ERILAN[j]
        #that was a lot, sorry! i love brute forcing <3
      }
    }
  }
  return(spec)
}

Colligra <- Match(Colligra, cendat)
Clapur <- Match(Clapur, cendat)

Colligra
#delte clarkia column
Clapur

#cenCol <- cendat[cendat$Sp=="COLLIGRA",]


#get COlligra seeds from complete samples


Clapur$strat <- "het"
Clapur$strat[Clapur$genSpec=="genOnly"] <- "homo"

model_rough_clapur <- lm(polldiff ~ strat + CLAAMO + CLAPUR + COLLIGRA + COLLOGRA +
                              GILCAP + MADGRA + MICGRA + EPIBRA + COLHET + PHAHET + 
                              EPIANG, data = Clapur)

plot(Clapur$ERILAN, Clapur$polldiff)
#took out ERILAN

#formula
#y~x 
#y : polldif ~ strat + CLAAMO + CLAPUR + COLLIGRA + ....... + 
summary(model_rough_clapur)

model_rough0_clapur <- lm(polldiff ~ strat, data = Clapur)
summary(model_rough0_clapur)

model_rough2_clapur <- lm(polldiff ~ strat + dens + rich, data = Clapur)
summary(model_rough2_clapur)

anova(model_rough0_clapur, model_rough_clapur)

Clapur$St <- as.character(Clapur$St)

model_mixed_rough_clapur <- lmer(polldiff ~ strat + CLAAMO + CLAPUR + COLLIGRA + COLLOGRA +
                           GILCAP + MADGRA + MICGRA + EPIBRA + PHAHET  
                         + (1 | St) + (1 | Fc), data = Clapur)

isSingular(model_mixed_rough_clapur, tol = 1e-4)

summary(model_mixed_rough_clapur)
str(summary(model_mixed_rough_clapur))
anova(model_mixed_rough_clapur)
colSums(Clapur[which(names(Clapur)=="CLAAMO"):which(names(Clapur)=="EPIANG")])

var(Clapur$St)
class(Clapur$S)


