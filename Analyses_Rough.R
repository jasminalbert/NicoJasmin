# rough draft of analyses---I'll clean it up later
library(ggplot2)
library(dplyr)
library(lme4)
#install.packages("lme4")
#install.packages("lmerTest")
library(lmerTest)
library(performance)
library(tidyverse)
library(modelr)

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
Colligra <- rbind(GLTC_Colligra, WW_Colligra)

Colligra$CseedsTotal <- Colligra$CseedsGood + Colligra$CseedsBad
Colligra_complete <- Colligra[Colligra$Ocomplete. == 1,]

Colligra_complete

Colligra_complete$OseedsTotal <- Colligra_complete$OseedsGood + Colligra_complete$OseedsBad

#open mass -> seed conversion
trend_1 = lm(Colligra_complete$OseedsTotal ~ Colligra_complete$OemptyInflorMass)
summary(trend_1)
str(trend_1)
trend_1$coefficients
class(trend_1$coefficients)
names(trend_1$coefficients)
m_1 <- trend_1$coefficients["Colligra_complete$OemptyInflorMass"]
b_1 <- trend_1$coefficients["(Intercept)"]
RangeOmass_1 <- range(Colligra_complete$OemptyInflorMass)
#y <- (m * x) + b

for (i in 1:nrow(Colligra)) {
  Colligra$OseedsTotal[i] <- (m_1 * (Colligra$OemptyInflorMass[i])) + b_1
}


Cleandata <- function(spec) {
  # Remove rows with NA values in Oseeds or Cseeds
  spec <- spec[!is.na(spec$OseedsTotal) & !is.na(spec$CseedsTotal), ]
  
  spec$polldiff <- ((spec$OseedsTotal / spec$Oblooms) - (spec$CseedsTotal / spec$Cblooms))
  spec <- spec[!is.na(spec$polldiff), ]
  return(spec)
}

GLTC_Clapur <- Cleandata(GLTC_Clapur)
WW_Clapur <- Cleandata(WW_Clapur)
Colligra <- Cleandata(Colligra)

#dropping a column for WW clapur so i can join df's
WW_Clapur <- subset(WW_Clapur, select = -notesClosed)
GLTC_Clapur <- subset(GLTC_Clapur, select = -X)
GLTC_Clapur <- subset(GLTC_Clapur, select = -X.1)

WW_Clapur
GLTC_Clapur

####
Clapur <- rbind(GLTC_Clapur, WW_Clapur)
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

Colligra$strat <- NA
for (i in 1:nrow(Colligra)) {
  if (Colligra$genSpec[i] == "specOnly" | Colligra$genSpec[i] == "ss") {
    Colligra$strat[i] <- "homo"
  }
  else {
    Colligra$strat[i] <- "het"
  }
}
####################################
#Analysis starts here more or less
####################################
#model_rough_clapur <- lm(polldiff ~ strat + CLAAMO + CLAPUR + COLLIGRA + COLLOGRA +
                              # GILCAP + MADGRA + MICGRA + EPIBRA + COLHET + PHAHET + 
                              # EPIANG, data = Clapur)

#plot(Clapur$ERILAN, Clapur$polldiff)
#took out ERILAN

#formula
#y~x 
#y : polldif ~ strat + CLAAMO + CLAPUR + COLLIGRA + ....... + 
#summary(model_rough_clapur)

# model_rough0_clapur <- lm(polldiff ~ strat, data = Clapur)
# summary(model_rough0_clapur)
# 
# model_rough2_clapur <- lm(polldiff ~ strat + dens + rich, data = Clapur)
# summary(model_rough2_clapur)
# 
# anova(model_rough0_clapur, model_rough_clapur)
#
#model_mixed_rough_clapur <- lmer(polldiff ~ strat + (1 | S),  data = Clapur)
#
#model_mixed_rough_colligra <- lmer(polldiff ~ strat + (1 | S),  data = Colligra)
#
# isSingular(model_mixed_rough_clapur, tol = 1e-4)
# 
# summary(model_mixed_rough_clapur)
# str(summary(model_mixed_rough_clapur))
# anova(model_mixed_rough_clapur)
# colSums(Clapur[which(names(Clapur)=="CLAAMO"):which(names(Clapur)=="EPIANG")])
# 
# var(Clapur$St)
# class(Clapur$S)
#  
# check_model(model_mixed_rough_clapur)
# 
# summary(model_mixed_rough_colligra)
# check_model(model_mixed_rough_colligra)

#TODO: change stand numbers to be (1,2,3,4) instead of (700, 702, 704 etc.)?
#gonna try and use mix ($m) as a nested random effect first...

Clapur_v1 <- Clapur
Colligra_v1 <- Colligra


Clapur_v1$St <- recode(Clapur_v1$St,'700' = 1,'604' = 1,'702' = 2,'605' = 2,'704' = 3,'607' = 3,'705' = 4,'611' = 4)
Colligra_v1$St <- recode(Colligra_v1$St,'700' = 1,'604' = 1,'702' = 2,'605' = 2,'704' = 3,'607' = 3,'705' = 4,'611' = 4)


Colligra_v1
Clapur_v1
Clapur_v1$St <- as.character(Clapur_v1$St)
Colligra_v1$St <- as.character(Colligra_v1$St)
Clapur_v1$St <- as.factor(Clapur_v1$St)
Colligra_v1$St <- as.factor(Colligra_v1$St)

Clapur$St <- as.character(Clapur$St)
Colligra$St <- as.character(Colligra$St)
Clapur$St <- as.factor(Clapur$St)
Colligra$St <- as.factor(Colligra$St)
##

################
#model_mixed_rough_clapur_v1 <- lmer(polldiff ~ strat +  (1 | S / St),  data = Clapur_v1)
# summary(model_mixed_rough_clapur_v1)
# str(summary(model_mixed_rough_clapur_v1))
# anova(model_mixed_rough_clapur_v1)
# check_model(model_mixed_rough_clapur_v1)
##
model_mixed_rough_clapur_v2 <- lmer(polldiff ~ strat +  (1 | St ),  data = Clapur)
check_model(model_mixed_rough_clapur_v2)
##
#model_mixed_rough_clapur_v3 <- lmer(polldiff ~ strat +  (1 | St ),  data = Clapur_v1)
#model_mixed_rough_clapur_v4 <- lmer(polldiff ~ strat +  (1 | S) + (1 | St),  data = Clapur_v1)
#check_model(model_mixed_rough_clapur_v4)
################

# rand(model_mixed_rough_clapur_v1)
# rand(model_mixed_rough_clapur_v4)

# hist(Clapur_v1$polldiff)
# hist(Colligra_v1$polldiff)
#######

#model_mixed_rough_colligra_v1 <- lmer(polldiff ~ strat + (1 | S / St),  data = Colligra_v1)
#check_model(model_mixed_rough_colligra_v1)
##
model_mixed_rough_colligra_v2 <- lmer(polldiff ~ strat + (1 | St),  data = Colligra)
##
check_model(model_mixed_rough_colligra_v2)
summary(model_mixed_rough_colligra_v2)

#model_mixed_rough_colligra_v3 <- lmer(polldiff ~ strat + (1 | St),  data = Colligra_v1)
#check_model(model_mixed_rough_colligra_v3)

#model_mixed_rough_colligra_v4 <- lmer(polldiff ~ strat +  (1 | S) + (1 | St),  data = Colligra_v1)

#graphs perchance?
################

predicted_vals_clapur_og <- predict(model_mixed_rough_clapur_v2, newdata = NULL, re.form = NULL)
predicted_vals_colligra_og <- predict(model_mixed_rough_colligra_v2, newdata = NULL, re.form = NULL)

newdata <- data.frame(strat=sample(as.factor(c("homo","het")), 138, replace=T), St=as.factor(sample(levels(Clapur$St),138, replace=T)))
predicted_vals_clapur <- predict(model_mixed_rough_clapur_v2, newdata = newdata, re.form = NULL)
pred <- cbind(predicted_vals_clapur,newdata)
boxplot(predicted_vals_clapur~strat, data=pred)
boxplot(predicted_vals_clapur_og~Clapur$strat)

newdata_colligra <- data.frame(strat=sample(as.factor(c("homo","het")), 138, replace=T), St=as.factor(sample(levels(Colligra$St),138, replace=T)))
predicted_vals_colligra <- predict(model_mixed_rough_clapur_v2, newdata = newdata, re.form = NULL)
pred_colligra <- cbind(predicted_vals_colligra, newdata_colligra)
boxplot(predicted_vals_colligra~strat, data=pred_colligra)
boxplot(predicted_vals_colligra_og~Colligra$strat)

predicted_vals_clapur
predicted_vals_colligra

#base graphs are here
Clapur_graph <- ggplot(Clapur, aes(x = strat, y = polldiff)) +
  geom_boxplot(fill="purple", alpha=0.2) +
  geom_jitter(width = 0.2, height = 0, color = "black", size = 1) +
  geom_line(aes(y = pred$predicted_vals_clapur, linewidth = 0.5, colour = "red")) +
  ggtitle("Clapur Seed Difference (Open-Closed)") +
  xlab("PollinationStrategy") +
  ylab("Seed Difference (OseedsTotal - CseedsTotal")

Clapur_graph

Colligra_graph <- ggplot(Colligra, aes(x = strat, y = polldiff)) +
  geom_boxplot(fill="skyblue", alpha=0.2) +
  geom_jitter(width = 0.2, height = 0, color = "black", size = 1) +
  geom_line(aes(y = pred_colligra$predicted_vals_colligra, linewidth = 0.5, colour = "red")) +
  ggtitle("Colligra Seed Difference (Open-Closed)") +
  xlab("PollinationStrategy") +
  ylab("Seed Difference (OseedsTotal - CseedsTotal")

Colligra_graph
#plot(Colligra$polldiff[Colligra$strat=="het"], Colligra$polldiff[Colligra$strat=="homo"])

meanperst <- Colligra %>% group_by(strat, St) %>% summarize(mean(polldiff, na.rm=T))
homomean <- unlist(meanperst[meanperst$strat=="homo",3])
hetmean <- unlist(meanperst[meanperst$strat=="het",3])
plot(homomean,hetmean, col=meanperst$St, pch=19, ylim=c(0,5), xlim=c(0,5)) 
cor(homomean,hetmean)
lines(0:5,0:5)
legend("topleft", legend=levels(meanperst$St), fill=meanperst$St)
