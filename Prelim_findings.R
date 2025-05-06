## WW & GLTC CLAPUR & COLLIGRA Prelim sweep
library(ggplot2)

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
GLTC_Clapur

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

Colligra <- rbind(GLTC_Colligra, WW_Colligra)
#seed conversion from mass (open samples)
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

Colligra

#dropping a column for WW clapur so i can join df's
WW_Clapur <- subset(WW_Clapur, select = -notesClosed)
GLTC_Clapur <- subset(GLTC_Clapur, select = -X)
GLTC_Clapur <- subset(GLTC_Clapur, select = -X.1)

WW_Clapur
GLTC_Clapur

Clapur <- rbind(GLTC_Clapur, WW_Clapur)

Clapur$PollinationStrat <- NA
for (i in 1:nrow(Clapur)) {
  if (Clapur$genSpec[i] == "genOnly") {
    Clapur$PollinationStrat[i] <- "Homogeneous"
  }
  else {
    Clapur$PollinationStrat[i] <- "Heterogeneous"
  }
}

Colligra$PollinationStrat <- NA
for (i in 1:nrow(Colligra)) {
  if (Colligra$genSpec[i] == "genOnly" | Colligra$genSpec[i] == "ss") {
    Colligra$PollinationStrat[i] <- "Homogeneous"
  }
  else {
    Colligra$PollinationStrat[i] <- "Heterogeneous"
  }
}
#graphs

Clapur_graph <- ggplot(Clapur, aes(x = PollinationStrat, y = polldiff )) +
  geom_boxplot(fill="purple", alpha=0.2) +
  geom_jitter(width = 0.2, height = 0, color = "black", size = 1) +
  ggtitle("Clapur Seed Difference (Open-Closed)") +
  xlab("PollinationStrategy") +
  ylab("Seed Difference (OseedsTotal - CseedsTotal")

Clapur_graph

Colligra_graph <- ggplot(Colligra, aes(x = PollinationStrat, y = polldiff )) +
  geom_boxplot(fill="skyblue", alpha=0.2) +
  geom_jitter(width = 0.2, height = 0, color = "black", size = 1) +
  ggtitle("Colligra Seed Difference (Open-Closed)") +
  xlab("PollinationStrategy") +
  ylab("Seed Difference (OseedsTotal - CseedsTotal")

Colligra_graph

t.test(Colligra$polldiff[Colligra$PollinationStrat == "Heterogeneous"], Colligra$polldiff[Colligra$PollinationStrat == "Homogeneous"])

