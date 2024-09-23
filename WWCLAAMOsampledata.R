### First little thingy in R ###
data_location <- "./data/"
WW_Claamo_fname <- "SeedsData - WWCLAAMO.csv"
WW_Claamo_loc <- paste0(data_location, WW_Claamo_fname)
data = read.csv(WW_Claamo_loc, skip = 1)
data$OseedsTotal == (data$OseedsGood + data$OseedsBad)
data$OseedsTotal <- (data$OseedsGood + data$OseedsBad)
#class(data)
omass = data$Omass
is.na(omass)
!is.na(omass)
is.na(omass) == !is.na(omass)
all(is.na(omass) == !is.na(omass))
omass_present = data[!is.na(omass)&!is.na(data$OseedsTotal),]
trendline = lm(omass_present$OseedsTotal ~ omass_present$Omass)
summary(trendline)
str(trendline)
trendline$coefficients
#y=mx+b
class(trendline$coefficients)
names(trendline$coefficients)
m <- trendline$coefficients["omass_present$Omass"]
b <- trendline$coefficients["(Intercept)"]
rangeomass <- range(omass_present$Omass)
x <- seq(0, 0.16, 0.01)
y <- (m * x) + b

str(summary(trendline))
p_val <- signif(summary(trendline)$coefficients[2,4], digits = 3)

plot(omass_present$Omass, omass_present$OseedsTotal, xlab = "Open Mass", ylab = "Open Seeds Total", main = "WWClaamo Mass vs Seeds")
lines(x, y)
text(0,max(omass_present$OseedsTotal), labels = paste("y =", round(m, digits = 2), "x +", round(b, digits = 2), "; p =", p_val), adj = 0)
