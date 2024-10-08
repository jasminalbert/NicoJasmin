## GLTC Clapur inital data sweep! ##
data_path <- "~/Documents/LabThings/NicoJasmin/Data/"
GLTC_Clapur_fname <- "SeedsData - GLTCCLAPUR.csv"
GLTC_Clapur_location <- paste0(data_path, GLTC_Clapur_fname)
the_data = read.csv(GLTC_Clapur_location, skip = 1)
the_data$OseedsTotal == (the_data$OseedsGood + the_data$OseedsBad)
the_data$OseedsTotal <- (the_data$OseedsGood + the_data$OseedsBad)
the_data$Omass <- as.numeric(the_data$Omass)
Omass = (the_data$Omass)
Omass_present = the_data[!is.na(Omass)&!is.na(the_data$OseedsTotal),]

trend = lm(Omass_present$OseedsTotal ~ Omass_present$Omass)
summary(trend)
str(trend)
trend$coefficients
class(trend$coefficients)
names(trend$coefficients)
m_1 <- trend$coefficients["Omass_present$Omass"]
b_1 <- trend$coefficients["(Intercept)"]
RangeOmass <- range(Omass_present$Omass)
x_1 <- seq(0, 0.16, 0.01)
y_1 <- (m_1 * x_1) + b_1

str(summary(trend))
p_value <- signif(summary(trend)$coefficients[2,4], digits = 3)
rsq <- signif(summary(trend)$r.squared, digits = 3)

plot(Omass_present$Omass, Omass_present$OseedsTotal, xlab = "Open Mass", ylab = "Open Seeds Total", main = "GLTCClapur Mass vs Seeds")
lines(x_1, y_1)
text(0,max(Omass_present$OseedsTotal) * .95, labels = paste("y =", round(m_1, digits = 2), "x +", round(b_1, digits = 2), "; p =", p_value, "; \nR^2 =", rsq), adj = 0)

