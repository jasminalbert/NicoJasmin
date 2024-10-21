## Nico Burns & Jasmin Albert
## Initial sweep of GLTC & WW Clapur & Claamo 
## 10/12/2024

data_path <- "~/Documents/LabThings/NicoJasmin/Data/"
f1 <- paste0(data_path, "SeedsData - GLTCCLAAMO.csv")
f2 <- paste0(data_path, "SeedsData - WWCLAAMO.csv")
f3 <- paste0(data_path, "SeedsData - GLTCCLAPUR.csv")

stats <- data.frame()

filenames <- c(f1, f2, f3)

names(filenames) <- c("GLTCCLAAMO", "WWCLAAMO", "GLTCCLAPUR")

##lets do this
for (f in names(filenames)) {
  
  values <- read.csv(filenames[f], skip = 1)
  values$OseedsTotal <- (values$OseedsGood + values$OseedsBad)
  values$Omass <- as.numeric(values$Omass)
  Omass_Present <- values[!is.na(values$Omass)&!is.na(values$OseedsTotal),]
  
  trend = lm(Omass_Present$OseedsTotal ~ Omass_Present$Omass)
  trend
  summary(trend)
  str(trend)
  trend$coefficients
  class(trend$coefficients)
  names(trend$coefficients)
  m <- trend$coefficients["Omass_Present$Omass"]
  m
  b <- trend$coefficients["(Intercept)"]
  RangeOmass <- range(Omass_Present$Omass)
  x <- seq(0, 0.16, 0.01)
  y <- (m * x) + b
  
  str(summary(trend))
  p_val <- signif(summary(trend)$coefficients[2,4], digits = 3)
  rsq <- signif(summary(trend)$r.squared, digits = 3)
  
  i <- which(filenames == f)
  
  #filling data frame

  stats[f,"m"] <- m
  stats[f,"b"] <- b
  stats[f,"p_value"] <- p_val
  stats[f,"R_squared"] <- rsq
  
  #constructing graph, need help on how to access individual graphs
  
  plot(Omass_Present$Omass, Omass_Present$OseedsTotal, xlab = "Open Mass", ylab = "Open Seeds Total", main ="Mass vs Seeds")
  lines(x, y)
  text(0,max(Omass_Present$OseedsTotal) * .95, labels = paste("y =", round(m, digits = 2), "x +", round(b, digits = 2), "; p =", p_val, "; \nR^2 =", rsq), adj = 0)
  
  
}

stats

