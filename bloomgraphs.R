#bloom count graphs here
library(ggplot2)
library(dplyr)
library(tidyverse)

bloom_counts <- read.csv("~/Documents/LabThings/NicoJasmin/Data/vegHoliday.csv")
bloom_counts$FlowerperPlant <- (bloom_counts$NumFlowers / bloom_counts$NumPlant)
bloom_counts$Date <- as.Date(bloom_counts$Date)
bloom_counts <- bloom_counts[-(1:24), ]
bloom_counts


bloom_graph <- ggplot(bloom_counts, aes(x = Date, y = FlowerperPlant, color = PlantGenusSpecies)) +
  geom_point() +
  geom_line(aes(group = PlantGenusSpecies)) 

bloom_graph
