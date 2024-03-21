# Sunflower Pericarp Trait Mapping Project
# Joey White
# July 2022 - 2024


# this is the complete script used for formatting and analyzing pericarp
# thickness data and strength data from 2022
# 

# 

# load packages
library(tidyverse)
library(tidyr)
library(dplyr)
library(lme4)

setwd("dir")


################################################################################
# Phenotype Data Analysis ------------------------------------------------------
################################################################################


####
# import from file
pericarp <- read.csv("pericarp_individual_combined.csv", header = TRUE)

par(mfrow = c(1,2))
# visualize distributions 
hist(pericarp$Mean_hardness, breaks = 15, main = "Pericarp Strength", xlab = "Force Required to Penetrate (N)", col = "darkslategray3")
# perent phenotype values
# abline(v = 2.62, col = "black", lty = 1, lwd = 3)
# abline(v = 4.49, col = "black", lty = 1, lwd = 3)

hist(pericarp$Avg_thickness, breaks = 15, main = "Pericarp Thickness", xlab = "Thickness (mm)", col = "coral3")
hist(pericarp$residuals, breaks = 30, main = "Pericarp Thickness", xlab = "Thickness (mm)", col = "coral3")

#HA 467 – mean 2.62
#PI 170415 – mean 4.49


summary(pericarp$Mean_hardness)
sd(pericarp$Mean_hardness)
IQR(pericarp$Mean_hardness)
var(pericarp$Mean_hardness)

### Plot hardness as a function of thickness
# linear model, hardness as a function of thickness
fit_ht <- lm(pericarp$Mean_hardness ~ pericarp$Avg_thickness)

#stats
modsum <- summary(fit_ht)
r2 = modsum$adj.r.squared
p = modsum$coefficients[2,4]

# labels
r2lab = bquote(italic(R)^2 == .(format(r2, digits = 3)))
plab = bquote(italic(p) == .(format(p, digits = 3)))

quartz()
# plot the two phenotypes
plot(1, type = "n", main = "Pericarp Strength as a Function of Thickness", 
     xlab = "Pericarp Thickness (mm)", 
     ylab = "Force to Penetrate (N)",
     xlim = c(0, .5),
     ylim = c(1.00, 5.50),)
points(pericarp$Mean_hardness ~ pericarp$Avg_thickness, pch = 20, cex = 0.65)

# add fit line
abline(fit_ht, lwd = 2.5, col = "firebrick")

#add stats
text(x = .03, y = 5.50, labels = r2lab)
text(x = .04, y = 5.25, labels = plab)




#### Do the same thing for RIL phenotype values instead of individuals

pericarp <- read.csv("pheno_data.csv", header = TRUE)
pericarp$hull_thickness <- pericarp$hull_thickness/100
pericarp$strength <- pericarp$strength/100

fit_ht <- lm(pericarp$strength ~ pericarp$hull_thickness)

#stats
modsum <- summary(fit_ht)
r2 = modsum$adj.r.squared
p = modsum$coefficients[2,4]

# labels
r2lab = bquote(italic(R)^2 == .(format(r2, digits = 3)))
plab = bquote(italic(p) == .(format(p, digits = 3)))

quartz()
plot(1, type = "n", main = "Pericarp Strength as a Function of Thickness", 
     xlab = "Pericarp Thickness (mm)", 
     ylab = "Force to Penetrate (N)",
     xlim = c(0, .4),
     ylim = c(1.00, 5.50),)
points(pericarp$strength ~ pericarp$hull_thickness, pch = 20, cex = 0.65)

# add fit line
abline(fit_ht, lwd = 2.5, col = "firebrick")

#add stats
text(x = .02, y = 5.50, labels = r2lab)
text(x = .025, y = 5.25, labels = plab)




