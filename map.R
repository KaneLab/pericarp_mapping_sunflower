#genetic map for pericarp RIL population


library(dplyr)
library(ggplot2)
library(qtl)
setwd("dir")

b <- read.cross("csv", "", "pericarp_subset.qtl.csv", estimate.map=FALSE, na.strings=c("H"), genotypes = c('A', 'B'))
b <- convert2riself(b)



# estimate geno probs
b <- calc.genoprob(b, step=10, error.prob=0.01)
summary(b)
summaryMap(b)



#### Composite interval mapping
#sit
b.c.perm <- cim(b, pheno.col = 4, n.perm = 100, map.function="kosambi")

b.c <- cim(b, pheno.col = 4, map.function="kosambi")
plot(b.c, col=c("darkgreen"), main = "a. Strength Independent of Thickness", ylab="LOD score")
abline(h = 6.59, col = "grey", lty = 2)

#strength
b.c_s <- cim(b, pheno.col = 3, map.function="kosambi")
plot(b.c_s, col=c("darkgreen"), main = "b. Strength", ylab="LOD score")
abline(h = 6.42, col = "grey", lty = 2)

#thick
b.c_t <- cim(b, pheno.col = 2, map.function="kosambi")
plot(b.c_t, col=c("darkgreen"), main = "c. Thickness", ylab="LOD score")
abline(h = 6.59, col = "grey", lty = 2)

summary(b.c)
summary(b.c_s)
summary(b.c_t)

b.c.perm <- cim(b, pheno.col = 4, n.perm = 100, map.function="kosambi")
b.c_s.perm <- cim(b, pheno.col = 3, n.perm = 100, map.function="kosambi")
b.c_t.perm <- cim(b, pheno.col = 2, n.perm = 100, map.function="kosambi")

### Get LOD significance values
summary(b.c.perm, alpha=0.05)
summary(b.c_s.perm, alpha=0.05)
summary(b.c_t.perm, alpha=0.05)

# pull out chromosomes
ch5 <- b.c[which(b.c$chr == 5),]
ch16 <- b.c_t[which(b.c_t$chr == 16),]

ch5_s <- b.c_s[which(b.c_s$chr == 5),]
ch16_s <- b.c_s[which(b.c_s$chr == 16),]



### stack chromosome manhattan plots
plot(1, type = "n", main = "Chromosome 5", 
     xlab = "Map position (cM)", 
     ylab = "LOD score",
     xlim = c(0, 185),
     ylim = c(0, 17),)
lines(ch5$pos, ch5$lod, col=c("purple4"), lwd = 2)
lines(ch5_s$pos, ch5_s$lod, col=c("firebrick"), lwd = 2)
# points(ch5$pos, ch5$lod, pch = 20, col = "purple4")
# points(ch5_s$pos, ch5_s$lod, pch = 20, col = "firebrick")
abline(h = 6.42, col = "grey", lty = 2)
legend(0, 17, legend = c("S     T", "Strength"), fill = c("purple4", "firebrick"))

plot(1, type = "n", main = "Chromosome 16", 
     xlab = "Map position (cM)", 
     ylab = "LOD score",
     xlim = c(0, 215),
     ylim = c(0, 20),)
lines(ch16$pos, ch16$lod, col=c("darkgreen"), lwd = 2)
lines(ch16_s$pos, ch16_s$lod, col=c("firebrick"), lwd = 2)
abline(h = 6.59, col = "grey", lty = 2)
legend(0, 20, legend = c("Thickness", "Strength"), fill = c("darkgreen", "firebrick"))





## pull out all markers above the significance thesholds

high_lods5 <- ch5[which(ch5$lod > 6.59),]
qtl_5 <- ch5[which(ch5$lod > 14),]

high_lods5_s <- ch5_s[which(ch5_s$lod > 6.42),]
qtl_5_s <- ch5_s[which(ch5_s$lod > 6),]


high_lods16 <- ch16[which(ch16$lod > 6.59),]
qtl_16 <- ch16[which(ch16$lod > 18),]

high_lods16_s <- ch16_s[which(ch16_s$lod > 6.42),]


summary(ch5)
summary(ch5_s)

summary(ch9)
summary(ch16)
summary(ch16_s)









## Interval mapping
#residuals phenotype
b.em <- scanone(b, pheno.col = 4, method="em")
plot(b.em, col=c("red"), main = "Residuals Phenotype (Strength Independent of Thickness)", ylab="LOD score")
#strength
b.em_s <- scanone(b, pheno.col = 3, method="em")
plot(b.em_s, col=c("red"), main = "Strength (Penetrometer)", ylab="LOD score")
#thick
b.em_t <- scanone(b, pheno.col = 2, method="em")
plot(b.em_t, col=c("red"), main = "Thickness", ylab="LOD score")


### Haley-Knott
#residuals phenotype
b.hk <- scanone(b, pheno.col = 4, method="hk")
plot(b.hk, col=c("blue"), main = "Residuals Phenotype (Strength Independent of Thickness)", ylab="LOD score")

#both together
plot(b.hk, b.em, col=c("blue","red"), main = "", ylab="LOD score")

#strength
b.hk_s <- scanone(b, pheno.col = 3, method="hk")
plot(b.hk_s, col=c("blue","red"), main = "Strength (Penetrometer)", ylab="LOD score")
#thick
b.hk_t <- scanone(b, pheno.col = 2, method="hk")
plot(b.hk_t, col=c("blue","red"), main = "Thickness", ylab="LOD score")

summary(b.hk)
summary(b.hk_s)
summary(b.hk_t)





### estimate map, genetic distance
em <- est.map(b, error.prob = 0.0001, map.function = "kosambi", verbose = TRUE)

b_em<-replace.map(b,em)

# write.cross(b_em, format = "csv", filestem = "estimated_map")
# est_map <- read.csv("estimated_map.csv")

em.c <- cim(em, pheno.col = 4, map.function="kosambi")

em_r <- est.rf(b_em)
plot.rf(hyper, alternate.chrid=TRUE)


em2 <- est.map()
summary(em)
plotMap(em)

plot(em, b$pheno)












