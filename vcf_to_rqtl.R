#### convert vcf to rqtl format
# code modified from: Rim Gubaev, 2022
# https://github.com/RimGubaev/vcf_to_qtl/blob/master/README.md 


setwd("dir")

library(dplyr)


#prepearing a vcf for setting the alleles 
tx  <- readLines("chi_squared_filtered.vcf")

tx2  <- gsub(pattern = "#CHROM", replace = "CHROM", x = tx)
writeLines(tx2, con="cross.vcf.wohash.vcf")
rm(tx,tx2)

vcf <- read.table("cross.vcf.wohash.vcf", header = T, stringsAsFactors = F)
vcf[,c(10:ncol(vcf))] <- sapply(vcf[,c(10:ncol(vcf))], function(x)(
  substr(x, 1,3)
))




### reformat my vcf
### make new data frame, assigning marker IDs
key <- (vcf[,c(1,2,3)])

key$CHROM_num <- sub("Ha412HOChr", "", key$CHROM)
key$CHROM_num <- as.numeric(key$CHROM_num)

### marker number per chromosome
key <- key %>% 
  group_by(CHROM_num) %>%
  mutate(marker_num = row_number())

# write marker numbers to id column
key$ID <- paste0("c", key$CHROM_num, "m", key$marker_num)

# assign vcf marker ID 
vcf$ID <- key$ID
vcf$CHROM <- key$CHROM_num

# drop temp columns
key <- (key[,c(1,2,3)])

# write marker key to .csv file
write.csv(key, "marker_ID_key.csv", row.names=FALSE)
rm(key)




#setting the maternal genotype which will be set to As
parent.A = "HA467"

#replacing vcf encoded alleles to As, Bs, and Hs
vcf[,10:ncol(vcf)] <-  t(apply(vcf[,10:ncol(vcf)], 1, function(x)(
  ifelse(x == x[names(x) == parent.A], "A",
         ifelse(x == "0|1" | x  == "1|0", "H", 
                ifelse(x == ".|.", "-", "B")))
)))





# creating qtl formatted file
# extract ID, chromosome, position, and genotypes and transpose the vcf 
Rqtl <- t(vcf[,c(3,1,2,10:ncol(vcf))])
colnames(Rqtl) <- Rqtl[1,]

rm(vcf)

#approximation of physical distance to genetic distance (1cM ~ 1Mbp)
Rqtl[3,] <- round(as.numeric(Rqtl[3,]) / 1000000)

#formatting to qtl object
Rqtl <-  as.data.frame(Rqtl)
Rqtl <- data.frame(id = rownames(Rqtl), Rqtl, row.names = NULL)

### import the phonotype data set
pheno_data <- read.csv("pheno_data.csv")


# temp columns for merging data and retaining current row order
Rqtl <- data.frame(ID = rownames(Rqtl), Rqtl, row.names = NULL)
Rqtl$ID <- as.numeric(Rqtl$ID)
pheno_data <- data.frame(ID = rownames(pheno_data), pheno_data, row.names = NULL)
pheno_data$ID <- as.numeric(pheno_data$ID) + 3


# merge in the two phenotype data 
Rqtl <- merge(pheno_data, Rqtl, by="ID", all.y = TRUE)
#clean up
Rqtl$ID <- Rqtl$id.y
Rqtl <- (Rqtl[,-c(3,4)])

# clear the first three cells of each pheno column
Rqtl[c(2,3),1:2] <- ""
Rqtl[c(1),1] <- "id"
Rqtl[c(1),2] <- "pheno"

colnames(Rqtl) <- Rqtl[1,]

# write to file
write.table(Rqtl, "pericarp.qtl.csv", quote = F, sep = ",", row.names = F, col.names = F, na = "")


