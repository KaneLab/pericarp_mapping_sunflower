### filter out markers

setwd("dir")

library(dplyr)


# import 
data <- read.csv("filtered_pericarp.qtl.csv")


#subset data
df <- (data[,c(3:ncol(data))])
new_data <- data.frame(data[,c(1:2)])

df[1, ] <- as.numeric(df[1, ])

#get chromosome number
# n <- t(df[1,])
# n <- as.numeric(unique(n))
# top_chrom <- max(n)


for (i in 1:17){
  temp_chrom <- df[, df[1, ] == i, drop = FALSE]
  # remove 9 out of every 10 columns
  # columns_to_keep <- seq(1, ncol(temp_chrom), by=10)
  columns_to_keep <- seq(1, ncol(temp_chrom), by = 10)
  subset_chrom <- temp_chrom[, columns_to_keep]
  new_data <- cbind(new_data, subset_chrom)
  rm(temp_chrom, subset_chrom, columns_to_keep)
}
# new_data <- cbind(new_data, chrom_17)

# print stats
print(c("Markers before:", ncol(data)-2))
# print stats
print(c("Markers after:", ncol(new_data)-2))

#write filtered csv
write.table(new_data, "pericarp_subset.qtl.csv", quote = F, sep = ",", row.names = F, col.names = T, na = "")




