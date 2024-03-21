

########
# filter_duplicate_markers.R
#######
## filter out duplicate markers


# import  
full_data <- read.csv("pericarp.qtl.csv")


geno <- (full_data[,-c(1,2)])
#temp assignment
geno[1,] = 1
geno[2,] = 1


#num markers before filter
num_cols <- ncol(geno)


# check for duplicates
columns_to_drop <- as.data.frame(geno[,1])


size_before <- ncol(geno)
size_after <- 0
while (size_after < size_before) {
  size_before <- ncol(geno)
  for (i in 1:(size_before - 1)) {
    if (all(geno[, i] == geno[, i+1])) {
      duplicate_col <- names(geno)[i+1]
      columns_to_drop[duplicate_col] <- geno[, i+1]
      i = i+1
    }
  }
  geno <- geno[, !(names(geno))  %in% (names(columns_to_drop))]
  # columns_to_drop <- as.data.frame(geno[,1])
  size_after <- ncol(geno)
  print(size_after)
}


#drop duplicates
full_data <- full_data[, !(names(full_data))  %in% (names(columns_to_drop))]
# new_data <- full_data[, !(names(full_data))  %in% (names(columns_to_drop))]


#print results
print(c("Num markers before filtering:", num_cols))
print(c("Num markers after filtering:", ncol(full_data)-2))

#write filtered csv
write.table(full_data, "filtered_pericarp.qtl.csv", quote = F, sep = ",", row.names = F, col.names = T, na = "")

