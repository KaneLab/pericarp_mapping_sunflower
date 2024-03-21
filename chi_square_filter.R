# modified from Ashley Barstow et al. 2022
# 10.3389/fpls.2022.1056278 

# install.packages("vcfR")
library(vcfR)

# setwd("dir")


# Load VCF file
vcf <- read.vcfR("PericarpHardness.vcf")


# Extract the genotype data into a matrix
gt_matrix <- extract.gt(vcf, element="GT")

# drop the parent HA467, the 160th index 
gt_matrix <- gt_matrix[, -160]




# Total number of markers before filtering
total_markers_before <- nrow(gt_matrix)
print(paste("Total markers before filtering:", total_markers_before))

# Define a chi-square test filter for expected proportions
chi_square_filter <- function(row) {
  # Count the number of homozygous and heterozygous genotypes
  homozygous_count0 <- sum(row %in% c("0|0", "0/0"))
  homozygous_count1 <- sum(row %in% c("1|1", "1/1"))
  heterozygous_count <- sum(row %in% c("1|0", "0|1", "1/0", "0/1"))
  # total, excluding any missingness
  total_count <- sum(homozygous_count0, homozygous_count1, heterozygous_count) # length(row) 
  # Calculate the expected count based on F6 inbred expectation
  expected_homozygous0 <- 0.485 * total_count
  expected_homozygous1 <- 0.485 * total_count
  expected_heterozygous <- 0.03 * total_count
  # Perform the chi-square test
  observed <- c(homozygous_count0, homozygous_count1, heterozygous_count)
  expected <- c(expected_homozygous0, expected_homozygous1, expected_heterozygous)
  chisq_test <- chisq.test(observed, p=expected/total_count)
  # Return TRUE if p-value > 0.10, indicating the marker should be retained
  return(chisq_test$p.value > 0.10)
}

# Apply the chi-square test to the genotype data
filter_results <- apply(gt_matrix, 1, chi_square_filter)
vcf_filtered_by_chisq <- vcf[filter_results, ]


# Total number of markers after filtering
total_markers_after <- sum(filter_results)
# Print summary statistics
cat("Total markers before filtering:", total_markers_before, "\n")
cat("Total markers after filtering: ", total_markers_after, "\n")
cat("Markers retained: ", total_markers_after, "/", total_markers_before, "(", round(total_markers_after / total_markers_before * 100, 2), "%)\n")



# Save filtered VCF
write.vcf(vcf_filtered_by_chisq, file = "chi_squared_filtered.vcf")



