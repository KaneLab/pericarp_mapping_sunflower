Repository for all R code used in the pericarp mapping project.
-
chi_square_filter.R : filtering for segregation distortion
map.R : create genetic map, perform QTL analysis (CIM), generate plots

VCF filtering and formatting are written to follow this pipeline
chi_square_filter.R --> vcf_to_rqtl.R --> filter_duplicate_markers.R --> subset_data.R
