# Install Bioconductor and GO.db if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GO.db")

# Load necessary libraries
library(GO.db)
library(AnnotationDbi)

# Set the GO ID for Biological Process
biological_process_id <- "GO:0008150"

# Get all GO IDs
all_goids <- keys(GO.db::GO.db, keytype = "GOID")

# Get parents of all GO terms
parents <- as.list(GO.db::GOBPPARENTS)

# Identify direct children of Biological Process
children_ids <- all_goids[sapply(all_goids, function(goid) {
  biological_process_id %in% parents[[goid]]
})]

# Retrieve the terms for these direct children
direct_bp_terms <- AnnotationDbi::select(GO.db::GO.db, 
                                         keys = children_ids, 
                                         columns = c("GOID", "TERM"), 
                                         keytype = "GOID")

# Display the direct Biological Process terms
print(direct_bp_terms)

# Optionally, save to a text file
write.table(direct_bp_terms, file = "/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/direct_biological_process_terms.txt", sep = "\t", row.names = FALSE, quote = FALSE)
