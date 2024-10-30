# Load necessary libraries
library(GO.db)
library(AnnotationDbi)

# Function to get all descendants of a specific GO term
get_all_descendants <- function(goid) {
  all_goids <- keys(GO.db::GO.db, keytype = "GOID")
  descendants <- c()
  queue <- c(goid)  # Initialize the queue with the specified GO term
  visited <- c()    # Track visited terms
  
  while (length(queue) > 0) {
    current_term <- queue[1]
    queue <- queue[-1]  # Dequeue the first term
    
    if (!(current_term %in% visited)) {
      visited <- c(visited, current_term)  # Mark the term as visited
      
      # Get all children of the current term
      children_ids <- all_goids[sapply(all_goids, function(id) {
        current_term %in% as.list(GO.db::GOBPPARENTS)[[id]]
      })]
      
      if (length(children_ids) > 0) {
        descendants <- c(descendants, children_ids)
        queue <- c(queue, children_ids)  # Enqueue the children
      }
    }
  }
  
  return(unique(descendants))  # Return unique descendants
}

# Specify the GO term you want to analyze
target_goid <- "GO:0000003"
all_descendants <- get_all_descendants(target_goid)

# Retrieve the terms for these descendants
if (length(all_descendants) > 0) {
  descendant_terms <- AnnotationDbi::select(GO.db::GO.db, 
                                            keys = all_descendants, 
                                            columns = c("GOID", "TERM"), 
                                            keytype = "GOID")
} else {
  descendant_terms <- data.frame(GOID = character(0), TERM = character(0))
}

# Save to a new text file named after the GO term
output_file <- paste0("/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/all_descendants_of_", gsub(":", "_", target_goid), ".txt")
write.table(descendant_terms, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("File created successfully for", target_goid)
