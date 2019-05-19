# Scaffold functional profiling

# Set working directory to folder containing IMG download
setwd("")

# Enter IMG taxon OID code
taxon_oid <- "3300032269"

# # Import IMG .map.txt file which relates MEGAHIT contid ID with IMG scaffold ID
# img_map <- read.delim(file = paste0(taxon_oid, ".a.map.txt"), header = FALSE,
#                       sep = "\t", quote = "", stringsAsFactors = FALSE)
# colnames(img_map) <- c("megahit_contig_id", "img_scaffold_id")

# Import IMG .depth.txt file containing mean depth of coverage per scaffold
img_depth <- read.delim(file = paste0(taxon_oid, ".a.depth.txt"), header = TRUE,
                        sep = "\t", quote = "", stringsAsFactors = FALSE)

# Import IMG .cog.txt file which contains gene product information
img_cog <- read.delim(file = paste0(taxon_oid, ".a.cog.txt"), header = FALSE,
                      sep = "\t", quote = "", stringsAsFactors = FALSE)
colnames(img_cog) <- c("gene_id", "cog_id", "percent_identity", "align_length",
                       "query_start", "query_end", "subj_start", "subj_end",
                       "evalue", "bit_score")

# Import IMG .ko.txt file which contains gene product information
img_ko <- read.delim(file = paste0(taxon_oid, ".a.ko.txt"), header = FALSE,
                     sep = "\t", quote = "", stringsAsFactors = FALSE)
colnames(img_ko) <- c("gene_id", "img_ko_flag", "ko_term", "percent_identity",
                      "query_start", "query_end", "subj_start", "subj_end",
                      "evalue", "bit_score", "align_length")
img_ko$ko_term <- sapply(X = img_ko[,"ko_term"], FUN = function(ko_term) {
  ko <- sub(pattern = "KO:", replacement = "", x = ko_term)
  return(ko)
})  # Remove string "KO:" from KO term

# Extract IMG scaffold IDs embedded in gene IDs
extract_scaffold_id <- function(gene_id) {
  scaffold_id <- substr(x = gene_id, start = 1, stop = nchar(img_depth$ID[1]))
  return(scaffold_id)
}
img_cog$scaffold_id <- sapply(X = img_cog[,"gene_id"], FUN = extract_scaffold_id)
img_ko$scaffold_id <- sapply(X = img_ko[,"gene_id"], FUN = extract_scaffold_id)

