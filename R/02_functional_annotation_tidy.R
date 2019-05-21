# Scaffold functional profiling using tidyverse

library(tidyverse)

# Set working directory to folder containing IMG download
setwd("")

# Enter IMG taxon OID code
taxon_oid <- "3300032269"

# Import IMG .depth.txt file containing mean depth of coverage per scaffold
img_depth <- read_tsv(file = paste0(taxon_oid, ".a.depth.txt"), quote = "",
                      col_names = TRUE, progress = TRUE)

# Import IMG .cog.txt file which contains gene product information
img_cog <- read_tsv(file = paste0(taxon_oid, ".a.cog.txt"),
                    col_names = c("gene_id", "cog_id", "percent_identity",
                                  "align_length", "query_start", "query_end",
                                  "subj_start", "subj_end", "evalue",
                                  "bit_score"), quote = "", progress = TRUE)

# Import IMG .ko.txt file which contains gene product information
img_ko <- read_tsv(file = paste0(taxon_oid, ".a.ko.txt"),
                   col_names = c("gene_id", "img_ko_flag", "ko_term",
                                 "percent_identity", "query_start", "query_end",
                                 "subj_start", "subj_end", "evalue", "bit_score",
                                 "align_length"), quote = "", progress = TRUE)

# Remove string "KO:" from KO term
img_ko <- img_ko %>%
  mutate(ko_term = gsub(pattern = "KO:", replacement = "", x = ko_term))

# Extract IMG scaffold IDs embedded in gene IDs
extract_scaffold_id <- function(gene_id) {
  scaffold_id <- substr(x = gene_id, start = 1, stop = nchar(img_depth$ID[1]))
  return(scaffold_id)
}
img_cog <- img_cog %>% mutate(scaffold_id = extract_scaffold_id(gene_id))
img_cog <- left_join(x = img_cog, y = img_depth, by = c("scaffold_id" = "ID"))
img_ko <- img_ko %>% mutate(scaffold_id = extract_scaffold_id(gene_id))
img_ko <- left_join(x = img_ko, y = img_depth, by = c("scaffold_id" = "ID"))

# Import COGs table
cogs <- read_csv(file = "C:/Users/Gandalf/Projects/metagenomicon/data/cogs.csv",
                 col_names = TRUE, trim_ws = TRUE)

# COG annotations
img_cog <- left_join(x = img_cog, y = cogs, by = "cog_id")

# Visualize COG profile
# img_cog %>%
#   ggplot(aes(x = functional_category, y = AvgFold, fill = cog_category)) +
#   geom_bar(stat = "identity")
