# Functional annotations for > 1 metagenome

library(tidyverse)

# Enter path to directory where IMG downloads are saved
img_file_path <- "C:/Users/Gandalf/Desktop/img_files/"

img_depth_files <- dir(path = img_file_path, pattern = "*.a.depth.txt")

# Import and combine all IMG depth files into 1 data frame
# Assign taxon OID and GOLD analysis project ID
img_depth <- img_depth_files %>%
  map_dfr(function(x) {
    read_tsv(file.path(img_file_path, x)) %>%
      mutate(taxon_oid = gsub(pattern = ".a.depth.txt", replacement = "", x = x))
  }) %>%
  mutate(analysis_project_id = gsub(pattern = "_.*", replacement = "", x = ID))


######### For scaffold ID extraction, make sure to include info that scaffold ids may be of differing length between analysis projects


# Import IMG .cog.txt file which contains gene product information
img_cog <- read_tsv(file = paste0(taxon_oid, ".a.cog.txt"),
                    col_names = c("gene_id", "cog_id", "percent_identity",
                                  "align_length", "query_start", "query_end",
                                  "subj_start", "subj_end", "evalue",
                                  "bit_score"), quote = "", progress = TRUE)

