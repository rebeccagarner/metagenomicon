# Functional annotations for > 1 metagenome

library(tidyverse)

# Enter path to metagenomicon repo
setwd("PATH/metagenomicon/")


# Enter path to directory where IMG downloads are saved
img_file_path <- "C:/Users/Gandalf/Desktop/img_files/"


# Import and combine IMG depth files
# Assign taxon OID and GOLD analysis project ID
(img_depth_files <- dir(path = img_file_path, pattern = "*.a.depth.txt"))

img_depth <- img_depth_files %>%
  map_dfr(function(x) {
    read_tsv(file.path(img_file_path, x)) %>%
      mutate(taxon_oid = gsub(pattern = ".a.depth.txt", replacement = "", x = x))
  }) %>%
  rename(scaffold_id = ID) %>%
  mutate(project_id = gsub(pattern = "_.*", replacement = "", x = scaffold_id))


# Extract IMG scaffold ID string lengths
img_scaffold_ids <- img_depth %>%
  distinct(taxon_oid, project_id, .keep_all = TRUE) %>%
  mutate(scaffold_id_length = str_length(scaffold_id)) %>%
  select(-c(scaffold_id, AvgFold))


# Import and combine IMG COG files
(img_cog_files <- dir(path = img_file_path, pattern = "*.a.cog.txt"))

img_cog <- img_cog_files %>%
  map_dfr(function(x) {
    read_tsv(file.path(img_file_path, x),
             col_names = c("gene_id", "cog_id", "percent_identity", "align_length",
                           "query_start", "query_end", "subj_start", "subj_end",
                           "evalue", "bit_score"),
             cols_only(gene_id = col_character(), cog_id = col_character())) %>%
      mutate(taxon_oid = gsub(pattern = ".a.cog.txt", replacement = "", x = x))
  }) %>%
  left_join(img_scaffold_ids, by = "taxon_oid") %>%
  mutate(scaffold_id = substr(gene_id, start = 1, stop = scaffold_id_length)) %>%
  select(-c(taxon_oid, project_id, scaffold_id_length)) %>%
  left_join(img_depth, by = "scaffold_id")


# Import GOLD analysis projects summary
img_projects <- read_tsv(file = "data/gold_analysis_projects.txt", col_names = TRUE) %>%
  select(-`Add Date`) %>%
  rename(project_id = `GOLD Analysis Project ID`,
         project_name = `Analysis Project Name`,
         project_type = `Analysis Project Type`)
