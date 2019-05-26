# Functional annotations for > 1 metagenome

library(tidyverse)

# Enter path to metagenomicon repo
setwd("PATH/metagenomicon/")


# Have the following IMG downloads for ALL samples saved in a single directory:
# 1) depth.txt
# 2) cog.txt

# Enter path to directory where IMG downloads are saved
img_file_path <- "PATH/img_files/"


# Import GOLD analysis projects summary
# Extract lake ID for projects containing ##-### patterns
# Note: comment out %>% mutate(lake_id ...) if working with non-LakePulse samples
# or with LakePulse co-assemblies
img_projects <- read_tsv(file = "data/gold_analysis_projects.txt", col_names = TRUE) %>%
  select(-`Add Date`) %>%
  rename(project_id = `GOLD Analysis Project ID`,
         project_name = `Analysis Project Name`,
         project_type = `Analysis Project Type`) %>%
  mutate(lake_id = str_extract(project_name, "\\d\\d\\-\\d\\d\\d"))


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
cog_table <- read_csv("data/cogs.csv", col_names = TRUE, trim_ws = TRUE) %>%
  select(-functional_category_abbr, -functional_category_abbr_char1)

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
  left_join(img_depth, by = "scaffold_id") %>%
  left_join(cog_table, by = "cog_id")


# You can now use img_cog for visualizing COG profiles
# Optionally, save img_cog as .rda file for use in future R scripts
#save(img_cog, file = "PATH/img_cog.rda")


# Create sample by COG matrix
sample_by_cog <- img_cog %>%
  group_by(project_id, cog_id) %>%
  summarize(avg_cov = sum(AvgFold)) %>%
  spread(cog_id, avg_cov, fill = 0)

sample_by_cog_mat <- sample_by_cog %>% as.data.frame
rownames(sample_by_cog_mat) <- sample_by_cog_mat$project_id
sample_by_cog_mat$project_id <- NULL


# You can now use sample_by_cog_mat for statistical analyses on COG profiles
# Optionally, save sample_by_cog_mat as .rda file for use in future R scripts
#save(sample_by_cog_mat, file = "PATH/sample_by_cog_mat.rda")


# Optional: hand-off to basic COG profile visualization
img_cog %>%
  group_by(project_id, cog_id, functional_category) %>%
  summarize(coverage = sum(AvgFold)) %>%
  ggplot(aes(x = project_id, y = coverage, fill = functional_category)) +
  geom_bar(stat = "identity", position = "fill")


# Optional: hand-off to basic COG profile ordination
# Transform data using method of choice (example here: Chord transformation)
library(vegan)
sample_by_cog_transformed <- decostand(x = sample_by_cog_mat, method = "normalize")
cogs_pca <- rda(sample_by_cog_transformed)
biplot(cogs_pca)
