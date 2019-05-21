
library(tidyverse)

img_file_path <- "C:/Users/Gandalf/Desktop/img_files/"
img_depth_files <- dir(path = img_file_path, pattern = "*.a.depth.txt")

img_depth <- img_depth_files %>%
  map_dfr(function(x) {
    read_tsv(file.path(img_file_path, x)) %>%
      mutate(taxon_oid = gsub(pattern = ".a.depth.txt", replacement = "", x = x))
  }) %>%
  mutate(analysis_project_id = gsub(pattern = "_.*", replacement = "", x = ID))

