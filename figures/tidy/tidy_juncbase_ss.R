source("utils_tidy.R")

truths <- readRDS("~/projects/AS/data/sim/rsem/truths.rds")

dir_out <- "~/projects/AS/analysis/novel_splice_site/"



##### JuncBASE #####
message("JuncBASE\n")

dir_juncbase <- glue("{dir_out}/JuncBASE/")

juncbase <- readRDS(glue("{dir_juncbase}/juncbase_c.rds"))

cols_coord <- c(
  "exclusion_junctions", "inclusion_junctions", "exclusion_exons",
  "inclusion_exons", "intron-exon_junctions"
)
juncbase_e <- markEvent(
  juncbase, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_7"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_8"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_9"
)
Sys.sleep(5)

saveRDS(juncbase_e, glue("{dir_out}/JuncBASE/juncbase2.rds"))
Sys.sleep(5)




