source("utils_tidy.R")


file_tx2gene <- "~/doc/tx2gene/tx2gene_v32.txt.gz"
tx2gene <- vroom(
  file_tx2gene, delim = "\t", col_types = cols(),
  col_names = c("chr", "gene_id", "gene_symbol", "transcript_id")
)
tx2gene_convert <- tx2gene |>
  dplyr::distinct(gene_id, .keep_all = T) |>
  dplyr::select(gene_id, gene_symbol)

truths <- readRDS("~/projects/AS/data/sim/rsem/truths.rds")

dir_out <- "~/projects/AS/analysis/novel_junction/"


##### rMATS #####
message("rMATS\n")

dir_rmats <- glue("{dir_out}/rMATS/")

list_rmats <- read_rmats(dir_rmats)
rmats <- bind_rows(list_rmats) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(as_type, ": ", ID))

cols_coord <- c(
  "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES",
  "downstreamEE", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
  "2ndExonEnd", "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
  "flankingES", "flankingEE", "riExonStart_0base", "riExonEnd"
)
rmats_e <- markEvent(
  result = rmats, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_7"
)
Sys.sleep(5)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_8"
)
Sys.sleep(5)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_9"
)
Sys.sleep(5)

saveRDS(rmats_e, glue("{dir_out}/rMATS/rmats.rds"))
Sys.sleep(5)


##### CASH #####
message("CASH\n")

dir_cash <- glue("{dir_out}/CASH/")

cash <- read_cash(
  glue("{dir_cash}/cash.s1vss2.alldiff.txt"),
  event_type = "all",
  "s1", "s2"
) |>
  left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(gene_symbol, ": ", Location))
cash_cc <- cash |>
  separate(Location, c("chr", "coord1", "coord2"), sep = "[:-]")

cols_coord <- c("coord1", "coord2")
cash_e <- markEvent(
  cash_cc, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
cash_e <- markEvent(
  cash_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
cash_e <- markEvent(
  cash_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

cash_e <- markEvent(
  cash_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
cash_e <- markEvent(
  cash_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
cash_e <- markEvent(
  cash_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

saveRDS(cash_e, glue("{dir_out}/CASH/cash.rds"))
Sys.sleep(5)


##### LeafCutter #####
message("LeafCutter\n")

dir_leafcutter <- glue("{dir_out}/LeafCutter/")

leafcutter <- read_leafcutter(
  glue("{dir_leafcutter}/ds/leafcutter_cluster_significance.txt"),
  glue("{dir_leafcutter}/ds/leafcutter_effect_sizes.txt"),
  "s1", "s2", flatten = F
) |>
  # dplyr::left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = cluster)
leafcutter$as_type <- "OTHER"
leafcutter2 <- read_leafcutter(
  glue("{dir_leafcutter}/ds/leafcutter_cluster_significance.txt"),
  glue("{dir_leafcutter}/ds/leafcutter_effect_sizes.txt"),
  "s1", "s2", flatten = T
) |>
  dplyr::left_join(tx2gene_convert, by = "gene_symbol") |>
  left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = cluster)
leafcutter2$as_type <- "OTHER"

cols_coord <- c("myintron")
leafcutter_e <- markEvent(
  leafcutter, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 3,
  error_distance = 1, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
leafcutter_e <- markEvent(
  leafcutter_e, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 3,
  error_distance = 5, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
leafcutter_e <- markEvent(
  leafcutter_e, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 3,
  error_distance = 10, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

leafcutter_e <- markEvent(
  leafcutter_e, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
leafcutter_e <- markEvent(
  leafcutter_e, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 4,
  error_distance = 5, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
leafcutter_e <- markEvent(
  leafcutter_e, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 4,
  error_distance = 10, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

saveRDS(leafcutter_e, glue("{dir_out}/LeafCutter/leafcutter.rds"))
Sys.sleep(5)


##### MAJIQ #####
message("MAJIQ\n")

dir_majiq <- glue("{dir_out}/MAJIQ/")

majiq <- read_majiq(
  glue("{dir_majiq}/deltapsi/s1_s2.deltapsi.tsv"),
  glue("{dir_majiq}/voila/s1_s2_0.2_showall.tsv"),
  condition_1 = "s1", condition_2 = "s2", flatten = F
) |>
  dplyr::mutate(padj = p) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = lsv_id)
majiq2 <- read_majiq(
  glue("{dir_majiq}/deltapsi/s1_s2.deltapsi.tsv"),
  glue("{dir_majiq}/voila/s1_s2_0.2.tsv"),
  condition_1 = "s1", condition_2 = "s2",
  flatten = T
) |>
  dplyr::mutate(padj = p) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = lsv_id)

cols_coord <- c(
  "exons_coords", "junctions_coords" # , "ir_coords"
)
majiq_e <- markEvent(
  majiq, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
majiq_e <- markEvent(
  majiq_e, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 5, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
majiq_e <- markEvent(
  majiq_e, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 10, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

majiq_e <- markEvent(
  majiq_e, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
majiq_e <- markEvent(
  majiq_e, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 5, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
majiq_e <- markEvent(
  majiq_e, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 10, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

saveRDS(majiq_e, glue("{dir_out}/MAJIQ/majiq.rds"))
Sys.sleep(5)


majiq_e2 <- markEvent(
  majiq2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
majiq_e2 <- markEvent(
  majiq_e2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 5, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
majiq_e2 <- markEvent(
  majiq_e2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 10, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

majiq_e2 <- markEvent(
  majiq_e2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
majiq_e2 <- markEvent(
  majiq_e2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 5, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
majiq_e2 <- markEvent(
  majiq_e2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 10, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

saveRDS(majiq_e2, glue("{dir_out}/MAJIQ/majiq2.rds"))
Sys.sleep(5)


##### Whippet #####
message("Whippet\n")

dir_whippet <- glue("{dir_out}/Whippet/")

whippet <- read_whippet(glue("{dir_whippet}/delta/s1_s2.diff2.gz")) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  dplyr::mutate(padj = as.character(1 - as.numeric(Probability))) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(Coord, Type, Node))

whippet_cc <- whippet |>
  separate(Coord, c("chr", "coord1", "coord2"), "[:-]")

cols_coord <- c("coord1", "coord2")
whippet_e <- markEvent(
  whippet_cc, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
whippet_e <- markEvent(
  whippet_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
whippet_e <- markEvent(
  whippet_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

whippet_e <- markEvent(
  whippet_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
whippet_e <- markEvent(
  whippet_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
whippet_e <- markEvent(
  whippet_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

saveRDS(whippet_e, glue("{dir_out}/Whippet/whippet.rds"))
Sys.sleep(5)


##### psichomics #####
message("psichomics\n")

dir_psichomics <- glue("{dir_out}/psichomics/")

psichomics <- read_psichomics(result_file = glue("{dir_psichomics}/results.txt")) |>
  left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = event)

psichomics_symbol <- psichomics[!grepl("ENSG", psichomics$gene_symbol), ]
psichomics_id <- psichomics[grepl("ENSG", psichomics$gene_symbol), ] |>
  dplyr::select(!gene_id) |>
  dplyr::rename(gene_id = gene_symbol) |>
  left_join(tx2gene_convert, by = "gene_id")

psichomics <- bind_rows(psichomics_symbol, psichomics_id)

psichomics_cc <- psichomics |>
  separate(event, c("event", "coord"), sep = "\\_[+-]\\_")

cols_coord <- c("coord")
psichomics_e <- markEvent(
  psichomics_cc, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_7"
)
Sys.sleep(5)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 5, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_8"
)
Sys.sleep(5)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 10, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_9"
)
Sys.sleep(5)

saveRDS(psichomics_e, glue("{dir_out}/psichomics/psichomics.rds"))
Sys.sleep(5)


##### JuncBASE #####
message("JuncBASE\n")

dir_juncbase <- glue("{dir_out}/JuncBASE/")

juncbase <- read_juncbase(result_file = glue("{dir_juncbase}/sample_set_comparison.txt")) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(
    myID = paste0(
      as_event_type, exclusion_junctions, inclusion_exons, exclusion_exons,
      inclusion_exons, `intron-exon_junctions`
    )
  )

cols_coord <- c(
  "exclusion_junctions", "inclusion_junctions", "exclusion_exons",
  "inclusion_exons", "intron-exon_junctions"
)
juncbase_e <- markEvent(
  juncbase, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 5, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 10, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_3"
)
Sys.sleep(5)

juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 5, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_5"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 10, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_6"
)
Sys.sleep(5)

juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_7"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 5, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_8"
)
Sys.sleep(5)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 10, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_9"
)
Sys.sleep(5)

saveRDS(juncbase_e, glue("{dir_out}/JuncBASE/juncbase.rds"))
Sys.sleep(5)


##### DARTS #####
dir_darts <- glue("{dir_out}/DARTS/")

list_darts_bht <- read_darts(glue("{dir_darts}/BHT/"))
darts_bht <- dplyr::bind_rows(list_darts_bht) |>
  dplyr::mutate(padj = as.character(1 - as.numeric(Posterior))) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(as_type, ": ", ID))

cols_coord <- c(
  "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES",
  "downstreamEE", "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
  "flankingES", "flankingEE", "riExonStart_0base", "riExonEnd"
)
darts_bht_e <- markEvent(
  result = darts_bht, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
)

saveRDS(darts_bht_e, glue("{dir_darts}/darts_bht.rds"))


list_darts_bht_dnn <- read_darts(glue("{dir_darts}/BHT_DNN/"))
darts_bht_dnn <- dplyr::bind_rows(list_darts_bht_dnn) |>
  dplyr::mutate(padj = as.character(1 - as.numeric(Posterior))) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(as_type, ": ", ID))

darts_bht_dnn_e <- markEvent(
  result = darts_bht_dnn, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
)

saveRDS(darts_bht_dnn_e, glue("{dir_darts}/darts_bht_dnn.rds"))


##### JUM #####
dir_jum <- glue("{dir_out}/JUM/")

list_jum <- read_jum(glue("{dir_jum}/JUM_diff/FINAL_JUM_OUTPUT_adjusted_pvalue_1/"))
jum <- bind_rows(list_jum) |>
  dplyr::left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = AS_event_ID)

# cols_coord <- c(
#   "upstream_exon_end_coor", "cassette_exon_start_coor", "cassette_exon_end_coor",
#   "downstream_exon_start_coor", "MXE_exon_coordinates", "common_5_SS_coor",
#   "A3SS_coordinates", "common_3_SS_coor", "A5SS_coordinates", "retained_intron_start",
#   "retained_intron_end", "Composite_coordinates"
# )

cols_coord <- c("AS_event_ID")
jum_e <- markEvent(
  result = jum, truths = truths, cols_coord = cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 3, "RI" = 2, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
jum_e <- markEvent(
  result = jum_e, truths = truths, cols_coord = cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 4),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths_2"
)
jum_e <- markEvent(
  result = jum_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 2, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)
jum_e <- markEvent(
  result = jum_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 4),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_4"
)

saveRDS(jum_e, glue("{dir_jum}/jum.rds"))




