library(rtracklayer)
source("../utils_tidy.R")


file_tx2gene <- "~/doc/tx2gene/tx2gene_v32.txt.gz"
tx2gene <- vroom(
  file_tx2gene, delim = "\t", col_types = cols(),
  col_names = c("chr", "gene_id", "gene_symbol", "transcript_id")
)
tx2gene_convert <- tx2gene |>
  dplyr::distinct(gene_id, .keep_all = T) |>
  dplyr::select(gene_id, gene_symbol)

gtf <- import("~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf")
gr <- gtf[gtf$type == "gene"]

truths <- readRDS("~/projects/AS/data/sim/rsem/truths.rds")

dir_out <- "~/projects/as/analysis/revision2/"


##### BRIE2 #####
list_brie2 <- read_brie2(glue("{dir_out}/BRIE2/"), filter = F)
brie2 <- bind_rows(list_brie2) |>
  dplyr::rename(myID = gene_id) |>
  myAnnot(gr, cols_loc = c("loc_1", "loc_2", "loc_3", "loc_4"), T) |>
  tidyr::drop_na(padj)

cols_coord <- c("loc_1", "loc_2", "loc_3", "loc_4")
brie2_e <- markEvent(
  brie2, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
)
brie2_e <- markEvent(
  brie2_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths2", verbose = T
) # best
brie2_e <- markEvent(
  brie2_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths3", verbose = T
)

saveRDS(brie2_e, glue("{dir_out}/BRIE2/brie2.rds"))


##### SingleSplice #####
singlesplice <- read_singlesplice(glue("{dir_out}/SingleSplice/"))

singlesplice_annot <- singlesplice |>
  myAnnot(gr, cols_loc = "loc", T) |>
  dplyr::mutate(
    padj = case_when(significant == "yes" ~ "0.01", T ~ "1"),
    myID = asm_id
  )

cols_coord <- c("loc")
singlesplice_annot_e <- markEvent(
  singlesplice_annot, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
) # best
singlesplice_annot_e <- markEvent(
  singlesplice_annot_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths2", verbose = T
)

saveRDS(singlesplice_annot_e, glue("{dir_out}/SingleSplice/singlesplice.rds"))


##### Psix #####
psix <- read_psix(
  glue("{dir_out}/Psix/psix_results.tab"),
  "~/doc/psix/psix_annotation.tab"
) |>
  dplyr::left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = event, as_type = "OTHER")

cols_coord <- c("myintron")
psix_e <- markEvent(
  psix, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = 3,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
)
psix_e <- markEvent(
  psix_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths2", verbose = T
) # best

saveRDS(psix_e, glue("{dir_out}/Psix/psix.rds"))


##### DESJ #####
desj <- read_desj(glue("{dir_out}/DESJ/first.res.xls")) |>
  tidyr::drop_na(padj) |>
  dplyr::left_join(tx2gene_convert, by = "gene_symbol") |>
  dplyr::mutate(myID = junction_id, as_type = "OTHER")

cols_coord <- c("loc")
desj_e <- markEvent(
  desj, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = 2,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
) # best
desj_e <- markEvent(
  desj_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = 2,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths2", verbose = T
)

saveRDS(desj_e, glue("{dir_out}/DESJ/desj.rds"))


##### Outrigger #####
outrigger <- read_outrigger(glue("{dir_out}/Outrigger/result.txt")) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = ID)

outrigger_annot <- outrigger |>
  myAnnot(gr, cols_loc = c("loc_1", "loc_2", "loc_3", "loc_4"), T)

cols_coord <- c("loc_1", "loc_2", "loc_3", "loc_4")
outrigger_annot_e <- markEvent(
  outrigger_annot, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
)
outrigger_annot_e <- markEvent(
  outrigger_annot_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths2", verbose = T
)
outrigger_annot_e <- markEvent(
  outrigger_annot_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths3", verbose = T
) # best

saveRDS(outrigger_annot_e, glue("{dir_out}/Outrigger/outrigger.rds"))


##### MicroExonator #####
microexonator <- read_microexonator(
  glue("{dir_out}/MicroExonator/A_vs_B.all_nodes.microexons.txt")
) |>
  dplyr::left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(Coord, Type, Node))

cols_coord <- c("Coord")
microexonator_e <- markEvent(
  microexonator, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths", verbose = T
) # best
microexonator_e <- markEvent(
  microexonator_e, truths, cols_coord, split = T, split_by = "[\\:\\-]",
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths2", verbose = T
)

saveRDS(microexonator_e, glue("{dir_out}/MicroExonator/microexonator.rds"))




