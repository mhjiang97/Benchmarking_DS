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

dir_out <- "~/projects/AS/analysis/sim_1/"


##### rMATS #####
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
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
rmats_e <- markEvent(
  result = rmats_e, truths = truths, cols_coord = cols_coord, split = F,
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)

saveRDS(rmats_e, glue("{dir_out}/rMATS/rmats.rds"))


##### CASH #####
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
cash_e <- markEvent(
  cash_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(cash_e, glue("{dir_out}/CASH/cash.rds"))


##### LeafCutter #####
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
leafcutter_e <- markEvent(
  leafcutter_e, truths, cols_coord, split = T, split_by = "[\\:\\,]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = F, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(leafcutter_e, glue("{dir_out}/LeafCutter/leafcutter.rds"))


##### SplAdder #####
dir_spladder <- glue("{dir_out}/SplAdder/")

list_spladder <- read_spladder(
  glue("{dir_spladder}/graphs/testing_s1_vs_s2/"),
  glue("{dir_spladder}/graphs/")
)
spladder <- bind_rows(list_spladder) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  mutate(myID = event_id)

## coordinate columns ##
cols_coord <- c(
  "exon_pre_start", "exon_pre_end", "exon_start", "exon_end", "exon_aft_start", "exon_aft_end",
  "exon1_start", "exon1_end", "exon2_start", "exon2_end",
  "exon_const_start", "exon_const_end", "exon_alt1_start", "exon_alt1_end", "exon_alt2_start", "exon_alt2_end",
  "intron_start", "intron_end",
  "exon_starts", "exon_ends"
)
spladder_e <- markEvent(
  spladder, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
spladder_e <- markEvent(
  spladder_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
spladder_e <- markEvent(
  spladder_e, truths, cols_coord, split = F,
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)
saveRDS(spladder_e, glue("{dir_out}/SplAdder/spladder.rds"))


##### MAJIQ #####
dir_majiq <- glue("{dir_out}/MAJIQ/")

majiq <- read_majiq(
  glue("{dir_majiq}/deltapsi/s1_s2.deltapsi.tsv"),
  glue("{dir_majiq}/voila/s1_s2_0.05_showall.tsv"),
  condition_1 = "s1", condition_2 = "s2", flatten = F
) |>
  dplyr::mutate(padj = p) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = lsv_id)
majiq2 <- read_majiq(
  glue("{dir_majiq}/deltapsi/s1_s2.deltapsi.tsv"),
  glue("{dir_majiq}/voila/s1_s2_0.05.tsv"),
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
majiq_e <- markEvent(
  majiq_e, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
majiq_e2 <- markEvent(
  majiq2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 3,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
majiq_e2 <- markEvent(
  majiq_e2, truths, cols_coord, split = T, split_by = "[;-]",
  numbers_equal = 4,
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(majiq_e, glue("{dir_out}/MAJIQ/majiq.rds"))
saveRDS(majiq_e2, glue("{dir_out}/MAJIQ/majiq2.rds"))


##### SGSeq #####
dir_sgseq <- glue("{dir_out}/SGSeq/")

sgseq <- read_sgseq(RData = glue("{dir_sgseq}/SGSeq.RData"), flatten = F) |>
  tidyr::drop_na(padj)
sgseq2 <- read_sgseq(RData = glue("{dir_sgseq}/SGSeq.RData"), flatten = T) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj)

sgseq_cc <- sgseq2 |>
  separate(from, c("from_e", "from_chr", "from_coord", "from_strand")) |>
  separate(to, c("to_e", "to_chr", "to_coord", "to_strand"))

cols_coord <- c("from_coord", "to_coord")


##### SUPPA #####
dir_suppa <- glue("{dir_out}/SUPPA/")

list_suppa <- read_suppa(dir = glue("{dir_suppa}/ds/"))
suppa <- bind_rows(list_suppa) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = coord)

cols_coord <- c("coord")
suppa_e <- markEvent(
  suppa, truths, cols_coord, split = T, split_by = "[:-]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
suppa_e <- markEvent(
  suppa_e, truths, cols_coord, split = T, split_by = "[:-]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
suppa_e <- markEvent(
  suppa_e, truths, cols_coord, split = T, split_by = "[:-]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)
saveRDS(suppa_e, glue("{dir_out}/SUPPA/suppa.rds"))


##### Whippet #####
dir_whippet <- glue("{dir_out}/Whippet/")

whippet <- read_whippet(glue("{dir_whippet}/delta/s1_s2.diff2.gz")) |>
  dplyr::left_join(tx2gene_convert, by = "gene_id") |>
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
whippet_e <- markEvent(
  whippet_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(whippet_e, glue("{dir_out}/Whippet/whippet.rds"))


##### ASpli #####
dir_aspli <- glue("{dir_out}/ASpli/")

aspli <- read_aspli(
  file_exon_du = glue("{dir_aspli}/s1-s2/exon.du.tab"),
  file_intron_du = glue("{dir_aspli}/s1-s2/intron.du.tab")
) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = ID)

saveRDS(aspli, glue("{dir_out}/ASpli/aspli.rds"))


##### BANDITS #####
dir_bandits <- glue("{dir_out}/BANDITS/")

bandits <- read_bandits(
  transcript_result_file = glue("{dir_bandits}/results_transcript.txt"),
  gene_result_file = glue("{dir_bandits}/results_gene.txt")
) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(gene_id, transcript_id))
bandits$as_type <- "OTHER"

saveRDS(bandits, glue("{dir_out}/BANDITS/bandits.rds"))


##### NBSplice #####
dir_nbsplice <- glue("{dir_out}/NBSplice/")

nbsplice <- read_nbsplice(result_file = glue("{dir_nbsplice}/results.txt")) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(gene_id, transcript_id))
nbsplice$as_type <- "OTHER"

saveRDS(nbsplice, glue("{dir_out}/NBSplice/nbsplice.rds"))


##### IsoformSwitchAnalyzeR #####
dir_isoformswitchanalyzer <- glue("{dir_out}/IsoformSwitchAnalyzeR/")

isoformswitchanalyzer <- read_isoformswitchanalyzer(
  feature_file = glue("{dir_isoformswitchanalyzer}/isoformFeatures.txt"),
  splicing_file = glue("{dir_isoformswitchanalyzer}/AlternativeSplicingAnalysis.txt")
) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = iso_ref)

isoformswitchanalyzer_tmp <- isoformswitchanalyzer[, c("SE", "MXE", "SME", "RI", "A5SS", "A3SS", "AFE")]

saveRDS(isoformswitchanalyzer, glue("{dir_out}/IsoformSwitchAnalyzeR/isoformswitchanalyzer.rds"))


##### DiffSplice #####
dir_diffsplice <- glue("{dir_out}/DiffSplice/")

diffsplice <- read_diffsplice(
  transcript_file = glue("{dir_diffsplice}/result/differential_transcription.txt")
)


##### psichomics #####
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
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
psichomics_e <- markEvent(
  psichomics_e, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)
saveRDS(psichomics_e, glue("{dir_out}/psichomics/psichomics.rds"))


##### DIEGO #####
dir_diego <- glue("{dir_out}/DIEGO/")

diego <- read_diego(result_file = glue("{dir_diego}/results.txt")) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = junction)
diego$as_type <- "OTHER"

diego_cc <- diego |>
  separate(junction, c("chr", "coord1", "coord2"), "[:-]")

cols_coord <- c("coord1", "coord2")
diego_e <- markEvent(
  diego_cc, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
diego_e <- markEvent(
  diego_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = F, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(diego_e, glue("{dir_out}/DIEGO/diego.rds"))


##### VAST-TOOLS #####
dir_vasttools <- glue("{dir_out}/VAST-TOOLS/")

vasttools <- read_vasttools(
  diff_file = "/Volumes/MyBook/projects/AS/analysis/sim_1/VAST-TOOLS/combine/diff.tab",
  event_annot_file = "/Volumes/MyBook/projects/AS/analysis/sim_1/VAST-TOOLS/EVENT_INFO-hg38.tab.gz"
) |>
  left_join(tx2gene_convert, by = "gene_symbol") |>
  dplyr::mutate(padj = as.character(0.1 - as.numeric(`MV[dPsi]_at_0.95`))) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = EVENT)

cols_coord <- c("COORD_o", "FULL_CO", "REF_CO", "CO_C1", "CO_A", "CO_C2")
vasttools_e <- markEvent(
  vasttools, truths, cols_coord, split = T, split_by = "[\\,\\-\\+\\:]",
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
vasttools_e <- markEvent(
  vasttools_e, truths, cols_coord, split = T, split_by = "[\\,\\-\\+\\:]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
vasttools_e <- markEvent(
  vasttools_e, truths, cols_coord, split = T, split_by = "[\\,\\-\\+\\:]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)
saveRDS(vasttools_e, glue("{dir_out}/VAST-TOOLS/vasttools.rds"))


##### RATS #####
dir_rats <- glue("{dir_out}/RATS/")

rats <- read_rats(transcript_result_file = glue("{dir_rats}/transcripts.txt")) |>
  left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(transcript_id, gene_id))
rats$as_type <- "OTHER"

saveRDS(rats, glue("{dir_out}/RATS/rats.rds"))


##### dSpliceType #####
dir_dsplicetype <- glue("{dir_out}/dSpliceType/")

list_dsplicetype <- read_dsplicetype(dir = dir_dsplicetype, event_type = "all")
dsplicetype <- dplyr::bind_rows(list_dsplicetype) |>
  dplyr::left_join(tx2gene_convert, by = "gene_id") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = paste0(as_type, ": ", ID))

cols_coord <- c("event")
# dsplicetype_e <- markEvent(
#   dsplicetype, truths, cols_coord, split = T, split_by = "\\_",
#   numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
#   error_distance = 1, as = T, gene = T, gene_symbol = T,
#   truths_col = c("X16", "X17"), mark_col = "idx_truths"
# )
dsplicetype_e <- markEvent(
  dsplicetype, truths, cols_coord, split = T, split_by = "\\_",
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(dsplicetype_e, glue("{dir_out}/dSpliceType/dsplicetype.rds"))


##### PSI-Sigma #####
dir_psisigma <- glue("{dir_out}/PSI-Sigma/")

psisigma <- read_psisigma(result_file = glue("{dir_psisigma}/s1_s2_r5_ir3.sorted.txt")) |>
  left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(myID = `Database ID`)
psisigma_cc <- psisigma |>
  separate(`Target Exon`, c("chr", "coord1", "coord2"), "[:-]")

cols_coord <- c("coord1", "coord2")
psisigma_e <- markEvent(
  psisigma_cc, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17"), mark_col = "idx_truths"
)
psisigma_e <- markEvent(
  psisigma_e, truths, cols_coord, split = F, split_by = NULL,
  numbers_equal = c("SE" = 2, "RI" = 2, "A3SS" = 2, "A5SS" = 2, "MXE" = 2),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
saveRDS(psisigma_e, glue("{dir_out}/PSI-Sigma/psisigma.rds"))


##### JuncBASE #####
dir_juncbase <- glue("{dir_out}/JuncBASE/")

juncbase <- read_juncbase(result_file = glue("{dir_juncbase}/sample_set_comparison.txt")) |>
  tidyr::drop_na(padj) |>
  dplyr::mutate(
    myID = paste0(
      as_event_type, exclusion_junctions, inclusion_exons, exclusion_exons,
      inclusion_exons, `intron-exon_junctions`
    )
  )
# juncbase$as_type[juncbase$as_type == "coordinate_cassette"] <- "OTHER"
juncbase2 <- read_juncbase(result_file = glue("{dir_juncbase}/sample_set_comparison.txt"), T) |>
  left_join(tx2gene_convert, by = "gene_id") |>
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
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 5, "RI" = 3, "A3SS" = 5, "A5SS" = 5, "MXE" = 7),
  error_distance = 1, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_2"
)
juncbase_e <- markEvent(
  juncbase_e, truths, cols_coord, split = T, split_by = "[\\:\\-\\;]",
  numbers_equal = c("SE" = 4, "RI" = 3, "A3SS" = 4, "A5SS" = 4, "MXE" = 6),
  error_distance = 1, as = T, gene = F, gene_symbol = F,
  truths_col = c("X16", "X17", "coord_append"), mark_col = "idx_truths_3"
)
saveRDS(juncbase_e, glue("{dir_out}/JuncBASE/juncbase.rds"))


##### SpliceSeq #####
dir_spliceseq <- glue("{dir_out}/SpliceSeq/")

spliceseq <- read_spliceseq(glue("{dir_spliceseq}/spliceseq.tsv")) |>
  dplyr::left_join(tx2gene_convert, by = "gene_symbol") |>
  tidyr::drop_na(padj) %>%
  dplyr::mutate(myID = paste0(gene_symbol, "_", 1:nrow(.)))

saveRDS(spliceseq, glue("{dir_out}/SpliceSeq/spliceseq.rds"))


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




