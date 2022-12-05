pkgs <- c(
  "gtools", "dplyr", "biomaRt", "glue", "vroom", "ggplot2", "scales",
  "GenomicRanges", "GenomicFeatures", "Hmisc", "rtracklayer", "this.path"
)
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = T))
setwd(dirname(this.path()))

rdirichlet2 <- function(n, alpha) {
  a <- rdirichlet(n, alpha)
  a[is.na(a)] <- 0
  a
}
mutate_IsoPct <- function(w, ref) {
  o <- order(ref, decreasing = T)[1:2]
  w[rev(o)] <- w[o]
  w
}
mutated_IsoPct <- function(ref) {
  o <- order(ref, decreasing = T)[1:2]
  w <- rep(0, length(ref))
  w[o] <- 1
  w
}

supAsEvents <- function(ioe_file, isoform_table, ds_genes) {
  library(vroom)
  library(stringr)
  library(dplyr)
  # ioe_file <- "~/doc/suppa/events_A3_strict.ioe.gz"
  ioe <- vroom(ioe_file, col_types = cols()) %>%
    as.data.frame()
  ioe$idct <- 0
  # should get ds_genes clear first
  for (i in 1:nrow(ioe)) {
    if (!ioe$gene_id[i] %in% ds_genes) next
    
    isos_1 <- str_split(ioe$alternative_transcripts[i], ",") %>%
      unlist()
    isos_all <- str_split(ioe$total_transcripts[i], ",") %>%
      unlist()
    isos_2 <- setdiff(isos_all, isos_1)
    
    iso_switched1 <- dplyr::filter(
      isoform_table, gene_id == ioe$gene_id[i], major == T
    )$transcript_id
    iso_switched2 <- dplyr::filter(
      isoform_table, gene_id == ioe$gene_id[i], major2 == T
    )$transcript_id
    
    if (iso_switched1 %in% isos_1 && iso_switched2 %in% isos_2) ioe$idct[i] <- 1
    if (iso_switched1 %in% isos_2 && iso_switched2 %in% isos_1) ioe$idct[i] <- 1
  }
  dplyr::filter(ioe, idct == 1)
}


isoform_results_file <- glue(
  "../../files/SRR493366.isoforms.mod.results"
)
meandisp.file <- glue(
  "../../files/Pickrell.Cheung.Mu.Phi.Estimates.rds"
)
meandisp <- readRDS(meandisp.file)
meanvect <- meandisp$pickrell.cheung.mu
dispvect <- meandisp$pickrell.cheung.phi
librarysize <- 5000000
nbr_diff_expr <- 0
nbr_per_group <- 20
nbr_gene_expr <- 5000
nrep <- 1
out_dir <- "../../profiles/single_cell/20/"
ensembl <- useMart(
  "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://sep2019.archive.ensembl.org"
)
select_as_events <- F

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
if (!file.exists(glue("{out_dir}/all_affected.txt"))) {
  file.copy("../../profiles/all_affected.txt", out_dir)
}

## Define the isoform results file to start from (from RSEM).
## We will modify the TPM column of this file to create the files that
## are used in the simulation (the simulator uses only the TPM column, so
## we will not change the others)

isoform.initial <- vroom(
  isoform_results_file, delim = "\t", col_types = cols()
) %>%
  as.data.frame()

gene_summary_tmp <- isoform.initial %>%
  group_by(gene_id) %>%
  summarise(
    expected_gene_count_gr1 = sum(expected_count),
    effective_gene_length = sum(effective_length * IsoPct / 100),
    nbr_isoforms = length(IsoPct),
    nbr_expr_isoforms = length(which(IsoPct > 0)),
    nbr_expr_isoforms10 = length(which(IsoPct > 10))
  ) %>%
  ungroup() %>%
  as.data.frame()

all_aff <- read.table(glue("{out_dir}/all_affected.txt"), sep = "\t")
set.seed(1)
gene_filter <- sample(
  setdiff(
    intersect(
      unique(isoform.initial$gene_id),
      unique(gene_summary_tmp$gene_id[gene_summary_tmp$expected_gene_count_gr1 > 0])
    ),
    unique(all_aff[, 9])
  ),
  length(unique(isoform.initial$gene_id[isoform.initial$TPM > 0])) - nbr_gene_expr
)
isoform.initial$expected_count[which(isoform.initial$gene_id %in% gene_filter)] <- 0
# isoform.initial$effective_length[which(isoform.initial$gene_id %in% gene_filter)] <- 0
isoform.initial$TPM[which(isoform.initial$gene_id %in% gene_filter)] <- 0
isoform.initial$FPKM[which(isoform.initial$gene_id %in% gene_filter)] <- 0
isoform.initial$IsoPct[which(isoform.initial$gene_id %in% gene_filter)] <- 0
rm("gene_summary_tmp")

gene_summary <- isoform.initial %>%
  group_by(gene_id) %>%
  summarise(
    expected_gene_count_gr1 = sum(expected_count),
    effective_gene_length = sum(effective_length * IsoPct / 100),
    nbr_isoforms = length(IsoPct),
    nbr_expr_isoforms = length(which(IsoPct > 0)),
    nbr_expr_isoforms10 = length(which(IsoPct > 10))
  ) %>%
  ungroup() %>%
  as.data.frame()

gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr1
gene_summary$gene_de_status <- 0

set.seed(1)
fold_changes <- (2 + rexp(nbr_diff_expr, rate = 1))^
  (c(-1, 1)[round(runif(nbr_diff_expr)) + 1])
set.seed(1)
diff_expr_genes <- sample(1:nrow(gene_summary), nbr_diff_expr, replace = F)
gene_summary$expected_gene_count_gr2[diff_expr_genes] <-
  gene_summary$expected_gene_count_gr1[diff_expr_genes] * fold_changes
gene_summary$gene_de_status[diff_expr_genes] <- 1

gene_summary$fold_change <- 1
gene_summary$fold_change[diff_expr_genes] <- fold_changes

## Adjust the expected gene count to the desired library size, to obtain
## the right dispersion estimates
gene_summary$expected_gene_count_gr1 <- gene_summary$expected_gene_count_gr1 /
  sum(gene_summary$expected_gene_count_gr1) * librarysize
gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr2 /
  sum(gene_summary$expected_gene_count_gr2) * librarysize

##### assign dispersion to each gene in each conditions #####
gene_summary$dispersion_gr1 <- sapply(
  gene_summary$expected_gene_count_gr1,
  function(w) dispvect[which.min(abs(meanvect - w))]
)
gene_summary$dispersion_gr2 <- sapply(
  gene_summary$expected_gene_count_gr2,
  function(w) dispvect[which.min(abs(meanvect - w))]
)

##### find out which chromosome do genes belong to #####
bm <- getBM(
  attributes = c("ensembl_gene_id_version", "chromosome_name"),
  filters = "ensembl_gene_id_version",
  values = gene_summary$gene_id, mart = ensembl
)
gene_summary <- merge(
  gene_summary, bm, by.x = "gene_id", by.y = "ensembl_gene_id_version", all = T
)
## PAR doesn't appear in package biomaRt ##
gene_summary$chromosome_name[is.na(gene_summary$chromosome_name)] <- "PAR"

##### simulate gene counts for the desired number of samples and calculate RPK for genes #####
for (i in 1:nbr_per_group) {
  gene_summary[, paste0("s", i, "_geneCount")] <- sapply(
    1:nrow(gene_summary), function(j) {
      rnbinom(
        n = 1, size = 1 / gene_summary$dispersion_gr1[j],
        mu = gene_summary$expected_gene_count_gr1[j]
      )
    })
  gene_summary[, paste0("s", i, "_geneRPK")] <-
    gene_summary[, paste0("s", i, "_geneCount")] /
    gene_summary$effective_gene_length * 1e3
}

for (i in (nbr_per_group + 1):(2 * nbr_per_group)) {
  gene_summary[, paste0("s", i, "_geneCount")] <- sapply(
    1:nrow(gene_summary), function(j) {
      rnbinom(
        n = 1, size = 1 / gene_summary$dispersion_gr2[j],
        mu = gene_summary$expected_gene_count_gr2[j]
      )
    })
  gene_summary[, paste0("s", i, "_geneRPK")] <-
    gene_summary[, paste0("s", i, "_geneCount")] /
    gene_summary$effective_gene_length * 1e3
}

isoform_summary <- isoform.initial
##### add columns indicating the most expressed and the second most expressed isoform for each genes #####
id_major <- isoform_summary %>%
  group_by(gene_id) %>%
  mutate(id_major = transcript_id[which.max(IsoPct)]) %>%
  ungroup() %>%
  as.data.frame() %>%
  .$id_major %>%
  unique()
isoform_summary$major <- isoform_summary$transcript_id %in% id_major

id_major2 <- isoform_summary[-which(isoform_summary$major == T), ] %>%
  group_by(gene_id) %>%
  mutate(id_major2 = transcript_id[which.max(IsoPct)]) %>%
  ungroup() %>%
  as.data.frame() %>%
  .$id_major2 %>%
  unique()
isoform_summary$major2 <- isoform_summary$transcript_id %in% id_major2

##### run only once to define genes harboring alternative splicing events #####
if (select_as_events) {
  index_gr1 <- which(
    gene_summary$expected_gene_count_gr1 > 500 &
      gene_summary$nbr_expr_isoforms10 >= 2
  )
  asta_common_file <- "../../files/asta_table_common_coor.txt.gz" # asta_common <- read.table(asta_common_file, sep = "\t", header = F)
  candi <- read.table(asta_common_file, sep = "\t", header = F) %>%
    dplyr::filter(
      V9 %in% gene_summary$gene_id[index_gr1]
    ) # candi <- dplyr::filter(asta_common, V9 %in% gene_summary$gene_id[index_gr1]) # candi <- asta_common[asta_common[, 9] %in% gene_summary$gene_id[index_gr1], ]
  ##### mark out common events which will be alternative splicing between two conditions if switched the relative expression percentage of two main isoforms #####
  candi$idct <- 0
  for (i in 1:nrow(candi)) {
    tmp_gene <- candi[i, 9]
    isos <- strsplit(candi[i, 10], ",") %>% unlist()
    isos_1 <- strsplit(isos[1], "/") %>% unlist()
    isos_2 <- strsplit(isos[2], "/") %>% unlist()
    
    iso_max <- isoform_summary$transcript_id[
      isoform_summary$gene_id == tmp_gene & isoform_summary$major == T
    ]
    iso_max2 <- isoform_summary$transcript_id[
      isoform_summary$gene_id == tmp_gene & isoform_summary$major2 == T
    ]
    
    if (iso_max %in% isos_1 && iso_max2 %in% isos_2) candi$idct[i] <- 1
    if (iso_max %in% isos_2 && iso_max2 %in% isos_1) candi$idct[i] <- 1
  }
  
  EE <- dplyr::filter(candi, V13 == "1-2^,3-4^", idct == 1)
  IR <- dplyr::filter(candi, V13 == "0,1^2-", idct == 1)
  AD <- dplyr::filter(candi, V13 == "1^,2^", idct == 1)
  AA <- dplyr::filter(candi, V13 == "1-,2-", idct == 1)
  ES <- dplyr::filter(candi, V13 == "0,1-2^", idct == 1)
  
  genes_remain <- candi[candi$idct == 1, 9] %>% unique()
  length(unique(dplyr::filter(EE, V9 %in% genes_remain)[, 9]))
  genes_ee <- dplyr::filter(EE, V9 %in% genes_remain)[, 9] %>%
    unique() %>%
    sample(44)
  ee <- dplyr::filter(EE, V9 %in% genes_ee)
  genes_remain <- setdiff(genes_remain, genes_ee)
  
  length(unique(dplyr::filter(IR, V9 %in% genes_remain)[, 9]))
  genes_ir <- dplyr::filter(IR, V9 %in% genes_remain)[, 9] %>%
    unique() %>%
    sample(110)
  ir <- dplyr::filter(IR, V9 %in% genes_ir)
  genes_remain <- setdiff(genes_remain, genes_ir)
  
  length(unique(dplyr::filter(AD, V9 %in% genes_remain)[, 9]))
  genes_ad <- dplyr::filter(AD, V9 %in% genes_remain)[, 9] %>%
    unique() %>%
    sample(282)
  ad <- dplyr::filter(AD, V9 %in% genes_ad)
  genes_remain <- setdiff(genes_remain, genes_ad)
  
  length(unique(dplyr::filter(AA, V9 %in% genes_remain)[, 9]))
  genes_aa <- dplyr::filter(AA, V9 %in% genes_remain)[, 9] %>%
    unique() %>%
    sample(282)
  aa <- dplyr::filter(AA, V9 %in% genes_aa)
  genes_remain <- setdiff(genes_remain, genes_aa)
  
  length(unique(dplyr::filter(ES, V9 %in% genes_remain)[, 9]))
  genes_es <- dplyr::filter(ES, V9 %in% genes_remain)[, 9] %>%
    unique() %>%
    sample(282)
  es <- dplyr::filter(ES, V9 %in% genes_es)
  
  genes_sel <- c(genes_ee, genes_ir, genes_ad, genes_aa, genes_es)
  all_aff <- dplyr::filter(candi, idct == 1, V9 %in% genes_sel)
  
  write.table(
    all_aff, glue("{out_dir}/all_affected.txt"),
    sep = "\t", quote = F, col.names = F, row.names = F
  )
} else {
  all_aff <- read.table(glue("{out_dir}/all_affected.txt"), sep = "\t")
}

for (nsim in 1:nrep) {
  mydir <- glue("{out_dir}/sim_{nsim}")
  if (!dir.exists(mydir)) dir.create(mydir)
  setwd(mydir)
  
  isoform_summary <- isoform_summary %>%
    dplyr::select(
      transcript_id, gene_id, length, effective_length, expected_count, TPM, FPKM,
      IsoPct, major, major2
    ) %>%
    as.data.frame()
  
  ##### mutate isoform expression percentage slightly for each samples #####
  for (i in 1:(2 * nbr_per_group)) {
    isoform_summary <- isoform_summary %>%
      group_by(gene_id) %>%
      mutate(IsoPctDirichlet = c(rdirichlet2(1, IsoPct / 100 * 100))) %>%
      setNames(c(colnames(isoform_summary), paste0("s", i, "_IsoPct"))) %>%
      ungroup() %>%
      as.data.frame()
  }
  
  ##### switch the top 2 expressed isoforms of selected differential alternative splicing genes #####
  ds_genes <- unique(all_aff[, 9])
  isoform_summary_nonds <- isoform_summary[!(isoform_summary$gene_id %in% ds_genes), ]
  isoform_summary_ds <- isoform_summary[isoform_summary$gene_id %in% ds_genes, ]
  
  isoform_summary_ds <- isoform_summary_ds %>%
    group_by(gene_id) %>%
    mutate(
      s21_IsoPct = mutate_IsoPct(s21_IsoPct, IsoPct),
      s22_IsoPct = mutate_IsoPct(s22_IsoPct, IsoPct),
      s23_IsoPct = mutate_IsoPct(s23_IsoPct, IsoPct),
      s24_IsoPct = mutate_IsoPct(s24_IsoPct, IsoPct),
      s25_IsoPct = mutate_IsoPct(s25_IsoPct, IsoPct),
      s26_IsoPct = mutate_IsoPct(s26_IsoPct, IsoPct),
      s27_IsoPct = mutate_IsoPct(s27_IsoPct, IsoPct),
      s28_IsoPct = mutate_IsoPct(s28_IsoPct, IsoPct),
      s29_IsoPct = mutate_IsoPct(s29_IsoPct, IsoPct),
      s30_IsoPct = mutate_IsoPct(s30_IsoPct, IsoPct),
      s31_IsoPct = mutate_IsoPct(s31_IsoPct, IsoPct),
      s32_IsoPct = mutate_IsoPct(s32_IsoPct, IsoPct),
      s33_IsoPct = mutate_IsoPct(s33_IsoPct, IsoPct),
      s34_IsoPct = mutate_IsoPct(s34_IsoPct, IsoPct),
      s35_IsoPct = mutate_IsoPct(s35_IsoPct, IsoPct),
      s36_IsoPct = mutate_IsoPct(s36_IsoPct, IsoPct),
      s37_IsoPct = mutate_IsoPct(s37_IsoPct, IsoPct),
      s38_IsoPct = mutate_IsoPct(s38_IsoPct, IsoPct),
      s39_IsoPct = mutate_IsoPct(s39_IsoPct, IsoPct),
      s40_IsoPct = mutate_IsoPct(s40_IsoPct, IsoPct),
      diff_IsoPct = -diff(sort(IsoPct, decreasing = T))[1] / 100,
      gene_ds_status = 1,
      transcript_ds_status = mutated_IsoPct(IsoPct)
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  isoform_summary_nonds <- isoform_summary_nonds %>%
    group_by(gene_id) %>%
    mutate(
      diff_IsoPct = -diff(sort(IsoPct, decreasing = T))[1] / 100,
      gene_ds_status = 0,
      transcript_ds_status = 0
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  isoform_summary <- rbind(isoform_summary_nonds, isoform_summary_ds)
  
  ##### merge the gene data frame and the isoform data frame #####
  final_summary <- merge(
    isoform_summary, gene_summary,
    by.x = "gene_id", by.y = "gene_id", all = T
  )
  
  ##### calculate RPK for each isoforms by dividing gene RPK according to IsoPct columns #####
  for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformRPK")] <-
      final_summary[, paste0("s", i, "_geneRPK")] *
      final_summary[, paste0("s", i, "_IsoPct")]
  }
  
  ##### calculate count for isoforms #####
  for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformCount")] <-
      final_summary[, paste0("s", i, "_isoformRPK")] *
      final_summary$effective_length / 1000
  }
  
  ##### calculate FPKM for isoforms #####
  for (i in 1:(2 * nbr_per_group)) {
    fpkm <- (final_summary[, paste0("s", i, "_isoformCount")] * 10^9) /
      (
        final_summary$effective_length *
          sum(final_summary[, paste0("s", i, "_isoformCount")], na.rm = T)
      )
    fpkm[is.na(fpkm) | fpkm == Inf] <- 0
    
    final_summary[, paste0("s", i, "_isoformFPKM")] <- fpkm
  }
  
  ##### calculate TPM for isoforms by zooming the original TPM according to the ratios between final simulated counts and original counts #####
  for (i in 1:(2 * nbr_per_group)) {
    final_summary[, paste0("s", i, "_isoformTPM")] <- round(
      final_summary[, paste0("s", i, "_isoformCount")] /
        final_summary$expected_count * final_summary$TPM, 2
    )
  }
  
  final_summary[is.na(final_summary)] <- 0
  
  ##### scale each isoformTPM column so that it sums to 1 million #####
  idx <- grep("isoformTPM", colnames(final_summary))
  for (i in idx) {
    final_summary[, i] <- final_summary[, i] / sum(final_summary[, i]) * 1e6
  }
  
  ##### output the final summary data frame #####
  write.table(
    final_summary, "simulation_details.txt",
    row.names = F, col.names = T, quote = F, sep = "\t"
  )
  
  ##### output files for each simulated samples #####
  for (i in 1:(2 * nbr_per_group)) {
    tmp <- isoform.initial
    
    tmp$TPM <- final_summary[
      match(tmp$transcript_id, final_summary$transcript_id),
      paste0("s", i, "_isoformTPM")
    ]
    tmp$TPM <- round(tmp$TPM, digits = 2)
    tmp$TPM <- as.character(tmp$TPM)
    tmp$TPM[tmp$TPM == "0"] <- "0.00"
    
    tmp$expected_count <- as.character(tmp$expected_count)
    tmp$expected_count[tmp$expected_count == "0"] <- "0.00"
    
    tmp$FPKM <- final_summary[
      match(tmp$transcript_id, final_summary$transcript_id),
      paste0("s", i, "_isoformFPKM")
    ]
    tmp$FPKM <- round(tmp$FPKM, digits = 2)
    tmp$FPKM <- as.character(tmp$FPKM)
    tmp$FPKM[tmp$FPKM == "0"] <- "0.00"
    
    tmp$IsoPct <- as.character(tmp$IsoPct)
    tmp$IsoPct[tmp$IsoPct == "0"] <- "0.00"
    
    tmp$effective_length <- as.character(tmp$effective_length)
    tmp$effective_length[tmp$effective_length == "0"] <- "0.00"
    
    write.table(
      tmp, file = glue("sample{i}.txt"),
      row.names = F, col.names = T, quote = F, sep = "\t"
    )
  }
  
  ##### save the whole working environment #####
  save.image(file = "sim.RData", compress = T)
}


##### plot the simulation scenario #####
for (nsim in 1:nrep) {
  ##### read the simulation details file generated by codes above #####
  details <- read.delim(
    glue("{out_dir}/sim_{nsim}/simulation_details.txt"),
    header = T, as.is = T
  )
  
  details$gr1_isoformTPM <- 1 / 20 * (
    details$s1_isoformTPM + details$s2_isoformTPM + details$s3_isoformTPM + details$s4_isoformTPM + details$s5_isoformTPM +
      details$s6_isoformTPM + details$s7_isoformTPM + details$s8_isoformTPM + details$s9_isoformTPM + details$s10_isoformTPM +
      details$s11_isoformTPM + details$s12_isoformTPM + details$s13_isoformTPM + details$s14_isoformTPM + details$s15_isoformTPM +
      details$s16_isoformTPM + details$s17_isoformTPM + details$s18_isoformTPM + details$s19_isoformTPM + details$s20_isoformTPM
  )
  details$gr2_isoformTPM <- 1 / 20 * (
    details$s21_isoformTPM + details$s22_isoformTPM + details$s23_isoformTPM + details$s24_isoformTPM + details$s25_isoformTPM +
      details$s26_isoformTPM + details$s27_isoformTPM + details$s28_isoformTPM + details$s29_isoformTPM + details$s30_isoformTPM +
      details$s31_isoformTPM + details$s32_isoformTPM + details$s33_isoformTPM + details$s34_isoformTPM + details$s35_isoformTPM +
      details$s36_isoformTPM + details$s37_isoformTPM + details$s38_isoformTPM + details$s39_isoformTPM + details$s40_isoformTPM
  )
  details$transcript_ds_cat <- ifelse(
    details$transcript_ds_status == 0, "non-differentially used", "differentially used"
  )
  x <- details %>%
    group_by(gene_id) %>%
    summarise(
      s1_geneTPM = sum(s1_isoformTPM),
      s2_geneTPM = sum(s2_isoformTPM),
      s3_geneTPM = sum(s3_isoformTPM),
      s4_geneTPM = sum(s4_isoformTPM),
      s5_geneTPM = sum(s5_isoformTPM),
      s6_geneTPM = sum(s6_isoformTPM),
      s7_geneTPM = sum(s7_isoformTPM),
      s8_geneTPM = sum(s8_isoformTPM),
      s9_geneTPM = sum(s9_isoformTPM),
      s10_geneTPM = sum(s10_isoformTPM),
      s11_geneTPM = sum(s11_isoformTPM),
      s12_geneTPM = sum(s12_isoformTPM),
      s13_geneTPM = sum(s13_isoformTPM),
      s14_geneTPM = sum(s14_isoformTPM),
      s15_geneTPM = sum(s15_isoformTPM),
      s16_geneTPM = sum(s16_isoformTPM),
      s17_geneTPM = sum(s17_isoformTPM),
      s18_geneTPM = sum(s18_isoformTPM),
      s19_geneTPM = sum(s19_isoformTPM),
      s20_geneTPM = sum(s20_isoformTPM),
      s21_geneTPM = sum(s21_isoformTPM),
      s22_geneTPM = sum(s22_isoformTPM),
      s23_geneTPM = sum(s23_isoformTPM),
      s24_geneTPM = sum(s24_isoformTPM),
      s25_geneTPM = sum(s25_isoformTPM),
      s26_geneTPM = sum(s26_isoformTPM),
      s27_geneTPM = sum(s27_isoformTPM),
      s28_geneTPM = sum(s28_isoformTPM),
      s29_geneTPM = sum(s29_isoformTPM),
      s30_geneTPM = sum(s30_isoformTPM),
      s31_geneTPM = sum(s31_isoformTPM),
      s32_geneTPM = sum(s32_isoformTPM),
      s33_geneTPM = sum(s33_isoformTPM),
      s34_geneTPM = sum(s34_isoformTPM),
      s35_geneTPM = sum(s35_isoformTPM),
      s36_geneTPM = sum(s36_isoformTPM),
      s37_geneTPM = sum(s37_isoformTPM),
      s38_geneTPM = sum(s38_isoformTPM),
      s39_geneTPM = sum(s39_isoformTPM),
      s40_geneTPM = sum(s40_isoformTPM),
      gr1_geneTPM = 1 / 20 * (
        s1_geneTPM + s2_geneTPM + s3_geneTPM + s4_geneTPM + s5_geneTPM +
          s6_geneTPM + s7_geneTPM + s8_geneTPM + s9_geneTPM + s10_geneTPM +
          s11_geneTPM + s12_geneTPM + s13_geneTPM + s14_geneTPM + s15_geneTPM +
          s16_geneTPM + s17_geneTPM + s18_geneTPM + s19_geneTPM + s20_geneTPM
      ),
      gr2_geneTPM = 1 / 20 * (
        s21_geneTPM + s22_geneTPM + s23_geneTPM + s24_geneTPM + s25_geneTPM +
          s26_geneTPM + s27_geneTPM + s28_geneTPM + s29_geneTPM + s30_geneTPM +
          s31_geneTPM + s32_geneTPM + s33_geneTPM + s34_geneTPM + s35_geneTPM +
          s36_geneTPM + s37_geneTPM + s38_geneTPM + s39_geneTPM + s40_geneTPM
      ),
      s1_geneCount = sum(s1_isoformCount),
      s2_geneCount = sum(s2_isoformCount),
      s3_geneCount = sum(s3_isoformCount),
      s4_geneCount = sum(s4_isoformCount),
      s5_geneCount = sum(s5_isoformCount),
      s6_geneCount = sum(s6_isoformCount),
      s7_geneCount = sum(s7_isoformCount),
      s8_geneCount = sum(s8_isoformCount),
      s9_geneCount = sum(s9_isoformCount),
      s10_geneCount = sum(s10_isoformCount),
      s11_geneCount = sum(s11_isoformCount),
      s12_geneCount = sum(s12_isoformCount),
      s13_geneCount = sum(s13_isoformCount),
      s14_geneCount = sum(s14_isoformCount),
      s15_geneCount = sum(s15_isoformCount),
      s16_geneCount = sum(s16_isoformCount),
      s17_geneCount = sum(s17_isoformCount),
      s18_geneCount = sum(s18_isoformCount),
      s19_geneCount = sum(s19_isoformCount),
      s20_geneCount = sum(s20_isoformCount),
      s21_geneCount = sum(s21_isoformCount),
      s22_geneCount = sum(s22_isoformCount),
      s23_geneCount = sum(s23_isoformCount),
      s24_geneCount = sum(s24_isoformCount),
      s25_geneCount = sum(s25_isoformCount),
      s26_geneCount = sum(s26_isoformCount),
      s27_geneCount = sum(s27_isoformCount),
      s28_geneCount = sum(s28_isoformCount),
      s29_geneCount = sum(s29_isoformCount),
      s30_geneCount = sum(s30_isoformCount),
      s31_geneCount = sum(s31_isoformCount),
      s32_geneCount = sum(s32_isoformCount),
      s33_geneCount = sum(s33_isoformCount),
      s34_geneCount = sum(s34_isoformCount),
      s35_geneCount = sum(s35_isoformCount),
      s36_geneCount = sum(s36_isoformCount),
      s37_geneCount = sum(s37_isoformCount),
      s38_geneCount = sum(s38_isoformCount),
      s39_geneCount = sum(s39_isoformCount),
      s40_geneCount = sum(s40_isoformCount),
      gr1_geneCount = 1 / 20 * (
        s1_geneCount + s2_geneCount + s3_geneCount + s4_geneCount + s5_geneCount +
          s6_geneCount + s7_geneCount + s8_geneCount + s9_geneCount + s10_geneCount +
          s11_geneCount + s12_geneCount + s13_geneCount + s14_geneCount + s15_geneCount +
          s16_geneCount + s17_geneCount + s18_geneCount + s19_geneCount + s20_geneCount
      ),
      gr2_geneCount = 1 / 20 * (
        s21_geneCount + s22_geneCount + s23_geneCount + s24_geneCount + s25_geneCount +
          s26_geneCount + s27_geneCount + s28_geneCount + s29_geneCount + s30_geneCount +
          s31_geneCount + s32_geneCount + s33_geneCount + s34_geneCount + s35_geneCount +
          s36_geneCount + s37_geneCount + s38_geneCount + s39_geneCount + s40_geneCount
      ),
      status_ds = gene_ds_status[1],
      status_de = gene_de_status[1]
    )
  
  x$status_ds <- ifelse(
    x$status_ds == 0, "non-differentially spliced", "differentially spliced"
  )
  x$status_de <- ifelse(
    x$status_de == 0, "non-differetially expressed", "differentially expressed"
  )
  pdf(glue("{out_dir}/sim_{nsim}/simulation_details.pdf"))
  print(
    ggplot(x, aes(x = gr1_geneCount, y = gr2_geneCount, col = status_ds)) +
      geom_abline(intercept = 0, slope = 1) +
      geom_point(alpha = 1 / 2, aes(shape = status_de)) +
      scale_shape_manual(values = c(18, 16)) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(values = c("red", "black"), name = "") +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      xlab("Average gene count, condition 1") +
      ylab("Average gene count, condition 2") +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = NA, colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(
          angle = 90, vjust = 0.5, hjust = 1, size = 15
        ),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = NA, colour = "black")
      ) +
      ggtitle("")
  )
  
  print(
    ggplot(x, aes(x = gr1_geneTPM, y = gr2_geneTPM, col = status_ds)) +
      geom_abline(intercept = 0, slope = 1) +
      geom_point(alpha = 1 / 2) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(values = c("red", "black"), name = "") +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      xlab("Average gene TPM, condition 1") +
      ylab("Average gene TPM, condition 2") +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = NA, colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(
          angle = 90, vjust = 0.5, hjust = 1, size = 15
        ),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = NA, colour = "black")
      ) +
      ggtitle("")
  )
  
  print(
    ggplot(
      details, aes(x = gr1_isoformTPM, y = gr2_isoformTPM, col = transcript_ds_cat)
    ) +
      geom_abline(intercept = 0, slope = 1) +
      geom_point(alpha = 1 / 2) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(values = c("red", "black"), name = "") +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      xlab("Average isoform TPM, condition 1") +
      ylab("Average isoform TPM, condition 2") +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = NA, colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(
          angle = 90, vjust = 0.5, hjust = 1, size = 15
        ),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = NA, colour = "black")
      ) +
      ggtitle("")
  )
  dev.off()
}


##### generate a truth table based on simulation summary #####
gtf.file <- glue(
  "~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf"
)

for (nsim in 1:nrep) {
  final_summary <- read.delim(
    glue("{out_dir}/sim_{nsim}/simulation_details.txt"), as.is = T
  )
  
  gene_table <- final_summary %>%
    group_by(gene_id) %>%
    summarise(
      ds_status = unique(gene_ds_status),
      de_status = unique(gene_de_status),
      ds_de_status = paste0(ds_status, ":", de_status),
      TPM = sum(TPM),
      nbr_isoforms = unique(nbr_isoforms),
      nbr_isoforms_above5 = length(which(IsoPct >= 5)),
      nbr_isoforms_above10 = length(which(IsoPct >= 10)),
      nbr_isoforms_above15 = length(which(IsoPct >= 15)),
      nbr_isoforms_above25 = length(which(IsoPct >= 25)),
      diff_IsoPct = unique(diff_IsoPct)
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  gene_table$TPM_2 <- cut2(gene_table$TPM, g = 2)
  gene_table$TPM_5 <- cut2(gene_table$TPM, g = 5)
  gene_table$TPM_10 <- cut2(gene_table$TPM, g = 10)
  gene_table$nbr_isoforms_2 <- cut2(gene_table$nbr_isoforms, g = 2)
  gene_table$nbr_isoforms_5 <- cut2(gene_table$nbr_isoforms, g = 5)
  gene_table$nbr_isoforms_10 <- cut2(gene_table$nbr_isoforms, g = 10)
  gene_table$diff_IsoPct_2 <- cut2(gene_table$diff_IsoPct, g = 2)
  gene_table$diff_IsoPct_5 <- cut2(gene_table$diff_IsoPct, g = 5)
  gene_table$diff_IsoPct_10 <- cut2(gene_table$diff_IsoPct, g = 10)
  gene_table$diff_IsoPct_3eq <- cut2(
    gene_table$diff_IsoPct, cuts = c(0, 1 / 3, 2 / 3, 1)
  )
  
  tmp <- final_summary %>%
    group_by(gene_id) %>%
    summarise(
      topisoform = transcript_id[which.max(IsoPct)],
      secisoform = ifelse(
        length(transcript_id) > 1,
        transcript_id[order(IsoPct, decreasing = T)[2]],
        "NA"
      )
    ) %>%
    ungroup() %>%
    as.data.frame()
  txdb <- makeTxDbFromGFF(gtf.file, format = "gtf")
  ebt <- exonsBy(txdb, "tx", use.names = T)
  ebt2 <- ebt[setdiff(c(tmp$topisoform, tmp$secisoform), "NA")]
  tmp$diffbp <- sapply(1:nrow(tmp), function(i) {
    if (any(tmp[i, ] == "NA")) {
      0
    } else {
      a <- ebt2[[unlist(tmp[i, "topisoform"])]]
      b <- ebt2[[unlist(tmp[i, "secisoform"])]]
      sum(width(GenomicRanges::union(a, b))) -
        sum(width(GenomicRanges::intersect(a, b)))
    }
  })
  tmp$diffbp_2 <- cut2(tmp$diffbp, g = 2)
  tmp$diffbp_5 <- cut2(tmp$diffbp, g = 5)
  tmp$diffbp_10 <- cut2(tmp$diffbp, g = 10)
  gene_table <- merge(
    gene_table,
    tmp[, c("gene_id", "diffbp", "diffbp_2", "diffbp_5", "diffbp_10")],
    by = "gene_id", all = T
  )
  
  idx <- which(colnames(gene_table) == "gene_id")
  colnames(gene_table)[idx] <- "gene"
  
  write.table(
    gene_table, file = glue("{out_dir}/sim_{nsim}/truth.txt"),
    row.names = F, col.names = T, sep = "\t", quote = F
  )
}


