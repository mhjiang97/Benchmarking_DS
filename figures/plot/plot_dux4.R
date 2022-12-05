##### utils #####
psiFilter <- function(dat, na_nbr = NULL, mean_rate = NULL, anno_field = NULL, mode = "NA", ...) {
  library(BiocParallel)
  if (is.null(anno_field)) {
    psi <- dat
  } else {
    psi <- dat[, !(1:ncol(dat) %in% anno_field), drop = F]
    anno <- dat[, anno_field, drop = F]
  }
  
  if (mode == "NA" && !is.null(na_nbr)) {
    rm_idx <- do.call(base::c, bplapply(1:nrow(psi), function(i) {
      if (sum(is.na(psi[i, ])) >= na_nbr) return(i)
    }, BPPARAM = MulticoreParam(...)))
  }
  
  if (mode == "mean" && !is.null(mean_rate)) {
    row_means <- rowMeans(psi, na.rm = T)
    if (mean_rate < 0.5) rm_idx <- which(row_means <= mean_rate) else rm_idx <- which(row_means >= mean_rate)
  }
  
  if (mode == "both"  && !is.null(na_nbr) && !is.null(mean_rate)) {
    rm_na_idx <- do.call(base::c, bplapply(1:nrow(psi), function(i) {
      if (sum(is.na(psi[i, ])) >= na_nbr) return(i)
    }, BPPARAM = MulticoreParam(...)))
    row_means <- rowMeans(psi, na.rm = T)
    if (mean_rate < 0.5) rm_mean_idx <- which(row_means <= mean_rate) else rm_mean_idx <- which(row_means >= mean_rate)
    rm_idx <- unique(c(rm_na_idx, rm_mean_idx))
  }
  
  if (!length(rm_idx) || is.null(rm_idx)) {
    return(dat)
  } else {
    if (is.null(anno_field)) {
      return(psi[-rm_idx, ])
    } else {
      return(cbind(anno[-rm_idx, ], psi[-rm_idx, ]))
    }
  }
}

psiFilter01 <- function(dat, nbr = NULL, anno_field = NULL, ...) {
  library(BiocParallel)
  if (is.null(anno_field)) {
    psi <- dat
  } else {
    psi <- dat[, !(1:ncol(dat) %in% anno_field), drop = F]
    anno <- dat[, anno_field, drop = F]
  }
  
  rm_idx <- do.call(c, bplapply(1:nrow(psi), function(i) {
    if (length(which(psi[i, ] == 0)) >= nbr || length(which(psi[i, ] == 1)) >= nbr) return(i)
  }, BPPARAM = MulticoreParam(...)))
  
  if (is.null(rm_idx)) {
    return(dat)
  } else {
    if (is.null(anno_field)) {
      return(psi[-rm_idx, ])
    } else {
      return(cbind(anno[-rm_idx, ], psi[-rm_idx, ]))
    }
  }
}

# library(conflicted)
library(RColorBrewer)
library(vroom)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(magrittr)
library(yyplot)
library(UpSetR)
library(ggimage)
library(ggplotify)
library(ggcorrplot2)
library(corrplot)
library(clusterProfiler)
library(ComplexHeatmap)

##### back-up colors #####
color_lib <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#F4B3BE",
  "#F4A11D", "#8DC8ED", "#4C6CB0", "#8A1C1B", "#CBCC2B", "#EA644C",
  "#634795", "#005B1D", "#26418A", "#CB8A93", "#F1E404", "#E22826"
)
## ration out colors in RColorBrewer for different as types
colors_as <- c(brewer.pal(11, "Set3"), "#D3D3D3")
names(colors_as) <- c(
  "a3ss", "a5ss", "mxe", "ri", "se", "altstart",
  "altend", "sme", "ce", "te", "ts", "others"
)
## swap the color for ce to make figs look better
colors_as["ce"] <- "#FFED6F"

##### change the work directory, read psi or result files #####
setwd("~/projects/dux4/analysis/rnaseq/as/")
##### CASH #####
## read the cash result file then add ds_type and gene columns
cash_result <- vroom(
  "CASH/final.alldiff.non_novel.txt.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
cash_result <- cash_result %>% mutate(
  ds_type = case_when(
    `SplicingType` == "A3SS" ~ "a3ss",
    `SplicingType` == "A5SS" ~ "a5ss",
    `SplicingType` == "MXE" ~ "mxe",
    `SplicingType` == "IR" ~ "ri",
    `SplicingType` == "Cassette" ~ "se",
    `SplicingType` == "AltStart" ~ "altstart",
    `SplicingType` == "AltEnd" ~ "altend",
    `SplicingType` == "Cassette_multi" ~ "sme"
  ), gene = AccID
)
## remove rows which don't meet the criterion
rm_idx <- vector("list", nrow(cash_result))
for (i in 1:nrow(cash_result)) {
  in_1 <- str_split(
    cash_result$`dux4_Junc_Inclusive::Exclusive`[i], "::"
  )[[1]][1] %>% as.numeric()
  ex_1 <- str_split(
    cash_result$`dux4_Junc_Inclusive::Exclusive`[i], "::"
  )[[1]][2] %>% as.numeric()
  psi_1 <- in_1 / (in_1 + ex_1)

  in_2 <- str_split(
    cash_result$`others_dux4_Junc_Inclusive::Exclusive`[i], "::"
  )[[1]][1] %>% as.numeric()
  ex_2 <- str_split(
    cash_result$`others_dux4_Junc_Inclusive::Exclusive`[i], "::"
  )[[1]][2] %>% as.numeric()
  psi_2 <- in_2 / (in_2 + ex_2)
  
  if (is.nan(psi_1) || is.nan(psi_2)) {rm_idx[[i]] <- i; next}
  if (mean(c(psi_1, psi_2)) <= 0.03) rm_idx[[i]] <- i
}
rm_idx <- purrr::flatten_int(rm_idx)
cash_result <- cash_result[-rm_idx, ]

##### leafcutter #####
## read the filtered per individual count file (filtered out rows containing over 30 NAs among all 100 samples)
leafcutter <- vroom(
  "leafcutter/cluster/leafcutter.fil.txt.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
## filter out rows containing over 70 zeros or ones
leafcutter <- psiFilter01(leafcutter, 70, 1, 20) # leafcutter <- psiFilter(leafcutter, mean_rate = 0.05, mode = "mean", anno_field = 1, 10)
## filter out rows whose mean psi value is less than 0.03
leafcutter <- psiFilter(
  leafcutter, mean_rate = 0.03, mode = "mean", anno_field = 1, 10
)
## count how many clusters leafcutter finds
clus <- vector(mode = "character", length = nrow(leafcutter))
for (i in 1:nrow(leafcutter)) {
  clus[[i]] <- str_split(leafcutter[i, 1], ":")[[1]][4]
}
leafcutter_cluster <- data.frame(cluster = unique(clus), ds_type = "others")
## read the leafcutter ds result to annotate clusters with gene symbols
leafcutter_result <- vroom(
  "leafcutter/ds/leafcutter_cluster_significance.txt.gz", delim = "\t", col_types = cols()
) %>% group_by(cluster) %>%
  mutate(clus = str_split(cluster, ":")[[1]][2]) %>% ungroup() %>% as.data.frame() %>%
  dplyr::filter(status == "Success", !is.na(p), !is.na(genes))

leafcutter_multi <- dplyr::filter(leafcutter_result, grepl(",", genes))
lf_ml_list <- list()
for (i in 1:nrow(leafcutter_multi)) {
  tmp_genes <- str_split(leafcutter_multi$genes[i], ",") %>% unlist()
  tmp_tbl <- data.frame()
  for (g in tmp_genes) {
    tmp_tbl <- rbind(tmp_tbl, leafcutter_multi[i, ])
  }
  tmp_tbl$genes <- tmp_genes
  lf_ml_list[[i]] <- tmp_tbl
}
lf_ml <- bind_rows(lf_ml_list)
leafcutter_result <- dplyr::filter(leafcutter_result, !grepl(",", genes))
leafcutter_result <- rbind(leafcutter_result, lf_ml)

leafcutter_cluster <- left_join(
  leafcutter_cluster, leafcutter_result[, c("genes", "clus")], by = c("cluster" = "clus")
)
colnames(leafcutter_cluster)[3] <- "gene"

##### leafcutter with regtools #####
leafcutter_regtools <- vroom(
  "leafcutter/cluster_regtools/leafcutter.fil.txt.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
leafcutter_regtools <- psiFilter01(leafcutter_regtools, 70, 1, 20) # leafcutter_regtools <- psiFilter(leafcutter_regtools, mean_rate = 0.05, mode = "mean", anno_field = 1, 10)
leafcutter_regtools <- psiFilter(
  leafcutter_regtools, mean_rate = 0.03, mode = "mean", anno_field = 1, 10
)
clus_regtools <- vector("character", nrow(leafcutter_regtools))
for (i in 1:nrow(leafcutter_regtools)) {
  clus_regtools[[i]] <- str_split(leafcutter_regtools[i, 1], ":")[[1]][4]
}
leafcutter_regtools_cluster <- data.frame(
  cluster = unique(clus_regtools), ds_type = "others"
)
leafcutter_regtools_result <- vroom(
  "leafcutter/ds_regtools/leafcutter_cluster_significance.txt.gz", delim = "\t", col_types = cols()
) %>% group_by(cluster) %>%
  mutate(clus = str_split(cluster, ":")[[1]][2]) %>% ungroup() %>% as.data.frame() %>%
  dplyr::filter(!is.na(p), !is.na(genes))

leafcutter_regtools_multi <- dplyr::filter(leafcutter_regtools_result, grepl(",", genes))
lf_rt_ml_list <- list()
for (i in 1:nrow(leafcutter_regtools_multi)) {
  tmp_genes <- str_split(leafcutter_regtools_multi$genes[i], ",") %>% unlist()
  tmp_tbl <- data.frame()
  for (g in tmp_genes) {
    tmp_tbl <- rbind(tmp_tbl, leafcutter_regtools_multi[i, ])
  }
  tmp_tbl$genes <- tmp_genes
  lf_rt_ml_list[[i]] <- tmp_tbl
}
lf_rt_ml <- bind_rows(lf_rt_ml_list)
leafcutter_regtools_result <- dplyr::filter(leafcutter_regtools_result, !grepl(",", genes))
leafcutter_regtools_result <- rbind(leafcutter_regtools_result, lf_rt_ml)

leafcutter_regtools_cluster <- left_join(
  leafcutter_regtools_cluster, leafcutter_regtools_result[, c("genes", "clus")],
  by = c("cluster" = "clus")
)
colnames(leafcutter_regtools_cluster)[3] <- "gene"

##### MAJIQ #####
## read the majiq psi file and annotate different lsv types and tag complex ones as others
majiq_psi <- vroom(
  "MAJIQ/non_novel/psi/both.psi.tsv.gz",
  delim = "\t", col_types = cols()
)
majiq_psi <- majiq_psi %>% mutate(
  ds_type = case_when(
  `Num. Junctions` == 1 & `ES` == "FALSE" & `A5SS` == "FALSE" & `A3SS` == "FALSE" & !is.na(`IR coords`) ~ "ri",
  `Num. Junctions` == 2 & `Num. Exons` == 3 ~ "se",
  `Num. Junctions` == 4 & `Num. Exons` == 4 ~ "mxe",
  `Num. Junctions` %in% c(2, 3) & `A3SS` == "TRUE" ~ "a3ss",
  `Num. Junctions` %in% c(2, 3) & `A5SS` == "TRUE" ~ "a5ss"
  ), gene = `Gene ID`
)
majiq_psi$ds_type[is.na(majiq_psi$ds_type)] <- "others"
## filter out rows that don't meet the criterion
rm_idx <- vector("list", nrow(majiq_psi))
for (i in 1:nrow(majiq_psi)) {
  psis <- str_split(majiq_psi$`E(PSI) per LSV junction`[i], ";")[[1]] %>%
    as.numeric()
  psi2 <- psis[head(order(psis, decreasing = T), 2)]
  if (min(psi2) <= 0.03) rm_idx[[i]] <- i
}
rm_idx <- flatten_int(rm_idx)
majiq_psi <- majiq_psi[-rm_idx, ]

majiq_gene_id <- vector("list", nrow(majiq_psi))
majiq_lsv_id <- vector("list", nrow(majiq_psi))
majiq_lsv_type <- vector("list", nrow(majiq_psi))
majiq_ds_type <- vector("list", nrow(majiq_psi))
majiq_psii <- vector("list", nrow(majiq_psi))
for (i in 1:nrow(majiq_psi)) {
  tmp_psi <- str_split(majiq_psi$`E(PSI) per LSV junction`[i], ";")[[1]] %>% as.numeric()
  tmp_gene_id <- rep(majiq_psi$`Gene ID`[i], length(tmp_psi))
  tmp_lsv_id <- rep(majiq_psi$`LSV ID`[i], length(tmp_psi))
  tmp_lsv_type <- rep(majiq_psi$`LSV Type`[i], length(tmp_psi))
  tmp_ds_type <- rep(majiq_psi$ds_type[i], length(tmp_psi))
  majiq_gene_id[[i]] <- tmp_gene_id
  majiq_lsv_id[[i]] <- tmp_lsv_id
  majiq_lsv_type[[i]] <- tmp_lsv_type
  majiq_psii[[i]] <- tmp_psi
  majiq_ds_type[[i]] <- tmp_ds_type
}
majiq_gene_id <- flatten_chr(majiq_gene_id)
majiq_lsv_id <- flatten_chr(majiq_lsv_id)
majiq_lsv_type <- flatten_chr(majiq_lsv_type)
majiq_ds_type <- flatten_chr(majiq_ds_type)
majiq_psii <- flatten_dbl(majiq_psii)
majiq_psi2 <- data.frame(
  gene_id = majiq_gene_id, psi_median = majiq_psii,
  lsv_type = majiq_lsv_type, lsv_id = majiq_lsv_id, ds_type = majiq_ds_type
)

##### MISO #####
## read filtered miso psi files (filtered out rows containing over 30 NAs among all 100 samples)
miso_psi_list <- vector("list", 5)
names(miso_psi_list) <- c("a3ss", "a5ss", "mxe", "ri", "se")
for (type in names(miso_psi_list)) {
  miso_psi_list[[type]] <- vroom(sprintf("MISO/psi/psi.%s.fil.txt.gz", type),
    delim = "\t", col_types = cols()
  ) %>% mutate(ds_type = type) %>% as.data.frame()
}
miso_psi <- bind_rows(miso_psi_list)
miso_psi <- psiFilter01(miso_psi, 70, c(1:2, 103), 10) # miso_psi <- psiFilter(miso_psi, mean_rate = 0.05, mode = "mean", anno_field = c(1:2, 103), 10)
miso_psi <- psiFilter(
  miso_psi, mean_rate = 0.03, mode = "mean", anno_field = c(1:3), 10
)
## annotate the miso psi file with gene symbol
miso_psi <- miso_psi %>% group_by(event_name) %>%
  mutate(
    chr = str_split(event_name, "[:|-]")[[1]][1],
    coor_1 = str_split(event_name, "[:|-]")[[1]][2],
    coor_2 = str_split(event_name, "[:|-]")[[1]][3]
  ) %>% ungroup() %>% as.data.frame()
for (i in 1:nrow(miso_psi)) {
  if (as.numeric(miso_psi$coor_1[i]) > as.numeric(miso_psi$coor_2[i])) {
    a <- miso_psi$coor_2[i]
    miso_psi$coor_2[i] <- miso_psi$coor_1[i]
    miso_psi$coor_1[i] <- a
  }
}
gtf <- vroom(
  "~/doc/tx2gene/gencode.v29.txt.gz", col_names = F, delim = "\t",
  col_types = paste0(rep("c", 11), collapse = "")
) %>% as.data.frame()
miso_psi <- left_join(miso_psi, gtf[, c(1, 4:5, 10)],
  by = c("chr" = "X1", "coor_1" = "X4", "coor_2" = "X5")
)
miso_psi <- miso_psi[match(unique(miso_psi$event_name), miso_psi$event_name), ]
colnames(miso_psi)[colnames(miso_psi) == "X10"] <- "gene"

##### rMATs #####
## read filtered rmats psi files (filtered out rows containing over 30 NAs among all 100 samples)
rmats_psi_list <- vector("list", 5)
names(rmats_psi_list) <- c("a3ss", "a5ss", "mxe", "ri", "se")
for (type in names(rmats_psi_list)) {
  rmats_psi_list[[type]] <- vroom(sprintf("rMATs/psi/psi.%s.fil.txt.gz", type),
    delim = "\t", col_types = cols(), col_names = F
  ) %>% mutate(ds_type = type) %>% as.data.frame()
}
## read rmats annotation files to determine whether a event is novel or not
rmats_anno_list <- vector("list", 5)
names(rmats_anno_list) <- c("a3ss", "a5ss", "mxe", "ri", "se")
for (type in names(rmats_anno_list)) {
  rmats_anno_list[[type]] <- rbind(
    vroom(sprintf(
      "rMATs/results/fromGTF.novelJunction.%s.txt.gz",
      toupper(type)
    ), delim = "\t", col_types = cols()) %>% as.data.frame(),
    vroom(sprintf(
      "rMATs/results/fromGTF.novelSpliceSite.%s.txt.gz",
      toupper(type)
    ), delim = "\t", col_types = cols()) %>% as.data.frame()
  ) %>% mutate(ds_typ = type)
}
## filter out novel ones
for (type in names(rmats_psi_list)) {
  novel_id <- rmats_anno_list[[type]][, 1]
  rmats_psi_list[[type]] <- rmats_psi_list[[type]][
    -which(rmats_psi_list[[type]][, 1] %in% novel_id),
  ]
}
rmats_psi <- bind_rows(rmats_psi_list)
rmats_psi <- psiFilter01(rmats_psi, 70, c(1:2, 103), 10) # rmats_psi <- psiFilter(rmats_psi, mean_rate = 0.05, mode = "mean", anno_field = c(1:2, 103), 10)
rmats_psi <- psiFilter(
  rmats_psi, mean_rate = 0.03, mode = "mean", anno_field = c(1:3), 10
)
rmats_psi$gene <- rmats_psi$X2

##### SplAdder #####
## read filtered spladder psi files (filtered out rows containing over 30 NAs among all 100 samples)
spladder_psi_list <- vector("list", 6)
names(spladder_psi_list) <- c(
  "alt_3prime", "alt_5prime", "mutex_exons",
  "intron_retention", "exon_skip", "mult_exon_skip"
)
for (type in names(spladder_psi_list)) {
  spladder_psi_list[[type]] <- vroom(
    sprintf("spladder/non_novel/psi/psi.%s.fil.txt.gz", type),
    delim = "\t", col_types = cols()
  ) %>% mutate(type = type) %>% as.data.frame()
}
spladder_psi <- bind_rows(spladder_psi_list)
spladder_psi <- psiFilter01(spladder_psi, 70, c(1:10, 111:125), 10) # spladder_psi <- psiFilter(spladder_psi, mean_rate = 0.05, mode = "mean", anno_field = c(1:10, 111:125), 10)
spladder_psi <- psiFilter(
  spladder_psi, mean_rate = 0.03, mode = "mean", anno_field = c(1:25), 10
)
spladder_psi <- spladder_psi %>% mutate(
  ds_type = case_when(
  `type` == "alt_3prime" ~ "a3ss",
  `type` == "alt_5prime" ~ "a5ss",
  `type` == "mutex_exons" ~ "mxe",
  `type` == "intron_retention" ~ "ri",
  `type` == "exon_skip" ~ "se",
  `type` == "mult_exon_skip" ~ "sme"
  ), gene = gene_name
)

##### SUPPA #####
## read filtered suppa psi files (filtered out rows containing over 30 NAs among all 100 samples)
suppa_psi_list <- vector("list", 7)
names(suppa_psi_list) <- c("a3", "a5", "mx", "ri", "se", "af", "al")
for (type in names(suppa_psi_list)) {
  suppa_psi_list[[type]] <- vroom(
    sprintf("SUPPA/psi/psi.%s.fil.txt.gz", type), delim = "\t", col_types = cols()
  ) %>% as.data.frame() %>%
    {
      colnames(.)[1] <- "ID"
      .
    } %>% mutate(type = type)
}
suppa_psi <- bind_rows(suppa_psi_list)
suppa_psi <- psiFilter01(suppa_psi, 70, c(1, 102), 10) # suppa_psi <- psiFilter(suppa_psi, mean_rate = 0.05, mode = "mean", anno_field = c(1, 102), 10)
suppa_psi <- psiFilter(
  suppa_psi, mean_rate = 0.03, mode = "mean", anno_field = c(1:2), 10
)
suppa_psi <- suppa_psi %>%
  mutate(ds_type = case_when(
    `type` == "a3" ~ "a3ss",
    `type` == "a5" ~ "a5ss",
    `type` == "mx" ~ "mxe",
    `type` == "ri" ~ "ri",
    `type` == "se" ~ "se",
    `type` == "af" ~ "altstart",
    `type` == "al" ~ "altend"
    )
  ) %>% group_by(ID) %>% mutate(gene = str_split(ID, ";")[[1]][1]) %>% ungroup() %>% as.data.frame()

##### Whippet #####
## read the filtered whippet psi file (filtered out rows containing over 30 NAs among all 100 samples)
whippet_psi <- vroom("Whippet/psi/psi.fil.txt.gz",
  delim = "\t", col_types = cols()
) %>% as.data.frame()
whippet_psi <- psiFilter01(whippet_psi, 70, c(1:5), 10) # whippet_psi <- psiFilter(whippet_psi, mean_rate = 0.05, mode = "mean", anno_field = c(1:5), 10)
whippet_psi <- psiFilter(
  whippet_psi, mean_rate = 0.03, mode = "mean", anno_field = c(1:5), 10
)
whippet_psi <- whippet_psi %>% mutate(ds_type = case_when(
  `Type` == "AA" ~ "a3ss",
  `Type` == "AD" ~ "a5ss",
  `Type` == "RI" ~ "ri",
  `Type` == "AF" ~ "altstart",
  `Type` == "AL" ~ "altend",
  `Type` == "CE" ~ "ce",
  `Type` == "TE" ~ "te",
  `Type` == "TS" ~ "ts"
  ), gene = Gene
)

if (F) {
  save.image("psi_non_novel.RData")
}
##### create descriptive data and plot #####
##### create a list incorporating filtered events from all softwares #####
softwares <- c(
  "CASH", "leafcutter", "leafcutter_regtools", "MAJIQ",
  "MISO", "rMATs", "SplAdder", "SUPPA", "Whippet"
)
sw_list <- list(
  cash_result, leafcutter_cluster, leafcutter_regtools_cluster,
  majiq_psi, miso_psi, rmats_psi, spladder_psi, suppa_psi, whippet_psi
)
names(sw_list) <- softwares
##### plot event numbers and event proportion for each software #####
## event number data
dat_event_nbr <- data.frame(
  matrix(0, nrow = length(sofrwares), ncol = length(colors_as)), row.names = softwares
)
colnames(dat_event_nbr) <- names(colors_as)
for (sw in softwares) {
  sw_tbl <- table(sw_list[[sw]]$ds_type) %>%
    {
      .[match(colnames(dat_event_nbr), names(.))] -> .
      .
    } %>%
    as.numeric()
  dat_event_nbr[sw, ] <- sw_tbl
}
## melt the data for plotting
dat_event_nbr_m <- melt(dat_event_nbr)
## add a software column
dat_event_nbr_m$softwares <- rep(softwares, 12)
## adjust the plot order manually
dat_event_nbr_m$variable <- factor(
  dat_event_nbr_m$variable,
  levels = c(
    "mxe", "sme", "altend", "altstart", "a5ss", "a3ss",
    "others", "ri", "se", "ts", "ce", "te"
  )
)
dat_event_nbr_m$softwares <- factor(
  dat_event_nbr_m$softwares,
  levels = c(
    "leafcutter_regtools", "leafcutter", "SUPPA", "SplAdder",
    "CASH", "MISO", "rMATs", "MAJIQ", "Whippet"
  )
)
dat_event_nbr_m <- dat_event_nbr_m[!is.na(dat_event_nbr_m$value), ]

##### mosaicplot #####
a <- c()
b <- c()
for(i in 1:nrow(dat_event_nbr_m)){
a <- c(a, rep(as.character(dat_event_nbr_m$variable[i]), dat_event_nbr_m$value[i]))
b <- c(b, rep(as.character(dat_event_nbr_m$softwares[i]), dat_event_nbr_m$value[i]))
}
dat_event_nbr_m2 <- data.frame(AS_types = a, softwares = b)
dat_event_nbr_m2$AS_types <- factor(dat_event_nbr_m2$AS_types, levels = names(colors_as))
pdf("figs/mosaic_events.pdf")
mosaicplot(softwares ~ AS_types, data = dat_event_nbr_m2, color = colors_as, main = "", shade = F)
dev.off()

## the event number fig
## plot p1 and p2 separately and joint them
ggplot(dat_event_nbr_m, mapping = aes(y = softwares, x = value, fill = variable)) +
  geom_bar(
    width = 1, size = 0.1, col = "black", stat = "identity",
    position = position_dodge2(
      width = 1.2, padding = 0.25,
      preserve = "single"
    )
  ) + xlab("events") + scale_fill_manual(values = colors_as) + # scale_fill_brewer(type = "qual", palette = "Set3")
  coord_cartesian(xlim = c(0, 11000)) +
  theme(panel.grid.major = element_line(colour = NA)) + theme_bw() -> p1
  
ggplot(dat_event_nbr_m, mapping = aes(y = softwares, x = value, fill = variable)) +
  geom_bar(
    width = 1, size = 0.1, col = "black", stat = "identity",
    position = position_dodge2(
      width = 1.2, padding = 0.25,
      preserve = "single"
    )
  ) + labs(x = NULL, y = NULL, fill = NULL) + scale_fill_manual(values = colors_as) +
  theme(
    axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_line(colour = NA)
  ) + theme_bw() +
  scale_x_continuous(labels = scales::scientific, breaks = c(18000, 25000, 50000)) +
  coord_cartesian(xlim = c(15000, 50000)) -> p2

ggarrange(
  p1, p2, widths = c(8 / 9, 1 / 9), ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "right", align = "h"
) -> p
ggsave("figs/events_nbr.pdf", p, width = 6, height = 7)

## the proportion fig
if (F) {
dat_event_nbr_m_fil_other <- dat_event_nbr_m[!dat_event_nbr_m$variable == "others",]
dat_event_nbr_m_fil_other$variable <- factor(
  dat_event_nbr_m_fil_other$variable,
  levels = c(
    "mxe", "sme", "altend", "altstart", "a3ss", "a5ss", "ri", "se", "ts", "ce", "te"
  )
)
}
dat_event_nbr_m$variable <- factor(
  dat_event_nbr_m$variable,
  levels = c(
    "mxe", "sme", "altend", "altstart", "a5ss",
    "a3ss", "ri", "others", "se", "ts", "ce", "te"
  )
)
ggplot(dat_event_nbr_m, mapping = aes(y = softwares, x = value, fill = variable)) +
  geom_bar(
    size = 0.1, col = "black", stat = "identity",
    position = position_fill(0.5), width = 0.7
  ) +
  scale_fill_manual(values = colors_as) +
  xlab("events proportion") +
  theme_set(theme_bw()) +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  ) -> p
ggsave("figs/proportion_events.pdf", p, width = 8, height = 7)

##### plot gene numbers and gene proportion for each software #####
## gene number data
dat_gene_nbr <- data.frame(
  matrix(0, nrow = length(softwares), ncol = length(colors_as)), row.names = softwares
)
colnames(dat_gene_nbr) <- names(colors_as)
for (sw in softwares) {
  for (type in unique(sw_list[[sw]]$ds_type)) {
    tmp_genes <- sw_list[[sw]]$gene[sw_list[[sw]]$ds_type == type]
    dat_gene_nbr[sw, type] <- length(unique(tmp_genes[!is.na(tmp_genes)]))
  }
}
dat_gene_nbr_m <- melt(dat_gene_nbr)
dat_gene_nbr_m$softwares <- rep(softwares, 12)
dat_gene_nbr_m$variable <- factor(
  dat_gene_nbr_m$variable,
  levels = c(
    "mxe", "sme", "altend", "altstart", "a5ss",
    "a3ss", "others", "ri", "se", "ts", "ce", "te"
  )
)
dat_gene_nbr_m <- dat_gene_nbr_m[order(dat_gene_nbr_m$variable), ]
dat_gene_nbr_m$softwares <- factor(dat_gene_nbr_m$softwares,
  levels = c(
    "leafcutter_regtools", "leafcutter",
    "SUPPA", "SplAdder", "CASH", "MISO",
    "rMATs", "MAJIQ", "Whippet"
  )
)
dat_gene_nbr_m <- dat_gene_nbr_m[-which(dat_gene_nbr_m$value == 0), ]

##### mosaicplot #####
a <- c()
b <- c()
for(i in 1:nrow(dat_gene_nbr_m)){
  a <- c(a, rep(as.character(dat_gene_nbr_m$variable[i]), dat_gene_nbr_m$value[i]))
  b <- c(b, rep(as.character(dat_gene_nbr_m$softwares[i]), dat_gene_nbr_m$value[i]))
}
dat_gene_nbr_m2 <- data.frame(AS_types = a, softwares = b)
dat_gene_nbr_m2$AS_types <- factor(dat_gene_nbr_m2$AS_types, levels = names(colors_as))
pdf("figs/mosaic_genes.pdf")
mosaicplot(softwares ~ AS_types, data = dat_gene_nbr_m2, color = colors_as, main = "", shade = F)
dev.off()
## the gene number fig
ggplot(dat_gene_nbr_m, mapping = aes(y = softwares, x = value, fill = variable)) +
  geom_bar(
    size = 0.1, col = "black", stat = "identity",
    position = position_dodge2(
      width = 0.95, padding = 0.25,
      preserve = "single"
    ), width = 0.9
  ) + xlab("genes") + scale_fill_manual(values = colors_as) + theme_bw() +
  theme(panel.grid.major = element_line(colour = NA)) -> p
ggsave("figs/genes_nbr.pdf", p, width = 6, height = 7)

## the proportion fig
if (F) {
dat_gene_nbr_m_fil_other <- dat_gene_nbr_m[dat_gene_nbr_m$variable %in% c("altstart", "altend", "sme", "se", "ri", "a3ss", "a5ss", "mxe"), ]
dat_gene_nbr_m_fil_other$variable <- factor(
  dat_gene_nbr_m_fil_other$variable, levels = c("mxe", "sme", "altend", "altstart", "a3ss", "a5ss", "ri", "se")
)
}
dat_gene_nbr_m$variable <- factor(
  dat_gene_nbr_m$variable,
  levels = c(
    "mxe", "sme", "altend", "altstart", "a5ss", "a3ss",
    "ri", "others", "se", "ts", "ce", "te"
  )
)
ggplot(dat_gene_nbr_m, mapping = aes(y = softwares, x = value, fill = variable)) +
  geom_bar(
    size = 0.1, col = "black", stat = "identity",
    position = position_fill(0.5), width = 0.7
  ) +
  scale_fill_manual(values = colors_as) +
  xlab("events proportion") +
  theme_set(theme_bw()) +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  ) -> p
ggsave("figs/proportion_genes.pdf", p, width = 8, height = 7)

##### plot ratios between the gene number and the event number #####
## add a total number column
dat_event_nbr <- dat_event_nbr %>% group_by_all() %>% mutate(
  total = sum(a3ss, a5ss, mxe, ri, se, others, sme, altend, altstart, na.rm = T)
) %>% ungroup() %>% as.data.frame() %>% set_rownames(softwares)
dat_gene_nbr <- dat_gene_nbr %>% group_by_all() %>% mutate(
  total = sum(a3ss, a5ss, mxe, ri, se, others, sme, altend, altstart, na.rm = T)
) %>% ungroup() %>% as.data.frame() %>% set_rownames(softwares)
dat_per_gene_events <- dat_event_nbr / dat_gene_nbr
dat_per_gene_events_m <- melt(dat_per_gene_events)
dat_per_gene_events_m$softwares <- rep(c(softwares), 13)
dat_per_gene_events_m$variable <- factor(
  dat_per_gene_events_m$variable, levels = c(
    "total", "mxe", "sme", "altend", "altstart", "a5ss",
    "a3ss", "others", "ri", "se", "ts", "ce", "te"
  )
)
dat_per_gene_events_m$softwares <- factor(
  dat_per_gene_events_m$softwares, levels = c(
    "leafcutter_regtools", "leafcutter", "SUPPA", "SplAdder",
    "CASH", "MISO", "rMATs", "MAJIQ", "Whippet"
  )
)
dat_per_gene_events_m <- dat_per_gene_events_m[!is.na(dat_per_gene_events_m$value), ]
ggplot(
  dat_per_gene_events_m,
  mapping = aes(x = softwares, y = value, color = "#000000", fill = variable)
) +
  geom_bar(
    size = 0.1, col = "black", stat = "identity", position = position_dodge2(
      width = 0.95, padding = 0.25, preserve = "single"
    ), width = 0.9
  ) + ylab("ratio") + scale_fill_manual(values = c(colors_as, total = "#808080")) +
  theme_bw() + theme(panel.grid.major = element_line(colour = NA)) -> p
ggsave("figs/ratio.pdf", p, width = 10, height = 5)

##### facet plot gene and event number for each as type #####
a <- dat_event_nbr[, match(c(
  "a3ss", "a5ss", "altend", "altstart", "ce", "mxe",
  "ri", "se", "sme", "te", "others", "ts", "total"
), colnames(dat_event_nbr))]
colnames(a) <- paste0(colnames(a), "_event")
b <- dat_gene_nbr[, match(c(
  "a3ss", "a5ss", "altend", "altstart", "ce", "mxe",
  "ri", "se", "sme", "te", "others", "ts", "total"
), colnames(dat_gene_nbr))]
colnames(b) <- paste0(colnames(b), "_gene")
dat_nbr <- cbind(a, b)
dat_nbr_m <- melt(dat_nbr)
dat_nbr_m$softwares <- rep(softwares, 26)
dat_nbr_m$ds_type <- gsub("_.*", "", dat_nbr_m$variable)
dat_nbr_m <- dat_nbr_m[!is.na(dat_nbr_m$value), ]
colors_1 <- colors_2 <- c(colors_as, total = "#808080")
dat_nbr_m$variable <- factor(
  dat_nbr_m$variable,
  levels = c(
    "a3ss_gene", "a3ss_event", "a5ss_gene", "a5ss_event", "altend_gene", "altend_event",
    "altstart_gene", "altstart_event", "ce_gene", "ce_event", "mxe_gene", "mxe_event",
    "others_gene", "others_event", "ri_gene", "ri_event", "se_gene", "se_event",
    "sme_gene", "sme_event", "te_gene", "te_event", "total_gene", "total_event",
    "ts_gene", "ts_event"
  )
)
dat_nbr_m$ds_type <- factor(
  dat_nbr_m$ds_type,
  levels = c(
    "a3ss", "a5ss", "ri", "se", "mxe", "altstart",
    "altend", "sme", "te", "ts", "ce", "others", "total"
  )
)
dat_nbr_m$softwares <- factor(
  dat_nbr_m$softwares,
  levels = c("CASH", "MAJIQ", "MISO", "rMATs", "SplAdder", "SUPPA",
             "Whippet", "leafcutter", "leafcutter_regtools")
)
names(colors_1) <- paste0(names(colors_1), "_event")
names(colors_2) <- paste0(names(colors_2), "_gene")
ggplot(dat_nbr_m, mapping = aes(x = softwares, y = value, fill = variable)) +
  geom_bar(
    width = 0.8, size = 0.1, col = "black", stat = "identity", alpha = 5 / 7,
    position = "identity"
  ) + ylab("Number") + scale_fill_manual(values = c(colors_1, colors_2)) +
  facet_grid(ds_type ~ ., scales = "free_y", shrink = F) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text.y = element_text(size = 11, face = "bold"),
    panel.border = element_rect(fill = NA)
  ) -> p
ggsave("figs/gene_vs_events_nbr.pdf", p, width = 10, height = 10)

##### plot average gene and event numbers and ratios for the whole cohort #####
dat_as_table <- data.frame(
  row.names = c(
    "a3ss", "a5ss", "mxe", "ri", "se", "altend", "altstart", "sme", "te", "ce", "ts", "total"
  ),
  aver_gene_nbr = rep(0, 12), aver_event_nbr = rep(0, 12)
)
for (type in rownames(dat_as_table)) {
  dat_as_table[type, "aver_gene_nbr"] <- mean(dat_gene_nbr[dat_gene_nbr[, type] != 0, type])
  dat_as_table[type, "aver_event_nbr"] <- mean(dat_event_nbr[!is.na(dat_event_nbr[, type]), type])
}
dat_as_table$ratio <- dat_as_table$aver_event_nbr / dat_as_table$aver_gene_nbr
dat_as_table$ds_type <- rownames(dat_as_table)
dat_as_table$ds_type <- factor(
  dat_as_table$ds_type, levels = c(
  "mxe", "sme", "a5ss", "altend", "a3ss", "ri",
  "altstart", "se", "ts", "ce", "te", "total"
))
dat_as_table$ratio <- dat_as_table$ratio * 10000
dat_as_table_m <- melt(dat_as_table)
ggplot(dat_as_table_m, aes(x = ds_type, y = value, group = variable, linetype = variable)) +
  geom_line(aes(color = variable)) +
  scale_color_manual(values = color_lib[c(1, 13, 3)]) +
  scale_linetype_manual(values = c("solid", "solid", "dotdash")) +
  geom_point(aes(color = variable, shape = variable, stroke = 0.5)) + xlab("AS type")
  geom_text_repel(
    data = dat_as_table_m[dat_as_table_m$variable == "ratio", ],
    aes(ds_type, value, label = round(value / 10000, 3)), size = 2.5,
    direction = "both"
  ) + scale_y_continuous("Number", sec.axis = sec_axis(~ . / 10000, name = "Ratio")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 7),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = rel(1)),
    legend.key = element_blank(),
    strip.background = element_rect(
      fill = "white", colour = "black",
      size = rel(2)
    ),
    complete = TRUE
  ) -> p
ggsave("figs/average_per_as.pdf", p, width = 4, height = 3)

##### read ds result files and plot figs about ds events and ds genes #####
##### MAJIQ #####
majiq_result <- vroom(
  "MAJIQ/non_novel/voila/final.0.05.tsv.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
majiq_result <- majiq_result %>% mutate(
  ds_type = case_when(
  `Num. Junctions` == 1 & `ES` == "FALSE" & `A5SS` == "FALSE" & `A3SS` == "FALSE" & !is.na(`IR coords`) ~ "ri",
  `Num. Junctions` == 2 & `Num. Exons` == 3 ~ "se",
  `Num. Junctions` == 4 & `Num. Exons` == 4 ~ "mxe",
  `Num. Junctions` %in% c(2, 3) & `A3SS` == "TRUE" ~ "a3ss",
  `Num. Junctions` %in% c(2, 3) & `A5SS` == "TRUE" ~ "a5ss"
  )
)
majiq_result$ds_type[is.na(majiq_result$ds_type)] <- "others"
## reconstruct a new majiq result data
majiq_gene_id <- vector("list", nrow(majiq_result))
majiq_gene_symbol <- vector("list", nrow(majiq_result))
majiq_lsv_id <- vector("list", nrow(majiq_result))
majiq_lsv_type <- vector("list", nrow(majiq_result))
majiq_p <- vector("list", nrow(majiq_result))
majiq_ds_type <- vector("list", nrow(majiq_result))
majiq_dpsi <- vector("list", nrow(majiq_result))
for (i in 1:nrow(majiq_result)) {
  tmp_p <- str_split(majiq_result$`P(|dPSI|>=0.05) per LSV junction`[i], ";")[[1]] %>% as.numeric()
  tmp_gene_id <- rep(majiq_result$`Gene ID`[i], length(tmp_p))
  tmp_gene_symbol <- rep(majiq_result$`Gene Name`[i], length(tmp_p))
  tmp_lsv_id <- rep(majiq_result$`LSV ID`[i], length(tmp_p))
  tmp_lsv_type <- rep(majiq_result$`LSV Type`[i], length(tmp_p))
  tmp_ds_type <- rep(majiq_result$ds_type[i], length(tmp_p))
  tmp_dpsi <- str_split(majiq_result$`E(dPSI) per LSV junction`[i], ";")[[1]] %>% as.numeric()
  majiq_gene_id[[i]] <- tmp_gene_id
  majiq_gene_symbol[[i]] <- tmp_gene_symbol
  majiq_lsv_id[[i]] <- tmp_lsv_id
  majiq_lsv_type[[i]] <- tmp_lsv_type
  majiq_p[[i]] <- tmp_p
  majiq_ds_type[[i]] <- tmp_ds_type
  majiq_dpsi[[i]] <- tmp_dpsi
}
majiq_gene_id <- flatten_chr(majiq_gene_id)
majiq_gene_symbol <- flatten_chr(majiq_gene_symbol)
majiq_lsv_id <- flatten_chr(majiq_lsv_id)
majiq_lsv_type <- flatten_chr(majiq_lsv_type)
majiq_p <- flatten_dbl(majiq_p)
majiq_ds_type <- flatten_chr(majiq_ds_type)
majiq_dpsi <- flatten_dbl(majiq_dpsi)
majiq_result2 <- data.frame(
  gene_id = majiq_gene_id, gene_symbol = majiq_gene_symbol, dpsi = majiq_dpsi,
  lsv_type = majiq_lsv_type, lsv_id = majiq_lsv_id, p = majiq_p, ds_type = majiq_ds_type
)
majiq_result2 <- majiq_result2 %>% group_by(lsv_id) %>%
  dplyr::mutate(tmp = dpsi[which.max(abs(dpsi))]) %>% ungroup() %>% as.data.frame()
##### MAJIQ showall #####
majiq_all_result <- vroom(
  "MAJIQ/non_novel/voila/final.0.05.showall.tsv.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
majiq_all_result <- majiq_all_result %>% mutate(
  ds_type = case_when(
    `Num. Junctions` == 1 & `ES` == "FALSE" & `A5SS` == "FALSE" & `A3SS` == "FALSE" & !is.na(`IR coords`) ~ "ri",
    `Num. Junctions` == 2 & `Num. Exons` == 3 ~ "se",
    `Num. Junctions` == 4 & `Num. Exons` == 4 ~ "mxe",
    `Num. Junctions` %in% c(2, 3) & `A3SS` == "TRUE" ~ "a3ss",
    `Num. Junctions` %in% c(2, 3) & `A5SS` == "TRUE" ~ "a5ss"
  )
)
majiq_all_result$ds_type[is.na(majiq_all_result$ds_type)] <- "others"
## reconstruct a new majiq result data
majiq_gene_id <- vector("list", nrow(majiq_all_result))
majiq_gene_symbol <- vector("list", nrow(majiq_all_result))
majiq_lsv_id <- vector("list", nrow(majiq_all_result))
majiq_lsv_type <- vector("list", nrow(majiq_all_result))
majiq_p <- vector("list", nrow(majiq_all_result))
majiq_ds_type <- vector("list", nrow(majiq_all_result))
majiq_dpsi <- vector("list", nrow(majiq_all_result))
for (i in 1:nrow(majiq_all_result)) {
  tmp_p <- str_split(majiq_all_result$`P(|dPSI|>=0.05) per LSV junction`[i], ";")[[1]] %>% as.numeric()
  tmp_gene_id <- rep(majiq_all_result$`Gene ID`[i], length(tmp_p))
  tmp_gene_symbol <- rep(majiq_all_result$`Gene Name`[i], length(tmp_p))
  tmp_lsv_id <- rep(majiq_all_result$`LSV ID`[i], length(tmp_p))
  tmp_lsv_type <- rep(majiq_all_result$`LSV Type`[i], length(tmp_p))
  tmp_ds_type <- rep(majiq_all_result$ds_type[i], length(tmp_p))
  tmp_dpsi <- str_split(majiq_all_result$`E(dPSI) per LSV junction`[i], ";")[[1]] %>% as.numeric()
  majiq_gene_id[[i]] <- tmp_gene_id
  majiq_gene_symbol[[i]] <- tmp_gene_symbol
  majiq_lsv_id[[i]] <- tmp_lsv_id
  majiq_lsv_type[[i]] <- tmp_lsv_type
  majiq_p[[i]] <- tmp_p
  majiq_ds_type[[i]] <- tmp_ds_type
  majiq_dpsi[[i]] <- tmp_dpsi
}
majiq_gene_id <- flatten_chr(majiq_gene_id)
majiq_gene_symbol <- flatten_chr(majiq_gene_symbol)
majiq_lsv_id <- flatten_chr(majiq_lsv_id)
majiq_lsv_type <- flatten_chr(majiq_lsv_type)
majiq_p <- flatten_dbl(majiq_p)
majiq_ds_type <- flatten_chr(majiq_ds_type)
majiq_dpsi <- flatten_dbl(majiq_dpsi)
majiq_all_result2 <- data.frame(
  gene_id = majiq_gene_id, gene_symbol = majiq_gene_symbol, dpsi = majiq_dpsi,
  lsv_type = majiq_lsv_type, lsv_id = majiq_lsv_id, p = majiq_p, ds_type = majiq_ds_type
)
majiq_all_result2 <- majiq_all_result2 %>% group_by(lsv_id) %>%
  dplyr::mutate(tmp = dpsi[which.max(abs(dpsi))]) %>% ungroup() %>% as.data.frame()

##### rMATs #####
rmats_result_list <- vector("list", 5)
names(rmats_result_list) <- c("a3ss", "a5ss", "mxe", "ri", "se")
for (type in names(rmats_result_list)) {
  rmats_result_list[[type]] <- vroom(
    sprintf("rMATs/results/%s.MATS.JC.txt.gz", toupper(type)),
    delim = "\t",
    col_types = cols()
  ) %>% as.data.frame()
  rmats_result_list[[type]] <- mutate(rmats_result_list[[type]], ds_type = type)
}
for (type in names(rmats_result_list)) {
  novel_id <- rmats_anno_list[[type]][, 1]
  rmats_result_list[[type]] <- rmats_result_list[[type]][-which(rmats_result_list[[type]][, 1] %in% novel_id), ]
}
rmats_result <- bind_rows(rmats_result_list)

##### SplAdder #####
spladder_result_list <- vector("list", 6)
names(spladder_result_list) <- c(
  "alt_3prime", "alt_5prime", "mutex_exons",
  "intron_retention", "exon_skip", "mult_exon_skip"
)
for (type in names(spladder_result_list)) {
  spladder_result_list[[type]] <- vroom(
    sprintf(
      "spladder/non_novel/single_graphs/testing_dux4_vs_others_dux4_no_cap_exp_outlier/test_results_C3_%s.tsv.gz",
      type
    ), delim = "\t", col_types = cols()
  ) %>% as.data.frame()
  spladder_result_list[[type]] <- mutate(spladder_result_list[[type]], ds_type = type)
}
spladder_result <- bind_rows(spladder_result_list) %>% filter(!is.na(p_val_adj))

##### SUPPA #####
##### suppa mean result #####
suppa_mean_result_list <- vector("list", 7)
names(suppa_mean_result_list) <- c("a3", "a5", "mx", "ri", "se", "al", "af")
for (type in names(suppa_mean_result_list)) {
  suppa_mean_result_list[[type]] <- vroom(
    sprintf("SUPPA/diff/diff.%s.mean.dpsi.gz", toupper(type)), delim = "\t", col_types = cols()
  ) %>% as.data.frame()
  suppa_mean_result_list[[type]] <- mutate(
    suppa_mean_result_list[[type]], ds_type = type,
  ) %>% group_by(...1) %>% mutate(gene = strsplit(...1, ";")[[1]][1]) %>% ungroup() %>% as.data.frame()
}
suppa_mean_result <- bind_rows(suppa_mean_result_list) %>%
  dplyr::filter(!is.na(`dux4-others_dux4_p-val`), !is.na(`dux4-others_dux4_dPSI`)) %>%
  mutate(
    ds_type = case_when(
    `ds_type` == "a3" ~ "a3ss",
    `ds_type` == "a5" ~ "a5ss",
    `ds_type` == "mx" ~ "mxe",
    `ds_type` == "ri" ~ "ri",
    `ds_type` == "se" ~ "se",
    `ds_type` == "af" ~ "altstart",
    `ds_type` == "al" ~ "altend"
    )
  )
##### suppa mean no gene correction result #####
suppa_mean_nogc_result_list <- vector("list", 7)
names(suppa_mean_nogc_result_list) <- c("a3", "a5", "mx", "ri", "se", "al", "af")
for (type in names(suppa_mean_nogc_result_list)) {
  suppa_mean_nogc_result_list[[type]] <- vroom(
    sprintf("SUPPA/diff/nogc/diff.%s.mean.dpsi.gz", toupper(type)), delim = "\t", col_types = cols()
  ) %>% as.data.frame()
  suppa_mean_nogc_result_list[[type]] <- mutate(
    suppa_mean_nogc_result_list[[type]], ds_type = type,
  ) %>% group_by(...1) %>% mutate(gene = strsplit(...1, ";")[[1]][1]) %>% ungroup() %>% as.data.frame()
}
suppa_mean_nogc_result <- bind_rows(suppa_mean_nogc_result_list) %>%
  dplyr::filter(!is.na(`dux4-others_dux4_p-val`), !is.na(`dux4-others_dux4_dPSI`)) %>%
  mutate(
    ds_type = case_when(
      `ds_type` == "a3" ~ "a3ss",
      `ds_type` == "a5" ~ "a5ss",
      `ds_type` == "mx" ~ "mxe",
      `ds_type` == "ri" ~ "ri",
      `ds_type` == "se" ~ "se",
      `ds_type` == "af" ~ "altstart",
      `ds_type` == "al" ~ "altend"
    )
  )
##### suppa median result #####
suppa_median_result_list <- vector("list", 7)
names(suppa_median_result_list) <- c("a3", "a5", "mx", "ri", "se", "al", "af")
for (type in names(suppa_median_result_list)) {
  suppa_median_result_list[[type]] <- vroom(
    sprintf("SUPPA/diff/diff.%s.median.dpsi.gz", toupper(type)), delim = "\t", col_types = cols()
  ) %>% as.data.frame()
  suppa_median_result_list[[type]] <- mutate(
    suppa_median_result_list[[type]], ds_type = type,
  ) %>% group_by(...1) %>% mutate(gene = strsplit(...1, ";")[[1]][1]) %>% ungroup() %>% as.data.frame()
}
suppa_median_result <- bind_rows(suppa_median_result_list) %>%
  dplyr::filter(!is.na(`dux4-others_dux4_p-val`), !is.na(`dux4-others_dux4_dPSI`)) %>%
  mutate(
    ds_type = case_when(
      `ds_type` == "a3" ~ "a3ss",
      `ds_type` == "a5" ~ "a5ss",
      `ds_type` == "mx" ~ "mxe",
      `ds_type` == "ri" ~ "ri",
      `ds_type` == "se" ~ "se",
      `ds_type` == "af" ~ "altstart",
      `ds_type` == "al" ~ "altend"
    )
  )
##### suppa median no gene correction result #####
suppa_median_nogc_result_list <- vector("list", 7)
names(suppa_median_nogc_result_list) <- c("a3", "a5", "mx", "ri", "se", "al", "af")
for (type in names(suppa_median_nogc_result_list)) {
  suppa_median_nogc_result_list[[type]] <- vroom(
    sprintf("SUPPA/diff/nogc/diff.%s.median.dpsi.gz", toupper(type)), delim = "\t", col_types = cols()
  ) %>% as.data.frame()
  suppa_median_nogc_result_list[[type]] <- mutate(
    suppa_median_nogc_result_list[[type]], ds_type = type,
  ) %>% group_by(...1) %>% mutate(gene = strsplit(...1, ";")[[1]][1]) %>% ungroup() %>% as.data.frame()
}
suppa_median_nogc_result <- bind_rows(suppa_median_nogc_result_list) %>%
  dplyr::filter(!is.na(`dux4-others_dux4_p-val`), !is.na(`dux4-others_dux4_dPSI`)) %>%
  mutate(
    ds_type = case_when(
      `ds_type` == "a3" ~ "a3ss",
      `ds_type` == "a5" ~ "a5ss",
      `ds_type` == "mx" ~ "mxe",
      `ds_type` == "ri" ~ "ri",
      `ds_type` == "se" ~ "se",
      `ds_type` == "af" ~ "altstart",
      `ds_type` == "al" ~ "altend"
    )
  )

##### Whippet #####
whippet_result <- vroom(
  "Whippet/results/output.diff", delim = "\t", col_types = cols()
) %>% as.data.frame()
whippet_result <- whippet_result %>% mutate(
  ds_type = case_when(
  `Type` == "AA" ~ "a3ss",
  `Type` == "AD" ~ "a5ss",
  `Type` == "RI" ~ "ri",
  `Type` == "AF" ~ "altstart",
  `Type` == "AL" ~ "altend",
  `Type` == "CE" ~ "ce",
  `Type` == "TE" ~ "te",
  `Type` == "TS" ~ "ts"
  )
)

##### leafcutter #####
## assign dpsi to result data
leafcutter_size <- vroom(
  "leafcutter/ds/leafcutter_effect_sizes.txt.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
leafcutter_size <- leafcutter_size %>% group_by(intron) %>% 
  mutate(clus = str_split(intron, ":")[[1]][4]) %>% ungroup() %>% as.data.frame()
leafcutter_size <- leafcutter_size %>% group_by(clus) %>% mutate(
  tmp = deltapsi[which.max(abs(deltapsi))]
) %>% ungroup() %>% as.data.frame()
lc_size <- leafcutter_size[match(
  unique(leafcutter_size$clus), leafcutter_size$clus
), c("clus", "deltapsi")]
leafcutter_result <- left_join(leafcutter_result, lc_size, by = "clus")

##### leafcutter_regtools #####
leafcutter_regtools_size <- vroom(
  "leafcutter/ds_regtools/leafcutter_effect_sizes.txt.gz", delim = "\t", col_types = cols()
) %>% as.data.frame()
leafcutter_regtools_size <- leafcutter_regtools_size %>% group_by(intron) %>%
  mutate(clus = str_split(intron, ":")[[1]][4]) %>% ungroup() %>% as.data.frame()
leafcutter_regtools_size <- leafcutter_regtools_size %>% group_by(clus) %>% mutate(
  tmp = deltapsi[which.max(abs(deltapsi))]
) %>% ungroup() %>% as.data.frame()
lc_rt_size <- leafcutter_regtools_size[match(
  unique(leafcutter_regtools_size$clus), leafcutter_regtools_size$clus
), c("clus", "deltapsi")]
leafcutter_regtools_result <- left_join(leafcutter_regtools_result, lc_rt_size, by = "clus")

##### check column names #####
tx2gene <- vroom(
  "~/doc/tx2gene/tx2gene_v29.txt.gz", delim = "\t", col_names = F, col_types = cols()
) %>% as.data.frame()
colnames(tx2gene) <- c("transcript_id", "gene_id", "gene_symbol")
tx2gene <- tx2gene[match(unique(tx2gene$gene_id), tx2gene$gene_id), c("gene_id", "gene_symbol")]
rename(cash_result, gene_symbol = gene, dpsi = delta_PSI, p.adj = FDR) -> cash_result
rename(leafcutter_result, gene_symbol = genes, dpsi = deltapsi, p.adj = p.adjust) -> leafcutter_result
rename(leafcutter_regtools_result, gene_symbol = genes, dpsi = deltapsi, p.adj = p.adjust) -> leafcutter_regtools_result
rename(rmats_result, gene_symbol = geneSymbol, dpsi = IncLevelDifference, p.adj = FDR) -> rmats_result
rename(spladder_result, gene_id = gene, p.adj = p_val_adj) -> spladder_result
left_join(spladder_result, tx2gene, by = "gene_id") -> spladder_result
rename(suppa_median_result, gene_id = gene, dpsi_median = `dux4-others_dux4_dPSI`, p.adj_median = `dux4-others_dux4_p-val`) -> suppa_median_result
left_join(suppa_median_result, tx2gene, by = "gene_id") -> suppa_median_result
rename(suppa_median_nogc_result, gene_id = gene, dpsi_median_nogc = `dux4-others_dux4_dPSI`, p_median_nogc = `dux4-others_dux4_p-val`) -> suppa_median_nogc_result
left_join(suppa_median_nogc_result, tx2gene, by = "gene_id") -> suppa_median_nogc_result
rename(suppa_mean_result, gene_id = gene, dpsi_mean = `dux4-others_dux4_dPSI`, p.adj_mean = `dux4-others_dux4_p-val`) -> suppa_mean_result
left_join(suppa_mean_result, tx2gene, by = "gene_id") -> suppa_mean_result
rename(suppa_mean_nogc_result, gene_id = gene, dpsi_mean_nogc = `dux4-others_dux4_dPSI`, p_mean_nogc = `dux4-others_dux4_p-val`) -> suppa_mean_nogc_result
left_join(suppa_mean_nogc_result, tx2gene, by = "gene_id") -> suppa_mean_nogc_result
suppa_result <- left_join(
  suppa_median_result, suppa_mean_result,
  by = c("...1" = "...1", "ds_type" = "ds_type", "gene_id" = "gene_id", "gene_symbol" = "gene_symbol")
)
suppa_nogc_result <- left_join(
  suppa_median_nogc_result, suppa_mean_nogc_result,
  by = c("...1" = "...1", "ds_type" = "ds_type", "gene_id" = "gene_id", "gene_symbol" = "gene_symbol")
)
rename(whippet_result, gene_id = Gene, dpsi = DeltaPsi) -> whippet_result
left_join(whippet_result, tx2gene, by = "gene_id") -> whippet_result

##### extract ds events whose adjusted p value or probability is below 0.05 #####
cash_sig <- filter(cash_result, p.adj <= 0.05, abs(dpsi) >= 0.05)
leafcutter_sig <- filter(leafcutter_result, p.adj <= 0.05, abs(dpsi) >= 0.05)
leafcutter_regtools_sig <- filter(leafcutter_regtools_result, p.adj <= 0.05, abs(dpsi) >= 0.05)
majiq_sig <- filter(majiq_result2, p <= 0.05)
majiq_sig <- majiq_sig[match(unique(majiq_sig$lsv_id), majiq_sig$lsv_id), ]
majiq_sig$ds_type[grep("i", majiq_sig$lsv_type)] <- "ri"
rmats_sig <- filter(rmats_result, p.adj <= 0.05, abs(dpsi) >= 0.05)
spladder_sig <- filter(spladder_result, p.adj <= 0.05, abs(log2FC_event_count) >= 1)
suppa_sig <- filter(suppa_result, (p.adj_median <= 0.05 & abs(dpsi_median) >= 0.05) | (p.adj_mean <= 0.05 & abs(dpsi_mean) >= 0.05))
whippet_sig <- filter(whippet_result, Probability >= 0.9, abs(dpsi) >= 0.05)

##### create a list incorporating ds results #####
ds_list <- list(
  cash_sig, leafcutter_sig, leafcutter_regtools_sig, majiq_sig,
  rmats_sig, spladder_sig, suppa_sig, whippet_sig
)
names(ds_list) <- c(
  "CASH", "leafcutter", "leafcutter_regtools",
  "MAJIQ", "rMATs", "SplAdder", "SUPPA", "Whippet"
)
##### create summary data #####
ds_gene_all <- unique(
  c(
  unique(cash_sig$gene_symbol), unique(leafcutter_sig$gene_symbol),
  unique(leafcutter_regtools_sig$gene_symbol), unique(majiq_sig$gene_symbol),
  unique(rmats_sig$gene_symbol), unique(spladder_sig$gene_symbol),
  unique(suppa_sig$gene_symbol), unique(whippet_sig$gene_symbol)
  )
)
dat_ds_gene <- data.frame(matrix(0, nrow = length(ds_gene_all), ncol = length(softwares)))
colnames(dat_ds_gene) <- softwares
dat_ds_gene$ds_gene <- ds_gene_all
dat_ds_gene <- dplyr::select(dat_ds_gene, ds_gene, everything())
for (sw in names(ds_list)) {
  ds_tbl <- table(ds_list[[sw]]$gene_symbol) %>%
    {
      .[match(dat_ds_gene$ds_gene, names(.))] -> .
      .
    } %>%
    as.numeric()
  dat_ds_gene[, sw] <- ds_tbl
}
for (i in 2:ncol(dat_ds_gene)) {
  dat_ds_gene[is.na(dat_ds_gene[, i]), i] <- 0
  dat_ds_gene[dat_ds_gene[, i] != 0, i] <- 1
}
##### venn plot #####
color_sw <- color_lib[c(6, 7, 1, 3, 4, 9, 13, 29)]
names(color_sw) <- c(
  "leafcutter", "leafcutter_regtools", "CASH", "MAJIQ", "rMATs", "SplAdder", "SUPPA", "Whippet"
)
ggvenn(dat_ds_gene[, 2:ncol(dat_ds_gene)], alpha = 0.7) + guides(fill = F) +
  scale_fill_manual(values = color_sw) +
  theme(
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) -> p1
ggsave("figs/venn_ds_genes.pdf", p1, width = 5, height = 5, scale = 2)
##### upset plot #####
upset(
  dat_ds_gene, order.by = "freq", matrix.color = "#87CEFA", nsets = 8, nintersects = NA,
  main.bar.color = "#424242", sets.bar.color = color_sw[c(2, 1, 8, 4, 5, 6, 7, 3)],
  mainbar.y.label = "AS gene intersection", sets.x.label = "number of AS gene",
  point.size = 1.25, line.size = 0.6, mb.ratio = c(0.7, 0.3),
  att.pos = "top", shade.color = "black",
  shade.alpha = 0.1, matrix.dot.alpha = 0.5, text.scale = 0.8
) -> p2
pdf("figs/upset_ds_genes.pdf", width = 10, height = 7)
p2
dev.off()
p3 <- as.ggplot(p2)
p4 <- p3 + geom_subview(
  subview = p1 + theme_void(), x = .7, y = .7, w = .5, h = .5,
)
ggsave("figs/upset_venn_ds_gene.pdf", p4, width = 10, height = 7, scale = 1.5)
##### correlation plot #####
cor_jaccrd <- 1 - as.matrix(ecodist::distance(t(dat_ds_gene[, -1]), method = "jaccard"))
pdf("figs/cor_jaccard_ds_gene.pdf", width = 7, height = 7)
corrplot(as.matrix(cor_jaccrd), tl.col = "black", order = "AOE", addCoef.col = "black", number.digits = 3)
dev.off()

##### enrich GO terms and plot heatmap of p values #####
genesForGo <- function(sig, nbr) {
  tmp_sig <- sig[!is.na(sig$p.adj), ]
  tmp_sig <- tmp_sig %>% group_by(gene_symbol) %>% mutate(
    tmp = p.adj[which.min(p.adj)]
  ) %>% ungroup() %>% as.data.frame()
  tmp_sig2 <- tmp_sig[match(unique(tmp_sig$gene_symbol), tmp_sig$gene_symbol), c("gene_symbol", "tmp")]
  genes_for_go <- tmp_sig2$gene_symbol[head(order(tmp_sig2$tmp), nbr)]
  genes_for_go
}
##### use adjusted or unadjusted p to extract top 500 genes #####
genes_go_rmats <- genesForGo(dplyr::rename(rmats_result, pp = p.adj, p.adj = PValue), 500)
genes_go_leafcutter <- genesForGo(dplyr::rename(leafcutter_result, pp = p.adj, p.adj = p), 500)
genes_go_leafcutter_regtools <- genesForGo(dplyr::rename(leafcutter_regtools_result, pp = p.adj, p.adj = p), 500)
genes_go_majiq <- genesForGo(dplyr::rename(majiq_sig, p.adj = p), 500)
genes_go_whippet <- genesForGo(mutate(whippet_sig, p.adj = 1 - Probability), 500)
genes_go_spladder <- genesForGo(dplyr::rename(spladder_result, pp = p.adj, p.adj = p_val), 500)
##### enrichment #####
enrichGO(genes_go_rmats, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL") -> go_rmats
enrichGO(genes_go_leafcutter, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL") -> go_leafcutter
enrichGO(genes_go_leafcutter_regtools, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL") -> go_leafcutter_regtools
enrichGO(genes_go_majiq, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL") -> go_majiq
enrichGO(genes_go_whippet, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL") -> go_whippet
enrichGO(genes_go_spladder, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL") -> go_spladder
##### take all significant go terms into account #####
go_term_sig_rmats <- dplyr::filter(go_rmats@result, p.adjust <= 0.05)#go_rmats@result[head(order(go_rmats@result$p.adjust), 25), ]
go_term_sig_rmats <- dplyr::select(go_term_sig_rmats, ID, Description, p.adjust) %>% 
  dplyr::rename(rMATs = p.adjust) %>% mutate(GO_term = paste0(ID, "-", Description)) %>%
  dplyr::select(GO_term, rMATs)
go_term_sig_leafcutter <- dplyr::filter(go_leafcutter@result, p.adjust <= 0.05)#go_leafcutter@result[head(order(go_leafcutter@result$p.adjust), 25), ]
go_term_sig_leafcutter <- dplyr::select(go_term_sig_leafcutter, ID, Description, p.adjust) %>% 
  dplyr::rename(leafcutter = p.adjust) %>% mutate(GO_term = paste0(ID, "-", Description)) %>%
  dplyr::select(GO_term, leafcutter)
go_term_sig_leafcutter_regtools <- dplyr::filter(go_leafcutter_regtools@result, p.adjust <= 0.05)#go_leafcutter_regtools@result[head(order(go_leafcutter_regtools@result$p.adjust), 25), ]
go_term_sig_leafcutter_regtools <- dplyr::select(go_term_sig_leafcutter_regtools, ID, Description, p.adjust) %>% 
  dplyr::rename(leafcutter_regtools = p.adjust) %>% mutate(GO_term = paste0(ID, "-", Description)) %>%
  dplyr::select(GO_term, leafcutter_regtools)
go_term_sig_majiq <- dplyr::filter(go_majiq@result, p.adjust <= 0.05)#go_majiq@result[head(order(go_majiq@result$p.adjust), 25), ]
go_term_sig_majiq <- dplyr::select(go_term_sig_majiq, ID, Description, p.adjust) %>% 
  dplyr::rename(MAJIQ = p.adjust) %>% mutate(GO_term = paste0(ID, "-", Description)) %>%
  dplyr::select(GO_term, MAJIQ)
go_term_sig_whippet <- dplyr::filter(go_whippet@result, p.adjust <= 0.05)#go_whippet@result[head(order(go_whippet@result$p.adjust), 25), ]
go_term_sig_whippet <- dplyr::select(go_term_sig_whippet, ID, Description, p.adjust) %>% 
  dplyr::rename(Whippet = p.adjust) %>% mutate(GO_term = paste0(ID, "-", Description)) %>%
  dplyr::select(GO_term, Whippet)
go_term_sig_spladder <- dplyr::filter(go_spladder@result, p.adjust <= 0.05)#go_spladder@result[head(order(go_spladder@result$p.adjust), 25), ]
go_term_sig_spladder <- dplyr::select(go_term_sig_spladder, ID, Description, p.adjust) %>% 
  dplyr::rename(SplAdder = p.adjust) %>% mutate(GO_term = paste0(ID, "-", Description)) %>%
  dplyr::select(GO_term, SplAdder)
##### merge p values #####
go_merge <- go_term_sig_leafcutter
for (go in list(
  go_term_sig_leafcutter_regtools, go_term_sig_majiq,
  go_term_sig_rmats, go_term_sig_spladder, go_term_sig_whippet
)) {
  go_merge <- full_join(go_merge, go, by = "GO_term")
}
rownames(go_merge) <- go_merge$GO_term
go_merge <- go_merge[, -1]
for (i in 1:ncol(go_merge)) {
  go_merge[is.na(go_merge[, i]), i] <- 1
  go_merge[, i] <- -log10(go_merge[, i])
}
##### complexheatmap #####
pdf("figs/heatmap_GO.pdf", width = 10, height = 10)
Heatmap(
  as.matrix(go_merge), name = "-log10 p",
  col = colorRampPalette(c(rep("white", 1), c(rep("blue", 1))))(500),
  row_names_gp = gpar(fontsize = 0.5),
  column_names_gp = gpar(fontsize = 10),
  row_dend_reorder = T
)
dev.off()

##### plot median psi value #####
psiMedian <- function(dat, anno_field = NULL, attach_ori = F) {
  library(matrixStats)
  library(magrittr)
  if (is.null(anno_field)) {
    psi <- dat
  } else {
    psi <- dat[, !(1:ncol(dat) %in% anno_field), drop = F]
    anno <- dat[, anno_field, drop = F]
  }
  
  tmp_median <- rowMedians(as.matrix(psi), na.rm = T) %>% as.data.frame() %>%
    magrittr::set_names("psi_median")
  
  if (attach_ori) {
    return(cbind(dat, tmp_median))
  } else {
    if (is.null(anno_field)) {
      return(tmp_median)
    } else {
      return(cbind(anno, tmp_median))
    }
  }
}

miso_psi_median <- psiMedian(miso_psi, c(1:3, 104:107))
rmats_psi_median <- psiMedian(rmats_psi, c(1:3, 104))
spladder_psi_median <- psiMedian(spladder_psi, c(1:25, 126:127))
suppa_psi_median <- psiMedian(suppa_psi, c(1:2, 103:104))
whippet_psi_median <- psiMedian(whippet_psi, c(1:5, 106:107))

psi_median_list <- list(
  miso_psi_median[, c("ds_type", "psi_median")], rmats_psi_median[, c("ds_type", "psi_median")],
  spladder_psi_median[, c("ds_type", "psi_median")], suppa_psi_median[, c("ds_type", "psi_median")],
  whippet_psi_median[, c("ds_type", "psi_median")], majiq_psi2[, c("ds_type", "psi_median")]
)
names(psi_median_list) <- c("MISO", "rMATs", "SplAdder", "SUPPA", "Whippet", "MAJIQ")

dat_psi_median <- data.frame(
  software = c(
    rep("MISO", length(unique(miso_psi_median$ds_type))),
    rep("rMATs", length(unique(rmats_psi_median$ds_type))),
    rep("SplAdder", length(unique(spladder_psi_median$ds_type))),
    rep("SUPPA", length(unique(suppa_psi_median$ds_type))),
    rep("Whippet", length(unique(whippet_psi_median$ds_type))),
    rep("MAJIQ", length(unique(majiq_psi$ds_type)))
  ), # length(c(rep("MAJIQ", length(unique(majiq_psi$ds_type))), rep("MISO", length(unique(miso_psi_median$ds_type))),rep("rMATs", length(unique(rmats_psi_median$ds_type))),rep("SplAdder", length(unique(spladder_psi_median$ds_type))),rep("SUPPA", length(unique(suppa_psi_median$ds_type))),rep("Whippet", length(unique(whippet_psi_median$ds_type)))))
  ds_type = c(
    unique(miso_psi_median$ds_type), unique(rmats_psi_median$ds_type),
    unique(spladder_psi_median$ds_type), unique(suppa_psi_median$ds_type),
    unique(whippet_psi_median$ds_type), unique(majiq_psi$ds_type)
  ), 
  low = rep(0, 37), variable = rep(0, 37), high = rep(0, 37)
)

for (sw in names(psi_median_list)) {
  tmp_sw <- sw
  tmp_psi_median <- psi_median_list[[sw]]
  tmp_types <- unique(tmp_psi_median$ds_type)
  for (type in tmp_types) {
    tmp_psi_median_type <- tmp_psi_median[tmp_psi_median$ds_type == type,]
    tmp_low <- length(which(tmp_psi_median_type$psi_median <= 0.3))
    tmp_variable <- length(which(
      tmp_psi_median_type$psi_median > 0.3 & tmp_psi_median_type$psi_median < 0.7
      ))
    tmp_high <- length(which(tmp_psi_median_type$psi_median >= 0.7))
    dat_psi_median$low[dat_psi_median$software == sw & dat_psi_median$ds_type == type] <- tmp_low
    dat_psi_median$variable[dat_psi_median$software == sw & dat_psi_median$ds_type == type] <- tmp_variable
    dat_psi_median$high[dat_psi_median$software == sw & dat_psi_median$ds_type == type] <- tmp_high
  }
}

dat_psi_median_m <- melt(dat_psi_median)
dat_psi_median_m$ds_type <- factor(
  dat_psi_median_m$ds_type, levels = c(
    "a3ss", "a5ss", "ri", "se", "mxe", "altstart", "altend", "te", "ts", "ce", "sme"
  )
)
ggplot(dat_psi_median_m, mapping = aes(x = ds_type, y = value, fill = variable)) +
  geom_bar(
    size = 0.1, col = "black", stat = "identity",
    position = position_fill(0.5), width = 0.7
  ) +
  scale_fill_manual(values = c(
    "high" = "#117632", "low" = "#DDCC77", "variable" = "#999932"
  )) + facet_grid(software ~ ., scales = "free") + 
  xlab("AS type") + theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  ) -> p
ggsave("figs/proportion_psi.pdf", p, width = 7, height = 8)


##### plot deltapsi for all genes #####
deltapsi_list <- list(
  cash_result[, c("dpsi", "ds_type")], leafcutter_size[, c("clus", "tmp")],
  leafcutter_regtools_size[, c("clus", "tmp")], majiq_result2_all[, c("lsv_id", "tmp", "ds_type")],
  rmats_result[, c("dpsi", "ds_type")], suppa_median_result[, c("dpsi_median", "ds_type")],
  whippet_result[, c("dpsi", "ds_type")]
)
names(deltapsi_list) <- c(
  "CASH", "leafcutter", "leafcutter_regtools", "MAJIQ", "rMATs", "SUPPA", "Whippet"
)

deltapsi_list[["leafcutter"]] <- dplyr::rename(deltapsi_list[["leafcutter"]], dpsi = tmp)
deltapsi_list[["leafcutter"]]$ds_type <- "others"
deltapsi_list[["leafcutter"]] <- deltapsi_list[["leafcutter"]][match(unique(deltapsi_list[["leafcutter"]]$clus), deltapsi_list[["leafcutter"]]$clus), c("dpsi", "ds_type")]

deltapsi_list[["leafcutter_regtools"]] <- dplyr::rename(deltapsi_list[["leafcutter_regtools"]], dpsi = tmp)
deltapsi_list[["leafcutter_regtools"]]$ds_type <- "others"
deltapsi_list[["leafcutter_regtools"]] <- deltapsi_list[["leafcutter_regtools"]][match(unique(deltapsi_list[["leafcutter_regtools"]]$clus), deltapsi_list[["leafcutter_regtools"]]$clus), c("dpsi", "ds_type")]

deltapsi_list[["SUPPA"]] <- dplyr::rename(deltapsi_list[["SUPPA"]], dpsi = dpsi_median)
deltapsi_list[["MAJIQ"]] <- dplyr::rename(deltapsi_list[["MAJIQ"]], dpsi = tmp)
deltapsi_list[["MAJIQ"]] <- deltapsi_list[["MAJIQ"]][match(unique(deltapsi_list[["MAJIQ"]]$lsv_id), deltapsi_list[["MAJIQ"]]$lsv_id), c("dpsi", "ds_type")]

for (sw in names(deltapsi_list)) {
  deltapsi_list[[sw]]$softwares <- sw
}

dat_deltapsi <- bind_rows(deltapsi_list)
dat_deltapsi$softwares <- factor(dat_deltapsi$softwares, levels = c(
  "Whippet", "CASH", "rMATs", "SUPPA", "MAJIQ", "leafcutter", "leafcutter_regtools"
))
dat_deltapsi$ds_type <- factor(dat_deltapsi$ds_type, levels = c(
  "a3ss", "a5ss", "ri", "se", "mxe", "altstart", "altend", "te", "ts", "ce", "sme", "others"
))
ggplot(dat_deltapsi, aes(x = ds_type, y = dpsi, fill = ds_type)) +
  geom_boxplot(outlier.size = .5, outlier.alpha = .2, outlier.stroke = .1) +
  scale_fill_manual(values = colors_as) +
  facet_grid(softwares ~ ., scales = "free") +
  theme_bw() +
  theme(
    #panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  ) -> p
ggsave("figs/boxplot_dpsi_all.pdf", p, width = 7, height = 8)

##### plot deltapsi boxplot for ds genes #####
suppa_sig <- dplyr::mutate(suppa_sig, tmp = case_when(
  p.adj_median <= p.adj_mean ~ dpsi_median,
  p.adj_median > p.adj_mean ~ dpsi_mean
))

deltapsi_ds_list <- list(
  cash_sig[, c("dpsi", "ds_type")], leafcutter_sig[, "dpsi", drop = F],
  leafcutter_regtools_sig[, "dpsi", drop = F], majiq_sig[, c("tmp", "ds_type")],
  rmats_sig[, c("dpsi", "ds_type")], suppa_sig[, c("tmp", "ds_type")],
  whippet_sig[, c("dpsi", "ds_type")]
)
names(deltapsi_ds_list) <- c(
  "CASH", "leafcutter", "leafcutter_regtools", "MAJIQ", "rMATs", "SUPPA", "Whippet"
)
deltapsi_ds_list[["leafcutter"]]$ds_type <- "others"
deltapsi_ds_list[["leafcutter_regtools"]]$ds_type <- "others"
deltapsi_ds_list[["MAJIQ"]] <- dplyr::rename(deltapsi_ds_list[["MAJIQ"]], dpsi = tmp)
deltapsi_ds_list[["SUPPA"]] <- dplyr::rename(deltapsi_ds_list[["SUPPA"]], dpsi = tmp)

for (sw in names(deltapsi_ds_list)) {
  deltapsi_ds_list[[sw]]$softwares <- sw
}


dat_deltapsi_ds <- bind_rows(deltapsi_ds_list)
dat_deltapsi_ds$softwares <- factor(dat_deltapsi_ds$softwares, levels = c(
  "Whippet", "CASH", "rMATs", "SUPPA", "MAJIQ", "leafcutter", "leafcutter_regtools"
))
dat_deltapsi_ds$ds_type <- factor(dat_deltapsi_ds$ds_type, levels = c(
  "a3ss", "a5ss", "ri", "se", "mxe", "altstart", "altend", "te", "ts", "ce", "sme", "others"
))
ggplot(dat_deltapsi_ds, aes(x = ds_type, y = dpsi, fill = ds_type)) +
  geom_boxplot(outlier.size = .5, outlier.alpha = .2, outlier.stroke = .1) +
  geom_jitter(size = .2, alpha = .5, position = position_jitterdodge(jitter.width = 2)) +
  scale_fill_manual(values = colors_as) +
  facet_grid(softwares ~ ., scales = "free") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  ) -> p
ggsave("figs/boxplot_dpsi_ds.pdf", p, width = 7, height = 8)


