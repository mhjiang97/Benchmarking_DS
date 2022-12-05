pkgs <- c(
  "dplyr", "vroom", "glue", "purrr", "stringr", "tibble", "tidyr",
  "ggplot2", "RColorBrewer"
)
for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = T))


myVroom <- function(file, na_append = NULL, delim = "\t", comment = "", col_names = T) {
  na_origin <- c("", "NA")
  myna <- c(na_origin, na_append)

  n_col <- read.table(
    file, header = T, sep = "\t", na.strings = c("NA", na_append), nrows = 10
  ) |>
    ncol()

  dat <- vroom::vroom(
    file, delim = delim, col_types = paste0(rep("c", n_col), collapse = ""),
    na = myna, comment = comment, col_names = col_names
  )

  dat
}

findReadFiles <- function(types_all, types_event, dir, patterns) {
  et <- intersect(types_event, types_all)
  if ("all" %in% types_event) et <- types_all

  files_read <- vector("character", length(et)) |> setNames(et)
  for (e in names(files_read)) {
    files_read[e] <- list.files(path = dir, pattern = patterns[e], full.names = T)
  }

  files_read
}

markPsi <- function(result_psi, truths_psi, as_type = T) {
  list_delta_psi <- vector("list", nrow(result_psi))
  for (i in 1:nrow(result_psi)) {
    if (is.na(result_psi$idx[i]) || is.na(result_psi$sample1[i])) next

    idxs <- as.numeric(unlist(str_split(result_psi$idx[i], "\\,")))
    tmp <- data.frame(row.names = idxs)
    for (idx in idxs) {
      psi_t <- c(truths_psi$psi[idx], truths_psi$psi2[idx])
      psi_d <- as.numeric(result_psi$sample1[i]) - psi_t
      psi_idx <- which.min(abs(psi_d))
      tmp[as.character(idx), 1] <- psi_d[psi_idx]
      tmp[as.character(idx), 2] <- psi_t[psi_idx]
      tmp[as.character(idx), 3] <- as.numeric(result_psi$sample1[i])
      tmp[as.character(idx), 4] <- idx
      if (as_type) {
        tmp[as.character(idx), 5] <- result_psi$as_type[i]
      }
    }
    idx_tmp <- which.min(abs(na.omit(tmp[, 1])))
    tmp2 <- tmp[idx_tmp, ]
    list_delta_psi[[i]] <- tmp2
  }
  delta_psi <- bind_rows(list_delta_psi)

  if (as_type) {
    delta_psi |> setNames(c("delta", "truth", "original", "truths_idx", "as_type"))
  } else {
    delta_psi |> setNames(c("delta", "truth", "original", "truths_idx"))
  }
}


psi_rmats <- function(dir, read_type = "JC", event_type = "all") {
  types_all <- c("SE", "MXE", "A3SS", "A5SS", "RI")
  patterns <- glue::glue("{types_all}.*.JC.txt") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )

  list_rmats <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_rmats)) {
    f <- files_read[e]
    tmp <- f |>
      myVroom(na_append = c("nan", "NaN")) |>
      dplyr::mutate(as_type = e) |>
      tibble::remove_rownames() |>
      tibble::column_to_rownames("ID...1") |>
      dplyr::select(!dplyr::starts_with("ID")) |>
      tibble::rownames_to_column("ID") |>
      dplyr::rename(
        padj = FDR, dpsi = IncLevelDifference,
        gene_id = GeneID, gene_symbol = geneSymbol,
        p = PValue, psis_1 = IncLevel1, psis_2 = IncLevel2
      ) |>
      tibble::as_tibble() |>
      dplyr::mutate(myID = paste0(as_type, ": ", ID))
    tmp1 <- tmp[, c("myID", "as_type", "psis_1"), drop = F] |>
      tidyr::separate(col = "psis_1", into = paste0("sample", 1:4), sep = "\\,")
    tmp2 <- tmp[, c("psis_2"), drop = F] |>
      tidyr::separate(col = "psis_2", into = paste0("sample", 5:8), sep = "\\,")
    list_rmats[[e]] <- dplyr::bind_cols(tmp1, tmp2)
  }

  list_rmats
}

psi_leafcutter <- function(file_perind) {
  f <- myVroom(file_perind, delim = " ")
  final <- apply(f[, -1], 2, function(x) {
    DOSE::parse_ratio(x)
  }) |>
    tibble::as_tibble()
  f <- f |>
    dplyr::rowwise() |>
    dplyr::mutate(
      myID = paste0(
        unlist(stringr::str_split(chrom, "\\:"))[1],
        ":",
        unlist(stringr::str_split(chrom, "\\:"))[4]
      )
    ) |>
    ungroup()
  final$myID <- f$myID

  final
}

psi_spladder <- function(dir_event, event_type = "all", confidence_level = 3) {
  types_all <- c("exon_skip", "mutex_exons", "alt_3prime", "alt_5prime", "intron_retention", "mult_exon_skip")

  patterns_event <- glue::glue(
    "merge_graphs_{types_all}_C{confidence_level}.confirmed.txt"
  ) |>
    setNames(types_all)
  files_read_event <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir_event, patterns = patterns_event
  )

  list_spladder_event <- vector("list", length(files_read_event)) |>
    setNames(names(files_read_event))
  for (e in names(list_spladder_event)) {
    f <- files_read_event[e] |>
      myVroom(na_append = c("nan", "NaN", "inf", "-inf")) |>
      dplyr::rowwise() |>
      dplyr::mutate(event_id = stringr::str_replace_all(event_id, "\\.", "\\_")) |>
      dplyr::ungroup() |>
      dplyr::mutate(myID = event_id) |>
      dplyr::mutate(as_type = e) |>
      dplyr::mutate(
        as_type = dplyr::case_when(
          as_type == "exon_skip" ~ "SE",
          as_type == "mutex_exons" ~ "MXE",
          as_type == "alt_3prime" ~ "A3SS",
          as_type == "alt_5prime" ~ "A5SS",
          as_type == "intron_retention" ~ "RI",
          as_type == "mult_exon_skip" ~ "SME"
        )
      ) |>
      dplyr::select(myID, as_type, dplyr::ends_with("psi"))

    list_spladder_event[[e]] <- f
  }

  list_spladder_event
}

psi_suppa <- function(dir, s1 = "s1", s2 = "s2", event_type = "all") {
  types_all <- c("SE", "MX", "A3", "A5", "RI", "AF", "AL")

  patterns_s1 <- glue::glue("{s1}\\.{types_all}\\.psi") |>
    as.character() |>
    setNames(types_all)
  files_read_s1 <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns_s1
  )
  patterns_s2 <- glue::glue("{s2}\\.{types_all}\\.psi") |>
    as.character() |>
    setNames(types_all)
  files_read_s2 <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns_s2
  )

  list_suppa <- vector("list", length(files_read_s1)) |>
    setNames(names(files_read_s1))
  for (e in names(list_suppa)) {
    f1 <- files_read_s1[e] |>
      read.delim() |>
      dplyr::mutate(as_type = e) |>
      tibble::rownames_to_column("event") |>
      tidyr::separate(event, c("gene_id", "coord"), sep = ";") |>
      dplyr::mutate(
        myID = coord,
        as_type = dplyr::case_when(
          as_type == "MX" ~ "MXE",
          as_type == "A3" ~ "A3SS",
          as_type == "A5" ~ "A5SS",

          as_type == "SE" ~ "SE",
          as_type == "RI" ~ "RI",
          as_type == "AF" ~ "AFE",
          as_type == "AL" ~ "ALE"
        )
      ) |>
      dplyr::select(!c(gene_id, coord)) |>
      tibble::as_tibble()

    f2 <- files_read_s2[e] |>
      read.delim() |>
      tibble::as_tibble()

    list_suppa[[e]] <- dplyr::bind_cols(f1, f2)
  }

  list_suppa
}

psi_whippet <- function(dir, samples, event_type = "all") {
  types_all <- c("CE", "AA", "AD", "RI", "TS", "TE", "AF", "AL")
  et <- intersect(event_type, types_all)
  if ("all" %in% event_type) et <- types_all

  files_psi <- glue("{dir}/{samples}.psi.gz")
  list_psi <- vector("list", length(files_psi))
  for (i in 1:length(files_psi)) {
    s <- samples[i]
    tmp <- myVroom(files_psi[i]) |>
      dplyr::rename(!!s := Psi)
    if (i == 1) {
      tmp <- tmp |>
        dplyr::filter(Type %in% et) |>
        dplyr::mutate(
          as_type = dplyr::case_when(
            Type == "CE" ~ "CE",
            Type == "AA" ~ "A3SS",
            Type == "AD" ~ "A5SS",

            Type == "RI" ~ "RI",
            Type == "TS" ~ "TS",
            Type == "TE" ~ "TE",
            Type == "AF" ~ "AFE",
            Type == "AL" ~ "ALE"
          ),
          myID = paste0(Coord, Type, Node)
        ) |>
        dplyr::select(myID, as_type, dplyr::all_of(samples[i]))
    } else {
      tmp <- tmp |>
        dplyr::filter(Type %in% et) |>
        dplyr::select(dplyr::all_of(samples[i]))
    }

    list_psi[[i]] <- tmp
  }

  bind_cols(list_psi)
}

psi_psichomics <- function(psi_file) {
  psi <- myVroom(psi_file) |>
    dplyr::rename(myID = `...1`) |>
    dplyr::mutate(
      as_type = dplyr::case_when(
        grepl("SE_", myID) ~ "SE",
        grepl("A5SS_", myID) ~ "A5SS",
        grepl("AFE_", myID) ~ "AFE",
        grepl("ALE_", myID) ~ "ALE",
        grepl("MXE_", myID) ~ "MXE",
        grepl("A3SS_", myID) ~ "A3SS"
      )
    )

  psi
}




