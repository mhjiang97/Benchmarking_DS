pkgs <- c("dplyr", "vroom", "glue", "purrr", "stringr", "tibble", "tidyr")
for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = T))

##### generic #####
querySig <- function(
  dat, p_adj = 0.05, p_value = NULL, delta_psi = 0,
  statistic_col = NULL, statistic_threshold, filter_direction = c(">", "<"),
  effect_size_col = NULL, effect_size_threshod = 0
) {
  if (!is.null(statistic_col)) {
    dat[[statistic_col]] <- as.numeric(dat[[statistic_col]])
    if (filter_direction[1] == ">") dat_sig <- dat[dat[[statistic_col]] >= statistic_threshold, ]
    if (filter_direction[1] == "<") dat_sig <- dat[dat[[statistic_col]] <= statistic_threshold, ]
  } else {
    if (is.null(p_value)) {
      dat$padj <- as.numeric(dat$padj)
      dat_sig <- dplyr::filter(dat, padj <= p_adj)
    } else {
      dat$p <- as.numeric(dat$p)
      dat_sig <- dplyr::filter(dat, p <= p_value)
    }
  }

  ## always filter according effect size after p value ##
  if (!is.null(effect_size_col)) {
    dat_sig[[effect_size_col]] <- as.numeric(dat_sig[[effect_size_col]])
    dat_sig <- dat_sig[abs(dat_sig[[effect_size_col]]) >= effect_size_threshod, ]
  } else {
    dat_sig$dpsi <- as.numeric(dat_sig$dpsi)
    dat_sig <- dplyr::filter(dat_sig, abs(dpsi) >= delta_psi)
  }

  dat_sig
}

myVroom <- function(file, na_append = NULL, delim = "\t", comment = "", col_names = T) {
  na_origin <- c("", "NA")
  myna <- c(na_origin, na_append)

  n_col <- read.table(
    file, header = T, sep = "\t", comment.char = comment, na.strings = myna, nrows = 10
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

getOs <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

o <- function(file = ".") {
  os <- Sys.info()[1]
  if (F) {
    if (dir.exists(file)) {
      stop("open directory in RStudio Server is not supported.")
    }
    rserver_ip <- getOption("rserver_ip")
    if (!is.null(rserver_ip)) {
      rserver_port <- getOption("rserver_port") %||% "8787"
      if (!startsWith(rserver_ip, "http")) {
        rserver_ip <- paste0("http://", rserver_ip)
      }
      utils::browseURL(
        paste0(
          paste(rserver_ip, rserver_port, sep = ":"),
          "/file_show?path=",
          file
        )
      )
    } else {
      file.edit <- get("file.edit")
      file.edit(file)
    }
  } else if (os == "Darwin") {
    cmd <- paste("open", file)
    system(cmd)
  } else if (os == "Linux") {
    cmd <- paste("xdg-open", file, "&")
    system(cmd)
  } else if (os == "Windows") {
    cmd <- paste("start", file)
    shell(cmd)
  }
}

show_in_excel <- function(.data) {
  f <- tempfile(fileext = '.csv')
  write.csv(.data, file = f)
  o(f)
  invisible(.data)
}

markEvent <- function(
  result, truths, cols_coord, split = T, split_by = NULL,
  numbers_equal = c("SE" = 3, "RI" = 3, "A3SS" = 3, "A5SS" = 3, "MXE" = 5),
  error_distance = 1, as = T, gene = T, gene_symbol = T,
  mark_col = "idx_truths", truths_col = c("X16", "X17"),
  verbose = F
) {
  result[[mark_col]] <- NA

  if (gene_symbol) gene_col <- "gene_symbol"
  if (!gene_symbol) gene_col <- "gene_id"

  for (i in 1:nrow(result)) {
    if (verbose) message(glue("number of result: {i}"))

    gene_result <- result[[gene_col]][i] |>
      stringr::str_split(";") |>
      unlist()
    if (as) {
      as_result <- result$as_type[i] |>
        stringr::str_split("/") |>
        unlist()
    } else {
      as_result <- "NO"
    }
    if (gene && !any(gene_result %in% truths[[gene_col]])) next
    if (as && !any(as_result %in% c("A3SS", "A5SS", "MXE", "RI", "SE"))) next

    if (!split) {
      tmp_coords_result <- result[i, cols_coord]
    } else {
      tmp_coords_result <- unlist(str_split(result[i, cols_coord], split_by))
    }
    tmp_coords_result <- tmp_coords_result |>
      as.numeric() |>
      suppressWarnings() |>
      na.omit()
    coords_result <- c(
      unlist(purrr::map(tmp_coords_result, function(x) {x - c(0:error_distance)})),
      unlist(purrr::map(tmp_coords_result, function(x) {x + c(0:error_distance)}))
      # tmp_coords_result,
      # tmp_coords_result - error_distance,
      # tmp_coords_result + error_distance
    ) |>
      as.character() |>
      unique()

    if (as) {
      number_equal <- numbers_equal[as_result[1]]
    } else {
      number_equal <- numbers_equal[1]
    }

    if (length(coords_result) < number_equal) {
      if (length(coords_result) <= 3) {
        number_equal <- length(coords_result)
      } else {
        number_equal <- length(coords_result) - 1
      }
    }


    idx <- vector("list", nrow(truths))
    for (j in 1:nrow(truths)) {
      as_truths <- truths$as_type[j]
      gene_truths <- truths[[gene_col]][j]
      # if (as && as_truths != as_result) next
      if (as && !any(as_truths %in% as_result)) next
      if (!any(gene_truths %in% gene_result)) next

      coords_truths <- truths[j, truths_col] |>
        str_split(",") |>
        unlist() |>
        unique() |>
        na.omit()

      if (length(intersect(coords_truths, coords_result)) >= number_equal) idx[[j]] <- j
    }

    if (!is.na(flatten_int(idx)[1])) result[[mark_col]][i] <- paste(flatten_int(idx), collapse = ",")
  }

  result
}

myAnnot <- function(result, gr, cols_loc, verbose = F) {
  result[["gene_id"]] <- result[["gene_symbol"]] <- NA
  
  for (i in 1:nrow(result)) {
    if (verbose) message(glue("number of result: {i}"))
    
    tmp_loc <- result[i, cols_loc] |>
      as.character() |>
      na.omit() |>
      GenomicRanges::GRanges() |>
      GenomicRanges::reduce()
    
    idx <- IRanges::findOverlaps(tmp_loc, gr, type = "within")@to |>
      unique()
    if (!is.na(idx[1])) {
      tmp_gr <- gr[idx]
      tmp_gr <- tmp_gr[!is.na(tmp_gr$hgnc_id)]
      result$gene_id[i] <- paste(tmp_gr$gene_id, collapse = ";")
      result$gene_symbol[i] <- paste(tmp_gr$gene_name, collapse = ";")
    }
  }
  
  result
}


##### rMATS #####
read_rmats <- function(dir, read_type = c("JC", "JCEC"), event_type = "all") {
  types_all <- c("SE", "MXE", "A3SS", "A5SS", "RI")
  patterns <- glue::glue("{types_all}.*.{read_type[1]}.txt") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )

  list_rmats <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_rmats)) {
    f <- files_read[e]
    list_rmats[[e]] <- f |>
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
      tibble::as_tibble() # |>
      # dplyr::filter(!is.na(PValue), !is.na(FDR))
  }

  list_rmats
}

##### CASH #####
read_cash <- function(file, event_type = "all", condition_1, condition_2) {
  types_all <- c("Cassette", "MXE", "A3SS", "A5SS", "IR", "AltStart", "AltEnd", "Cassette_multi")
  et <- intersect(event_type, types_all)
  if ("all" %in% event_type) et <- types_all

  cash <- myVroom(file, na_append = c("nan", "NaN")) |>
    dplyr::rename(
      gene_symbol = AccID, dpsi = delta_PSI, p = `P-Value`, padj = FDR
    ) |>
    dplyr::filter(SplicingType %in% et) |>
    dplyr::mutate(
      as_type = dplyr::case_when(
        SplicingType == "IR" ~ "RI",
        SplicingType == "Cassette" ~ "SE",
        SplicingType == "Cassette_multi" ~ "SME",

        SplicingType == "MXE" ~ "MXE",
        SplicingType == "A3SS" ~ "A3SS",
        SplicingType == "A5SS" ~ "A5SS",
        SplicingType == "AltStart" ~ "AFE",
        SplicingType == "AltEnd" ~ "ALE"
      )
    ) |>
    tidyr::separate(
      dplyr::all_of(glue::glue("{condition_1}_Junc_Inclusive::Exclusive")),
      c("inclusive_1", "exlclusive_1"),
      "::"
    ) |>
    tidyr::separate(
      dplyr::all_of(glue::glue("{condition_2}_Junc_Inclusive::Exclusive")),
      c("inclusive_2", "exlclusive_2"),
      "::"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      average_psi_1 = as.numeric(inclusive_1) / (as.numeric(inclusive_1) + as.numeric(exlclusive_1)),
      average_psi_2 = as.numeric(inclusive_2) / (as.numeric(inclusive_2) + as.numeric(exlclusive_2))
    ) |>
    dplyr::ungroup() # |>
    # dplyr::filter(!is.na(`P-Value`), !is.na(FDR)) |>
    # dplyr::select(!SplicingType) |>

  cash
}

##### LeafCutter #####
read_leafcutter <- function(
  file_significance, file_effect_sizes, condition_1, condition_2, flatten = F
) {
  leafcutter_size <- myVroom(file_effect_sizes) |>
    dplyr::rowwise() |>
    dplyr::mutate(clus = stringr::str_split(intron, ":")[[1]][4]) |>
    dplyr::ungroup() |>
    dplyr::group_by(clus) |>
    dplyr::mutate(
      average_psi_1 = paste(.data[[condition_1]], collapse = ","),
      average_psi_2 = paste(.data[[condition_2]], collapse = ","),
      dpsi = paste(deltapsi, collapse = ","),
      myintron = paste(intron, collapse = ",")
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct(clus, .keep_all = T) |>
    dplyr::select(
      !c(logef, intron, deltapsi) &
        !dplyr::all_of(glue::glue("{condition_1}")) &
        !dplyr::all_of(glue::glue("{condition_2}"))
    )

  leafcutter_result <- myVroom(file_significance, na_append = c("nan", "NaN")) |>
    dplyr::rowwise() |>
    dplyr::mutate(clus = stringr::str_split(cluster, ":")[[1]][2]) |>
    dplyr::ungroup() # |>
    # tidyr::separate_rows(genes, sep = ",") # |>
    # dplyr::filter(status == "Success", !is.na(p), !is.na(p.adjust), !is.na(genes))

  if (flatten) {
    leafcutter <- leafcutter_result |>
      tidyr::separate_rows(genes, sep = ",")
  }

  leafcutter <- dplyr::left_join(leafcutter_result, leafcutter_size, by = "clus") |>
    dplyr::rename(padj = p.adjust, gene_symbol = genes)

  leafcutter
}

##### LeafCutter supplement #####
read_leafcutter2 <- function(
    file_significance, file_effect_sizes, condition_1, condition_2,
    dpsi = 0.05, flatten = F
) {
  mydpsi <- dpsi
  leafcutter_size <- myVroom(file_effect_sizes) |>
    dplyr::filter(abs(as.numeric(deltapsi)) > mydpsi) |>
    dplyr::rowwise() |>
    dplyr::mutate(clus = stringr::str_split(intron, ":")[[1]][4]) |>
    dplyr::ungroup() |>
    dplyr::group_by(clus) |>
    dplyr::mutate(
      average_psi_1 = paste(.data[[condition_1]], collapse = ","),
      average_psi_2 = paste(.data[[condition_2]], collapse = ","),
      dpsi = paste(deltapsi, collapse = ","),
      myintron = paste(intron, collapse = ",")
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct(clus, .keep_all = T) |>
    dplyr::select(
      !c(logef, intron, deltapsi) &
        !dplyr::all_of(glue::glue("{condition_1}")) &
        !dplyr::all_of(glue::glue("{condition_2}"))
    )
  
  leafcutter_result <- myVroom(file_significance, na_append = c("nan", "NaN")) |>
    dplyr::rowwise() |>
    dplyr::mutate(clus = stringr::str_split(cluster, ":")[[1]][2]) |>
    dplyr::ungroup() # |>
  # tidyr::separate_rows(genes, sep = ",") # |>
  # dplyr::filter(status == "Success", !is.na(p), !is.na(p.adjust), !is.na(genes))
  
  if (flatten) {
    leafcutter <- leafcutter_result |>
      tidyr::separate_rows(genes, sep = ",")
  }
  
  leafcutter <- dplyr::inner_join(leafcutter_result, leafcutter_size, by = "clus") |>
    dplyr::rename(padj = p.adjust, gene_symbol = genes)
  
  leafcutter
}

##### SplAdder #####
read_spladder <- function(dir_test, dir_event, event_type = "all", confidence_level = 3) {
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
    f <- files_read_event[e]
    list_spladder_event[[e]] <- myVroom(f, na_append = c("nan", "NaN", "inf", "-inf")) |>
      dplyr::rowwise() |>
      dplyr::mutate(event_id = stringr::str_replace_all(event_id, "\\.", "\\_")) |>
      dplyr::ungroup()
  }


  patterns <- glue::glue(".*extended_C{confidence_level}_{types_all}.tsv") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir_test, patterns = patterns
  )

  list_spladder <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_spladder)) {
    f <- files_read[e]
    list_spladder[[e]] <- myVroom(f, na_append = c("nan", "NaN", "inf", "-inf")) |>
      dplyr::mutate(as_type = e) |>
      dplyr::rename(
        padj = p_val_adj, gene_id = gene, p = p_val
      ) |>
      dplyr::mutate(
        as_type = dplyr::case_when(
          as_type == "exon_skip" ~ "SE",
          as_type == "mutex_exons" ~ "MXE",
          as_type == "alt_3prime" ~ "A3SS",
          as_type == "alt_5prime" ~ "A5SS",
          as_type == "intron_retention" ~ "RI",
          as_type == "mult_exon_skip" ~ "SME"
        )
      ) # |>
      # dplyr::filter(!is.na(p_val), !is.na(p_val_adj))
  }

  for (e in names(list_spladder)) {
    list_spladder[[e]] <- list_spladder[[e]] |>
      dplyr::left_join(
        list_spladder_event[[e]],
        by = c("gene_id" = "gene_name", "event_id" = "event_id")
      )
  }

  list_spladder
}

##### SplAdder new #####
read_spladder_new <- function(dir_test, event_type = "all", confidence_level = 3) {
  types_all <- c("exon_skip", "mutex_exons", "alt_3prime", "alt_5prime", "intron_retention", "mult_exon_skip")

  patterns <- glue::glue(".*extended_C{confidence_level}_{types_all}.tsv") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir_test, patterns = patterns
  )

  list_spladder <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_spladder)) {
    f <- files_read[e]
    list_spladder[[e]] <- myVroom(f, na_append = c("nan", "NaN", "inf", "-inf")) |>
      dplyr::mutate(as_type = e) |>
      dplyr::rename(chr = chrm, padj = p_val_adj, p = p_val, dpsi = dPSI) |>
      dplyr::mutate(
        as_type = dplyr::case_when(
          as_type == "exon_skip" ~ "SE",
          as_type == "mutex_exons" ~ "MXE",
          as_type == "alt_3prime" ~ "A3SS",
          as_type == "alt_5prime" ~ "A5SS",
          as_type == "intron_retention" ~ "RI",
          as_type == "mult_exon_skip" ~ "SME"
        )
      ) # |>
    # dplyr::filter(!is.na(p_val), !is.na(p_val_adj))
  }

  list_spladder
}

##### MAJIQ #####
read_majiq <- function(delta_psi_file, voila_file, flatten = T, condition_1, condition_2) {
  delta_psi <- myVroom(delta_psi_file) |>
    dplyr::select(`LSV ID`, A5SS, A3SS, ES, `Num. Junctions`, `Num. Exons`)

  voila <- myVroom(voila_file, comment = "#") |>
    dplyr::left_join(delta_psi, by = c("lsv_id" = "LSV ID")) |>
    dplyr::mutate(
      as_type = dplyr::case_when(
        as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
        as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
        as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
        as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
        as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
        T ~ "OTHER"
      )
    ) |>
    dplyr::rename(
      gene_symbol = gene_name,
      dpsi = mean_dpsi_per_lsv_junction,
      p = probability_changing,
      average_psi_1 = dplyr::all_of(glue::glue("{condition_1}_mean_psi")),
      average_psi_2 = dplyr::all_of(glue::glue("{condition_2}_mean_psi"))
    )

  if (flatten) {
    voila <- voila |>
      tidyr::separate_rows(
        dpsi, p, probability_non_changing, average_psi_1, average_psi_2,
        sep = ";"
      )
  }

  voila

  # if (flatten) {
  #   majiq_gene_id <- majiq_gene_symbol <- majiq_lsv_id <- majiq_lsv_type <-
  #     majiq_p <- majiq_as_type <- majiq_dpsi <- majiq_average_psi_1 <- majiq_average_psi_2 <-
  #     majiq_junctions_coords <- majiq_exons_coords <- majiq_ir_coords <-
  #     vector("list", nrow(voila))
  #   for (i in 1:nrow(voila)) {
  #     tmp_p <- str_split(voila$probability_changing[i], ";")[[1]]
  #     tmp_gene_id <- rep(voila$gene_id[i], length(tmp_p))
  #     tmp_gene_symbol <- rep(voila$gene_name[i], length(tmp_p))
  #     tmp_lsv_id <- rep(voila$lsv_id[i], length(tmp_p))
  #     tmp_lsv_type <- rep(voila$lsv_type[i], length(tmp_p))
  #     tmp_as_type <- rep(voila$as_type[i], length(tmp_p))
  #     tmp_dpsi <- str_split(voila$mean_dpsi_per_lsv_junction[i], ";")[[1]]
  #     tmp_average_psi_1 <- str_split(voila[i, paste0(condition_1, "_mean_psi")], ";")[[1]]
  #     tmp_average_psi_2 <- str_split(voila[i, paste0(condition_2, "_mean_psi")], ";")[[1]]
  #     tmp_junctions_coords <- rep(voila$junctions_coords[i], length(tmp_p)) # str_split(voila$junctions_coords[i], ";")[[1]]
  #     tmp_exons_coords <- rep(voila$exons_coords[i], length(tmp_p)) # str_split(voila$exons_coords[i], ";")[[1]]
  #     tmp_ir_coords <- rep(voila$ir_coords[i], length(tmp_p)) # str_split(voila$ir_coords[i], ";")[[1]]
  #
  #     # if (all(is.na(tmp_ir_coords))) tmp_ir_coords <- rep(NA, length(tmp_p))
  #     # if (all(is.na(tmp_exons_coords))) tmp_exons_coords <- rep(NA, length(tmp_p))
  #
  #     majiq_gene_id[[i]] <- tmp_gene_id
  #     majiq_gene_symbol[[i]] <- tmp_gene_symbol
  #     majiq_lsv_id[[i]] <- tmp_lsv_id
  #     majiq_lsv_type[[i]] <- tmp_lsv_type
  #     majiq_p[[i]] <- tmp_p
  #     majiq_as_type[[i]] <- tmp_as_type
  #     majiq_dpsi[[i]] <- tmp_dpsi
  #     majiq_average_psi_1[[i]] <- tmp_average_psi_1
  #     majiq_average_psi_2[[i]] <- tmp_average_psi_2
  #     majiq_junctions_coords[[i]] <- tmp_junctions_coords
  #     majiq_exons_coords[[i]] <- tmp_exons_coords
  #     majiq_ir_coords[[i]] <- tmp_ir_coords
  #   }
  #   gene_id <- flatten_chr(majiq_gene_id)
  #   gene_symbol <- flatten_chr(majiq_gene_symbol)
  #   lsv_id <- flatten_chr(majiq_lsv_id)
  #   lsv_type <- flatten_chr(majiq_lsv_type)
  #   p <- flatten_chr(majiq_p)
  #   as_type <- flatten_chr(majiq_as_type)
  #   dpsi <- flatten_chr(majiq_dpsi)
  #   average_psi_1 <- flatten_chr(majiq_average_psi_1)
  #   average_psi_2 <- flatten_chr(majiq_average_psi_2)
  #   junctions_coords <- flatten_chr(majiq_junctions_coords)
  #   exons_coords <- flatten_chr(majiq_exons_coords)
  #   ir_coords <- flatten_chr(majiq_ir_coords)
  #
  #   voila2 <- data.frame(
  #     gene_id = gene_id, gene_symbol = gene_symbol, dpsi = dpsi,
  #     lsv_type = lsv_type, lsv_id = lsv_id, p = p, as_type = as_type,
  #     average_psi_1 = average_psi_1, average_psi_2 = average_psi_2,
  #     junctions_coords = junctions_coords, exons_coords = exons_coords,
  #     ir_coords = ir_coords
  #   )
  #   # voila2 <- voila2 |> group_by(lsv_id) |>
  #   #   dplyr::mutate(tmp = dpsi[which.max(abs(dpsi))]) |>
  #   #   ungroup()
  #
  # } else {voila2 <- voila}
  #
  # voila2
}

##### SGSeq #####
read_sgseq <- function(dxd = NULL, sgvc = NULL, RData = NULL, flatten = F) {
  if (!is.null(RData)) {
    suppressPackageStartupMessages(suppressMessages(load(RData)))
    table_splice <- tibble::as_tibble(S4Vectors::mcols(sgvc))
    table_dexseq <- tibble::as_tibble(dxr1)
  } else {
    table_splice <- tibble::as_tibble(S4Vectors::mcols(sgvc))
    table_dexseq <- tibble::as_tibble(DEXSeq::DEXSeqResults(dxd))
  }

  table_splice <- table_splice |>
    dplyr::mutate(
      eventID = as.character(as.numeric(eventID)),
      variantID = as.character(as.numeric(variantID))
    )
  table_dexseq <- table_dexseq |>
    dplyr::mutate(
      groupID = as.character(as.numeric(groupID)),
      featureID = as.character(as.numeric(featureID))
    )

  sgseq <- dplyr::left_join(
    table_dexseq, table_splice,
    by = c("groupID" = "eventID", "featureID" = "variantID")
  ) |>
    dplyr::rename(p = pvalue, gene_id = geneName)

  sgseq <- sgseq |>
    dplyr::rowwise() |>
    dplyr::mutate(
      featureID5p_c = paste(purrr::flatten_chr(featureID5p), collapse = ","),
      featureID3p_c = paste(purrr::flatten_chr(featureID3p), collapse = ","),
      featureID5pEvent_c = paste(purrr::flatten_chr(featureID5pEvent), collapse = ","),
      featureID3pEvent_c = paste(purrr::flatten_chr(featureID3pEvent), collapse = ","),
      txName_c = paste(purrr::flatten_chr(txName), collapse = ","),
      gene_id_c = paste(purrr::flatten_chr(gene_id), collapse = ","),
      variantType_c = paste(purrr::flatten_chr(variantType), collapse = ","),
      as_type = tail(unlist(stringr::str_split(variantName, "_")), 1)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      !c(
        txName, variantType, genomicData, featureID5p, featureID3p,
        featureID5pEvent, featureID3pEvent
      )
    ) # |>
  # dplyr::filter(!is.na(padj))

  if (flatten) {
    sgseq <- sgseq |>
      tidyr::unnest_longer(gene_id)
  }

  sgseq
}

##### SUPPA #####
read_suppa <- function(dir, event_type = "all", gene_correction = T) {
  types_all <- c("SE", "MX", "A3", "A5", "RI", "AF", "AL")
  colnames_suppa <- c("dpsi", "p")
  if (gene_correction) colnames_suppa <- c("dpsi", "padj")

  patterns <- glue::glue(".*\\.{types_all}\\.dpsi") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )

  list_suppa <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_suppa)) {
    f <- files_read[e]
    list_suppa[[e]] <- read.delim(f) |>
      setNames(colnames_suppa) |>
      dplyr::mutate(as_type = e) |>
      tibble::rownames_to_column("event") |>
      tidyr::separate(event, c("gene_id", "coord"), sep = ";") |>
      dplyr::mutate(
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
      tibble::as_tibble()
  }

  list_suppa
}

##### Whippet #####
read_whippet <- function(file, event_type = "all") {
  types_all <- c("CE", "AA", "AD", "RI", "TS", "TE", "AF", "AL")
  et <- intersect(event_type, types_all)
  if ("all" %in% event_type) et <- types_all

  whippet <- myVroom(file) |>
    dplyr::filter(Type %in% et) |>
    dplyr::rename(
      gene_id = Gene, average_psi_1 = Psi_A,
      average_psi_2 = Psi_B, dpsi = DeltaPsi
    ) |>
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
      )
    )

  whippet
}

##### ASpli #####
read_aspli <- function(file_exon_du = NULL, file_intron_du = NULL, file_is = NULL) {
  ## non-novel ##
  if (!is.null(file_exon_du) && !is.null(file_intron_du)) {
    du_exon <- myVroom(file_exon_du)
    du_intron <- myVroom(file_intron_du)
    du <- bind_rows(du_exon, du_intron) |>
      remove_rownames() |>
      mutate(
        as_type = case_when(
          grepl("ES", event) ~ "SE",
          grepl("Alt3ss", event) ~ "A3SS",
          grepl("Alt5ss", event) ~ "A5SS",
          grepl("IR", event) ~ "RI",
          T ~ "OTHER"
        )
      ) |>
      dplyr::rename(gene_id = symbol, p = pvalue, padj = bin.fdr, ID = `...1`) # |>
      # dplyr::select(!event)
  }
  ## inclued novel ##
  if (!is.null(file_is)) {
    du <- myVroom(file_is) |>
      mutate(
        as_type = case_when(
          grepl("ES", bin.event) ~ "SE",
          grepl("Alt3ss", bin.event) ~ "A3SS",
          grepl("Alt5ss", bin.event) ~ "A5SS",
          grepl("IR", bin.event) ~ "RI",
          bin.event == "Alt 5'/3'" ~ "Alt 5'/3'",
          bin.event == "ASCE" ~ "ASCE (alternative aplicing affecting a consensus exon)",
          bin.event == "CSP" ~ "CSP (complex splicing pattern)",
          bin.event == "IoR" ~ "IoR",
          bin.event == "Novel ASP" ~ "Novel ASP",
          bin.event == "Novel CSP",
          T ~ "OTHER"
        )
      ) # |>
      # dplyr::select(!bin.event)
  }

  du
}

##### BANDITS #####
read_bandits <- function(
  transcript_result_file = NULL, gene_result_file = NULL, RData = NULL
) {
  result <- myVroom(transcript_result_file) |>
    dplyr::rename(
      gene_id = Gene_id, transcript_id = Transcript_id,
      p = p.values, padj = adj.p.values
    )

  if (!is.null(gene_result_file)) {
    result_gene <- myVroom(gene_result_file) |>
      select(Gene_id, DTU_measure)

    result <- result |>
      left_join(
        result_gene, by = c("gene_id" = "Gene_id")
      )
  }

  result
}

##### NBSplice #####
read_nbsplice <- function(result_file = NULL, RData = NULL) {
  result <- myVroom(result_file) |>
    dplyr::rename(transcript_id = iso, gene_id = gene, p = pval, padj = FDR)

  result
}

##### IsoformSwitchAnalyzeR #####
read_isoformswitchanalyzer <- function(
  feature_file = NULL, splicing_file = NULL, RData = NULL
) {
  feature <- myVroom(feature_file) |>
    dplyr::select(!IR)
  splicing <- myVroom(splicing_file)

  result <- left_join(feature, splicing, by = "isoform_id") |>
    dplyr::rename(
      transcript_id = isoform_id, gene_symbol = gene_name, average_psi_1 = IF1,
      average_psi_2 = IF2, dpsi = dIF, padj = isoform_switch_q_value,
      SE = ES, SE_genomic_start = ES_genomic_start, SE_genomic_end = ES_genomic_end,
      RI = IR, RI_genomic_start = IR_genomic_start, RI_genomic_end = IR_genomic_end,
      MXE = MEE, MXE_genomic_start = MEE_genomic_start, MXE_genomic_end = MEE_genomic_end,
      SME = MES, SME_genomic_start = MES_genomic_start, SME_genomic_end = MES_genomic_end,
      A5SS = A5, A5SS_genomic_start = A5_genomic_start, A5SS_genomic_end = A5_genomic_end,
      A3SS = A3, A3SS_genomic_start = A3_genomic_start, A3SS_genomic_end = A3_genomic_end,
      AFE = ATSS, AFE_genomic_start = ATSS_genomic_start, AFE_genomic_end = ATSS_genomic_end,
      ALE = ATTS, ALE_genomic_start = ATTS_genomic_start, ALE_genomic_end = ATTS_genomic_end
    )

  result
}

##### DiffSplice #####
read_diffsplice <- function(transcript_file) {
  result <- myVroom(transcript_file) |>
    dplyr::rename(chr = chromosome, as_type = category) |>
    dplyr::mutate(
      as_type = case_when(
        as_type == "alter_trans_start/end" ~ "AFE/ALE",
        as_type == "exon_skipping" ~ "SE",
        as_type == "intron_retention" ~ "RI",
        as_type == "mutual_exclusive" ~ "MXE",
        T ~ "OTHER"
      )
    )

  result
}

##### psichomics #####
read_psichomics <- function(result_file, condition_1, condition_2, RData = NULL) {
  result <- myVroom(result_file) |>
    dplyr::rename(
      event = `...1`, chr = Chromosome, strand = Strand, gene_symbol = Gene,
      p = `T-test p-value`, padj = `T-test p-value (BH adjusted)`, dpsi = `∆ Median`,
      average_psi_1 = dplyr::all_of(glue::glue("Median ({condition_1})")),
      average_psi_2 = dplyr::all_of(glue::glue("Median ({condition_2})"))
    ) |>
    dplyr::mutate(
      as_type = case_when(
        `Event type` == "Alternative 3' splice site (A3SS)" ~ "A3SS",
        `Event type` == "Alternative 5' splice site (A5SS)" ~ "A5SS",
        `Event type` == "Alternative first exon (AFE)" ~ "AFE",
        `Event type` == "Alternative last exon (ALE)" ~ "ALE",
        `Event type` == "Mutually exclusive exon (MXE)" ~ "MXE",
        `Event type` == "Skipped exon (SE)" ~ "SE"
      )
    ) # |>
    # dplyr::select(!`Event type`)

  result
}

##### DIEGO #####
read_diego <- function(result_file) {
  result <- myVroom(result_file) |>
    dplyr::rename(p = p_val, padj = q_val, gene_id = geneID, gene_symbol = geneName)

  result
}

##### VAST-TOOLS #####
read_vasttools <- function(
  diff_file = NULL, compare_file = NULL, event_annot_file,
  condition_1 = "s1", condition_2 = "s2"
) {
  zcat <- "zcat"
  if (getOs() == "osx") zcat <- "gzcat"

  f <- tempfile()
  cmd <- glue("{zcat} {event_annot_file} | cut -f 1,2,3,4,5,6,7,8,9,10,11 > {f}")
  system(cmd)
  annot <- myVroom(f)

  if (!is.null(diff_file) && is.null(compare_file)) {
    result <- myVroom(diff_file) |>
      dplyr::rename(
        average_psi_1 = all_of(condition_1),
        average_psi_2 = all_of(condition_2),
        dpsi = `E[dPsi]`
      ) |>
      left_join(annot, by = c("EVENT", "GENE")) |>
      dplyr::rename(gene_symbol = GENE)
  }

  if (!is.null(compare_file) && is.null(diff_file)) {
    result <- myVroom(compare_file) |>
      left_join(
        annot,
        by = c(
          "GENE", "EVENT", "COMPLEX", "COORD" = "COORD_o", "LENGTH" = "LE_o",
          "FullCO" = "FULL_CO"
        )
      )
  }

  if (!is.null(diff_file) && !is.null(compare_file)) {
    stop("either diff or compare, not both!\n")
  }

  result <- result |>
    dplyr::mutate(
      as_type = case_when(
        grepl(".*EX.*", EVENT) ~ "SE",
        grepl(".*INT.*", EVENT) ~ "RI",
        grepl(".*ALTD.*", EVENT) ~ "A5SS",
        grepl(".*ALTA.*", EVENT) ~ "A3SS",
        T ~ "OTHER"
      )
    )

  result
}

##### RATS #####
read_rats <- function(gene_result_file = NULL, transcript_result_file = NULL) {
  if (!is.null(gene_result_file)) {
    result <- myVroom(gene_result_file) |>
      dplyr::rename(gene_id = parent_id, p = pval, padj = pval_corr)
  }

  if (!is.null(transcript_result_file)) {
    result <- myVroom(transcript_result_file) |>
      dplyr::rename(
        transcript_id = target_id, gene_id = parent_id, p = pval, padj = pval_corr
      )
  }

  result
}

##### dSpliceType #####
read_dsplicetype <- function(dir, event_type = "all") {
  types_all <- c("SE", "MXE", "A3SS", "A5SS", "RI")
  patterns <- glue("dSpliceType_{types_all}_tidy.txt") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )

  list_dsplicetype <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_dsplicetype)) {
    f <- files_read[e]
    list_dsplicetype[[e]] <- myVroom(f, na_append = c("nan", "NaN")) |>
      dplyr::rename_with(~ gsub("[^'\\'']*ID", "ID", .x), !dplyr::starts_with("gene")) |> # [^'\\'']*ID
      dplyr::rename_with(~ gsub("[^'\\'']*\\_Event", "event", .x)) |>
      mutate(as_type = e) |>
      dplyr::rename(p = DS_p_value, padj = DS_adj_p_value, gene_id = geneName)
  }

  list_dsplicetype
}

##### PSI-Sigma #####
read_psisigma <- function(result_file) {
  result <- myVroom(result_file) |>
    dplyr::rename(
      gene_symbol = `Gene Symbol`, transcript_id = `Reference Transcript`,
      dpsi = `ΔPSI (%)`, p = `T-test p-value`, padj = `FDR (BH)`
    ) |>
    dplyr::mutate(
      as_type = case_when(
        `Event Type` == "A3SS" ~ "A3SS",
        `Event Type` == "A5SS" ~ "A5SS",
        `Event Type` == "IR" ~ "RI",
        `Event Type` == "MES" ~ "SME",
        `Event Type` == "SES" ~ "SE",
        `Event Type` == "TSS|A3SS" ~ "A3SS",
        `Event Type` == "TSS|A5SS" ~ "A5SS",
      )
    ) # |>
    # dplyr::select(!`Event Type`)
}

##### JuncBASE #####
read_juncbase <- function(result_file, flatten = F) {
  result <- myVroom(result_file) |>
    dplyr::mutate(
      as_type = case_when(
        as_event_type == "cassette" ~ "SE",
        as_event_type == "alternative_donor" ~ "A5SS",
        as_event_type == "alternative_acceptor" ~ "A3SS",
        as_event_type == "mutually_exclusive" ~ "MXE",
        as_event_type == "coord_cassette" ~ "coordinate_cassette",
        as_event_type == "alternative_first_exon" ~ "AFE",
        as_event_type == "alternative_last_exon" ~ "ALE",
        as_event_type == "jcn_only_AD" ~ "A5SS",
        as_event_type == "jcn_only_AA" ~ "A3SS",
        as_event_type == "intron_retention" ~ "RI"
      )
    ) |>
    dplyr::rename(
      gene_id = gene_name, average_psi_1 = set1_med, average_psi_2 = set2_med,
      dpsi = delta_val, p = raw_pval, padj = corrected_pval
    ) # |>
    # dplyr::select(!as_event_type)

  if (flatten) {
    result <- result |>
      tidyr::separate_rows(gene_id, sep = "\\,")
  }

  result
}

##### SpliceSeq #####
read_spliceseq <- function(result_file) {
  result <- myVroom(result_file, comment = "#") |>
    dplyr::mutate(
      as_type = case_when(
        `Splice Type` == "ES" ~ "SE",
        `Splice Type` == "RI" ~ "RI",
        `Splice Type` == "AA" ~ "A3SS",
        `Splice Type` == "AT" ~ "ALE",
        `Splice Type` == "AD" ~ "A5SS",
        `Splice Type` == "AP" ~ "AFE",
        `Splice Type` == "ME" ~ "MXE"
      )
    ) |>
    dplyr::rename(dpsi = dPSI, padj = `p-value`, gene_symbol = `Gene Symbol`)
    
  result
}

##### DARTS #####
read_darts <- function(dir, event_type = "all") {
  types_all <- c("SE", "A3SS", "A5SS", "RI")
  patterns <- glue::glue("{types_all}.Darts_BHT.results.xlsx") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )
  
  list_darts <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_darts)) {
    f <- files_read[e]
    list_darts[[e]] <- f |>
      openxlsx::read.xlsx(na.strings = c("nan", "NaN")) |>
      # myVroom(na_append = c("nan", "NaN")) |>
      dplyr::mutate(as_type = e) |>
      # tibble::remove_rownames() |>
      # tibble::column_to_rownames("ID...1") |>
      # dplyr::select(!dplyr::starts_with("ID")) |>
      # tibble::rownames_to_column("ID") |>
      dplyr::rename(
        dpsi = IncLevelDiff, gene_id = GeneID, gene_symbol = geneSymbol,
        psis_1 = PSI1, psis_2 = PSI2
      ) |>
      tibble::as_tibble()
  }
  
  list_darts
}

##### DARTS supplement #####
read_darts2 <- function(file, event_type = "all", read_type = c("flat", "info")) {
  types_all <- c("SE", "A3SS", "A5SS", "RI")
  if (event_type == "all") {
    et <- types_all
  } else {
    et <- intersect(types_all, event_type)
  }
  
  list_darts <- vector("list", length(et)) |>
    setNames(et)
  for (e in names(list_darts)) {
    list_darts[[e]] <- file |>
      openxlsx::read.xlsx(
        na.strings = c("nan", "NaN"), sheet = glue("{e}-{read_type[1]}")
      ) |>
      dplyr::mutate(as_type = e) |>
      dplyr::rename(
        dpsi = IncLevelDiff, gene_id = GeneID, gene_symbol = geneSymbol,
        psis_1 = PSI1, psis_2 = PSI2
      ) |>
      tibble::as_tibble()
  }
  
  list_darts
}

##### BRIE2 #####
read_brie2 <- function(dir, event_type = "all", filter = F) {
  if (filter) {
    types_all <- c("SE")
    patterns <- glue::glue("brie_quant_.*\\.brie_ident\\.{types_all}\\.filtered\\.tsv") |>
      as.character() |>
      setNames(types_all)
  } else {
    types_all <- c("SE", "MXE")
    patterns <- glue::glue("brie_quant_.*\\.brie_ident\\.{types_all}\\.tsv") |>
      as.character() |>
      setNames(types_all)
  }
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )
  
  list_brie2 <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_brie2)) {
    f <- files_read[e]
    list_brie2[[e]] <- f |>
      myVroom() |>
      dplyr::mutate(as_type = e) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        loc_1 = paste(
          unlist(stringr::str_split(unlist(stringr::str_split(GeneID, "@"))[1], ":"))[1],
          paste(unlist(stringr::str_split(unlist(stringr::str_split(GeneID, "@"))[1], ":"))[2:3], collapse = "-"),
          sep = ":"
        ),
        loc_2 = paste(
          unlist(stringr::str_split(unlist(stringr::str_split(GeneID, "@"))[2], ":"))[1],
          paste(unlist(stringr::str_split(unlist(stringr::str_split(GeneID, "@"))[2], ":"))[2:3], collapse = "-"),
          sep = ":"
        ),
        loc_3 = paste(
          unlist(stringr::str_split(unlist(stringr::str_split(GeneID, "@"))[3], ":"))[1],
          paste(unlist(stringr::str_split(unlist(stringr::str_split(GeneID, "@"))[3], ":"))[2:3], collapse = "-"),
          sep = ":"
        )
      ) |>
      dplyr::ungroup() |>
      tibble::remove_rownames() |>
      dplyr::rename_with(~ gsub("cell_type_.*_", "", .x)) |>
      dplyr::rename(gene_id = GeneID, p = pval, padj = FDR) |>
      tibble::as_tibble()
    
    if (e == "MXE") {
      list_brie2[[e]] <- list_brie2[[e]] |>
        dplyr::rowwise() |>
        dplyr::mutate(
          loc_4 = paste(
            unlist(stringr::str_split(unlist(stringr::str_split(gene_id, "@"))[4], ":"))[1],
            paste(unlist(stringr::str_split(unlist(stringr::str_split(gene_id, "@"))[4], ":"))[2:3], collapse = "-"),
            sep = ":"
          )
        )
    }
  }
  
  list_brie2
}

##### SingleSplice #####
read_singlesplice <- function(dir) {
  # results_exp <- glue::glue("{dir}/differential_expression.txt") |>
  #   myVroom()
  result_transcript <- glue::glue("{dir}/differential_transcription.txt") |>
    myVroom() |>
    dplyr::mutate(
      as_type = case_when(
        category == "exon_skipping" ~ "SE",
        category == "alter_splice_site" ~ "A5SS/A3SS",
        category == "intron_retention" ~ "RI",
        category == "alter_trans_start/end" ~ "AFE/ALE",
        T ~ "OTHER"
      )
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      loc = paste(
        unlist(stringr::str_split(location, " "))[1],
        paste(unlist(stringr::str_split(location, " "))[2:3], collapse = "-"),
        sep = ":"
      )
    ) |>
    dplyr::ungroup()
  
  result_transcript
}

##### Psix #####
read_psix <- function(result_file, annot_file) {
  result <- result_file |>
    myVroom() |>
    dplyr::rename(event = `...1`, p = pvals, padj = qvals)
  
  annot <- annot_file |>
    myVroom() |>
    dplyr::rename(ID = `...1`, gene_symbol = gene) |>
    dplyr::group_by(event) |>
    dplyr::mutate(
      myID = paste(ID, collapse = ";"),
      myintron = paste(intron, collapse = ";")
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct(event, .keep_all = T) |>
    dplyr::select(!c(ID, intron))
  
  result_final <- left_join(result, annot, by = "event")
  
  result_final
}

##### DESJ #####
read_desj <- function(file) {
  result <- file |>
    myVroom() |>
    dplyr::rename(p = P.Value, padj = adj.P.Val, gene_symbol = genename) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      loc = paste(
        unlist(stringr::str_split(junction_id, "_"))[1],
        paste(unlist(stringr::str_split(junction_id, "_"))[2:3], collapse = "-"),
        sep = ":"
      )
    )
  
  result
}

##### Outrigger #####
read_outrigger <- function(file) {
  myLoc <- function(string) {
    tmp <- string |>
      stringr::str_split("\\:") |>
      unlist()
    loc <- paste(tmp[2:3], collapse = ":")
    
    loc
  }
  
  result <- file |>
    myVroom() |>
    dplyr::mutate(
      loc_1 = NA, loc_2 = NA, loc_3 = NA, loc_4 = NA
    )
    
  for (i in 1:nrow(result)) {
    tmp <- result$ID[i] |>
      stringr::str_split("[\\|\\@]") |>
      unlist()
    
    result[["loc_1"]][i] <- tmp[1] |>
      myLoc()
    
    result[["loc_2"]][i] <- tmp[2] |>
      myLoc()
    
    result[["loc_3"]][i] <- tmp[3] |>
      myLoc()
    
    result[["loc_4"]][i] <- tmp[4] |>
      myLoc()
  }
  
  result
}

##### MicroExonator #####
read_microexonator <- function(file) {
  result <- file |>
    myVroom() |>
    dplyr::rename(
      gene_id = Gene, average_psi_1 = Psi_A.mean,
      average_psi_2 = Psi_B.mean, dpsi = DeltaPsi.mean,
      padj = cdf.beta
    ) |>
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
      )
    )
  
  result
}

##### JUM #####
read_jum <- function(
  dir, event_type = "all", condition_1 = "s1", condition_2 = "s2"
) {
  types_all <- c("cassette_exon", "MXE", "A3SS", "A5SS", "intron_retention", "composite")
  patterns <- glue::glue("{types_all}.*_final_simplified.txt") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )
  
  list_jum <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_jum)) {
    f <- files_read[e]
    list_jum[[e]] <- f |>
      myVroom(na_append = c("NaN", "Inf", "-Inf")) |>
      dplyr::mutate(as_type = e) |>
      dplyr::mutate(
        as_type = case_when(
          as_type == "cassette_exon" ~ "SE",
          as_type == "intron_retention" ~ "RI",
          as_type == "composite" ~ "OTHER",
          as_type == "MXE" ~ "MXE",
          as_type == "A5SS" ~ "A5SS",
          as_type == "A3SS" ~ "A3SS"
        )
      ) |>
      dplyr::rename(
        padj = qvalue, p = pvalue, gene_symbol = Gene, 
        dpsi = dplyr::all_of(glue::glue("deltaPSI_{condition_1}-{condition_2}"))
      )
  }
  
  list_jum
}




