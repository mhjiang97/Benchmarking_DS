run_cmds <- function(cmds, run = T) {
  if (run) myflog <- futile.logger::flog.info
  if (!run) myflog <- futile.logger::flog.debug
  futile.logger::flog.threshold(futile.logger::DEBUG)

  myflog("The whole workflow:")
  for (name_cmd in names(cmds)) {
    myflog(glue::glue("{name_cmd}: {cmds[name_cmd]}"))
  }
  for (name_cmd in names(cmds)) {
    cmds_now <- cmds[[name_cmd]]
    for (cn in cmds_now) {
      myflog(glue::glue("The running command is: {cn}."))
      if (run) {
        fail <- system(cn, wait = T)
        if (fail) {
          return(futile.logger::flog.error(glue::glue("Error: step {name_cmd} failed.")))
          # return(fail)
        } else {futile.logger::flog.info("Finished.")}
      }
    }
    if (run) futile.logger::flog.info(glue::glue("Step {name_cmd} finished successfully."))
  }
  return(0)
}

myCreatedir <- function(dirs) {
  futile.logger::flog.info("Creating directories:")
  for(d in dirs) {
    futile.logger::flog.info(d)
    if (!dir.exists(d)) {
      dir.create(d, recursive = T)
    } else {
      futile.logger::flog.warn(glue::glue("{d} exits already."))
    }
  }
  futile.logger::flog.info("Creating directories was done.")
}

myParallel <- function(cmd) {
  cmd_new <- cmd |>
    paste(collapse = " & ") |>
    paste("& wait")
  cmd_new
}

myParallel2 <- function(cmd, n = 5) {
  list_cmd <- split(cmd, ceiling(seq_along(cmd) / n))
  cmd_new <- sapply(list_cmd, myParallel) |>
    paste(collapse = " ; ")
  cmd_new
}

myParallel3 <- function(cmd, n = 5) {
  list_cmd <- split(cmd, ceiling(seq_along(cmd) / n))
  cmd_new <- sapply(list_cmd, myParallel)
  cmd_new
}

`%>%` <- dplyr::`%>%`
getBasics <- function(sampletable) {
  condition_1 <- unique(sampletable$conditions)[1]
  condition_2 <- unique(sampletable$conditions)[2]

  samples_1 <- sampletable$samples[sampletable$conditions == condition_1]
  samples_2 <- sampletable$samples[sampletable$conditions == condition_2]

  bams_1 <- sampletable$files_bam[sampletable$conditions == condition_1] |>
    stats::setNames(samples_1)
  bams_2 <- sampletable$files_bam[sampletable$conditions == condition_2] |>
    stats::setNames(samples_2)

  if ("dirs_salmon" %in% colnames(sampletable)) {
    salmons_1 <- sampletable$dirs_salmon[sampletable$conditions == condition_1] |>
      stats::setNames(samples_1)
    salmons_2 <- sampletable$dirs_salmon[sampletable$conditions == condition_2] |>
      stats::setNames(samples_2)
  } else {
    salmons_1 <- salmons_2 <- NULL
  }

  if ("files_fq" %in% colnames(sampletable)) {
    fqs_R1_1 <- sampletable$files_fq[sampletable$conditions == condition_1] |>
      strsplit(";") |>
      stats::setNames(function(x) {return(x[1])}) |>
      stats::setNames(samples_1)
    fqs_R2_1 <- sampletable$files_fq[sampletable$conditions == condition_1] |>
      strsplit(";") |>
      sapply(function(x) {return(x[2])}) |>
      stats::setNames(samples_1)

    fqs_R1_2 <- sampletable$files_fq[sampletable$conditions == condition_2] |>
      strsplit(";") |>
      sapply(function(x) {return(x[1])}) |>
      stats::setNames(samples_2)
    fqs_R2_2 <- sampletable$files_fq[sampletable$conditions == condition_2] |>
      strsplit(";") |>
      sapply(function(x) {return(x[2])}) |>
      stats::setNames(samples_2)
  } else {
    fqs_R1_1 <- fqs_R2_1 <- fqs_R1_2 <- fqs_R2_2 <- NULL
  }


  files_star_junc <- vector("character", length = nrow(sampletable)) |>
    stats::setNames(sampletable$samples)
  for (s in names(files_star_junc)) {
    if (is.null(sampletable$files_bam[sampletable$samples == s])) {
      files_star_junc[s] <- ""
    } else {
      tmp_star_junc <- list.files(
        dirname(sampletable$files_bam[sampletable$samples == s]),
        pattern = "SJ.out.tab", full.names = T
      )
      if (all(is.na(tmp_star_junc))) {
        files_star_junc[s] <- ""
      } else {
        files_star_junc[s] <- tmp_star_junc
      }
    }
  }
  juncs_star_1 <- files_star_junc[samples_1] |>
    stats::setNames(samples_1)
  juncs_star_2 <- files_star_junc[samples_2] |>
    stats::setNames(samples_2)

  files_star_quant <- vector("character", length = nrow(sampletable)) |>
    stats::setNames(sampletable$samples)
  for (s in names(files_star_quant)) {
    if (is.null(sampletable$files_bam[sampletable$samples == s])) {
      files_star_quant[s] <- ""
    } else {
      tmp_star_quant <- list.files(
        dirname(sampletable$files_bam[sampletable$samples == s]),
        pattern = "ReadsPerGene.out.tab", full.names = T
      )
      if (all(is.na(tmp_star_quant))) {
        files_star_quant[s] <- ""
      } else {
        files_star_quant[s] <- tmp_star_quant
      }
    }
  }
  quants_star_1 <- files_star_quant[samples_1] |>
    stats::setNames(samples_1)
  quants_star_2 <- files_star_quant[samples_2] |>
    stats::setNames(samples_2)

  # files_tophat_junc <- vector("character", length = nrow(sampletable)) |>
  #   stats::setNames(sampletable$samples)
  # for (s in names(files_tophat_junc)) {
  #   if (is.null(sampletable$files_bam_tophat[sampletable$samples == s])) {
  #     files_tophat_junc[s] <- ""
  #   } else {
  #     tmp_tophat_junc <- list.files(
  #       dirname(sampletable$files_bam_tophat[sampletable$samples == s]),
  #       pattern = "junctions.bed", full.names = T
  #     )
  #     if (all(is.na(tmp_tophat_junc))) {
  #       files_tophat_junc[s] <- ""
  #     } else {
  #       files_tophat_junc[s] <- tmp_tophat_junc
  #     }
  #   }
  # }
  # juncs_tophat_1 <- files_tophat_junc[samples_1] |>
  #   stats::setNames(samples_1)
  # juncs_tophat_2 <- files_tophat_junc[samples_2] |>
  #   stats::setNames(samples_2)

  library_type_most <- table(sampletable$library_types) %>% .[. == max(.)] |> names()
  read_length_most <- table(sampletable$read_lengths) %>% .[. == max(.)] |> names() |> as.numeric()
  strandedness_most <- table(sampletable$strandedness) %>% .[. == max(.)] |> names()

  list(
    condition_1, condition_2,
    samples_1, samples_2,
    bams_1, bams_2,
    salmons_1, salmons_2,
    fqs_R1_1, fqs_R2_1,
    fqs_R1_2, fqs_R2_2,
    juncs_star_1, juncs_star_2,
    quants_star_1, quants_star_2,
    # juncs_tophat_1, juncs_tophat_2,
    library_type_most, read_length_most, strandedness_most
  ) |>
    stats::setNames(
      c(
        "condition_1", "condition_2", "samples_1", "samples_2", "bams_1", "bams_2",
        "salmons_1", "salmons_2", "fqs_R1_1", "fqs_R2_1", "fqs_R1_2", "fqs_R2_2",
        "juncs_star_1", "juncs_star_2", "quants_star_1", "quants_star_2",
        # "juncs_tophat_1", "juncs_tophat_2",
        "library_type_most", "read_length_most", "strandedness_most"
      )
    )
}




