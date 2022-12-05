pkgs <- c(
  "stringr", "glue", "tibble", "dplyr", "tidyr", "ggplot2", "RColorBrewer",
  "vroom", "ggbreak", "ggthemes", "ggsci", "scales", "openxlsx", "forcats",
  "reshape2", "Hmisc", "purrr", "corrplot", "DEXSeq"
)
for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = T))

##### function #####
calc_FDR <- function(
  results, tools, truths_m, colors_tool, g_or_e = "g",
  thresholds = c(0.01, 0.05, 0.1), pch_tool = NULL, filter_as = F
) {
  if (g_or_e == "g") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (filter_as) {
        keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
        if (length(keep_as) > 0) {
          tmp <- tmp |> filter(as_type %in% keep_as)
          truths_m <- truths_m |> filter(as_type %in% keep_as)
        }
      }

      sapply(thresholds, function(i) {
        tp <- length(
          intersect(
            truths_m$gene_id[truths_m$ds_status == 1],
            unlist(stringr::str_split(tmp$gene_id[as.numeric(tmp$padj) <= i], ";"))
          )
        )
        p <- length(unique(tmp$gene_id[as.numeric(tmp$padj) <= i]))
        1 - tp / p
      }) |>
        as.character() |>
        union_all(t) |>
        setNames(c(as.character(thresholds), "tool"))
    })
  }

  if (g_or_e == "e") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (!"idx" %in% colnames(tmp)) return(NULL)

      keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
      if (length(keep_as) > 0) {
        tmp <- tmp |> filter(as_type %in% keep_as)
      }

      sapply(thresholds, function(i) {
        tp <- length(unique(tmp$myID[as.numeric(tmp$padj) <= i & !is.na(tmp$idx)]))
        p <- length(unique(tmp$myID[as.numeric(tmp$padj) <= i]))
        1 - tp / p
      }) |>
        as.character() |>
        union_all(t) |>
        setNames(c(as.character(thresholds), "tool"))
    })
  }

  out <- out |>
    bind_rows() |>
    mutate(
      tool = factor(tool, levels = unique(tool)),
      colors_tool = colors_tool[as.character(tool)],
      g_or_e = g_or_e # ,
      # `0.1` = as.double(`0.1`),
      # `0.05` = as.double(`0.05`),
      # `0.01` = as.double(`0.01`)
    ) |>
    mutate(
      g_or_e = case_when(
        g_or_e == "g" ~ "Gene Level",
        g_or_e == "e" ~ "Event Level"
      )
    )
  
  for (threshold in thresholds) {
    out[[as.character(threshold)]] <- as.double(out[[as.character(threshold)]])
  }

  # rownames(out) <- out$tool
  if (!is.null(pch_tool)) out$pch_tool <- pch_tool[out$tool]

  out
}

calc_TPR <- function(
  results, tools, truths_m, colors_tool, g_or_e = "g",
  thresholds = c(0.01, 0.05, 0.1), pch_tool = NULL, filter_as = F
) {
  if (g_or_e == "g") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (filter_as) {
        keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
        if (length(keep_as) > 0) {
          tmp <- tmp |> filter(as_type %in% keep_as)
          truths_m <- truths_m |> filter(as_type %in% keep_as)
        }
      }

      sapply(thresholds, function(i) {
        tp <- length(
          intersect(
            truths_m$gene_id[truths_m$ds_status == 1],
            stringr::str_split(tmp$gene_id[as.numeric(tmp$padj) <= i], ";")
          )
        )
        pp <- length(unique(truths_m$gene_id[truths_m$ds_status == 1]))
        tp / pp
      }) |>
        as.character() |>
        union_all(t) |>
        setNames(c(as.character(thresholds), "tool"))
    })
  }

  if (g_or_e == "e") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (!"idx" %in% colnames(tmp)) return(NULL)

      keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
      if (length(keep_as) > 0) {
        tmp <- tmp |> filter(as_type %in% keep_as)
        truths_m <- truths_m |> filter(as_type %in% keep_as)
      } else {
        truths_m <- truths_m |> filter(!is.na(as_type))
      }

      sapply(thresholds, function(i) {
        tp <- tmp$idx[as.numeric(tmp$padj) <= i & !is.na(tmp$idx)] |>
          str_split("\\,") |>
          unlist() |>
          unique() |>
          length()
        # tp <- length(unique(tmp$myID[as.numeric(tmp$padj) <= i & !is.na(tmp$idx)]))
        pp <- nrow(truths_m[truths_m$ds_status == 1, ])
        tp / pp
      }) |>
        as.character() |>
        union_all(t) |>
        setNames(c(as.character(thresholds), "tool"))
    })
  }

  out <- out |>
    bind_rows() |>
    mutate(
      tool = factor(tool, levels = unique(tool)),
      colors_tool = colors_tool[as.character(tool)],
      g_or_e = g_or_e # ,
      # `0.1` = as.double(`0.1`),
      # `0.05` = as.double(`0.05`),
      # `0.01` = as.double(`0.01`)
    ) |>
    mutate(
      g_or_e = case_when(
        g_or_e == "g" ~ "Gene Level",
        g_or_e == "e" ~ "Event Level"
      )
    )

  for (threshold in thresholds) {
    out[[as.character(threshold)]] <- as.double(out[[as.character(threshold)]])
  }
  
  # rownames(out) <- out$tool
  if (!is.null(pch_tool)) out$pch_tool <- pch_tool[out$tool]

  out
}

calc_FDR_as <- function(
  results, tools, truths_m, colors_tool, g_or_e = "g",
  thresholds = c(0.01, 0.05, 0.1), pch_tool = NULL, ...
) {
  if (g_or_e == "g") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (!"as_type" %in% colnames(tmp)) return(NULL)
      keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
      if (length(keep_as) > 0) {
        tmp <- tmp |> filter(as_type %in% keep_as)
        truths_m <- truths_m |> filter(as_type %in% keep_as)
      } else {
        return(NULL)
      }

      as.data.frame(do.call(rbind, lapply(keep_as, function(as) {
        sapply(thresholds, function(i) {
          tp <- length(intersect(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == as], tmp$gene_id[as.numeric(tmp$padj) <= i & tmp$as_type == as]))
          p <- length(unique(tmp$gene_id[as.numeric(tmp$padj) <= i & tmp$as_type == as]))
          1 - tp / p
        }) |>
          as.character() |>
          union_all(c(t, as)) |>
          setNames(c(as.character(thresholds), "tool", "as_type"))
      })))
    })
  }

  if (g_or_e == "e") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (!"idx" %in% colnames(tmp) || !"as_type" %in% colnames(tmp)) return(NULL)
      keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
      if (length(keep_as) > 0) {
        tmp <- tmp |> filter(as_type %in% keep_as)
      } else {
        return(NULL)
      }

      as.data.frame(do.call(rbind, lapply(keep_as, function(as) {
        sapply(thresholds, function(i) {
          tp <- length(unique(tmp$myID[as.numeric(tmp$padj) <= i & !is.na(tmp$idx) & tmp$as_type == as]))
          p <- length(unique(tmp$myID[as.numeric(tmp$padj) <= i & tmp$as_type == as]))
          1 - tp / p
        }) |>
          as.character() |>
          union_all(c(t, as)) |>
          setNames(c(as.character(thresholds), "tool", "as_type"))
      })))
    })
  }

  out <- out |>
    bind_rows() |>
    mutate(
      tool = factor(tool, levels = unique(tool)),
      colors_tool = colors_tool[as.character(tool)],
      g_or_e = g_or_e # ,
      # `0.1` = as.double(`0.1`),
      # `0.05` = as.double(`0.05`),
      # `0.01` = as.double(`0.01`)
    ) |>
    mutate(
      g_or_e = case_when(
        g_or_e == "g" ~ "Gene Level",
        g_or_e == "e" ~ "Event Level"
      )
    )

  for (threshold in thresholds) {
    out[[as.character(threshold)]] <- as.double(out[[as.character(threshold)]])
  }
  
  # rownames(out) <- out$tool
  if (!is.null(pch_tool)) out$pch_tool <- pch_tool[out$tool]

  out
}

calc_TPR_as <- function(
  results, tools, truths_m, colors_tool, g_or_e = "g",
  thresholds = c(0.01, 0.05, 0.1), pch_tool = NULL, ...
) {
  if (g_or_e == "g") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (!"as_type" %in% colnames(tmp)) return(NULL)
      keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
      if (length(keep_as) > 0) {
        tmp <- tmp |> filter(as_type %in% keep_as)
        truths_m <- truths_m |> filter(as_type %in% keep_as)
      } else {
        return(NULL)
      }

      as.data.frame(do.call(rbind, lapply(keep_as, function(as) {
        sapply(thresholds, function(i) {
          tp <- length(intersect(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == as], tmp$gene_id[as.numeric(tmp$padj) <= i & tmp$as_type == as]))
          pp <- length(unique(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == as]))
          tp / pp
        }) |>
          as.character() |>
          union_all(c(t, as)) |>
          setNames(c(as.character(thresholds), "tool", "as_type"))
      })))
    })
  }

  if (g_or_e == "e") {
    out <- lapply(tools, function(t) {
      if (t == "VAST-TOOLS") thresholds <- 0.05
      
      tmp <- results[[t]] |>
        filter(gene_id %in% truths_m$gene_id)

      if (!"idx" %in% colnames(tmp) || !"as_type" %in% colnames(tmp)) return(NULL)
      keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI"), unique(tmp$as_type))
      if (length(keep_as) > 0) {
        tmp <- tmp |> filter(as_type %in% keep_as)
        truths_m <- truths_m |> filter(as_type %in% keep_as)
      } else {
        return(NULL)
      }

      as.data.frame(do.call(rbind, lapply(keep_as, function(as) {
        sapply(thresholds, function(i) {
          tp <- tmp$idx[as.numeric(tmp$padj) <= i & !is.na(tmp$idx) & tmp$as_type == as] |>
            str_split("\\,") |>
            unlist() |>
            unique() |>
            length()
          # tp <- length(unique(tmp$myID[as.numeric(tmp$padj) <= i & !is.na(tmp$idx) & tmp$as_type == as]))
          pp <- nrow(truths_m[truths_m$ds_status == 1 & truths_m$as_type == as, ])
          tp / pp
        }) |>
          as.character() |>
          union_all(c(t, as)) |>
          setNames(c(as.character(thresholds), "tool", "as_type"))
      })))
    })
  }

  out <- out |>
    bind_rows() |>
    mutate(
      tool = factor(tool, levels = unique(tool)),
      colors_tool = colors_tool[as.character(tool)],
      g_or_e = g_or_e # ,
      # `0.1` = as.double(`0.1`),
      # `0.05` = as.double(`0.05`),
      # `0.01` = as.double(`0.01`)
    ) |>
    mutate(
      g_or_e = case_when(
        g_or_e == "g" ~ "Gene Level",
        g_or_e == "e" ~ "Event Level"
      )
    )

  for (threshold in thresholds) {
    out[[as.character(threshold)]] <- as.double(out[[as.character(threshold)]])
  }
  
  # rownames(out) <- out$tool
  if (!is.null(pch_tool)) out$pch_tool <- pch_tool[out$tool]

  out
}

plot_fdr_tpr_paper <- function(
  results, truths_m, tools = NULL, colors_tool, thresholds = c(0.01, 0.05, 0.1),
  split_variable = NULL, pch_tool = NULL, filter_as = F, as = F, new_plot_type = F,
  output_filename, pathsize = 1, pointsize = 3, stripsize = 20, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 2, fig_height = 7, fig_width = 10, long = F, legendtextsize = 10
) {
  
  funcs_fdr <- list(
    "as" = calc_FDR_as, "no_as" = calc_FDR
  )
  funcs_tpr <- list(
    "as" = calc_TPR_as, "no_as" = calc_TPR
  )
  
  if (as) {
    fdr <- funcs_fdr[["as"]]
    tpr <- funcs_tpr[["as"]]
  } else {
    fdr <- funcs_fdr[["no_as"]]
    tpr <- funcs_tpr[["no_as"]]
  }
  ## Calculate FDR
  if (is.null(split_variable)) {
    FDR_g <- fdr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    TPR_g <- tpr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )

    FDR_e <- fdr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    TPR_e <- tpr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )

  } else {
    FDR_g <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) fdr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
    TPR_g <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) tpr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )

    FDR_e <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) fdr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
    TPR_e <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) tpr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
  }

  FDR_g <- melt(FDR_g, variable.name = "threshold", value.name = "FDR")
  TPR_g <- melt(TPR_g, variable.name = "threshold", value.name = "TPR")

  FDR_e <- melt(FDR_e, variable.name = "threshold", value.name = "FDR")
  TPR_e <- melt(TPR_e, variable.name = "threshold", value.name = "TPR")

  ## Put results together
  FDR_TPR_g <- merge(FDR_g, TPR_g)

  FDR_TPR_e <- merge(FDR_e, TPR_e)

  FDR_TPR <- rbind(FDR_TPR_g, FDR_TPR_e)
  FDR_TPR$threshold <- as.numeric(as.character(FDR_TPR$threshold))
  FDR_TPR$fill_color <- FDR_TPR$colors_tool
  FDR_TPR$fill_color[FDR_TPR$FDR > FDR_TPR$threshold] <- "white"
  FDR_TPR$g_or_e <- factor(FDR_TPR$g_or_e, levels = c("Gene Level", "Event Level"))
  fill_color <- unique(FDR_TPR$fill_color)
  names(fill_color) <- fill_color
  uniq <- !duplicated(FDR_TPR$tool)
  col_color <- FDR_TPR$colors_tool[uniq]
  names(col_color) <- FDR_TPR$tool[uniq]
  # uniq2 <- !duplicated(pch_tool)
  # pchs <- pch_tool[uniq2]
  # names(pchs) <- pch_tool[uniq2]

  write.table(
    FDR_TPR, file = gsub("\\.pdf$", ".txt", output_filename),
    row.names = F, col.names = T, sep = "\t", quote = F
  )
  
  if (!as) {
    if (is.null(split_variable)) {
      ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = tool)) +
        geom_vline(xintercept = seq(0, 1, 0.1),
                   colour = "lightgrey", linetype = "dashed") +
        geom_vline(xintercept = thresholds, linetype = "dashed") -> p
      
      if (length(thresholds) > 1) p <- p + geom_path(size = pathsize, aes(colour = tool))
        # geom_path(size = pathsize, aes(colour = tool)) +
      p <- p + facet_wrap(~ g_or_e) +
        geom_point(size = pointsize + 1,
                   aes(colour = tool), shape = 19) +
        geom_point(size = pointsize,
                   aes(fill = fill_color, colour = tool), shape = 21) +
        scale_fill_manual(values = fill_color, guide = "none") +
        scale_color_manual(values = col_color, name = "") +
        guides(col = guide_legend(nrow = legend_nrow)) +
        ylim(0, 1) +
        scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                           labels = c("", thresholds, seq(0.1, 1, 0.1)),
                           limits = c(0, 1)) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) +
        ggtitle("")
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    } else {
      ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = tool)) +
        geom_vline(xintercept = seq(0, 1, 0.1),
                   colour = "lightgrey", linetype = "dashed") +
        geom_vline(xintercept = thresholds, linetype = "dashed") -> p
      
      if (length(thresholds) > 1) p <- p + geom_path(size = pathsize, aes(colour = tool))
        # geom_path(size = pathsize, aes(colour = tool)) +
      p <- p + facet_wrap(g_or_e ~ L1) +
        geom_point(size = pointsize + 1,
                   aes(colour = tool), shape = 19) +
        geom_point(size = pointsize,
                   aes(fill = fill_color, colour = tool), shape = 21) +
        scale_fill_manual(values = fill_color, guide = "none") +
        scale_color_manual(values = col_color, name = "") +
        guides(col = guide_legend(nrow = legend_nrow)) +
        ylim(0, 1) +
        scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                           labels = c("", thresholds, seq(0.1, 1, 0.1)),
                           limits = c(0, 1)) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) +
        ggtitle("")
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    }
  } else {
      if (is.null(split_variable)) {
        ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = tool)) +
          geom_vline(xintercept = seq(0, 1, 0.1),
                     colour = "lightgrey", linetype = "dashed") +
          geom_vline(xintercept = thresholds, linetype = "dashed") -> p
        
        if (length(thresholds) > 1) p <- p + geom_path(size = pathsize, aes(colour = tool))
          # geom_path(size = pathsize, aes(colour = tool)) +
        p <- p + geom_point(size = pointsize + 1,
                            aes(colour = tool), shape = 19) +
          geom_point(size = pointsize,
                     aes(fill = fill_color, colour = tool), shape = 21) +
          scale_fill_manual(values = fill_color, guide = "none") +
          scale_color_manual(values = col_color, name = "") +
          guides(col = guide_legend(nrow = legend_nrow)) +
          ylim(0, 1) +
          scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                             labels = c("", thresholds, seq(0.1, 1, 0.1)),
                             limits = c(0, 1)) +
          theme(legend.position = "bottom",
                legend.text = element_text(size = legendtextsize),
                panel.background = element_rect(fill = NA, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black", size = 1),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                           size = axistextsize),
                axis.text.y = element_text(size = axistextsize),
                axis.title.x = element_text(size = axistitlesize),
                axis.title.y = element_text(size = axistitlesize),
                strip.text = element_text(size = stripsize),
                strip.background = element_rect(fill = NA, colour = "black")) +
          ggtitle("")
        if (long) p <- p + facet_grid(as_type ~ g_or_e)
        if (!long) p <- p + facet_grid(g_or_e ~ as_type)
        pdf(output_filename, width = fig_width, height = fig_height)
        print(p)
        dev.off()
      } else {
        if (new_plot_type) {
          FDR_TPR_new <- FDR_TPR |>
            mutate(
              tool_col = tool,
              tool = paste(tool_col, g_or_e)
            )
          ggplot(FDR_TPR_new, aes(x = FDR, y = TPR, group = tool)) + 
            geom_vline(xintercept = seq(0, 1, 0.1), 
                       colour = "lightgrey", linetype = "dashed") + 
            geom_vline(xintercept = thresholds, linetype = "dashed") -> p
          
          if (length(thresholds) > 1) p <- p + geom_path(size = pathsize, aes(colour = tool_col))
            # geom_path(size = pathsize, aes(colour = tool_col)) + 
            # facet_grid(as_type ~ L1) + 
          p <- p + geom_point(size = pointsize + 1, 
                              aes(colour = tool_col, shape = g_or_e)) + 
            geom_point(size = pointsize + 0.5, 
                       aes(colour = tool_col, shape = g_or_e)) + 
            geom_point(size = pointsize, 
                       aes(fill = fill_color, colour = tool_col, shape = g_or_e)) + 
            scale_fill_manual(values = fill_color, guide = "none") + 
            scale_color_manual(values = col_color, name = "") +
            scale_shape_manual(values = c("Gene Level" = 21, "Event Level" = 24), name = "") + 
            guides(col = guide_legend(nrow = legend_nrow)) + 
            ylim(0, 1) +
            scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                               labels = c("", thresholds, seq(0.1, 1, 0.1)),
                               limits = c(0, 1)) + 
            theme(legend.position = "bottom", 
                  legend.text = element_text(size = legendtextsize),
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle("")
          if (long) p <- p + facet_grid(as_type ~ L1)
          if (!long) p <- p + facet_grid(L1 ~ as_type)
          pdf(output_filename, width = fig_width, height = fig_height)
          print(p)
          dev.off()
        } else {
          ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = tool)) +
            geom_vline(xintercept = seq(0, 1, 0.1),
                       colour = "lightgrey", linetype = "dashed") +
            geom_vline(xintercept = thresholds, linetype = "dashed") -> p
          
          if (length(thresholds) > 1) p <- p + geom_path(size = pathsize, aes(colour = tool))
            # geom_path(size = pathsize, aes(colour = tool)) +
            # facet_grid(L1 + as_type ~ g_or_e) +
          p <- p + geom_point(size = pointsize + 1,
                              aes(colour = tool), shape = 19) +
            geom_point(size = pointsize,
                       aes(fill = fill_color, colour = tool), shape = 21) +
            scale_fill_manual(values = fill_color, guide = "none") +
            scale_color_manual(values = col_color, name = "") +
            guides(col = guide_legend(nrow = legend_nrow)) +
            ylim(0, 1) +
            scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                               labels = c("", thresholds, seq(0.1, 1, 0.1)),
                               limits = c(0, 1)) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = legendtextsize),
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) +
            ggtitle("")
          if (long) p <- p + facet_grid(L1 + as_type ~ g_or_e)
          if (!long) p <- p + facet_grid(L1 + g_or_e ~ as_type)
          pdf(output_filename, width = fig_width, height = fig_height)
          print(p)
          dev.off()
        }
      }
  }
}

plot_tpr_paper <- function(
  results, truths_m, tools = NULL, colors_tool, thresholds = c(0.01, 0.05, 0.1),
  split_variable = NULL, pch_tool = NULL, filter_as = F, as = F, new_plot_type = F,
  output_filename, pathsize = 1, pointsize = 3, stripsize = 20, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 2, fig_height = 7, fig_width = 10, long = F, legendtextsize = 10
) {
  funcs_tpr <- list(
    "as" = calc_TPR_as, "no_as" = calc_TPR
  )
  
  if (as) {
    tpr <- funcs_tpr[["as"]]
  } else {
    tpr <- funcs_tpr[["no_as"]]
  }
  ## Calculate FDR
  if (is.null(split_variable)) {
    TPR_g <- tpr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    
    TPR_e <- tpr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    
  } else {
    TPR_g <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) tpr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
    
    TPR_e <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) tpr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
  }


  TPR_g <- melt(TPR_g, variable.name = "threshold", value.name = "TPR")
  
  TPR_e <- melt(TPR_e, variable.name = "threshold", value.name = "TPR")
  
  ## Put results together
  TPR <- rbind(TPR_g, TPR_e)
  TPR$threshold <- as.numeric(as.character(TPR$threshold))
  uniq <- !duplicated(TPR$tool)
  col_color <- TPR$colors_tool[uniq]
  names(col_color) <- TPR$tool[uniq]
  
  write.table(TPR, file = gsub("\\.pdf$", ".txt", output_filename), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  if (!as) {
    if (is.null(split_variable)) {
      ggplot(TPR, aes(x = TPR, y = tool, group = tool)) -> p
      if (length(thresholds) > 1) p <- p + geom_path(size = 1, aes(colour = tool))
        # geom_path(size = 1, aes(colour = tool)) + 
        
      p <- p + facet_wrap(~ g_or_e) + 
        geom_point(size = pointsize + 1, 
                   aes(colour = tool), shape = 19) + 
        scale_color_manual(values = col_color, name = "") +
        guides(col = guide_legend(nrow = legend_nrow)) + 
        xlim(0, 1) + ylab("") +
        theme(legend.position = "bottom", 
              legend.text = element_text(size = legendtextsize),
              panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dotted"),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) + 
        ggtitle("")
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    } else {
      ggplot(TPR, aes(x = TPR, y = tool, group = tool)) -> p
      if (length(thresholds) > 1) p <- p + geom_path(size = 1, aes(colour = tool))
        # geom_path(size = 1, aes(colour = tool)) + 
        
      p <- p + facet_wrap(g_or_e ~ L1) + 
        geom_point(size = pointsize + 1, 
                   aes(colour = tool), shape = 19) + 
        scale_color_manual(values = col_color, name = "") +
        guides(col = guide_legend(nrow = legend_nrow)) + 
        xlim(0, 1) + ylab("") + 
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dotted"),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) + 
        ggtitle("")
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    }
  } else {
    if (is.null(split_variable)) {
      ggplot(TPR, aes(x = TPR, y = tool, group = tool)) -> p
      if (length(thresholds) > 1) p <- p + geom_path(size = 1, aes(colour = tool))
        # geom_path(size = 1, aes(colour = tool)) + 
      
      p <- p + # facet_grid(g_or_e ~ as_type) + 
        geom_point(size = pointsize + 1, 
                   aes(colour = tool), shape = 19) + 
        scale_color_manual(values = col_color, name = "") +
        guides(col = guide_legend(nrow = legend_nrow)) + 
        xlim(0, 1) + ylab("") +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dotted"),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) + 
        ggtitle("")
      if (long) p <- p + facet_grid(as_type ~ g_or_e)
      if (!long) p <- p + facet_grid(g_or_e ~ as_type)
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    } else {
      if (new_plot_type) {
        TPR_new <- TPR |>
          mutate(
            tool_col = tool,
            tool = paste(tool_col, g_or_e)
          )
        ggplot(TPR_new, aes(x = TPR, y = tool, group = tool)) -> p
        if (length(thresholds) > 1) p <- p + geom_path(size = 1, aes(colour = tool_col))
          # geom_path(size = 1, aes(colour = tool_col)) + 
        
        p <- p + # facet_grid(as_type ~ L1) + 
          geom_point(size = pointsize + 1, 
                     aes(colour = tool_col, shape = g_or_e)) + 
          scale_color_manual(values = col_color, name = "") +
          scale_shape_manual(values = c("Gene Level" = 19, "Event Level" = 17), name = "") + 
          guides(col = guide_legend(nrow = legend_nrow)) + 
          xlim(0, 1) + ylab("") + 
          theme(legend.position = "bottom",
                legend.text = element_text(size = legendtextsize),
                panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dotted"),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = NA, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black", size = 1),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                           size = axistextsize),
                axis.text.y = element_text(size = axistextsize),
                axis.title.x = element_text(size = axistitlesize),
                axis.title.y = element_text(size = axistitlesize),
                strip.text = element_text(size = stripsize),
                strip.background = element_rect(fill = NA, colour = "black")) + 
          ggtitle("")
        if (long) p <- p + facet_grid(L1 ~ as_type)
        if (!long) p <- p + facet_grid(as_type ~ L1)
        pdf(output_filename, width = fig_width, height = fig_height)
        print(p)
        dev.off()
      } else {
        ggplot(TPR, aes(x = TPR, y = tool, group = tool)) -> p
        if (length(thresholds) > 1) p <- p + geom_path(size = 1, aes(colour = tool))
          # geom_path(size = 1, aes(colour = tool)) + 
          
        p <- p + # facet_grid(L1 + g_or_e ~ as_type) + 
          geom_point(size = pointsize + 1, 
                     aes(colour = tool), shape = 19) + 
          scale_color_manual(values = col_color, name = "") +
          guides(col = guide_legend(nrow = legend_nrow)) + 
          xlim(0, 1) + ylab("") + 
          theme(legend.position = "bottom",
                legend.text = element_text(size = legendtextsize),
                panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dotted"),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = NA, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black", size = 1),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                           size = axistextsize),
                axis.text.y = element_text(size = axistextsize),
                axis.title.x = element_text(size = axistitlesize),
                axis.title.y = element_text(size = axistitlesize),
                strip.text = element_text(size = stripsize),
                strip.background = element_rect(fill = NA, colour = "black")) + 
          ggtitle("")
        if (long) p <- p + facet_grid(L1 + as_type ~ g_or_e)
        if (!long) p <- p + facet_grid(L1 + g_or_e ~ as_type)
        pdf(output_filename, width = fig_width, height = fig_height)
        print(p)
        dev.off()
      }
    }
  }
}

plot_f_score <- function(
  results, truths_m, tools = NULL, colors_tool, thresholds = c(0.05),
  split_variable = NULL, pch_tool = NULL, filter_as = F, as = F,
  output_filename, stripsize = 20, axistitlesize = 20, axistextsize = 15,
  fig_height = 7, fig_width = 10, long = F, legendtextsize = 10
) {
  
  funcs_fdr <- list(
    "as" = calc_FDR_as, "no_as" = calc_FDR
  )
  funcs_tpr <- list(
    "as" = calc_TPR_as, "no_as" = calc_TPR
  )
  
  if (as) {
    fdr <- funcs_fdr[["as"]]
    tpr <- funcs_tpr[["as"]]
  } else {
    fdr <- funcs_fdr[["no_as"]]
    tpr <- funcs_tpr[["no_as"]]
  }
  ## Calculate FDR
  if (is.null(split_variable)) {
    FDR_g <- fdr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    TPR_g <- tpr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    
    FDR_e <- fdr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    TPR_e <- tpr(
      results = results, tools = tools, truths_m = truths_m, colors_tool = colors_tool,
      g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
    )
    
  } else {
    FDR_g <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) fdr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
    TPR_g <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) tpr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "g", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
    
    FDR_e <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) fdr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
    TPR_e <- lapply(
      split(truths_m, truths_m[, split_variable]),
      function(w) tpr(
        results = results, tools = tools, truths_m = w, colors_tool = colors_tool,
        g_or_e = "e", thresholds = thresholds, pch_tool = pch_tool, filter_as = filter_as
      )
    )
  }
  
  FDR_g <- melt(FDR_g, variable.name = "threshold", value.name = "FDR")
  TPR_g <- melt(TPR_g, variable.name = "threshold", value.name = "TPR")
  
  FDR_e <- melt(FDR_e, variable.name = "threshold", value.name = "FDR")
  TPR_e <- melt(TPR_e, variable.name = "threshold", value.name = "TPR")
  
  ## Put results together
  FDR_TPR_g <- merge(FDR_g, TPR_g)
  
  FDR_TPR_e <- merge(FDR_e, TPR_e)
  
  FDR_TPR <- rbind(FDR_TPR_g, FDR_TPR_e)
  FDR_TPR$threshold <- as.numeric(as.character(FDR_TPR$threshold))
  FDR_TPR$fill_color <- FDR_TPR$colors_tool
  FDR_TPR$fill_color[FDR_TPR$FDR > FDR_TPR$threshold] <- "white"
  FDR_TPR$g_or_e <- factor(FDR_TPR$g_or_e, levels = c("Gene Level", "Event Level"))
  fill_color <- unique(FDR_TPR$fill_color)
  names(fill_color) <- fill_color
  uniq <- !duplicated(FDR_TPR$tool)
  col_color <- FDR_TPR$colors_tool[uniq]
  names(col_color) <- FDR_TPR$tool[uniq]
  # uniq2 <- !duplicated(pch_tool)
  # pchs <- pch_tool[uniq2]
  # names(pchs) <- pch_tool[uniq2]
  
  F_score <- FDR_TPR |>
    mutate(`F-score` = (2 * (1 - FDR) * TPR) / ((1 - FDR) + TPR)) |>
    mutate(tool = fct_reorder(tool, `F-score`))
  if (!as && length(unique(F_score$L1)) == 1) {
    F_score <- F_score |>
      group_by(tool, g_or_e) |>
      mutate(`F-score` = mean(`F-score`)) |>
      ungroup() |>
      arrange(g_or_e, `F-score`) |>
      mutate(order = row_number())
    # mynbr <- seq(length(unique(F_score$tool)))
    # names(mynbr) <- unique(F_score$tool)
    # F_score$order <- mynbr[F_score$tool]
  }
  
  
  if (!as) {
    if (is.null(split_variable)) {
      if (length(unique(F_score$L1)) == 1) {
        p <- F_score |>
          ggplot(aes(x = order, y = `F-score`, fill = tool)) +
          geom_col() +
          facet_wrap(.~g_or_e, scale = 'free_x') +
          scale_x_continuous(breaks = F_score$order, labels = F_score$tool)
      } else {
        p <- F_score |>
          ggplot(aes(x = tool, y = `F-score`, fill = tool)) +
          geom_col() +
          facet_wrap(.~g_or_e, scale = 'free_x')
      }
      p <- p +
        scale_fill_manual(values = col_color, guide = "none") +
        # guides(col = guide_legend(nrow = legend_nrow)) +
        ylim(0, 1) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) +
        ggtitle("") +
        xlab(label = NULL)
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    } else {
      if (length(unique(F_score$L1)) == 1) {
        p <- F_score |>
          ggplot(aes(x = order, y = `F-score`, fill = tool)) +
          geom_col() +
          facet_wrap(g_or_e ~ L1, scale = 'free_x') +
          scale_x_continuous(breaks = F_score$order, labels = F_score$tool)
      } else {
        p <- F_score |>
          ggplot(aes(x = tool, y = `F-score`, fill = tool)) +
          geom_col() +
          facet_wrap(g_or_e ~ L1, scale = 'free_x')
      }
      p <- p +
        scale_fill_manual(values = col_color, guide = "none") +
        # guides(col = guide_legend(nrow = legend_nrow)) +
        ylim(0, 1) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) +
        ggtitle("") +
        xlab(label = NULL)
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    }
  } else {
    if (is.null(split_variable)) {
      F_score |>
        ggplot(aes(x = tool, y = `F-score`, fill = tool)) +
        geom_col() +
        # facet_wrap(g_or_e ~ L1, scale = 'free_x') +
        # scale_x_continuous(breaks = F_score$order, labels = F_score$tool) +
        scale_fill_manual(values = col_color, guide = "none") +
        # guides(col = guide_legend(nrow = legend_nrow)) +
        ylim(0, 1) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) +
        ggtitle("") +
        xlab(label = NULL) -> p
      if (long) p <- p + facet_grid(as_type ~ g_or_e)
      if (!long) p <- p + facet_grid(g_or_e ~ as_type)
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    } else {
      F_score |>
        ggplot(aes(x = tool, y = `F-score`, fill = tool)) +
        geom_col() +
        # facet_wrap(g_or_e ~ L1, scale = 'free_x') +
        # scale_x_continuous(breaks = F_score$order, labels = F_score$tool) +
        scale_fill_manual(values = col_color, guide = "none") +
        # guides(col = guide_legend(nrow = legend_nrow)) +
        ylim(0, 1) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = legendtextsize),
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                         size = axistextsize),
              axis.text.y = element_text(size = axistextsize),
              axis.title.x = element_text(size = axistitlesize),
              axis.title.y = element_text(size = axistitlesize),
              strip.text = element_text(size = stripsize),
              strip.background = element_rect(fill = NA, colour = "black")) +
        ggtitle("") +
        xlab(label = NULL) -> p
      if (long) p <- p + facet_grid(L1 + as_type ~ g_or_e)
      if (!long) p <- p + facet_grid(L1 + g_or_e ~ as_type)
      pdf(output_filename, width = fig_width, height = fig_height)
      print(p)
      dev.off()
    }
  }
}

##### extend #####
perGeneQValueExact <- function(pGene, theta, geneSplit) {
  stopifnot(length(pGene) == length(geneSplit))
  numExons <- listLen(geneSplit)
  tab <- tabulate(numExons)
  notZero <- (tab > 0)
  numerator <- mapply(function(m, n) m * (1 - (1 - theta)^n),
                      m = tab[notZero],
                      n = which(notZero)
  )
  numerator <- rowSums(numerator)
  bins <- cut(pGene, breaks = c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts <- tabulate(bins, nbins = nlevels(bins))
  denom <- cumsum(counts)
  stopifnot(denom[length(denom)] == length(pGene))
  return(numerator / denom)
}

qValue <- function(result_tbl, gene_col_name, p_col_name) {
  geneSplit <- split(seq(along = result_tbl[[gene_col_name]]), result_tbl[[gene_col_name]])
  pGene <- sapply(geneSplit, function(i) min(as.numeric(result_tbl[[p_col_name]])[i]))
  theta <- unique(sort(pGene))
  q <- perGeneQValueExact(pGene, theta, geneSplit)
  qres_P <- rep(NA_real_, length(pGene))
  qres_P <- q[match(pGene, theta)]
  qres_P <- pmin(1, qres_P)
  names(qres_P) <- names(geneSplit)
  qres_P
}


