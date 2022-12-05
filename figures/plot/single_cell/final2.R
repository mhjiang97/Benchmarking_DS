source("utils_final.R")

##### colors #####
colors_as <- c(brewer.pal(11, "Set3"), "#D3D3D3")
names(colors_as) <- c(
  "A3SS", "A5SS", "MXE", "RI", "SE", "AFE", "ALE", "SME", "CE", "TE", "TS", "OTHER"
)
colors_as["CE"] <- "#FFED6F"
colors_as["A5SS/A3SS"] <- "#E22826"

colors_lib <- ggsci::pal_npg("nrc")(10)
colors_tool <- colors_lib[c(1:6)]
names(colors_tool) <- c(
  "BRIE2", "DESJ", "MicorExonator", "Outrigger", "Psix", "SingleSplice"
)

##### tx2gene #####
file_tx2gene <- "~/doc/tx2gene/tx2gene_v32.txt.gz"
tx2gene <- vroom(
  file_tx2gene,
  delim = "\t", col_types = cols(),
  col_names = c("chr", "gene_id", "gene_symbol", "transcript_id")
)
tx2gene_convert <- tx2gene |>
  dplyr::distinct(gene_id, .keep_all = T) |>
  dplyr::select(gene_id, gene_symbol)

##### truths #####
if (F) {
  truths <- readRDS("~/projects/AS/data/sim/rsem/truths.rds")
  truths <- truths |>
    dplyr::mutate(
      as_type_former = as_type,
      as_type = case_when(
        as_type_former == "SE" ~ "SE",
        as_type_former == "RI" ~ "RI",
        as_type_former == "MXE" ~ "MXE",
        as_type_former == "A3SS" ~ "A5SS/A3SS",
        as_type_former == "A5SS" ~ "A5SS/A3SS"
      )
    )
  file_truths_generated <- "~/projects/as/data/sim/rsem/revision2/sim_1/truth.txt"
  truths_generated <- vroom(
    file_truths_generated,
    delim = "\t", col_types = paste(rep("c", 25), collapse = "")
  )
  truths_m <- dplyr::full_join(truths_generated, truths, by = c("gene" = "gene_id")) |>
    dplyr::rename(gene_id = gene)
}

truths <- readRDS("~/projects/AS/data/sim/rsem/truths.rds")
truths_m <- readRDS("~/projects/AS/data/sim/rsem/revision2/truths_m.rds")


##### plot #####
dir_result <- "~/projects/as/analysis/revision2/"

tools <- c(
  "BRIE2", "DESJ", "MicroExonator", "Outrigger", "Psix", "SingleSplice"
)

results <- lapply(
  tools, function(x) {
    t <- str_replace_all(x, "\\-", "") |>
      tolower()
    r <- readRDS(glue("{dir_result}/{x}/{t}.rds"))
    r
  }
) |>
  setNames(tools)

results$BRIE2$idx <- results$BRIE2$idx_truths3
results$DESJ$idx <- results$DESJ$idx_truths2
results$MicroExonator$idx <- results$MicroExonator$idx_truths2
results$Outrigger$idx <- results$Outrigger$idx_truths
results$Psix$idx <- results$Psix$idx_truths
results$SingleSplice$idx <- results$SingleSplice$idx_truths2


number_detect <- lapply(
  tools,
  function(x) {
    results[[x]] |>
      dplyr::group_by(as_type) |>
      dplyr::summarise(
        g = dplyr::n_distinct(gene_id), e = dplyr::n_distinct(myID)
      ) |>
      dplyr::mutate(tool = x)
  }
)

number_detect_l <- number_detect |>
  bind_rows() |>
  arrange(desc(e)) |>
  mutate(
    tool = fct_reorder(tool, e),
    as_type = fct_reorder(as_type, e)
  ) |>
  pivot_longer(cols = 2:3, names_to = "g_or_e", values_to = "nbr")

number_detect_l |>
  ggplot(aes(reorder(tool, nbr), nbr, fill = as_type)) +
  geom_col(
    data = . %>% filter(g_or_e == "g"),
    position = position_dodge2(width = 2.2, padding = 0.2, preserve = "single"),
    alpha = 1, color = "black", size = 0.15, width = 1
  ) +
  geom_col(
    data = . %>% filter(g_or_e == "e"),
    position = position_dodge2(width = 2.2, padding = 0.2, preserve = "single"),
    alpha = 0.3, color = "black", size = 0.15, width = 1
  ) +
  scale_alpha_manual(values = c(1, 0.3)) +
  scale_fill_manual(values = colors_as) +
  scale_y_break(c(5000, 15000), scales = 0.2) +
  # scale_y_break(c(50000, 60000), scales = 0.1) +
  theme_light() +
  theme(
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  ) +
  labs(fill = "Types of alternative splicing") +
  xlab("Tools") +
  ylab("Number of genes / events") -> p
pdf(glue("{dir_result}/figs/nbr.pdf"), width = 8, height = 5, onefile = F)
print(p)
dev.off()

write.xlsx(bind_rows(number_detect), glue("{dir_result}/tables/nbr.xlsx"))


percent_detection <- lapply(
  names(results),
  function(x) {
    results[[x]] |>
      # group_by(as_type) |>
      summarise(
        g = length(unique(gene_id[as.numeric(padj) < 0.05])) / n_distinct(gene_id),
        e = length(unique(myID[as.numeric(padj) < 0.05])) / n_distinct(myID)
      ) |>
      mutate(tool = x)
  }
)

percent_detection_l <- percent_detection |>
  bind_rows() |>
  pivot_longer(cols = 1:2, names_to = "g_or_e", values_to = "per")

percent_detection_l |>
  ggplot(aes(reorder(tool, per, FUN = max), per, fill = g_or_e)) +
  geom_col(position = position_dodge2(), color = "black", size = 0.2) +
  scale_fill_brewer(type = "qual", palette = "Paired", direction = -1) +
  theme_light() +
  scale_y_continuous(labels = percent) +
  # scale_y_break(c(0.1287, 0.3035), scales = 0.2) + # 0.1003
  # scale_y_break(c(50000, 60000), scales = 0.1) +
  theme(
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  ) +
  labs(fill = "gene or event") +
  xlab("Tools") +
  ylab("Percent") -> p
pdf(glue("{dir_result}/figs/per_detection.pdf"), width = 6, height = 5, onefile = F)
print(p)
dev.off()

write.xlsx(bind_rows(percent_detection), glue("{dir_result}/tables/per_detection.xlsx"))


percent_detection_as <- lapply(
  names(results),
  function(x) {
    results[[x]] |>
      group_by(as_type) |>
      summarise(
        e = length(unique(myID[as.numeric(padj) < 0.05])) / n_distinct(myID),
        g = length(unique(gene_id[as.numeric(padj) < 0.05])) / n_distinct(gene_id)
      ) |>
      mutate(tool = x)
  }
)

percent_detection_as_l <- percent_detection_as |>
  bind_rows() |>
  mutate(
    tool = fct_reorder(tool, g),
    as_type = fct_reorder(as_type, g)
  ) |>
  pivot_longer(cols = 2:3, names_to = "g_or_e", values_to = "per")

percent_detection_as_l |>
  ggplot(aes(reorder(tool, per), per, fill = as_type)) +
  geom_col(
    data = . %>% filter(g_or_e == "e"),
    position = position_dodge2(width = 2.2, padding = 0.2, preserve = "single"),
    alpha = 1, color = "black", size = 0.15, width = 1
  ) +
  geom_col(
    data = . %>% filter(g_or_e == "g"),
    position = position_dodge2(width = 2.2, padding = 0.2, preserve = "single"),
    alpha = 0.3, color = "black", size = 0.15, width = 1
  ) +
  scale_alpha_manual(values = c(1, 0.3)) +
  scale_fill_manual(values = colors_as) +
  # scale_y_break(c(0.125, 0.25), scales = 0.2) +
  # scale_y_break(c(0.5, 0.6), scales = 0.1) +
  theme_light() +
  scale_y_continuous(labels = percent) +
  theme(
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  ) +
  labs(fill = "Types of alternative splicing") +
  xlab("Tools") +
  ylab("Percent of genes / events") -> p
pdf(glue("{dir_result}/figs/per_detection_as.pdf"), width = 8, height = 5, onefile = F)
print(p)
dev.off()

write.xlsx(bind_rows(percent_detection_as), glue("{dir_result}/tables/per_detection_as.xlsx"))


percent_detection_tf <- lapply(
  names(results),
  function(x) {
    results[[x]] |>
      # group_by(as_type) |>
      summarise(
        t = n_distinct(intersect(gene_id[as.numeric(padj) < 0.05], truths$gene_id)) / n_distinct(gene_id[as.numeric(padj) < 0.05]),
        f = n_distinct(setdiff(gene_id[as.numeric(padj) < 0.05], truths$gene_id)) / n_distinct(gene_id[as.numeric(padj) < 0.05])
      ) |>
      mutate(
        tool = x,
        g_or_e = "g"
      )
  }
)

percent_detection_tf_l <- percent_detection_tf |>
  bind_rows() |>
  mutate(tool = fct_reorder(tool, t)) |>
  pivot_longer(cols = 1:2, names_to = "t_or_f", values_to = "per")

percent_detection_tf_l |>
  ggplot(aes(tool, per, fill = t_or_f)) +
  geom_col(position = position_dodge2(padding = 0.05), color = "black", size = 0.2, width = 0.65) +
  scale_fill_manual(values = c("t" = "#00468BFF", "f" = "#ED0000FF")) +
  theme_light() +
  scale_y_continuous(labels = percent) +
  theme(
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  ) +
  labs(fill = "Ture or false") +
  xlab("Tools") +
  ylab("Percent (gene)") -> p
pdf(glue("{dir_result}/figs/per_detection_tf.pdf"), width = 7, height = 5, onefile = F)
print(p)
dev.off()

percent_detection_tf_e <- lapply(
  tools,
  function(x) {
    tmp <- results[[x]]
    keep_as <- intersect(c("SE", "MXE", "A3SS", "A5SS", "RI", "A5SS/A3SS"), unique(results[[x]]$as_type))
    if (length(keep_as) > 0) tmp <- tmp |> filter(as_type %in% keep_as)

    tmp |>
      # group_by(as_type) |>
      summarise(
        t = n_distinct(myID[as.numeric(padj) < 0.05 & !is.na(idx)]) / n_distinct(myID[as.numeric(padj) < 0.05]),
        f = n_distinct(myID[as.numeric(padj) < 0.05 & is.na(idx)]) / n_distinct(myID[as.numeric(padj) < 0.05])
      ) |>
      mutate(
        tool = x,
        g_or_e = "e"
      )
  }
)

percent_detection_tf_e_l <- percent_detection_tf_e |>
  bind_rows() |>
  mutate(tool = fct_reorder(tool, t)) |>
  pivot_longer(cols = 1:2, names_to = "t_or_f", values_to = "per")

percent_detection_tf_e_l |>
  ggplot(aes(tool, per, fill = t_or_f)) +
  geom_col(position = position_dodge2(padding = 0.05), color = "black", size = 0.2, width = 0.65) +
  scale_fill_manual(values = c("t" = "#00468BFF", "f" = "#ED0000FF")) +
  theme_light() +
  scale_y_continuous(labels = percent) +
  theme(
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  ) +
  labs(fill = "Ture or false") +
  xlab("Tools") +
  ylab("Percent (event)") -> p
pdf(glue("{dir_result}/figs/per_detection_tf_e.pdf"), width = 5.5, height = 5, onefile = F)
print(p)
dev.off()


##### tables #####
lapply(
  names(results),
  function(x) {
    results[[x]] |>
      # group_by(as_type) |>
      summarise(
        g = length(unique(gene_id[as.numeric(padj) < 0.05])),
        e = length(unique(myID[as.numeric(padj) < 0.05]))
      ) |>
      mutate(tool = x)
  }
) |> setNames(names(results)) -> diff_nbr

tbl_diff_nbr <- diff_nbr |> bind_rows()
write.xlsx(tbl_diff_nbr, glue("{dir_result}/tables/nbr_diff.xlsx"))


##### as_type #####
results$MicroExonator <- results$MicroExonator |>
  mutate(
    as_type_former = as_type,
    as_type = case_when(
      as_type_former == "AFE" ~ "AFE",
      as_type_former == "ALE" ~ "ALE",
      as_type_former == "CE" ~ "CE",
      as_type_former == "RI" ~ "RI",
      as_type_former == "TE" ~ "TE",
      as_type_former == "TS" ~ "TS",
      as_type_former == "MXE" ~ "MXE",
      as_type_former == "A3SS" ~ "A5SS/A3SS",
      as_type_former == "A5SS" ~ "A5SS/A3SS"
    )
  )


##### all #####
if (F) {
  truths_m$all <- "(n = 5723, n.ds = 1000, n.dse = 1723)"
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "all", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/all.pdf"),
  pathsize = 1, pointsize = 2.5, stripsize = 12, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 6, fig_width = 10, legendtextsize = 8
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "all", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/all_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 12,
  legend_nrow = 1, fig_height = 6, fig_width = 15, long = F
)


plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "all", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_all.pdf"),
  stripsize = 12, axistitlesize = 20, axistextsize = 10,
  fig_height = 8, fig_width = 10
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "all", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_all_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 8, fig_width = 13
)


##### isoform percent difference #####
if (F) {
  truths_m$diff_IsoPct_3 <- cut2(as.double(truths_m$diff_IsoPct), cuts = c(0, 1 / 3, 2 / 3, 1))
  truths_m$diff_IsoPct3 <- sapply(truths_m$diff_IsoPct_3, function(i) {
    paste0(i, " (n = ", length(which(truths_m$diff_IsoPct_3 == i)), ")")
  })
  truths_m$diff_IsoPct3_2 <- sapply(truths_m$diff_IsoPct_3, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$diff_IsoPct_3 == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$diff_IsoPct_3 == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$diff_IsoPct_3 == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "diff_IsoPct3_2", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_per_diff.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "diff_IsoPct3", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_per_diff_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)


plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "diff_IsoPct3", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_iso_per_diff.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "diff_IsoPct3", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_iso_per_diff_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### gene expression level #####
if (F) {
  truths_m$exprlevel <- "low"
  truths_m$exprlevel[as.double(truths_m$TPM) > median(as.double(truths_m$TPM)[truths_m$ds_status == 1])] <- "high"
  truths_m$exprlevel <- factor(truths_m$exprlevel, levels = c("low", "high"))
  truths_m$exprlevel2 <- sapply(truths_m$exprlevel, function(i) {
    paste0(i, " (n = ", length(which(truths_m$exprlevel == i)), ")")
  })
  truths_m$exprlevel2_2 <- sapply(truths_m$exprlevel, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$exprlevel == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$exprlevel == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$exprlevel == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "exprlevel2_2", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/gene_exp.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 9, fig_width = 10, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "exprlevel2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/gene_exp_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 9, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 9, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "exprlevel2_2", pch_tool = NULL, filter_as = F, as = F,
  output_filename = glue("{dir_result}/figs/bar_gene_exp.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "exprlevel2_2", pch_tool = NULL, filter_as = F, as = T,
  output_filename = glue("{dir_result}/figs/bar_gene_exp_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 5% #####
if (F) {
  truths_m$nbr_isoforms_a5 <- cut2(
    as.numeric(truths_m$nbr_isoforms_above5),
    cuts = c(2, 3, max(as.numeric(truths_m$nbr_isoforms_above5))),
    oneval = F
  )
  truths_m$nbr_isoforms_a5_2 <- sapply(truths_m$nbr_isoforms_a5, function(i) {
    paste0(i, " (n = ", length(which(truths_m$nbr_isoforms_a5 == i)), ")")
  })
  truths_m$nbr_isoforms_a5_3 <- sapply(truths_m$nbr_isoforms_a5, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$nbr_isoforms_a5 == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$nbr_isoforms_a5 == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$nbr_isoforms_a5 == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a5_3", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_5.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a5_2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_5_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a5_2", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_iso_over_5.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a5_2", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_iso_over_5_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 10% #####
if (F) {
  truths_m$nbr_isoforms_a10 <- cut2(
    as.numeric(truths_m$nbr_isoforms_above10),
    cuts = c(2, 3, max(as.numeric(truths_m$nbr_isoforms_above10))),
    oneval = F
  )
  truths_m$nbr_isoforms_a10_2 <- sapply(truths_m$nbr_isoforms_a10, function(i) {
    paste0(i, " (n = ", length(which(truths_m$nbr_isoforms_a10 == i)), ")")
  })
  truths_m$nbr_isoforms_a10_3 <- sapply(truths_m$nbr_isoforms_a10, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$nbr_isoforms_a10 == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$nbr_isoforms_a10 == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$nbr_isoforms_a10 == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a10_3", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_10.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a10_2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_10_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a10_2", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_iso_over_10.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a10_2", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_iso_over_10_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 15% #####
if (F) {
  truths_m$nbr_isoforms_a15 <- cut2(
    as.numeric(truths_m$nbr_isoforms_above15),
    cuts = c(2, 3, max(as.numeric(truths_m$nbr_isoforms_above15))),
    oneval = F
  )
  truths_m$nbr_isoforms_a15_2 <- sapply(truths_m$nbr_isoforms_a15, function(i) {
    paste0(i, " (n = ", length(which(truths_m$nbr_isoforms_a15 == i)), ")")
  })
  truths_m$nbr_isoforms_a15_3 <- sapply(truths_m$nbr_isoforms_a15, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$nbr_isoforms_a15 == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$nbr_isoforms_a15 == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$nbr_isoforms_a15 == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a15_3", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_15.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a15_2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_15_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a15_2", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_iso_over_15.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a15_2", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_iso_over_15_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 25% #####
if (F) {
  truths_m$nbr_isoforms_a25 <- cut2(
    as.numeric(truths_m$nbr_isoforms_above25),
    cuts = c(2, 3, max(as.numeric(truths_m$nbr_isoforms_above25))),
    oneval = F
  )
  truths_m$nbr_isoforms_a25_2 <- sapply(truths_m$nbr_isoforms_a25, function(i) {
    paste0(i, " (n = ", length(which(truths_m$nbr_isoforms_a25 == i)), ")")
  })
  truths_m$nbr_isoforms_a25_3 <- sapply(truths_m$nbr_isoforms_a25, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$nbr_isoforms_a25 == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$nbr_isoforms_a25 == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$nbr_isoforms_a25 == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a25_3", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_25.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_a25_2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_over_25_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a25_2", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_iso_over_25.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_a25_2", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_iso_over_25_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### isoform number #####
if (F) {
  truths_m$nbr_isoforms_cat <- cut2(as.numeric(truths_m$nbr_isoforms), cuts = c(1, 4, 10, max(as.numeric(truths_m$nbr_isoforms))))
  truths_m$nbr_isoforms_cat2 <- sapply(truths_m$nbr_isoforms_cat, function(i) {
    paste0(i, " (n = ", length(which(truths_m$nbr_isoforms_cat == i)), ")")
  })
  truths_m$nbr_isoforms_cat3 <- sapply(truths_m$nbr_isoforms_cat, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$nbr_isoforms_cat == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$nbr_isoforms_cat == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$nbr_isoforms_cat == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_cat3", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_nbr.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 10, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "nbr_isoforms_cat2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_nbr_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_cat2", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_iso_nbr.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "nbr_isoforms_cat2", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_iso_nbr_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### number of differing bp between diff used isoforms #####
if (F) {
  truths_m$diffbp_3 <- cut2(as.numeric(truths_m$diffbp), cuts = c(100, 1000, max(as.numeric(truths_m$diffbp), na.rm = T)))
  truths_m$diffbp3_2 <- sapply(truths_m$diffbp_3, function(i) {
    paste0(i, " (n = ", length(which(truths_m$diffbp_3 == i)), ")")
  })
  truths_m$diffbp4 <- sapply(truths_m$diffbp3_2, function(i) {
    paste0(
      i, " (n = ", length(which(truths_m$diffbp3_2 == i)),
      ", n.ds = ", length(intersect(truths_m$gene_id[which(truths_m$diffbp3_2 == i)], truths_m$gene_id[which(truths_m$ds_status == 1)])),
      ", n.dse = ", length(intersect(which(truths_m$diffbp3_2 == i), which(truths_m$ds_status == 1))), ")"
    )
  })
}

plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "diffbp4", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/diffbp.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 9, axistitlesize = 20, axistextsize = 15,
  legend_nrow = 1, fig_height = 8, fig_width = 12, long = F
)
plot_fdr_tpr_paper(
  results, truths_m, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "diffbp3_2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/diffbp_as.pdf"),
  pathsize = 1, pointsize = 2, stripsize = 8, axistitlesize = 20, axistextsize = 10,
  legend_nrow = 1, fig_height = 11, fig_width = 15, long = F
)

plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "diffbp3_2", pch_tool = NULL, filter_as = T, as = F,
  output_filename = glue("{dir_result}/figs/bar_diffbp.pdf"),
  stripsize = 10, axistitlesize = 10, axistextsize = 10,
  fig_height = 10, fig_width = 10, long = F
)
plot_f_score(
  results, truths_m, tools, colors_tool,
  thresholds = 0.05,
  split_variable = "diffbp3_2", pch_tool = NULL, filter_as = T, as = T,
  output_filename = glue("{dir_result}/figs/bar_diffbp_as.pdf"),
  stripsize = 8, axistitlesize = 20, axistextsize = 10,
  fig_height = 10, fig_width = 13
)


##### number of isoforms envolved in AS #####
if (F) {
  truths_m <- truths_m |>
    rowwise() |>
    mutate(
      iso_env = length(unlist(str_split(X10, "[\\,\\/]")))
    ) |>
    ungroup()
}

truths_m_fil <- truths_m |>
  filter(iso_env != 1)
truths_m_fil$iso_env2 <- cut2(truths_m_fil$iso_env, cuts = c(3, 5, max(truths_m_fil$iso_env, na.rm = T)))
truths_m_fil$iso_env3 <- sapply(truths_m_fil$iso_env2, function(i) {
  paste0(
    i, " (n = ", length(which(truths_m_fil$iso_env2 == i)),
    ", n.ds = ", length(intersect(truths_m_fil$gene_id[which(truths_m_fil$iso_env2 == i)], truths_m_fil$gene_id[which(truths_m_fil$ds_status == 1)])),
    ", n.dse = ", length(intersect(which(truths_m_fil$iso_env2 == i), which(truths_m_fil$ds_status == 1))), ")"
  )
})

plot_tpr_paper(
  results, truths_m_fil, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "iso_env3", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_env.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 9, axistitlesize = 15, axistextsize = 12,
  legend_nrow = 1, fig_height = 9, fig_width = 12, long = F
)
plot_tpr_paper(
  results, truths_m_fil, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "iso_env2", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_env_as.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 12, axistitlesize = 12, axistextsize = 12,
  legend_nrow = 1, fig_height = 13, fig_width = 15, long = F
)


##### psi difference #####
if (F) {
  mutate_IsoPct <- function(w, ref) {
    o <- order(ref, decreasing = T)[1:2]
    w[rev(o)] <- w[o]
    w
  }
  tmp <- read.delim(
    "~/projects/AS/data/sim/rsem/revision2/sim_1/simulation_details.txt",
    header = T, as.is = T
  )
  tmp_ds <- subset(tmp, gene_ds_status == 1)
  tmp_nonds <- subset(tmp, gene_ds_status == 0)

  tmp_ds <- tmp_ds %>%
    group_by(gene_id) %>%
    mutate(IsoPct2 = mutate_IsoPct(IsoPct, IsoPct))
  tmp_nonds$IsoPct2 <- tmp_nonds$IsoPct

  tmp <- rbind(tmp_ds, tmp_nonds)
  tmp <- tmp[order(tmp$gene_id), ]
  tmp <- ungroup(tmp)

  truths_m$diff_psi <- NA
  for (i in 1:nrow(truths_m)) {
    if (is.na(truths_m$X10[i])) next
    isos1 <- str_split(truths_m$X10[i], ",")[[1]][1] |>
      str_split("/") |>
      unlist()
    isos2 <- str_split(truths_m$X10[i], ",")[[1]][2] |>
      str_split("/") |>
      unlist()
    psi1 <- sum(tmp$IsoPct[tmp$transcript_id %in% isos1]) / sum(c(tmp$IsoPct[tmp$transcript_id %in% isos1], tmp$IsoPct[tmp$transcript_id %in% isos2]))
    psi2 <- sum(tmp$IsoPct2[tmp$transcript_id %in% isos1]) / sum(c(tmp$IsoPct2[tmp$transcript_id %in% isos1], tmp$IsoPct2[tmp$transcript_id %in% isos2]))
    truths_m$diff_psi[i] <- abs(psi1 - psi2)
  }
  # truths$ds_status <- 1
}

truths_m_fil <- truths_m[!is.na(truths_m$diff_psi), ]

truths_m_fil$diff_psi_cut <- cut2(truths_m_fil$diff_psi, g = 3)
truths_m_fil$diff_psi_cut2 <- sapply(truths_m_fil$diff_psi_cut, function(i) {
  paste0(
    i, " (n = ", length(which(truths_m_fil$diff_psi_cut == i)),
    ", n.ds = ", length(intersect(truths_m_fil$gene_id[which(truths_m_fil$diff_psi_cut == i)], truths_m_fil$gene_id[which(truths_m_fil$ds_status == 1)])),
    ", n.dse = ", length(intersect(which(truths_m_fil$diff_psi_cut == i), which(truths_m_fil$ds_status == 1))), ")"
  )
})

plot_tpr_paper(
  results, truths_m_fil, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "diff_psi_cut2", pch_tool = NULL, filter_as = T, as = F, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_diff_psi.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 9, axistitlesize = 15, axistextsize = 12,
  legend_nrow = 1, fig_height = 9, fig_width = 12, long = F
)
plot_tpr_paper(
  results, truths_m_fil, tools, colors_tool,
  thresholds = c(0.01, 0.05, 0.1),
  split_variable = "diff_psi_cut", pch_tool = NULL, filter_as = T, as = T, new_plot_type = F,
  output_filename = glue("{dir_result}/figs/iso_diff_psi_as.pdf"),
  pathsize = 1, pointsize = 3, stripsize = 12, axistitlesize = 12, axistextsize = 12,
  legend_nrow = 1, fig_height = 13, fig_width = 15, long = F
)


##### combine each two tools #####
pairs_tool <- matrix(
  nrow = length(results), ncol = length(results)
) |>
  as.data.frame()
colnames(pairs_tool) <- rownames(pairs_tool) <- names(results)

for (t1 in rownames(pairs_tool)) {
  for (t2 in colnames(pairs_tool)) {
    r1 <- results[[t1]]
    r2 <- results[[t2]]
    g <- c(r1$gene_id[as.numeric(r1$padj) <= 0.05], r2$gene_id[as.numeric(r2$padj) <= 0.05]) |>
      unique()
    
    tp <- length(intersect(truths_m$gene_id[truths_m$ds_status == 1], g))
    p <- length(unique(g))
    pp <- length(unique(truths_m$gene_id[truths_m$ds_status == 1]))
    
    fdr <- round(1 - tp / p, 5)
    tpr <- round(tp / pp, 5)
    f <- (1 + 0.5^2) * (((1 - fdr) * tpr) / (0.5^2 * (1 - fdr) + tpr))
    
    tmp <- glue(
      "FDR: {fdr}\nTPR: {tpr}\nF-socre: {round(f, 5)}"
    )
    pairs_tool[t1, t2] <- tmp
  }
}

pairs_tool <- pairs_tool |>
  rownames_to_column("tool")
write.xlsx(pairs_tool, glue("{dir_result}/tables/pairs.xlsx"))


pairs_tool <- matrix(
  nrow = length(results), ncol = length(results)
) |>
  as.data.frame()
colnames(pairs_tool) <- rownames(pairs_tool) <- names(results)

for (t1 in rownames(pairs_tool)) {
  for (t2 in colnames(pairs_tool)) {
    r1 <- results[[t1]]
    r2 <- results[[t2]]
    g <- intersect(r1$gene_id[as.numeric(r1$padj) <= 0.05], r2$gene_id[as.numeric(r2$padj) <= 0.05]) |>
      unique() |>
      na.omit()
    
    tp <- length(intersect(truths_m$gene_id[truths_m$ds_status == 1], g))
    p <- length(unique(g))
    pp <- length(unique(truths_m$gene_id[truths_m$ds_status == 1]))
    
    fdr <- round(1 - tp / p, 5)
    tpr <- round(tp / pp, 5)
    f <- (1 + 0.5^2) * (((1 - fdr) * tpr) / (0.5^2 * (1 - fdr) + tpr))
    
    tmp <- glue(
      "FDR: {fdr}\nTPR: {tpr}\nF-score: {round(f, 5)}"
    )
    pairs_tool[t1, t2] <- tmp
  }
}

pairs_tool <- pairs_tool |>
  rownames_to_column("tool")
write.xlsx(pairs_tool, glue("{dir_result}/tables/pairs.intersect.xlsx"))


##### save truths #####
if (F) saveRDS(truths_m, "~/projects/AS/data/sim/rsem/revision2/truths_m.rds")


##### extend #####
ds_gene_all <- sapply(tools[tools != "SingleSplice"], function(x) {
  unique(results[[x]]$gene_id[as.numeric(results[[x]]$padj) <= 0.05]) # |> head(500)
}) |>
  flatten_chr() |>
  unique()
dat_ds_gene <- data.frame(matrix(0, nrow = length(ds_gene_all), ncol = length(tools[tools != "SingleSplice"])))
colnames(dat_ds_gene) <- tools[tools != "SingleSplice"]
dat_ds_gene$ds_gene <- ds_gene_all
dat_ds_gene <- dplyr::select(dat_ds_gene, ds_gene, everything())
for (t in tools[tools != "SingleSplice"]) {
  ds_tbl <- table(results[[t]]$gene_id) %>%
    {
      .[match(dat_ds_gene$ds_gene, names(.))] -> .
      .
    } %>%
    as.numeric()
  dat_ds_gene[, t] <- ds_tbl
}
for (i in 2:ncol(dat_ds_gene)) {
  dat_ds_gene[is.na(dat_ds_gene[, i]), i] <- 0
  dat_ds_gene[dat_ds_gene[, i] != 0, i] <- 1
}

cor_jaccrd <- 1 - as.matrix(ecodist::distance(t(dat_ds_gene[, -1]), method = "jaccard"))
pdf(glue("{dir_result}/figs/cor_jaccard_ds_gene.pdf"), width = 10, height = 10)
corrplot(as.matrix(cor_jaccrd), method = "ellipse", type = "lower", tl.col = "black", order = "AOE", addCoef.col = "black", number.digits = 3, tl.cex = 0.8)
dev.off()

q_inte <- lapply(tools[!tools %in% c("SingleSplice", "Outrigger")], function(x) {
  tmp <- results[[x]] |> filter(as.numeric(padj) <= 0.05)
  qValue(tmp, "gene_id", "padj")
}) |>
  setNames(tools[!tools %in% c("SingleSplice", "Outrigger")])

ds_genes <- sapply(q_inte, function(x) {
  names(x)
}) |>
  flatten_chr() |>
  unique()
dat_q <- data.frame(matrix(0, nrow = length(ds_genes), ncol = length(q_inte)), row.names = ds_genes)
colnames(dat_q) <- names(q_inte)

for (t in names(q_inte)) {
  tmp_vec <- q_inte[[t]][match(rownames(dat_q), names(q_inte[[t]]))]
  dat_q[, t] <- tmp_vec
}

cor_spearman <- cor(dat_q, method = "spearman", use = "pairwise.complete.obs") %>% as.data.frame()
for (i in 1:ncol(cor_spearman)) {
  cor_spearman[is.na(cor_spearman[, i]), i] <- 0
}
pdf(glue("{dir_result}/figs/cor_spearman_all_gene.pdf"), width = 10, height = 10)
corrplot(as.matrix(cor_spearman), method = "ellipse", type = "lower", tl.col = "black", order = "AOE", addCoef.col = "black", number.digits = 3, tl.cex = 0.8)
dev.off()


for (at in c("SE", "RI", "A5SS/A3SS", "MXE")) {
  pairs_tool <- matrix(
    nrow = length(tools[!tools %in% c("DESJ", "Psix")]), ncol = length(tools[!tools %in% c("DESJ", "Psix")])
  ) |>
    as.data.frame()
  colnames(pairs_tool) <- rownames(pairs_tool) <- tools[!tools %in% c("DESJ", "Psix")]
  
  for (t1 in rownames(pairs_tool)) {
    for (t2 in colnames(pairs_tool)) {
      r1 <- results[[t1]] |> filter(as_type == at)
      r2 <- results[[t2]] |> filter(as_type == at)
      g <- c(r1$gene_id[as.numeric(r1$padj) <= 0.05], r2$gene_id[as.numeric(r2$padj) <= 0.05]) |>
        unique() |>
        na.omit()
      
      tp <- length(intersect(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == at], g))
      p <- length(unique(g))
      pp <- length(unique(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == at]))
      
      fdr <- round(1 - tp / p, 5)
      tpr <- round(tp / pp, 5)
      f <- (1 + 0.5^2) * (((1 - fdr) * tpr) / (0.5^2 * (1 - fdr) + tpr))
      
      tmp <- glue(
        "FDR: {fdr}\nTPR: {tpr}\nF-socre: {round(f, 5)}"
      )
      pairs_tool[t1, t2] <- tmp
    }
  }
  
  pairs_tool <- pairs_tool |>
    rownames_to_column("tool")
  write.xlsx(pairs_tool, glue('{dir_result}/tables/pairs.{str_replace_all(at, "/", "_")}.xlsx'))
}


for (at in c("SE", "RI", "A5SS/A3SS", "MXE")) {
  pairs_tool <- pairs_tool_f <- matrix(
    nrow = length(tools[!tools %in% c("DESJ", "Psix")]), ncol = length(tools[!tools %in% c("DESJ", "Psix")])
  ) |>
    as.data.frame()
  colnames(pairs_tool) <- rownames(pairs_tool) <-
    colnames(pairs_tool_f) <- rownames(pairs_tool_f) <- tools[!tools %in% c("DESJ", "Psix")]
  
  for (t1 in rownames(pairs_tool)) {
    for (t2 in colnames(pairs_tool)) {
      r1 <- results[[t1]] |> filter(as_type == at)
      r2 <- results[[t2]] |> filter(as_type == at)
      g <- intersect(r1$gene_id[as.numeric(r1$padj) <= 0.05], r2$gene_id[as.numeric(r2$padj) <= 0.05]) |>
        unique() |>
        na.omit()
      
      tp <- length(intersect(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == at], g))
      p <- length(unique(g))
      pp <- length(unique(truths_m$gene_id[truths_m$ds_status == 1 & truths_m$as_type == at]))
      
      fdr <- round(1 - tp / p, 5)
      tpr <- round(tp / pp, 5)
      f <- (1 + 0.5^2) * (((1 - fdr) * tpr) / (0.5^2 * (1 - fdr) + tpr))
      
      tmp <- glue(
        "FDR: {fdr}\nTPR: {tpr}\nF-socre: {round(f, 5)}"
      )
      pairs_tool[t1, t2] <- tmp
    }
  }
  
  pairs_tool <- pairs_tool |>
    rownames_to_column("tool")
  write.xlsx(pairs_tool, glue('{dir_result}/tables/pairs.intersect.{str_replace_all(at, "/", "_")}.xlsx'))
}


# for (at in c("SE", "RI", "A5SS/A3SS", "MXE")) {
#   ds_gene_all <- sapply(tools[!tools %in% c("DESJ", "Psix", "SingleSplice")], function(x) {
#     unique(results[[x]]$gene_id[as.numeric(results[[x]]$padj) <= 0.05 & results[[x]]$as_type == at])
#   }) |>
#     flatten_chr() |>
#     unique() |>
#     na.omit()
#
#   dat_ds_gene <- data.frame(matrix(0, nrow = length(ds_gene_all), ncol = length(tools[!tools %in% c("DESJ", "Psix", "SingleSplice")])))
#   colnames(dat_ds_gene) <- tools[!tools %in% c("DESJ", "Psix", "SingleSplice")]
#   dat_ds_gene$ds_gene <- ds_gene_all
#   dat_ds_gene <- dplyr::select(dat_ds_gene, ds_gene, everything())
#
#   for (t in tools[!tools %in% c("DESJ", "Psix", "SingleSplice")]) {
#     ds_tbl <- table(results[[t]]$gene_id[which(results[[t]]$as_type == at)]) %>%
#       {
#         .[match(dat_ds_gene$ds_gene, names(.))] -> .
#         .
#       } %>%
#       as.numeric()
#     dat_ds_gene[, t] <- ds_tbl
#   }
#   for (i in 2:ncol(dat_ds_gene)) {
#     dat_ds_gene[is.na(dat_ds_gene[, i]), i] <- 0
#     dat_ds_gene[dat_ds_gene[, i] != 0, i] <- 1
#   }
#
#   idx_rm <- vector("list", ncol(dat_ds_gene))
#   for (i in 2:ncol(dat_ds_gene)) {
#     if (all(dat_ds_gene[, i] == 0)) idx_rm[[i]] <- i
#   }
#   if (length(flatten_int(idx_rm)) != 0) dat_ds_gene <- dat_ds_gene[, -(flatten_int(idx_rm))]
#
#   cor_jaccrd <- 1 - as.matrix(ecodist::distance(t(dat_ds_gene[, -1]), method = "jaccard"))
#   cor_jaccrd[is.nan(cor_jaccrd)] <- 0
#
#   pdf(glue('{dir_result}/figs/cor_jaccard_ds_gene.{str_replace_all(at, "/", "_")}.pdf'), width = 10, height = 10)
#   corrplot(as.matrix(cor_jaccrd), method = "ellipse", type = "lower", tl.col = "black", order = "AOE", addCoef.col = "black", number.digits = 3, tl.cex = 0.8)
#   dev.off()
# }
#
#
# for (at in c("SE", "RI", "A5SS/A3SS", "A5SS", "MXE")) {
#   q_inte <- lapply(tools[!tools %in% c("DESJ", "Psix", "SingleSplice")], function(x) {
#     tmp <- results[[x]] |> filter(as.numeric(padj) <= 0.05, as_type == at)
#     if (nrow(tmp) == 0) return(NULL)
#     qValue(tmp, "gene_id", "padj")
#   }) |>
#     setNames(tools[!tools %in% c("DESJ", "Psix", "SingleSplice")])
#
#   ds_genes <- sapply(q_inte, function(x) {names(x)}) |> flatten_chr() |> unique()
#   dat_q <- data.frame(matrix(0, nrow = length(ds_genes), ncol = length(q_inte)), row.names = ds_genes)
#   colnames(dat_q) <- names(q_inte)
#
#   for (t in names(q_inte)) {
#     tmp_vec <- q_inte[[t]][match(rownames(dat_q), names(q_inte[[t]]))]
#     dat_q[, t] <- tmp_vec
#   }
#
#   idx_rm <- vector("list", ncol(dat_q))
#   for (i in 2:ncol(dat_q)) {
#     if (all(is.na(dat_q[, i]))) idx_rm[[i]] <- i
#   }
#   if (length(flatten_int(idx_rm)) != 0) dat_q <- dat_q[, -(flatten_int(idx_rm))]
#
#
#   cor_spearman <- cor(dat_q, method = "spearman", use = "pairwise.complete.obs") %>% as.data.frame()
#   for (i in 1:ncol(cor_spearman)) {
#     cor_spearman[is.na(cor_spearman[, i]), i] <- 0
#   }
#   pdf(glue("{dir_result}/figs/cor_spearman_all_gene.{at}.pdf"), width = 10, height = 10)
#   corrplot(as.matrix(cor_spearman), method = "ellipse", type = "lower", tl.col = "black", order = "AOE", addCoef.col = "black", number.digits = 3, tl.cex = 0.8)
#   dev.off()
# }




