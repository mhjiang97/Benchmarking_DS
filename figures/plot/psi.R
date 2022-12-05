source("~/projects/pkgs/R/asp/R/preliminary/utils_tidy_psi.R")

dir_result <- "~/projects/AS/analysis/sim_1/"

truths <- readRDS("~/projects/as/data/sim/rsem/truths.rds")

sample1 <- read.delim("~/projects/as/data/sim/rsem/sim_1/sample1.txt")

colors_lib <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#F4B3BE",
  "#F4A11D", "#8DC8ED", "#4C6CB0", "#8A1C1B", "#CBCC2B", "#EA644C",
  "#634795", "#005B1D", "#26418A", "#CB8A93", "#F1E404", "#E22826"
)
colors_tool <- colors_lib[c(6, 7, 1, 3, 5, 9, 27, 29, 2, 4, 12, 15, 17, 22, 24, 26, 13, 28)]
names(colors_tool) <- c(
  "rMATS", "CASH", "LeafCutter", "SplAdder", "MAJIQ", "SUPPA", "Whippet", "ASpli",
  "BANDITS", "NBSplice", "IsoformSwitchAnalyzeR", "psichomics", "DIEGO", "VAST-TOOLS",
  "RATS", "dSpliceType", "PSI-Sigma", "JuncBASE"
)

colors_as <- c(brewer.pal(11, "Set3"), "#D3D3D3")
names(colors_as) <- c(
  "A3SS", "A5SS", "MXE", "RI", "SE", "AFE", "ALE", "SME", "CE", "TE", "TS", "OTHER"
)
colors_as["CE"] <- "#FFED6F"
colors_as <- colors_as[c("A3SS", "A5SS", "MXE", "RI", "SE", "OTHER")]

## A3: 1/1+2
## A5: 2/1+2
truths_psi <- truths |>
  rowwise() |>
  mutate(
    psi = sum(sample1$TPM[
      match(unlist(str_split(
        unlist(str_split(X10, ","))[1], "/"
      )), sample1$transcript_id)
    ]) / sum(
      sample1$TPM[
        match(unlist(str_split(
          unlist(str_split(X10, ","))[1], "/"
        )), sample1$transcript_id)
      ],
      sample1$TPM[
        match(unlist(str_split(
          unlist(str_split(X10, ","))[2], "/"
        )), sample1$transcript_id)
      ]
    ),

    psi2 = 1 - psi
  ) |>
  ungroup()

tools <- c(
  "rMATS", "LeafCutter", "SplAdder", "SUPPA", "Whippet", "psichomics"
)

colors_tool <- colors_tool[tools]

psi <- lapply(
  tools, function(x) {
    tmp <- str_replace_all(x, "\\-", "") |>
      tolower()
    readRDS(glue("{dir_result}/{x}/{tmp}_psi.rds"))
  }
) |>
  setNames(tools)

delta_psi <- lapply(tools, function(x) {
  if (x != "LeafCutter") {
    markPsi(psi[[x]], truths_psi)
  } else {
    markPsi(psi[[x]], truths_psi, F)
  }
}) |>
  setNames(tools)

delta_psi_dedu <- lapply(tools, function(x) {
  delta_psi[[x]] |>
    group_by(truths_idx) |>
    mutate(final_delta = min(abs(delta))) |>
    ungroup() |>
    distinct(truths_idx, .keep_all = T) |>
    mutate(tool = x)
}) |>
  setNames(tools)
delta_psi_dedu$LeafCutter$as_type <- "OTHER"

bind_rows(delta_psi_dedu) |>
  ggplot(aes(x = final_delta, color = tool)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Error Rate (simulated psi - ground truth psi)") +
  ylab('Cumulative Frequency F(x)') +
  theme(axis.text.y = element_text(size = 12), aspect.ratio = 1) +
  scale_color_manual(values = colors_tool) -> p
pdf(glue("{dir_result}/figs/error_curve.pdf"), width = 6, height = 6)
print(p)
dev.off()


data_summary <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = T),
      sd = sd(x[[col]], na.rm = T)
    )
  }
  data_sum <- ddply(data, groupnames,
    .fun = summary_func,
    varname
  )
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df <- data_summary(
  bind_rows(delta_psi_dedu), varname = "final_delta", groupnames = c("tool")
)


ggplot(df, aes(x = tool, y = final_delta, fill = tool)) +
  geom_bar(
    stat = "identity", color = "black",
    position = position_dodge()
  ) +
  geom_errorbar(aes(ymin = final_delta, ymax = final_delta + sd),
    width = .2,
    position = position_dodge(.9)
  ) +
  labs(x = NULL, y = "Error Rate (simulated psi - ground truth psi)") +
  theme_classic() +
  scale_fill_manual(values = colors_tool) -> p
pdf(glue("{dir_result}/figs/error_bar.pdf"), width = 6, height = 6)
print(p)
dev.off()



df_as <- data_summary(
  bind_rows(delta_psi_dedu), varname = "final_delta", groupnames = c("tool", "as_type")
)


ggplot(df_as, aes(x = tool, y = final_delta, fill = as_type)) +
  geom_bar(
    stat = "identity", color = "black",
    position = position_dodge() 
  ) +
  geom_errorbar(
    aes(ymin = final_delta, ymax = final_delta + sd),
    width = .2,
    position = position_dodge(.9)
  ) +
  labs(x = NULL, y = "Error Rate (simulated psi - ground truth psi)") +
  theme_classic() +
  scale_fill_manual(values = colors_as) -> p
pdf(glue("{dir_result}/figs/error_bar_as.pdf"), width = 6, height = 6)
print(p)
dev.off()





