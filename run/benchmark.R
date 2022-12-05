library(glue)
library(dplyr)
library(purrr)
library(vroom)
source("utils.R")
##### samples #####
samples_1 <- paste0("sample", 1:4)
samples_2 <- paste0("sample", 5:8)

suffix_bam <- "SortedByCoord.bam"
dir_bam <- "~/projects/AS/data/sim/rsem/sim_1/star/"
bams_1 <- sprintf("%s/%s/%s.%s", dir_bam, samples_1, samples_1, suffix_bam)
bams_2 <- sprintf("%s/%s/%s.%s", dir_bam, samples_2, samples_2, suffix_bam)

dir_bam_tophat <- "~/projects/AS/data/sim/rsem/sim_1/tophat2/"
bams_tophat_1 <- glue("{dir_bam_tophat}/{samples_1}/accepted_hits.bam")
bams_tophat_2 <- glue("{dir_bam_tophat}/{samples_2}/accepted_hits.bam")

suffix_fq <- "fq.gz"
dir_fq <- "~/projects/AS/data/sim/rsem/sim_1/fq/"
fqs_1 <- paste(
  sprintf("%s/%s_R1.%s", dir_fq, samples_1, suffix_fq),
  sprintf("%s/%s_R2.%s", dir_fq, samples_1, suffix_fq),
  sep = ";"
)
fqs_2 <- paste(
  sprintf("%s/%s_R1.%s", dir_fq, samples_2, suffix_fq),
  sprintf("%s/%s_R2.%s", dir_fq, samples_2, suffix_fq),
  sep = ";"
)
dir_salmon <- "~/projects/AS/data/sim/rsem/sim_1/salmon/"

##### generate a sample table #####
sampleTable <- data.frame(
  samples = c(samples_1, samples_2),
  conditions = c(rep("s1", 4), rep("s2", 4)),
  files_bam = c(bams_1, bams_2),
  files_bam_tophat = c(bams_tophat_1, bams_tophat_2),
  files_fq = c(fqs_1, fqs_2),
  dirs_salmon = glue("{dir_salmon}/{c(samples_1, samples_2)}/"),
  read_lengths = rep(101, 8),
  library_types = rep("paired-end", 8),
  strandedness = rep("no", 8)
)

##### required reference and annotation files #####
gtf <- "~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf"
gff <- "~/doc/reference/gtf/gff/gencode.v32.annotation.protein_coding2.gff3"
fa <- "~/doc/reference/fa/GRCh38.primary_assembly.genome.fa"
fa_transcript <- "~/doc/reference/rsem_pa_bowtie2/index.idx.fa"
genome_version <- "hg38"

##### other parameters #####
basics <- getBasics(sampleTable)

dir_out <- "~/projects/AS/analysis/sim_1/"
nproc <- 15
novel <- F
write_log <- T
parallel <- T
np <- 5

##### rMATS #####
myCreatedir(glue("{dir_out}/rMATS"))

variable_read_length <- F
paired_stats <- F
cstat <- 0.05

rL <- table(sampleTable$read_lengths) %>% .[. == max(.)] %>% names() %>% as.numeric()
library_type <- table(sampleTable$library_types) %>% .[. == max(.)] %>% names()
tt <- ifelse(library_type == "paired-end", "paired", "single")
strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") lT <- "fr-unstranded"
if (strand == "yes") lT <- "fr-secondstrand"
if (strand == "reverse") lT <- "fr-firststrand"
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]]
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]]
b1 <- paste(path.expand(bams_1), collapse = ",")
b2 <- paste(path.expand(bams_2), collapse = ",")

cmds <- list()

cmds[["generating input files"]] <- c(
  glue('echo -e "{b1}\\n" >{dir_out}/rMATS/{conditions[1]}.txt'),
  glue('echo -e "{b2}\\n" >{dir_out}/rMATS/{conditions[2]}.txt')
)

if (dir.exists(glue("{dir_out}/rMATS/tmp/"))) cmds[["remove files in tmp dir"]] <- glue("rm -rf {dir_out}/rMATS/tmp/*")

cmds[["running rMATS"]] <- glue(
  "rmats.py \\
  --b1 {dir_out}/rMATS/{conditions[1]}.txt \\
  --b2 {dir_out}/rMATS/{conditions[2]}.txt \\
  --od {dir_out}/rMATS/ \\
  --gtf {gtf} \\
  --readLength {rL} \\
  -t {tt} \\
  --nthread {nproc} \\
  --cstat {cstat} \\
  --tstat {nproc} \\
  --libType {lT} \\
  --tmp {dir_out}/rMATS/tmp/"
)

if (novel) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--novelSS")
if (variable_read_length) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--variable-read-length")
if (paired_stats) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--paired-stats")
if (write_log) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], glue("1>{dir_out}/rMATS/_log 2>&1"))

run_cmds(cmds, run = T)


##### CASH #####
myCreatedir(glue("{dir_out}/CASH/"))

cash <- "~/projects/as/software/cash_v2.2.0/cash.jar"
MergePval <- "G"
ram <- "100g"

strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") SS <- "NONE"
if (strand == "yes") SS <- "F"
if (strand == "reverse") SS <- "R"
SC <- ifelse(novel, "True", "False")
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]]
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]]
b1 <- paste(path.expand(bams_1), collapse = ",")
b2 <- paste(path.expand(bams_2), collapse = ",")

cmds <- list()

cmds[["running CASH"]] <- glue(
  "java -jar -Xmx{ram} -Xms{ram} {cash} \\
  --Case:{conditions[1]} {b1} \\
  --Control:{conditions[2]} {b2} \\
  --MergePval {MergePval} \\
  --GTF {gtf} \\
  --Output {dir_out}/CASH/cash \\
  --SpliceCons {SC} \\
  --StrandSpecific {SS}"
)

if (write_log) cmds[["running CASH"]] <- paste(cmds[["running CASH"]], glue("1>{dir_out}/CASH/_log 2>&1"))

run_cmds(cmds, run = T)


##### LeafCutter #####
myCreatedir(glue('{dir_out}/LeafCutter/{c("ds", "juncs", "cluster")}/'))

method_junc <- "regtools"
exon_file <- "~/doc/leafcutter/gencode_v32_exons.txt.gz"

strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") ss <- "0"
if (strand == "yes") ss <- "2"
if (strand == "reverse") ss <- "1"
script_dir <- "~/projects/as/software/leafcutter"

message("generating groups file used when analyzing differential intron excision")
write.table(
  sampleTable[, c("samples", "conditions")],
  glue("{dir_out}/LeafCutter/ds/groups_file.txt"),
  row.names = F, col.names = F, sep = "\t", quote = F
)

cmds <- list()

cmds[["converting bam to junc"]] <- glue(
  "regtools junctions extract -s {ss} {sampleTable$files_bam} \\
  >{dir_out}/LeafCutter/juncs/{sampleTable$samples}.junc && \\
  echo {dir_out}/LeafCutter/juncs/{sampleTable$samples}.junc \\
  >>{dir_out}/LeafCutter/juncs/juncfiles.txt"
)
if (parallel) cmds[["converting bam to junc"]] <- myParallel(cmds[["converting bam to junc"]])

cmds[["intron clustering"]] <- glue(
  "python {script_dir}/clustering/leafcutter_cluster_regtools.py \\
  -j {dir_out}/LeafCutter/juncs/juncfiles.txt \\
  -r {dir_out}/LeafCutter/cluster/ \\
  -o leafcutter \\
  -s \\
  -k"
)

cmds[["analyzing differential intron excision"]] <- glue(
  "{script_dir}/scripts/leafcutter_ds.R --num_threads {nproc} \\
  --exon_file {exon_file} \\
  {dir_out}/LeafCutter/cluster/leafcutter_perind_numers.counts.gz \\
  {dir_out}/LeafCutter/ds/groups_file.txt \\
  -o {dir_out}/LeafCutter/ds/leafcutter"
)

if (min(table(sampleTable$conditions)) < 5) cmds[["analyzing differential intron excision"]] <- paste(cmds[["analyzing differential intron excision"]], glue("-i {min(table(sampleTable$conditions))}"))
if (min(table(sampleTable$conditions)) < 3) cmds[["analyzing differential intron excision"]] <- paste(cmds[["analyzing differential intron excision"]], glue("-g {min(table(sampleTable$conditions))}"))

if (write_log) {
  cmds[["intron clustering"]] <- paste(
    cmds[["intron clustering"]],
    glue("1>{dir_out}/LeafCutter/cluster/_log 2>&1")
  )
  cmds[["analyzing differential intron excision"]] <- paste(
    cmds[["analyzing differential intron excision"]],
    glue("1>{dir_out}/LeafCutter/ds/_log 2>&1")
  )
}

run_cmds(cmds, run = T)


##### SplAdder #####
myCreatedir(glue("{dir_out}/SplAdder/graphs/"))

rl <- table(sampleTable$read_lengths) %>% .[. == max(.)] %>% names() %>% as.numeric()
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]]
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]]
b1 <- paste(path.expand(bams_1), collapse = ",")
b2 <- paste(path.expand(bams_2), collapse = ",")

cmds <- list()

cmds[["building graphs"]] <- glue(
  "spladder build \\
  -o {dir_out}/SplAdder/graphs/ \\
  --readlen {rl} \\
  -a {gtf} \\
  -b {paste(b1, b2, sep = ',')}"
)

if (!novel) cmds[["building graphs"]] <- paste(cmds[["building graphs"]], "--no-insert-ir --no-insert-es --no-insert-ni --no-quantify-graph")

cmds[["testing alternative splicing events"]] <- glue(
  "spladder test \\
  --readlen {rl} \\
  --conditionA {b1} \\
  --conditionB {b2} \\
  --labelA {conditions[1]} \\
  --labelB {conditions[2]} \\
  --no-cap-exp-outliers \\
  --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons \\
  --outdir {dir_out}/SplAdder/graphs"
)

if (write_log) {
  cmds[["building graphs"]] <- paste(
    cmds[["building graphs"]],
    glue("1>{dir_out}/SplAdder/graphs/_log.build 2>&1")
  )
  cmds[["testing alternative splicing events"]] <- paste(
    cmds[["testing alternative splicing events"]],
    glue("1>{dir_out}/SplAdder/graphs/_log.test 2>&1")
  )
}

run_cmds(cmds, run = T)


##### MAJIQ #####
myCreatedir(glue('{dir_out}/MAJIQ/{c("psi", "deltapsi", "voila")}/'))

show_all <- F
threshold <- 0.05

rl <- max(sampleTable$read_lengths)
strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") ss <- "None"
if (strand == "yes") ss <- "forward"
if (strand == "reverse") ss <- "reverse"
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]] %>% gsub("\\.bam", "", .)
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]] %>% gsub("\\.bam", "", .)
b1 <- paste(basename(bams_1), collapse = ",")
b2 <- paste(basename(bams_2), collapse = ",")
bd <- paste(path.expand(dirname(sampleTable$files_bam)), collapse = ",")

cmds <- list()

cmds[["generating config"]] <- glue(
  'echo -e "[info]\\nreadlen={rl}\\nbamdirs={bd}\\
  \\ngenome={fa}\\nstrandness={ss}\\n\\
  [experiments]\\n{conditions[1]}={b1}\\n{conditions[2]}={b2}" > {dir_out}/MAJIQ/config'
)

if (!dir.exists(glue("{dir_out}/MAJIQ/build/"))) dir.create(glue("{dir_out}/MAJIQ/build/"), recursive = T)

cmds[["building majiq"]] <- glue(
  "majiq build {gff} \\
  -c {dir_out}/MAJIQ/config \\
  -j {nproc} \\
  -o {dir_out}/MAJIQ/build/"
)

if (!novel) cmds[["building majiq"]] <- paste(cmds[["building majiq"]], "--disable-denovo --disable-denovo-ir")

m1 <- basename(sampleTable$files_bam[sampleTable$conditions == conditions[1]]) %>%
  gsub("\\.bam", "\\.majiq", .)
m2 <- basename(sampleTable$files_bam[sampleTable$conditions == conditions[2]]) %>%
  gsub("\\.bam", "\\.majiq", .)
files_majiq1 <- glue("{dir_out}/MAJIQ/build/{m1}") %>% paste(collapse = " ")
files_majiq2 <- glue("{dir_out}/MAJIQ/build/{m2}") %>% paste(collapse = " ")
cmds[["calculating psi"]] <- glue(
  "majiq psi \\
  {c(files_majiq1, files_majiq2)} \\
  -j {nproc} \\
  -o {dir_out}/MAJIQ/psi \\
  -n {conditions}"
)

cmds[["calculating deltapsi"]] <- glue(
  "majiq deltapsi \\
  -grp1 {files_majiq1} \\
  -grp2 {files_majiq2} \\
  -j {nproc} \\
  -o {dir_out}/MAJIQ/deltapsi \\
  -n {conditions[1]} {conditions[2]}"
)

cmds[["generating voila tsv"]] <- glue(
  "voila tsv \\
  --threshold {threshold} \\
  {dir_out}/MAJIQ/build/splicegraph.sql \\
  {dir_out}/MAJIQ/deltapsi/{conditions[1]}_{conditions[2]}.deltapsi.voila \\
  -f {dir_out}/MAJIQ/voila/{conditions[1]}_{conditions[2]}_{threshold}.tsv"
)

if (show_all) {
  cmds[["generating voila tsv"]] <- glue(
    "voila tsv \\
    --threshold {threshold} \\
    {dir_out}/MAJIQ/build/splicegraph.sql \\
    {dir_out}/MAJIQ/deltapsi/{conditions[1]}_{conditions[2]}.deltapsi.voila \\
    -f {dir_out}/MAJIQ/voila/{conditions[1]}_{conditions[2]}_{threshold}_showall.tsv \\
    --show-all"
  )
}

if (write_log) {
  cmds[["building majiq"]] <- paste(
    cmds[["building majiq"]],
    glue("1>{dir_out}/MAJIQ/build/_log 2>&1")
  )
  cmds[["calculating psi"]] <- paste(
    cmds[["calculating psi"]],
    glue("1>{dir_out}/MAJIQ/psi/_log.{conditions} 2>&1")
  )
  cmds[["calculating deltapsi"]] <- paste(
    cmds[["calculating deltapsi"]],
    glue("1>{dir_out}/MAJIQ/deltapsi/_log 2>&1")
  )
  cmds[["generating voila tsv"]] <- paste(
    cmds[["generating voila tsv"]],
    glue("1>{dir_out}/MAJIQ/voila/_log 2>&1")
  )
}


##### SGSeq #####
myCreatedir(glue("{dir_out}/SGSeq/"))

library(SGSeq)
library(DEXSeq)

si <- select(sampleTable, samples, files_bam) %>%
  dplyr::rename(sample_name = samples, file_bam = files_bam)
si <- getBamInfo(sample_info = si, core = nproc)
txdb <- importTranscripts(gtf)
txf <- convertToTxFeatures(txdb)
sgf <- convertToSGFeatures(txf)
sgfc <- analyzeFeatures(sample_info = si, features = txf, predict = novel, cores = nproc)
sgvc <- analyzeVariants(sgfc, cores = nproc)
sgvq <- variantFreq(sgvc)
sgv <- rowRanges(sgvc)
sgvc_ds <- getSGVariantCounts(sgv, sample_info = si, cores = nproc)
x <- counts(sgvc_ds)
vid <- variantID(sgvc_ds)
eid <- eventID(sgvc_ds)

sampleData <- data.frame(
  row.names = sampleTable$samples,
  condition = sampleTable$conditions
)
design <- formula(~sample + exon + condition:exon)
dxd <- DEXSeqDataSet(
  countData = x, sampleData = sampleData, design = design,
  featureID = as.character(vid), groupID = as.character(eid)
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxr1 <- DEXSeqResults(dxd)
q_value <- perGeneQValue(dxr1)

save.image(glue("{dir_out}/SGSeq/SGSeq.RData"), compress = T)

##### MISO #####
ds <- T
events <- c("A3SS", "A5SS", "MXE", "RI", "SE")
dir_indexed_events <- "~/doc/miso/indexed_events/v32/pa/"
gff_miso <- "~/projects/as/benchmark/miso/hg38/exons/gencode.v32.annotation.min_1000.const_exons.gff"

strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") ss <- "fr-unstranded"
if (strand == "yes") ss <- "fr-secondstrand"
if (strand == "reverse") ss <- "fr-firststrand"
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]]
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]]
b1 <- paste(bams_1, collapse = " ")
b2 <- paste(bams_2, collapse = " ")
samples_pe <- sampleTable$samples[sampleTable$library_types == "paired-end"]
samples_se <- sampleTable$samples[sampleTable$library_types == "single-end"]
bams_pe <- sampleTable$files_bam[sampleTable$library_types == "paired-end"]
lt <- table(sampleTable$library_types) %>% .[. == max(.)] %>% names()
rl <- table(sampleTable$read_lengths) %>% .[. == max(.)] %>% names()

indexed_events <- vector("character", length = 5)
names(indexed_events) <- events
for (e in events) {
  indexed_events[e] <- list.files(dir_indexed_events, pattern = glue("(?i)({e})"), full.names = T)
}

if (length(bams_pe)) {
  cmds <- list()
  cmds[["computing insert length"]] <- glue(
    "pe_utils --compute-insert-len \\
    {bams_pe} \\
    {gff_miso} \\
    --output-dir {dir_out}/MISO/insert_len"
  )

  run_cmds(cmds)
  cmds <- list()

  dat_il <- data.frame(
    row.names = samples_pe,
    insert_mean = rep(0, length(samples_pe)),
    insert_sdev = rep(0, length(samples_pe))
  )
  for (s in samples_pe) {
    bam <- basename(sampleTable$files_bam[sampleTable$samples == s])
    tmp_chr <- read.table(
      glue("{dir_out}/MISO/insert_len/{bam}.insert_len"), nrows = 1, sep = ",", comment.char = "!"
    ) %>% as.character()
    dat_il[s, "insert_mean"] <- strsplit(tmp_chr[1], "=")[[1]][2]
    dat_il[s, "insert_sdev"] <- strsplit(tmp_chr[2], "=")[[1]][2]
  }
}

cmds[["generating setting file"]] <- glue(
  "echo -e \"[data]\\nfilter_results = True\\nmin_event_reads = 20\\nstrand = {ss}\\n\\n\\
  [cluster]\\ncluster_command = qsub\\n\\n[sampler]\\nburn_in = 500\\nlag = 10\\n\\
  num_iters = 5000\\nnum_processors = 4\\n\" > {dir_out}/MISO/miso_settings.txt"
)

if (length(samples_pe)) {
  run_miso_pe <- list()
  for (e in names(indexed_events)) {
    run_miso_pe[[e]] <- glue(
      "miso --run \\
      -p {nproc} \\
      --event-type {e} \\
      {indexed_events[e]} \\
      {sampleTable$files_bam[match(samples_pe, sampleTable$samples)]} \\
      --output-dir {dir_out}/MISO/{e}/{samples_pe} \\
      --read-len {sampleTable$read_lengths[match(samples_pe, sampleTable$samples)]} \\
      --paired-end {dat_il[samples_pe, 1]} {dat_il[samples_pe, 2]} \\
      --settings-filename {dir_out}/MISO/miso_settings.txt"
    )
  }
  run_miso_pe <- flatten_chr(run_miso_pe)
} else {run_miso_pe <- NULL}
if (length(samples_se)) {
  run_miso_se <- list()
  for (e in names(indexed_events)) {
    run_miso_se[[e]] <- glue(
      "miso --run \\
      -p {nproc} \\
      --event-type {e} \\
      {indexed_events[e]} \\
      {sampleTable$files_bam[match(samples_pe, sampleTable$samples)]} \\
      --output-dir {dir_out}/MISO/{e}/{samples_pe} \\
      --read-len {sampleTable$read_lengths[match(samples_pe, sampleTable$samples)]} \\
      --settings-filename {dir_out}/MISO/miso_settings.txt"
    )
  }
  run_miso_se <- flatten_chr(run_miso_se)
} else {run_miso_se <- NULL}
cmds[["running miso"]] <- c(run_miso_pe, run_miso_se)

summarise_miso <- list()
for (e in names(indexed_events)) {
  summarise_miso[[e]] <- glue(
    "summarize_miso \\
    --summarize-samples \\
    {dir_out}/MISO/{e}/{sampleTable$samples} \\
    {dir_out}/MISO/summary/{e}/{sampleTable$samples}"
  )
}
cmds[["summarising miso"]] <- flatten_chr(summarise_miso)

run_cmds(cmds)
cmds <- list()

if (ds) {
  if (nrow(sampleTable) > 2) {
    cmds[["merging bam"]] <- glue(
      "samtools merge -@ {nproc} {dir_out}/MISO/merge/{conditions}.merge.bam {c(b1, b2)}"
    )
    cmds[["indexing merged bam"]] <- glue(
      "samtools index {dir_out}/MISO/merge/{conditions}.merge.bam"
    )
    if (lt == "paired-end") {
      cmds[["computing merged insert length"]] <- glue(
        "pe_utils --compute-insert-len \\
        --no-bam-filter \\
        {dir_out}/MISO/merge/{conditions}.merge.bam \\
        {gff_miso} \\
        --output-dir {dir_out}/MISO/insert_len"
      )
      run_cmds(cmds)
      cmds <- list()

      dat_il <- data.frame(
        row.names = glue("{conditions}.merge.bam"),
        insert_mean = rep(0, 2),
        insert_sdev = rep(0, 2)
      )
      for (i in 1:2) {
        tmp_chr <- read.table(
          glue("{dir_out}/MISO/insert_len/{conditions[i].merge.bam}.insert_len"), nrows = 1, sep = ",", comment.char = "!"
        ) %>% as.character()
        dat_il[i, "insert_mean"] <- strsplit(tmp_chr[1], "=")[[1]][2]
        dat_il[i, "insert_sdev"] <- strsplit(tmp_chr[2], "=")[[1]][2]
      }
    }

    if (lt == "paired-end") {
      run_miso <- list()
      for (e in names(indexed_events)) {
        run_miso[[e]] <- glue(
          "miso --run \\
          -p {nproc} \\
          --no-bam-filter \\
          --event-type {e} \\
          {indexed_events[e]} \\
          {dir_out}/MISO/merge/{conditions}.merge.bam \\
          --output-dir {dir_out}/MISO/{e}/{conditions}_merge \\
          --read-len {rl} \\
          --paired-end {dat_il[, 1]} {dat_il[, 2]} \\
          --settings-filename {dir_out}/MISO/miso_settings.txt"
        )
      }
    } else {
      run_miso <- list()
      for (e in names(indexed_events)) {
        run_miso[[e]] <- glue(
          "miso --run \\
          -p {nproc} \\
          --no-bam-filter \\
          --event-type {e} \\
          {indexed_events[e]} \\
          {dir_out}/MISO/merge/{conditions}.merge.bam \\
          --output-dir {dir_out}/MISO/{e}/{conditions}_merge \\
          --read-len {rl} \\
          --settings-filename {dir_out}/MISO/miso_settings.txt"
        )
      }
    }
    cmds[["running miso on merged bam"]] <- flatten_chr(run_miso)

    summarise_miso <- list()
    for (e in names(indexed_events)) {
      summarise_miso[[e]] <- glue(
        "summarize_miso \\
        --summarize-samples \\
        {dir_out}/MISO/{e}/{conditions}_merge \\
        {dir_out}/MISO/summary/{e}/{conditions}_merge"
      )
    }
    cmds[["summarising miso on merged bam"]] <- flatten_chr(summarise_miso)

    compare_miso <- list()
    for (e in names(indexed_events)) {
      compare_miso[[e]] <- glue(
        "compare_miso --compare-samples \\
        {dir_out}/MISO/{e}/{conditions[1]}_merge \\
        {dir_out}/MISO/{e}/{conditions[2]}_merge \\
        {dir_out}/MISO/comparison/{conditions[1]}_{conditions[2]}_merged/${e}"
      )
    }
    cmds[["comparing miso on merged bam"]] <- flatten_chr(compare_miso)

  } else {
    compare_miso <- list()
    for (e in names(indexed_events)) {
      compare_miso[[e]] <- glue(
        "compare_miso --compare-samples \\
        {dir_out}/MISO/{e}/{sampleTable$samples[sampleTable$conditions == conditions[1]]} \\
        {dir_out}/MISO/{e}/{sampleTable$samples[sampleTable$conditions == conditions[2]]} \\
        {dir_out}/MISO/comparison/{conditions[1]}_{conditions[2]}/${e}"
      )
    }
    cmds[["comparing miso"]] <- flatten_chr(compare_miso)
  }
  run_cmds(cmds, T)
}


##### SUPPA #####
myCreatedir(glue('{dir_out}/SUPPA/{c("", "tpm", "psi", "ds")}/'))

library(tximport)
tpm_type <- "salmon"
events <- c("A3", "A5", "AF", "AL", "MX", "RI", "SE")
dir_events <- "~/doc/suppa/gencode_v32/pa/"

tpm_isoform_dir <- glue("~/projects/AS/data/sim/rsem/sim_1/{tpm_type}/")
conditions <- unique(sampleTable$conditions)

cmds <- list()
if (length(list.files(dir_events, pattern = glue(".*\\.ioe"), full.names = T)) < length(events)) {
  cmds[["generating events"]] <- glue(
    "suppa.py generateEvents -i {gtf} -o {dir_events} -f ioe -e SE SS MX RI FL"
  )
  run_cmds(cmds)
  cmds <- list()
}

ioe_suppa <- vector("character", length = length(events))
names(ioe_suppa) <- events
for (e in events) ioe_suppa[e] <- list.files(dir_events, pattern = glue(".*{e}.*\\.ioe"), full.names = T)

st <- sampleTable
for (s in st$samples) {
  st$files_tpm[st$samples == s] <- list.files(
    glue("{tpm_isoform_dir}/{s}/"),
    pattern = "quant.sf*",
    full.names = T
  )
}
files_tpm <- st$files_tpm %>% setNames(st$samples)
txi <- tximport(files_tpm, type = tpm_type, txOut = T)
tpm_1 <- txi$abundance[, st$samples[st$conditions == conditions[1]]]
tpm_2 <- txi$abundance[, st$samples[st$conditions == conditions[2]]]
write.table(tpm_1, glue("{dir_out}/SUPPA/tpm/{conditions[1]}.tpm"), quote = F, sep = "\t")
write.table(tpm_2, glue("{dir_out}/SUPPA/tpm/{conditions[2]}.tpm"), quote = F, sep = "\t")

calculate_psi <- vector("list", length = length(events))
for (e in events) {
  calculate_psi[[e]] <- glue(
    "suppa.py psiPerEvent \\
    -i {ioe_suppa[e]} \\
    -e {dir_out}/SUPPA/tpm/{conditions}.tpm \\
    -o {dir_out}/SUPPA/psi/{conditions}.{e}"
  )
  if (write_log) calculate_psi[[e]] <- paste(calculate_psi[[e]], glue("1>{dir_out}/SUPPA/psi/_log.{conditions}.{e} 2>&1"))
}
cmds[["calculating psi"]] <- flatten_chr(calculate_psi)

perform_ds <- vector("list", length = length(events))
for (e in events) {
  perform_ds[[e]] <- glue(
    "suppa.py diffSplice \\
    -m empirical \\
    -p {dir_out}/SUPPA/psi/{conditions[1]}.{e}.psi {dir_out}/SUPPA/psi/{conditions[2]}.{e}.psi \\
    -e {dir_out}/SUPPA/tpm/{conditions[1]}.tpm {dir_out}/SUPPA/tpm/{conditions[2]}.tpm \\
    -i {ioe_suppa[e]} \\
    --gene-correction \\
    -o {dir_out}/SUPPA/ds/{conditions[1]}_{conditions[2]}.{e}"
  )
  if (write_log) perform_ds[[e]] <- paste(perform_ds[[e]], glue("1>{dir_out}/SUPPA/ds/_log.{e} 2>&1"))
}
cmds[["performing differential splicing analysis"]] <- flatten_chr(perform_ds)

run_cmds(cmds, T)


##### Whippet #####
myCreatedir(glue('{dir_out}/Whippet/{c("index", "quant", "delta", "merge")}/'))

conditions <- unique(sampleTable$conditions)

prefix_index <- "~/doc/whippet/hg38/pa/index"
bam_min_reads <- 2

if (length(list.files(dirname(prefix_index), pattern = basename(prefix_index))) != 2) {
  myCreatedir(dirname(prefix_index))
  cmds <- list()
  if (novel) {
    cmds[["merging bam"]] <- glue(
      "samtools merge -@ {nproc} -l 9 {dir_out}/Whippet/merge/merged.bam {paste(sampleTable$files_bam, collapse = ' ')}"
    )
    cmds[["removing duplicates"]] <- glue(
      "samtools rmdup -S {dir_out}/Whippet/merge/merged.bam {dir_out}/Whippet/merge/merged.rmdup.bam"
    )
    cmds[["indexing rmdup bam"]] <- glue(
      "samtools index -@ {nproc} {dir_out}/Whippet/merge/merged.rmdup.bam"
    )
    cmds[["generating index"]] <- glue(
      "whippet-index.jl \\
      --fasta {fa}.gz \\
      --gtf {gtf}.gz \\
      --suppress-low-tsl \\
      --bam-both-novel \\
      --bam-min-reads {bam_min_reads} \\
      --bam {dir_out}/Whippet/merge/merged.rmdup.bam \\
      -x {prefix_index}"
    )
  } else {
    cmds[["generating index"]] <- glue(
      "whippet-index.jl \\
      --fasta {fa}.gz \\
      --gtf {gtf}.gz \\
      --suppress-low-tsl \\
      -x {prefix_index}"
    )
  }
  run_cmds(cmds)
}

cmds <- list()
quantify_fq <- vector("list", length = nrow(sampleTable))
names(quantify_fq) <- sampleTable$samples
for (s in names(quantify_fq)) {
  tmp_fq <- strsplit(sampleTable$files_fq[sampleTable$samples == s], ";") %>% unlist() %>% paste(collapse = " ")
  quantify_fq[[s]] <- glue(
    "whippet-quant.jl \\
    {tmp_fq} \\
    -x {prefix_index} \\
    -o {dir_out}/Whippet/quant/{s}"
  )
  if (write_log) quantify_fq[[s]] <- paste(quantify_fq[[s]], glue("1>{dir_out}/Whippet/quant/_log.{s} 2>&1"))
}
cmds[["quantifying fastq"]] <- flatten_chr(quantify_fq)

cmds[["quantifying fastq"]] <- paste(paste(cmds[["quantifying fastq"]], collapse = " & "), "& wait")

files_quant1 <- vector("character", length = length(which(sampleTable$conditions == conditions[1]))) %>%
  setNames(sampleTable$samples[sampleTable$conditions == conditions[1]])
for (s in names(files_quant1)) {
  files_quant1[s] <- list.files(glue("{dir_out}/Whippet/quant/"), pattern = glue("{s}.psi.gz"), full.names = T)
}
files_quant2 <- vector("character", length = length(which(sampleTable$conditions == conditions[2]))) %>%
  setNames(sampleTable$samples[sampleTable$conditions == conditions[2]])
for (s in names(files_quant2)) {
  files_quant2[s] <- list.files(glue("{dir_out}/Whippet/quant/"), pattern = glue("{s}.psi.gz"), full.names = T)
}

cmds[["comparing psi"]] <- glue(
  "whippet-delta.jl \\
  -a {paste(files_quant1, collapse = ',')} \\
  -b {paste(files_quant2, collapse = ',')} \\
  -o {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}"
)
if (write_log) cmds[["comparing psi"]] <- paste(cmds[["comparing psi"]], glue("1>{dir_out}/Whippet/delta/_log 2>&1"))

cmds[["tidy output"]] <- glue(
  "zcat {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}.diff.gz | \\
  cut -f 1,2,3,4,5,6,7,8,9,10,11 > {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}.diff2 | \\
  gzip {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}.diff2"
)

##### ASpli #####
myCreatedir(glue("{dir_out}/ASpli/"))

# method_ds <- NULL

rl <- table(sampleTable$read_lengths) %>% .[. == max(.)] %>% names() %>% as.numeric()
strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") strandMode <- 0
if (strand == "yes") strandMode <- 1
if (strand == "reverse") strandMode <- 2
lt <- table(sampleTable$library_types) %>% .[. == max(.)] %>% names()
if (lt == "paired-end") libType <- "PE"
if (lt == "single-end") libType <- "SE"

library(ASpli)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF(gtf, format = "gtf")
symbols <- data.frame(
  row.names = genes(txdb),
  symbol = paste("This is symbol of gene:", genes(txdb))
)
features <- binGenome(txdb, geneSymbols = symbols, logTo = glue("{dir_out}/ASpli/ASpli_binFeatures.log"))
if (F) {
  geneCoord <- featuresg(features)
  binCoord <- featuresb(features)
  junctionCoord <- featuresj(features)
  binMetadata <- mcols(binCoord)
}
targets <- dplyr::select(sampleTable, files_bam, conditions) %>%
  dplyr::rename(bam = files_bam)
rownames(targets) <- sampleTable$samples
# bam <- loadBAM(targets, cores = nproc)
# counts <- readCounts(
#   features,
#   bam,
#   targets,
#   cores = nproc,
#   readLength = rl,
#   maxISize = 50000
# )
counts <- gbCounts(
  features = features, targets = targets, libType = libType, strandMode = strandMode,
  minReadLength = 101, maxISize = 50000
)
if (F) {
  GeneCounts <- countsg(counts)
  GeneRd <- rdsg(counts)
  BinCounts <- countsb(counts)
  BinRd <- rdsb(counts)
  JunctionCounts <- countsj(counts)
  # Accessing counts for intron flanking regions
  e1iCounts <- countse1i(counts)
  ie2Counts <- countsie2(counts)
}
writeCounts(counts = counts, output.dir = glue("{dir_out}/ASpli/"))
writeRds(counts = counts, output.dir = glue("{dir_out}/ASpli/"))
# as <- AsDiscover(
#   counts, targets, features, bam, readLength = rl, cores = nproc
# )
as <- jCounts(
  counts = counts, features = features, libType = libType, strandMode = strandMode,
  minReadLength = 101
) # Accesors: irPIR, altPSI, esPSI, junctionsPIR, junctionsPJU
if (F) {
  # Access result from an ASpliAS object
  irPIR <- irPIR(as)
  altPSI <- altPSI(as)
  esPSI <- esPSI(as)
  junctionsPIR <- junctionsPIR(as)
  junctionsPJU <- junctionsPJU(as)
  allBins <- joint(as)
}
writeAS(as = as, output.dir = glue("{dir_out}/ASpli/"))
counts@condition.order <- c(basics$condition_1, basics$condition_2)
contrast <- c(1, -1)
# if (method_ds == "edgeR") {
du <- gbDUreport( # DUreportBinSplice
  counts,
  contrast = contrast,
  minGenReads = 0,
  minBinReads = 0
)
writeDU(du, output.dir = glue("{dir_out}/ASpli/"))
# } else {
# }
# du <- junctionDUreport(
#   counts,
#   targets,
#   appendTo = du,
#   contrast = contrast
# )
jdu <- jDUreport(asd = as, contrast = contrast, maxFDRForParticipation = 1)
writeJDU(jdu, output.dir = glue("{dir_out}/ASpli/"))
# results <- mergeBinDUAS(du, as, targets)
sr <- splicingReport(du, jdu, counts)
writeSplicingReport(sr, glue("{dir_out}/ASpli/"))
is <- integrateSignals(sr, as, bin.FC = 0, bin.fdr = 1, bjs.fdr = 1, a.fdr = 1, l.fdr = 1)
write.table(signals(is), glue("{dir_out}/ASpli/is.txt"), quote = F, sep = "\t", row.names = F)

write.table(binsDU(du), glue("{dir_out}/ASpli/du.txt"), quote = F, sep = "\t", row.names = F)
save.image(glue("{dir_out}/ASpli/ASpli.RData"), compress = T)


##### BANDITS #####
library(BANDITS)
library(tximport)
myCreatedir(glue("{dir_out}/BANDITS/"))

infer_prior <- F
salmon_dir <- "~/projects/AS/data/sim/rsem/sim_1/salmon/"

file_tx2gene <- "~/doc/tx2gene/tx2gene_v32.txt"
tx2gene <- read.delim(file_tx2gene, header = F) %>%
  dplyr::select(V2, V4) %>%
  setNames(c("gene_id", "transcript_id"))

samples_design <- sampleTable %>%
  dplyr::rename(sample_id = samples, group = conditions) %>%
  dplyr::select(sample_id, group)

quant_files <- vector("character", length = nrow(samples_design)) %>%
  setNames(samples_design$sample_id)
for (s in names(quant_files)) {
  quant_files[s] <- list.files(
    glue("{salmon_dir}/{s}/"),
    pattern = "quant.sf*",
    full.names = T
  )
}

txi <- tximport(quant_files, type = "salmon", txOut = T)
# counts <- txi$counts
eff_len = eff_len_compute(x_eff_len = txi$length)

equiv_classes_files <- vector("character", length = nrow(samples_design)) %>%
  setNames(samples_design$sample_id)
for (s in names(equiv_classes_files)) {
  equiv_classes_files[s] <- list.files(
    glue("{salmon_dir}/{s}/aux_info/"),
    pattern = "eq_classes.txt*",
    full.names = T
  )
}

input_data <- create_data(
  salmon_or_kallisto = "salmon",
  gene_to_transcript = tx2gene,
  salmon_path_to_eq_classes = equiv_classes_files,
  eff_len = eff_len,
  max_genes_per_group = 100,
  n_cores = nrow(samples_design)
)

if (infer_prior) {
  precision <- prior_precision(
    gene_to_transcript = tx2gene,
    transcript_counts = txi$counts,
    n_cores = nrow(samples_design),
    transcripts_to_keep = NULL
  )
} else {
  precision <- NULL
}


set.seed(1)
results <- test_DTU(
  BANDITS_data = input_data,
  precision = precision$prior,
  samples_design = samples_design,
  group_col_name = "group",
  R = 10^4,
  burn_in = 2*10^3,
  n_cores = nproc,
  gene_to_transcript = tx2gene
)

results_gene <- top_genes(results)
results_transcript <- top_transcripts(results)

write.table(
  results_gene, glue("{dir_out}/BANDITS/results_gene.txt"),
  row.names = F, sep = "\t", quote = F
)
write.table(
  results_transcript, glue("{dir_out}/BANDITS/results_transcript.txt"),
  row.names = F, sep = "\t", quote = F
)
save.image(glue("{dir_out}/BANDITS/bandits.RData"), compress = T)


##### NBSplice #####
myCreatedir(glue("{dir_out}/NBSplice/"))

library(NBSplice)
library(tximport)
library(BiocParallel)

tpm_isoform_dir <- "~/projects/AS/data/sim/rsem/sim_1/salmon/"
colName <- "condition"
file_tx2gene <- "~/doc/tx2gene/tx2gene_v32.txt"
test <- "F"

BPPARAM <- MulticoreParam(nproc)
designMatrix <- sampleTable %>%
  dplyr::select(samples, conditions) %>%
  setNames(c("sample", "condition"))
designMatrix$condition <- factor(designMatrix$condition, levels = c("s1", "s2"))
rownames(designMatrix) <- designMatrix$sample

geneIso <- read.delim(file_tx2gene, header = F) %>%
  dplyr::select(V2, V4) %>%
  setNames(c("gene_id", "isoform_id"))
rownames(geneIso) <- geneIso$isoform_id

files_quant <- vector("character", length = nrow(designMatrix)) %>%
  setNames(designMatrix$sample)
for (s in names(files_quant)) {
  files_quant[s] <- list.files(
    glue("{tpm_isoform_dir}/{s}/"),
    pattern = "quant.sf*",
    full.names = T
  )
}
txi <- tximport(files_quant, type = "salmon", txOut = T)
isoCounts <- as.data.frame(txi$counts)
geneIso <- geneIso[geneIso$isoform_id %in% rownames(isoCounts),]

myIsoDataSet <- IsoDataSet(isoCounts, designMatrix, colName, geneIso, BPPARAM = BPPARAM)
myDSResults <- NBTest(myIsoDataSet, colName, test = test, BPPARAM = BPPARAM)
results <- results(myDSResults, filter = F)
write.table(
  results, glue("{dir_out}/NBSplice/results.txt"),
  sep = "\t", row.names = F, quote = F
)
save.image(glue("{dir_out}/NBSplice/nbsplice.RData"), compress = T)


##### IsoformSwitchAnalyzeR #####
library(IsoformSwitchAnalyzeR)
myCreatedir(glue("{dir_out}/IsoformSwitchAnalyzeR/"))

tpm_isoform_dir <- "~/projects/AS/data/sim/rsem/sim_1/salmon/"

sampleVector <- vector("character", length = nrow(sampleTable)) %>%
  setNames(sampleTable$samples)
for (s in names(sampleVector)) {
  sampleVector[s] <- list.files(
    glue("{tpm_isoform_dir}/{s}/"),
    pattern = "quant.sf*",
    full.names = T
  )
}
salmonQuant <- importIsoformExpression(
  sampleVector = sampleVector,
  addIsofomIdAsColumn = TRUE
)

myDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = sampleTable$conditions
)

aSwitchList <- importRdata(
  isoformCountMatrix = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix = myDesign,
  isoformExonAnnoation = gtf,
  isoformNtFasta = fa_transcript,
  removeTECgenes = F,
  estimateDifferentialGeneRange = F,
  comparisonsToMake = data.frame(
    condition_1 = unique(sampleTable$conditions)[1],
    condition_2 = unique(sampleTable$conditions)[2]
  )
)

aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList,
  dIFcutoff = 0,
  reduceToSwitchingGenes = F,
  reduceFurtherToGenesWithConsequencePotential = F
)

aSwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = aSwitchListAnalyzed,
  onlySwitchingGenes = F,
  alpha = 1,
  dIFcutoff = 0
)
write.table(
  aSwitchListAnalyzed$isoformFeatures,
  glue("{dir_out}/IsoformSwitchAnalyzeR/isoformFeatures.txt"),
  quote = F, sep = "\t", row.names = F
)
write.table(
  aSwitchListAnalyzed$AlternativeSplicingAnalysis,
  glue("{dir_out}/IsoformSwitchAnalyzeR/AlternativeSplicingAnalysis.txt"),
  quote = F, sep = "\t", row.names = F
)
save.image(glue("{dir_out}/IsoformSwitchAnalyzeR/isoformswitchanalyzer.RData"), compress = T)


##### DiffSplice #####
myCreatedir(glue("{dir_out}/DiffSplice/"))

dir_diffsplice <- "~/bin/diffsplice_0.1.2beta"
bam2sam <- T

settings.cfg <- t(data.frame(
  thresh_junction_filter_max_read_support = 5,
  thresh_junction_filter_mean_read_support = 5,
  thresh_junction_filter_num_samples_presence = 2,
  ignore_minor_alternative_splicing_variants = "no",
  thresh_average_read_coverage_exon = 0,
  thresh_average_read_coverage_intron = 0,
  balanced_design_for_permutation_test = "no",
  false_discovery_rate = 0.05,
  thresh_foldchange_up = 2,
  thresh_foldchange_down = 0.5,
  thresh_sqrtJSD = 0.0
))
write.table(
  settings.cfg, glue("{dir_out}/DiffSplice/settings.cfg"),
  quote = F, sep = "\t", col.names = F
)

datafile.cfg <- data.frame(
  factor(sampleTable$conditions, levels = unique(sampleTable$conditions), labels = c("g1", "g2")),
  rep("id1", nrow(sampleTable)),
  c(
    paste0("s", 1:table(sampleTable$conditions)[unique(sampleTable$conditions)[1]]),
    paste0("s", 1:table(sampleTable$conditions)[unique(sampleTable$conditions)[2]])
  ),
  path.expand(gsub(".bam",".sam",sampleTable$files_bam))
)
write.table(
  datafile.cfg, glue("{dir_out}/DiffSplice/datafile.cfg"),
  quote = F, sep = "\t", col.names = F, row.names = F
)

cmds <- list()
if (bam2sam) {
  cmds[["converting bams to sams"]] <- glue(
    'samtools view -h -@ {nproc} {sampleTable$files_bam} -O SAM \\
    -o {gsub(".bam",".sam",sampleTable$files_bam)}'
  )
  if (parallel) {
    cmds[["converting bams to sams"]] <- myParallel2(cmds[["converting bams to sams"]], np)
  }
}
run_cmds(cmds, T)

cmds <- list()
cmds[["running diffsplice"]] <- glue(
  "cd {dir_diffsplice} && \\
  diffsplice {dir_out}/DiffSplice/settings.cfg \\
  {dir_out}/DiffSplice/datafile.cfg \\
  {dir_out}/DiffSplice/"
)
if (write_log) cmds[["running diffsplice"]] <- paste(cmds[["running diffsplice"]], glue("1>{dir_out}/DiffSplice/_log 2>&1"))

run_cmds(cmds, T)


##### psichomics #####
myCreatedir(glue("{dir_out}/psichomics/"))

library(psichomics)
strandedness <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strandedness == "no") sn <- "unstranded"
if (strandedness == "yes") sn <- "stranded"
if (strandedness == "reverse") sn <- "stranded (reverse)"

files_genequant <- glue(
  "~/projects/AS/data/sim/rsem/sim_1/star/sample{1:8}/ReadsPerGene.out.tab"
) %>% setNames(sampleTable$samples)
files_junctionquant <- glue(
  "~/projects/AS/data/sim/rsem/sim_1/star/sample{1:8}/SJ.out.tab"
) %>% setNames(sampleTable$samples)
anno_suppa <- "~/doc/suppa/gencode_v32/pa/"
anno_rmats <- glue("{dir_out}/rMATS/")
anno_miso <- "~/doc/miso/gff/gff_v32/pa/commonshortest/"
anno_vast <- "~/doc/vast_tools/Hs2/TEMPLATES/"
prefix_suppa <- "events"
prefix_miso <- "hg38"
prefix_vast <- "Hs2"
fil <- F

suppa <- parseSuppaAnnotation(anno_suppa, genome = prefix_suppa)
mats <- parseMatsAnnotation(anno_rmats, novelEvents = novel)
miso <- parseMisoAnnotation(anno_miso, genome = prefix_miso, types = c("SE", "MXE", "A5SS", "A3SS", "RI"))
vast <- parseVastToolsAnnotation(anno_vast, genome = prefix_vast)
annot <- prepareAnnotationFromEvents(suppa, vast, mats, miso)
saveRDS(annot, file = glue("{dir_out}/psichomics/annot.rds"))

prepareGeneQuant(
  files_genequant[sampleTable$samples], strandedness = sn,
  output = glue("{dir_out}/psichomics/psichomics_gene_counts.txt")
)
prepareJunctionQuant(
  files_junctionquant, output = glue("{dir_out}/psichomics/psichomics_junctions.txt")
)
write.table(
  dplyr::select(sampleTable, samples, conditions) %>% dplyr::rename(`Sample ID` = samples),
  glue("{dir_out}/psichomics/psichomics_metadata.txt"),
  sep = "\t", quote = F, row.names = F
)

data <- loadLocalFiles(glue("{dir_out}/psichomics/"))
sampleInfo <- data[[1]]$`Sample metadata`
geneExpr <- data[[1]]$`Gene expression`
junctionQuant <- data[[1]]$`Junction quantification`
if (fil) {
  filter <- filterGeneExpr(geneExpr)
} else {
  filter <- NULL
}
geneExprNorm <- normaliseGeneExpression(geneExpr, geneFilter = filter)
psi <- quantifySplicing(annot, junctionQuant)

condition_1 <- unique(sampleTable$conditions)[1]
condition_2 <- unique(sampleTable$conditions)[2]
s1VSs2 <- list(
  sampleTable$samples[sampleTable$conditions == condition_1],
  sampleTable$samples[sampleTable$conditions == condition_2]
) %>% setNames(unique(sampleTable$conditions))
diffSplicing <- diffAnalyses(psi, s1VSs2)
write.table(
  diffSplicing, glue("{dir_out}/psichomics/results.txt"),
  quote = F, sep = "\t"
)
write.table(
  psi, glue("{dir_out}/psichomics/psi.txt"),
  quote = F, sep = "\t"
)
save.image(glue("{dir_out}/psichomics/psichomics.RData"), compress = T)


##### DIEGO #####
myCreatedir(glue("{dir_out}/DIEGO/"))

files_star_junc <- glue(
  "~/projects/AS/data/sim/rsem/sim_1/star/sample{1:8}/SJ.out.tab"
) %>% setNames(sampleTable$samples)
dir_diego <- "~/bin/DIEGO/"

star_junc_table <- data.frame(
  samples = sampleTable$samples,
  files_star_junc = path.expand(files_star_junc[sampleTable$samples])
)
write.table(
  star_junc_table, glue("{dir_out}/DIEGO/star_junc_info.txt"),
  quote = F, sep = "\t", col.names = F, row.names = F
)

write.table(
  dplyr::select(sampleTable, conditions, samples), glue("{dir_out}/DIEGO/condition.txt"),
  quote = F, col.names = F, row.names = F, sep = "\t"
)

cmds <- list()
cmds[["generating bed file"]] <- glue(
  "perl {dir_diego}/gfftoDIEGObed.pl -g {gtf} -o {dir_out}/DIEGO/diego.bed"
)
cmds[["generating junction table"]] <- glue(
  "python {dir_diego}/pre_STAR.py \\
  -l {dir_out}/DIEGO/star_junc_info.txt \\
  -d {dir_out}/DIEGO/diego.bed \\
  -o {dir_out}/DIEGO/"
)
cmds[["running DIEGO"]] <- glue(
  "python {dir_diego}/diego.py \\
  -a {dir_out}/DIEGO/junction_table.txt \\
  -b {dir_out}/DIEGO/condition.txt \\
  -x {sampleTable$condition[1]} \\
  -q 0.05 \\
  -z 0.0 1>{dir_out}/DIEGO/results.txt"
)

run_cmds(cmds, T)


##### VAST-TOOLS #####
myCreatedir(glue('{dir_out}/VAST-TOOLS/{c("align", "combine", "diff")}'))

sp <- "hg38"
min_Fr <- 0.5
method_das <- "diff"

cmds <- list()
files_read1 <- sapply(strsplit(sampleTable$files_fq, ";"), function(x) {return(x[1])})
files_read2 <- sapply(strsplit(sampleTable$files_fq, ";"), function(x) {return(x[2])})
myCreatedir(glue("{dir_out}/VAST-TOOLS/align/{sampleTable$samples}"))
cmds[["aligning reads"]] <- glue(
  "vast-tools align \\
  --sp {sp} \\
  --cores {nproc} \\
  --name {sampleTable$samples} \\
  --output {dir_out}/VAST-TOOLS/align/{sampleTable$samples}/ \\
  {files_read1} {files_read2}"
)
if (write_log) cmds[["aligning reads"]] <- paste(cmds[["aligning reads"]], glue("1>{dir_out}/VAST-TOOLS/align/{sampleTable$samples}/_log 2>&1"))
if (parallel) cmds[["aligning reads"]] <- myParallel2(cmds[["aligning reads"]], np)

myCreatedir(glue("{dir_out}/VAST-TOOLS/combine/to_combine/"))
cmds[["joining align results to a same directory"]] <- glue(
  "cp {dir_out}/VAST-TOOLS/align/{sampleTable$samples}/to_combine/* \\
  {dir_out}/VAST-TOOLS/combine/to_combine/"
)

cmds[["combing align outputs"]] <- glue(
  "vast-tools combine \\
  -sp {sp} \\
  --cores 5 \\
  --output {dir_out}/VAST-TOOLS/combine/"
)
if (write_log) cmds[["combing align outputs"]] <- paste(cmds[["combing align outputs"]], glue("1>{dir_out}/VAST-TOOLS/combine/_log.combine 2>&1"))

table_groups <- dplyr::select(sampleTable, samples, conditions)
write.table(
  table_groups, glue("{dir_out}/VAST-TOOLS/groups"),
  sep = "\t", quote = F, col.names = F, row.names = F
)

file_inclusion <- list.files(glue("{dir_out}/VAST-TOOLS/combine/"), pattern = "*.tab", full.names = T)
cmds[["tidying combine outputs"]] <- glue(
  "vast-tools tidy \\
  {file_inclusion} \\
  -min_Fr {min_Fr} \\
  -groups {dir_out}/VAST-TOOLS/groups \\
  -outFile {dir_out}/VAST-TOOLS/combine/INCLUSION_LEVELS_FULL_tidy.tab"
)
condition_1 <- unique(sampleTable$conditions)[1]
condition_2 <- unique(sampleTable$conditions)[2]
samples_1 <- sampleTable$samples[sampleTable$conditions == condition_1]
samples_2 <- sampleTable$samples[sampleTable$conditions == condition_2]
if ((max(table(sampleTable$conditions)) > 3 && max(table(sampleTable$conditions)) <= 5) || "diff" %in% method_das) {
  cmds[["running vast-tools diff"]] <- glue(
    "vast-tools diff \\
    -i {file_inclusion} \\
    -a {paste(samples_1, collapse = ',')} \\
    -b {paste(samples_2, collapse = ',')} \\
    --sampleNameA {condition_1} \\
    --sampleNameB {condition_2} \\
    -c {nproc} \\
    -o {dir_out}/VAST-TOOLS/combine/ \\
    -d diff"
  )
  if (write_log) cmds[["running vast-tools diff"]] <- paste(cmds[["running vast-tools diff"]], glue("1>{dir_out}/VAST-TOOLS/combine/_log.diff 2>&1"))
}
if (max(table(sampleTable$conditions)) <= 3 || "compare" %in% method_das) {
  cmds[["running vast-tools compare"]] <- glue(
    "vast-tools compare \\
    {file_inclusion} \\
    -a {paste(samples_1, collapse = ',')} \\
    -b {paste(samples_2, collapse = ',')} \\
    -name_A {condition_1} \\
    -name_B {condition_2} \\
    --print_dPSI \\
    --outFile compare.tab"
  )
  if (write_log) cmds[["running vast-tools compare"]] <- paste(cmds[["running vast-tools compare"]], glue("1>{dir_out}/VAST-TOOLS/combine/_log.compare 2>&1"))
}
run_cmds(cmds)


##### RATS #####
## 0.05 ##
myCreatedir(glue("{dir_out}/RATS/"))

library(rats)

condition_1 <- unique(sampleTable$conditions)[1]
condition_2 <- unique(sampleTable$conditions)[2]
samples_1 <- sampleTable$samples[sampleTable$conditions == condition_1]
samples_2 <- sampleTable$samples[sampleTable$conditions == condition_2]

samples_A <- glue("~/projects/AS/data/sim/rsem/sim_1/salmon/{samples_1}/")
samples_B <- glue("~/projects/AS/data/sim/rsem/sim_1/salmon/{samples_2}/")
scales <- rep(40, 8)

myannot <- annot2ids(gtf)
mydata <- fish4rodents(
  A_paths = samples_A, B_paths = samples_B,
  annot = myannot, scaleto = 1e+06, threads = nproc
)
mydtu <- call_DTU(
  annot = myannot,
  boot_data_A = mydata$boot_data_A,
  boot_data_B = mydata$boot_data_B,
  name_A = condition_1,
  name_B = condition_2,
  verbose = T,
  scaling = scales,
  abund_thresh = 0,
  dprop_thresh = 0.05,
  threads = nproc,
  seed = 1,
  qbootnum = 200
)
write.table(
  mydtu$Genes, glue("{dir_out}/RATS/genes.txt"),
  quote = F, row.names = F, sep = "\t"
)
write.table(
  mydtu$Transcripts, glue("{dir_out}/RATS/transcripts.txt"),
  quote = F, row.names = F, sep = "\t"
)
save.image(glue("{dir_out}/RATS/rats.RData"), compress = T)


##### dSpliceType #####
## 0.05 ##
myCreatedir(glue("{dir_out}/dSpliceType/"))

dir_dsplicetype <- "~/bin/dSpliceType-2.0.0/"
ram <- "100g"
type_gtf <- "gencode"
delta <- 0.05

if (type_gtf %in% c("gencode", "ucsc")) file_gtf2gff <- "UCSC_gtf_to_gff.pl"
if (type_gtf %in% c("ensembl_hg38")) file_gtf2gff <- "ensembl_gtf_to_gff_v2.pl"
if (type_gtf %in% c("ensembl_hg19")) file_gtf2gff <- "ensembl_gtf_to_gff.pl"

cmds <- list()
cmds[["converting gtf to gff"]] <- glue(
  'perl {dir_dsplicetype}/gtfTogff/{file_gtf2gff} {gtf} \\
  > {dir_out}/dSpliceType/{gsub(".gtf", "", basename(gtf))}.gff'
)
cmds[["converting bams to bedgraphs"]] <- glue(
  "genomeCoverageBed -bg -ibam {sampleTable$files_bam} \\
  > {dir_out}/dSpliceType/{sampleTable$samples}.bedgraph"
)
if (parallel) cmds[["converting bams to bedgraphs"]] <- myParallel2(cmds[["converting bams to bedgraphs"]], np)

b1 <- paste(path.expand(glue("{dir_out}/dSpliceType/{basics$samples_1}.bedgraph")), collapse = ",")
b2 <- paste(path.expand(glue("{dir_out}/dSpliceType/{basics$samples_2}.bedgraph")), collapse = ",")
j1 <- paste(path.expand(glue("{basics$juncs_star_1}.sorted.bed")), collapse = ",")
j2 <- paste(path.expand(glue("{basics$juncs_star_2}.sorted.bed")), collapse = ",")

cmds[["running dsplicetype"]] <- glue(
  'java -jar -Xmx{ram} -Xms{ram} {dir_dsplicetype}/dSpliceType.jar \\
  -g {dir_out}/dSpliceType/{gsub(".gtf", "", basename(gtf))}.gff \\
  -b1 {b1} \\
  -b2 {b2} \\
  -j1 {j1} \\
  -j2 {j2} \\
  -C 10 \\
  -a 0.05 \\
  -L {1 - delta} \\
  -U {1 + delta} \\
  -o {dir_out}/dSpliceType/'
)
if (write_log) cmds[["running dsplicetype"]] <- paste(cmds[["running dsplicetype"]], glue("1>{dir_out}/dSpliceType/_log 2>&1"))

cmds[["tidying dsplicetype outputs"]] <- glue(
  "awk -F'\\t' 'BEGIN { OFS = FS }; NF { NF -= 1 }; 1' \\
  [dir_out]/dSpliceType/dSpliceType_[c('SE', 'MXE', 'A3SS', 'A5SS', 'RI')].txt \\
  > [dir_out]/dSpliceType/dSpliceType_[c('SE', 'MXE', 'A3SS', 'A5SS', 'RI')]_tidy.txt",
  .open = "[", .close = "]",
)


##### PSI-Sigma #####
myCreatedir(glue("{dir_out}/PSI-Sigma/"))

nread <- 5
dir_psi_sigma <- "/usr/local/bin/PSI-Sigma-1.9k/"

if (novel) {irmode <- 1} else {irmode <- 0}

mygtf <- vroom(gtf, col_types = "ccccccccc", col_names = F) %>% as.data.frame()
mygtf_reduced <- filter(mygtf, !X3 %in% c("gene"))
vroom_write(
  mygtf_reduced, glue('{dir_out}/PSI-Sigma/afolder/{gsub(".gtf", "", basename(gtf))}_reduced.gtf'),
  col_names = F, quote = "none", escape = "none"
)
gtf_reduced <- glue('{dir_out}/PSI-Sigma/afolder/{gsub(".gtf", "", basename(gtf))}_reduced.gtf')

cmds <- list()
cmds[["linking bams"]] <- c(
  glue("ln -s {sampleTable$files_bam} {dir_out}/PSI-Sigma/{sampleTable$samples}.Aligned.sortedByCoord.out.bam"),
  glue("ln -s {sampleTable$files_bam}.bai {dir_out}/PSI-Sigma/{sampleTable$samples}.Aligned.sortedByCoord.out.bam.bai")
)
cmds[["linking junction files"]] <- glue(
  "ln -s {dirname(sampleTable$files_bam)}/*SJ.out.tab {dir_out}/PSI-Sigma/{sampleTable$samples}.SJ.out.tab"
)
gtf_sorted <- gsub("\\.gtf", "_sorted\\.gtf", basename(gtf_reduced))
cmds[["sorting gtf"]] <- glue(
  '(grep "^#" {gtf_reduced}; grep -v "^#" {gtf_reduced} | sort -k1,1 -k4,4n) > {dir_out}/PSI-Sigma/afolder/{gtf_sorted}'
)
write.table(
  as.data.frame(glue("{basics$samples_1}.Aligned.sortedByCoord.out.bam")),
  glue("{dir_out}/PSI-Sigma/groupa.txt"),
  col.names = F, row.names = F, quote = F
)
write.table(
  as.data.frame(glue("{basics$samples_2}.Aligned.sortedByCoord.out.bam")),
  glue("{dir_out}/PSI-Sigma/groupb.txt"),
  col.names = F, row.names = F, quote = F
)
dummyai <- paste(c(dir_psi_sigma, "dummyai.pl"), collapse = "/")
cmds[["running psi-sigma"]] <- glue(
  "cd {dir_out}/PSI-Sigma/ && \\
  perl {dummyai} \\
  --gtf {gtf_sorted} \\
  --name {basics$condition_1}_{basics$condition_2} \\
  --type 1 \\
  --irmode {irmode} \\
  --nread {nread} \\
  --fmode 3 \\
  --adjp 2 \\
  --denominator 1"
)


##### JuncBASE #####
## 0.05 ##
myCreatedir(glue("{dir_out}/JuncBASE/"))

dir_juncbase <- "~/bin/juncBASE-1.2-beta/"
db <- "juncbase"

write.table(
  dplyr::select(sampleTable, samples, files_bam), glue("{dir_out}/JuncBASE/bams.txt"),
  sep = "\t", col.names = F, row.names = F, quote = F
)

library_type <- c()
for (i in 1:nrow(sampleTable)) {
  if (sampleTable$strandedness[i] == "no") library_type[i] <- "fr-unstranded"
  if (sampleTable$strandedness[i] == "yes") library_type[i] <- "fr-secondstrand"
  if (sampleTable$strandedness[i] == "reverse") library_type[i] <- "fr-firststrand"
}

cmds <- list()
if (novel) {
  cmds[["building a de-novo transcript database (1): cufflinks"]] <- glue(
    "cufflinks \\
    -p {nproc} \\
    -o {dir_out}/JuncBASE/cufflinks/{sampleTable$samples} \\
    -g {gtf} \\
    -u \\
    --library-type {library_type} \\
    {sampleTable$files_bam}"
  )
  cmds[["building a de-novo transcript database (2): cuffmerge"]] <- glue(
    "find {dir_out}/JuncBASE/cufflinks/ | grep transcripts.gtf > {dir_out}/JuncBASE/cufflinks/cufflinks_denovo.txt && \\
    mkdir -p {dir_out}/JuncBASE/tmp/ && \\
    python {dir_juncbase}/other_scripts/runCuffmerge.py \\
    -f {dir_out}/JuncBASE/cufflinks/cufflinks_denovo.txt \\
    --no_single_exons
    -g {gtf} \\
    -o {dir_out}/JuncBASE/cuffmerge/ \\
    -s {fa} \\
    --tmp_dir {dir_out}/JuncBASE/tmp/"
  )
  cmds[["building a de-novo transcript database (3): build DB"]] <- glue(
    "python {dir_juncbase}/build_DB_FromGTF.py \\
    -g {dir_out}/JuncBASE/cuffmerge/merged.gtf \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/ \\
    --initialize && \\
    python {dir_juncbase}/build_DB_FromGTF.py \\
    -g {dir_out}/JuncBASE/cuffmerge/merged.gtf \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/"
  )
} else {
  cmds[["building annotated AS DB in mysql"]] <- glue(
    "python {dir_juncbase}/build_DB_FromGTF.py \\
    -g {gtf} \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/ \\
    --initialize && \\
    python {dir_juncbase}/build_DB_FromGTF.py \\
    -g {gtf} \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/"
  )
}

cmds[["processing bam files"]] <- glue(
  "python {dir_juncbase}/run_preProcess_by_chr_step1.py \\
  -i {dir_out}/JuncBASE/bams.txt \\
  -o {dir_out}/JuncBASE/input_files/ \\
  -p {nrow(sampleTable)} \\
  --preProcess_options \"\" \\
  --nice"
)

cmds[["disambiguating splice junction orietations"]] <- glue(
  "python {dir_juncbase}/disambiguate_junctions.py \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -g {fa} \\
  --by_chr"
)

cmds[["identifying all junctions"]] <- glue(
  "python {dir_juncbase}/preProcess_getASEventReadCounts_by_chr_step2.py \\
  -i {dir_out}/JuncBASE/input_files/ \\
  --by_chr"
)

cmds[["creating exon-intron junction count files"]] <- glue(
  "python {dir_juncbase}/run_preProcess_step3_by_chr.py \\
  --input_dir {dir_out}/JuncBASE/input_files/ \\
  --num_processes {nrow(sampleTable)}"
)

cmds[["creating a pseudo all junction sample"]] <- glue(
  "python {dir_juncbase}/createPseudoSample.py \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -s {sample(sampleTable$samples, 1)} \\
  --by_chr"
)

cmds[["identifying alternative splicing events and quantify events from each sample"]] <- glue(
  "python {dir_juncbase}/run_getASEventReadCounts_multiSample.py \\
  -s {paste(sampleTable$samples, collapse = ',')} \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -o {dir_out}/JuncBASE/aseventreadcounts/ \\
  --txt_db1 {db} \\
  --txt_db2 {db} \\
  --txt_db3 {db} \\
  --sqlite_db_dir {dir_out}/JuncBASE/ \\
  --jcn_seq_len {(basics$read_length_most - 6) * 2} \\
  -p {nrow(sampleTable)} \\
  --nice \\
  --by_chr"
)

myCreatedir(glue("{dir_out}/JuncBASE/tables/"))

cmds[["creating tables of raw and length-normalized read counts of exclusion and inclusion isoforms"]] <- glue(
  "python {dir_juncbase}/run_createAS_CountTables.py \\
  -d {dir_out}/JuncBASE/aseventreadcounts/ \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -s {paste(sampleTable$samples, collapse = ',')} \\
  --num_processes {nrow(sampleTable)} \\
  --jcn_seq_len {(basics$read_length_most - 6) * 2} && \\
  python {dir_juncbase}/combine_createAS_CountTables_by_chr.py \\
  -d {dir_out}/JuncBASE/aseventreadcounts/ \\
  -o {dir_out}/JuncBASE/tables/juncbase"
)

cmds[["differential splicing between two groups of samples"]] <- glue(
  "python {dir_juncbase}/compareSampleSets.py \\
  --in_prefix {dir_out}/JuncBASE/tables/juncbase \\
  --all_psi_output {dir_out}/JuncBASE/sample_set_comparison.txt \\
  --mt_correction BH \\
  --thresh 10 \\
  --delta_thresh 0.0 \\
  --sample_set1 {paste(basics$samples_1, collapse = ',')} \\
  --sample_set2 {paste(basics$samples_2, collapse = ',')} \\
  --html_dir {dir_out}/JuncBASE/html/ \\
  --html_out_sign_thresh 1.0"
)


##### JUM #####
myCreatedir(glue("{dir_out}/JUM/"))

dir_jum <- "~/bin/JUM_2.0.2/"

for (i in 1:nrow(sampleTable)) {
  file.symlink(
    sampleTable$files_bam[i],
    glue("{dir_out}/JUM/{sampleTable$samples[i]}Aligned.out_sorted.bam")
  )
  
  file.symlink(
    glue("{dirname(sampleTable$files_bam[i])}/SJ.out.tab"),
    glue("{dir_out}/JUM/{sampleTable$samples[i]}SJ.out.tab")
  )
  
  cmd <- glue(
    "~/bin/anaconda3/envs/assw/bin/samtools view \\
    -h \\
    -o {dir_out}/JUM/{sampleTable$samples[i]}Aligned.out.sam \\
    {dir_out}/JUM/{sampleTable$samples[i]}Aligned.out_sorted.bam"
  )
  system(cmd)
}

cmds <- list()

cmds[["JUM A"]] <- glue(
  "cd {dir_out}/JUM/ && \\
  bash {dir_jum}/JUM_A.sh \\
  --Folder {dir_jum} \\
  --JuncThreshold 5 \\
  --Condition1_fileNum_threshold 2 \\
  --Condition2_fileNum_threshold 2 \\
  --IRthreshold 5 \\
  --Readlength 101 \\
  --Thread {floor(nproc/nrow(sampleTable))} \\
  --Condition1SampleName {paste(sampleTable$samples[sampleTable$conditions == 's1'], collapse = ',')} \\
  --Condition2SampleName {paste(sampleTable$samples[sampleTable$conditions == 's2'], collapse = ',')}"
)

ed <- data.frame(
  row.names = sampleTable$samples, condition = sampleTable$conditions
)
write.table(ed, glue("{dir_out}/JUM/JUM_diff/experiment_design.txt"), sep = "\t", quote = F)

cmds[["JUM R"]] <- glue(
  "cd {dir_out}/JUM/JUM_diff/ && \\
  Rscript {dir_jum}/R_script_JUM.R experiment_design.txt \\
  > outputFile.Rout 2> errorFile.Rout"
)

cmds[["JUM B"]] <- glue(
  "cd {dir_out}/JUM/JUM_diff/ && \\
  bash {dir_jum}/JUM_B.sh \\
  --Folder {dir_jum} \\
  --Test adjusted_pvalue \\
  --Cutoff 1 \\
  --TotalFileNum {nrow(sampleTable)} \\
  --Condition1_fileNum_threshold 2 \\
  --Condition2_fileNum_threshold 2 \\
  --Condition1SampleName {paste(sampleTable$samples[sampleTable$conditions == 's1'], collapse = ',')} \\
  --Condition2SampleName {paste(sampleTable$samples[sampleTable$conditions == 's2'], collapse = ',')}"
)

cmds[["JUM C"]] <- glue(
  "cd {dir_out}/JUM/JUM_diff/FINAL_JUM_OUTPUT_adjusted_pvalue_1/ && \\
  bash {dir_jum}/JUM_C.sh \\
  --Folder {dir_jum} \\
  --Test adjusted_pvalue \\
  --Cutoff 1 \\
  --TotalCondition1FileNum 4 \\
  --TotalCondition2FileNum 4 \\
  --REF ~/doc/reference/genepred/gencode_v32.primary_assembly.protein_coding2.reordered.gpd"
)


##### DARTS #####
myCreatedir(glue("{dir_out}/DARTS/{c('BHT', 'DNN', 'BHT_DNN')}"))

dir_rmats <- glue("{dir_out}/rMATS/")

chain <- "~/bin/liftOver/hg38ToHg19.over.chain"
chain_hg19 <- "~/bin/liftOver/hg19ToHg38.over.chain"

rmats_count <- vector("character", length = 5) |>
  setNames(c("SE", "RI", "A3SS", "A5SS", "MXE"))
for (n in names(rmats_count)) {
  rmats_count[[n]] <- list.files(
    dir_rmats, pattern = glue("JC.raw.input.{n}.txt"), full.names = T
  )
}
rmats_anno <- vector("character", length = 5) |>
  setNames(c("SE", "RI", "A3SS", "A5SS", "MXE"))
for (n in names(rmats_anno)) {
  rmats_anno[[n]] <- list.files(
    dir_rmats, pattern = glue("fromGTF.{n}.txt"), full.names = T
  )
}

cmds <- list()

for (s in c("SE", "RI", "A3SS", "A5SS")) {
  cmds[[glue("darts_bht {s}")]] <- glue(
    "Darts_BHT bayes_infer \\
    --rmats-count {rmats_count[s]} \\
    --od {dir_out}/DARTS/BHT/ \\
    --annot {rmats_anno[s]} \\
    -c 0.05 \\
    --replicate-model unpaired \\
    --nthread {nproc} \\
    -t {s} \\
    1>{dir_out}/DARTS/BHT/_log.bht.{s} 2>&1 && \\
    mv {dir_out}/DARTS/BHT/Darts_BHT.results.xlsx \\
    {dir_out}/DARTS/BHT/{s}.Darts_BHT.results.xlsx"
  )
}

for (s in c("SE", "RI", "A3SS", "A5SS")) {
  cmds[[glue("cut the first row to liftover {s}")]] <- glue(
    "cut -f 1 {dir_out}/DARTS/BHT/{s}.darts_bht.flat.txt > {dir_out}/DARTS/coord.{s}"
  )
}

for (s in c("SE", "RI", "A3SS", "A5SS")) {
  if (s == "SE") {
    tmp <- read.table(glue("{dir_out}/DARTS/coord.{s}"), sep = ":", skip = 1)
    tmp_all <- read.delim(glue("{dir_out}/DARTS/BHT/{s}.darts_bht.flat.txt"))
    
    write.table(
      bind_cols(tmp[, c(1, 3, 4)], tmp_all$ID, tmp[, 2]),
      glue("{dir_out}/DARTS/coord.{s}.1.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    write.table(
      bind_cols(tmp[, c(1, 5, 6)], tmp_all$ID),
      glue("{dir_out}/DARTS/coord.{s}.2.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    cmds[[glue("liftover {s}")]] <- glue(
      "liftOver \\
      {dir_out}/DARTS/coord.{s}.1.bed \\
      {chain} \\
      {dir_out}/DARTS/coord.{s}.1.hg19.bed \\
      {dir_out}/DARTS/coord.{s}.1.hg19.unmapped.bed && \\
      liftOver \\
      {dir_out}/DARTS/coord.{s}.2.bed \\
      {chain} \\
      {dir_out}/DARTS/coord.{s}.2.hg19.bed \\
      {dir_out}/DARTS/coord.{s}.2.hg19.unmapped.bed"
    )
    bed_1 <- read.delim(glue("{dir_out}/DARTS/coord.{s}.1.hg19.bed"), header = F)
    bed_2 <- read.delim(glue("{dir_out}/DARTS/coord.{s}.2.hg19.bed"), header = F)
    bed_f <- inner_join(bed_1, bed_2, by = c("V1", "V4"))
    bed_f$newID <- paste(
      bed_f[, 1], bed_f[, 5], bed_f[, 2], bed_f[, 3], bed_f[, 6], bed_f[, 7],
      sep = ":"
    )
    all_hg19 <- inner_join(tmp_all, bed_f[, c("V4", "newID")], by = c("ID" = "V4"))
    all_hg19 <- all_hg19 |>
      dplyr::mutate(ID = newID) |>
      dplyr::select(!newID)
    write.table(
      all_hg19, glue("{dir_out}/DARTS/BHT/{s}.darts_bht.flat.hg19.txt"),
      quote = F, row.names = F, col.names = T, sep = "\t"
    )
  }
  
  if (s %in% c("RI", "A3SS", "A5SS")) {
    tmp <- read.table(glue("{dir_out}/DARTS/coord.{s}"), sep = ":", skip = 1)
    tmp_all <- read.delim(glue("{dir_out}/DARTS/BHT/{s}.darts_bht.flat.txt"))
    
    write.table(
      bind_cols(tmp[, c(1, 3, 4)], tmp_all$ID, tmp[, 2]),
      glue("{dir_out}/DARTS/coord.{s}.1.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    write.table(
      bind_cols(tmp[, c(1, 5, 6)], tmp_all$ID),
      glue("{dir_out}/DARTS/coord.{s}.2.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    write.table(
      bind_cols(tmp[, c(1, 7, 8)], tmp_all$ID),
      glue("{dir_out}/DARTS/coord.{s}.3.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    cmds[[glue("liftover {s}")]] <- glue(
      "liftOver \\
      {dir_out}/DARTS/coord.{s}.1.bed \\
      {chain} \\
      {dir_out}/DARTS/coord.{s}.1.hg19.bed \\
      {dir_out}/DARTS/coord.{s}.1.hg19.unmapped.bed && \\
      liftOver \\
      {dir_out}/DARTS/coord.{s}.2.bed \\
      {chain} \\
      {dir_out}/DARTS/coord.{s}.2.hg19.bed \\
      {dir_out}/DARTS/coord.{s}.2.hg19.unmapped.bed && \\
      liftOver \\
      {dir_out}/DARTS/coord.{s}.3.bed \\
      {chain} \\
      {dir_out}/DARTS/coord.{s}.3.hg19.bed \\
      {dir_out}/DARTS/coord.{s}.3.hg19.unmapped.bed"
    )
    bed_1 <- read.delim(glue("{dir_out}/DARTS/coord.{s}.1.hg19.bed"), header = F)
    bed_2 <- read.delim(glue("{dir_out}/DARTS/coord.{s}.2.hg19.bed"), header = F)
    bed_3 <- read.delim(glue("{dir_out}/DARTS/coord.{s}.3.hg19.bed"), header = F)
    bed_f <- inner_join(bed_1, bed_2, by = c("V1", "V4"))
    bed_f <- inner_join(bed_f, bed_3, by = c("V1", "V4"))
    bed_f$newID <- paste(
      bed_f[, 1], bed_f[, 5], bed_f[, 2], bed_f[, 3], bed_f[, 6], bed_f[, 7], bed_f[, 8], bed_f[, 9],
      sep = ":"
    )
    all_hg19 <- inner_join(tmp_all, bed_f[, c("V4", "newID")], by = c("ID" = "V4"))
    all_hg19 <- all_hg19 |>
      dplyr::mutate(ID = newID) |>
      dplyr::select(!newID)
    write.table(
      all_hg19, glue("{dir_out}/DARTS/BHT/{s}.darts_bht.flat.hg19.txt"),
      quote = F, row.names = F, col.names = T, sep = "\t"
    )
  }
}

library(tximport)

tx2gene <- vroom(
  "~/doc/tx2gene/tx2gene_v32.txt", delim = "\t", col_types = cols(),
  col_names = c("chr", "gene_id", "gene_symbol", "transcript_id")
)

txi <- tximport(
  glue("{sampleTable$dirs_salmon}/quant.sf.gz") |> setNames(sampleTable$samples),
  type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_symbol")]
)
tpm_1 <- txi$abundance[, sampleTable$samples[sampleTable$conditions == "s1"]]
tpm_2 <- txi$abundance[, sampleTable$samples[sampleTable$conditions == "s2"]]

tpm_1 <- as.data.frame(tpm_1) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    mean = mean(c(sample1, sample2, sample3, sample4), na.rm = T)
  ) |>
  dplyr::ungroup() |>
  as.data.frame()
rownames(tpm_1) <- rownames(txi$abundance)
tpm_2 <- as.data.frame(tpm_2) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    mean = mean(c(sample5, sample6, sample7, sample8), na.rm = T)
  ) |>
  dplyr::ungroup() |>
  as.data.frame()
rownames(tpm_2) <- rownames(txi$abundance)
tpm_mean <- data.frame(row.names = rownames(tpm_1), condition_1 = tpm_1$mean, condition_2 = tpm_2$mean)

rbp <- read.delim("~/bin/DARTS/Darts_DNN/Darts_DNN/resources/rbp_gene_list.txt", header = F)
tpm_final <- data.frame(row.names = rbp[, 2])
tpm_final$s1 <- tpm_mean$condition_1[match(rownames(tpm_final), rownames(tpm_mean))]
tpm_final$s2 <- tpm_mean$condition_2[match(rownames(tpm_final), rownames(tpm_mean))]

write.table(
  tpm_final, glue("{dir_out}/DARTS/DNN/RBP_tpm.txt"),
  row.names = T, col.names = T, quote = F, sep = "\t"
)

for (s in c("SE", "RI", "A3SS", "A5SS")) {
  cmds[[glue("darts_dnn {s}")]] <- glue(
    "Darts_DNN predict \\
    -i {dir_out}/DARTS/BHT/{s}.darts_bht.flat.hg19.txt \\
    -e {dir_out}/DARTS/DNN/RBP_tpm.txt \\
    -o {dir_out}/DARTS/DNN/pred.{s}.txt \\
    -t {s} \\
    1>{dir_out}/DARTS/DNN/_log.dnn.{s} 2>&1"
  )
}


for (s in c("SE", "RI", "A3SS", "A5SS")) {
  cmds[[glue("cut the first row to liftover {s}")]] <- glue(
    "cut -f 1 {dir_out}/DARTS/DNN/pred.{s}.txt > {dir_out}/DARTS/coord.pred.{s}"
  )
}


for (s in c("SE", "RI", "A3SS", "A5SS")) {
  if (s == "SE") {
    tmp <- read.table(glue("{dir_out}/DARTS/coord.pred.{s}"), sep = ":", skip = 1)
    tmp_all <- read.delim(glue("{dir_out}/DARTS/DNN/pred.{s}.txt"))
    
    write.table(
      bind_cols(tmp[, c(1, 3, 4)], tmp_all$ID, tmp[, 2]),
      glue("{dir_out}/DARTS/coord.pred.{s}.1.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    write.table(
      bind_cols(tmp[, c(1, 5, 6)], tmp_all$ID),
      glue("{dir_out}/DARTS/coord.pred.{s}.2.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    cmds[[glue("liftover {s}")]] <- glue(
      "liftOver \\
      {dir_out}/DARTS/coord.pred.{s}.1.bed \\
      {chain_hg19} \\
      {dir_out}/DARTS/coord.pred.{s}.1.hg38.bed \\
      {dir_out}/DARTS/coord.pred.{s}.1.hg38.unmapped.bed && \\
      liftOver \\
      {dir_out}/DARTS/coord.pred.{s}.2.bed \\
      {chain_hg19} \\
      {dir_out}/DARTS/coord.pred.{s}.2.hg38.bed \\
      {dir_out}/DARTS/coord.pred.{s}.2.hg38.unmapped.bed"
    )
    bed_1 <- read.delim(glue("{dir_out}/DARTS/coord.pred.{s}.1.hg38.bed"), header = F)
    bed_2 <- read.delim(glue("{dir_out}/DARTS/coord.pred.{s}.2.hg38.bed"), header = F)
    bed_f <- inner_join(bed_1, bed_2, by = c("V1", "V4"))
    bed_f$newID <- paste(
      bed_f[, 1], bed_f[, 5], bed_f[, 2], bed_f[, 3], bed_f[, 6], bed_f[, 7],
      sep = ":"
    )
    all_hg38 <- inner_join(tmp_all, bed_f[, c("V4", "newID")], by = c("ID" = "V4"))
    all_hg38 <- all_hg38 |>
      dplyr::mutate(ID = newID) |>
      dplyr::select(!newID)
    write.table(
      all_hg38, glue("{dir_out}/DARTS/DNN/pred.{s}.hg38.txt"),
      quote = F, row.names = F, col.names = T, sep = "\t"
    )
  }
  
  if (s %in% c("RI", "A3SS", "A5SS")) {
    tmp <- read.table(glue("{dir_out}/DARTS/coord.pred.{s}"), sep = ":", skip = 1)
    tmp_all <- read.delim(glue("{dir_out}/DARTS/DNN/pred.{s}.txt"))
    
    write.table(
      bind_cols(tmp[, c(1, 3, 4)], tmp_all$ID, tmp[, 2]),
      glue("{dir_out}/DARTS/coord.pred.{s}.1.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    write.table(
      bind_cols(tmp[, c(1, 5, 6)], tmp_all$ID),
      glue("{dir_out}/DARTS/coord.pred.{s}.2.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    write.table(
      bind_cols(tmp[, c(1, 7, 8)], tmp_all$ID),
      glue("{dir_out}/DARTS/coord.pred.{s}.3.bed"),
      quote = F, row.names = F, col.names = F, sep = "\t"
    )
    cmds[[glue("liftover {s}")]] <- glue(
      "liftOver \\
      {dir_out}/DARTS/coord.pred.{s}.1.bed \\
      {chain_hg19} \\
      {dir_out}/DARTS/coord.pred.{s}.1.hg38.bed \\
      {dir_out}/DARTS/coord.pred.{s}.1.hg38.unmapped.bed && \\
      liftOver \\
      {dir_out}/DARTS/coord.pred.{s}.2.bed \\
      {chain_hg19} \\
      {dir_out}/DARTS/coord.pred.{s}.2.hg38.bed \\
      {dir_out}/DARTS/coord.pred.{s}.2.hg38.unmapped.bed && \\
      liftOver \\
      {dir_out}/DARTS/coord.pred.{s}.3.bed \\
      {chain_hg19} \\
      {dir_out}/DARTS/coord.pred.{s}.3.hg38.bed \\
      {dir_out}/DARTS/coord.pred.{s}.3.hg38.unmapped.bed"
    )
    bed_1 <- read.delim(glue("{dir_out}/DARTS/coord.pred.{s}.1.hg38.bed"), header = F)
    bed_2 <- read.delim(glue("{dir_out}/DARTS/coord.pred.{s}.2.hg38.bed"), header = F)
    bed_3 <- read.delim(glue("{dir_out}/DARTS/coord.pred.{s}.3.hg38.bed"), header = F)
    bed_f <- inner_join(bed_1, bed_2, by = c("V1", "V4"))
    bed_f <- inner_join(bed_f, bed_3, by = c("V1", "V4"))
    bed_f$newID <- paste(
      bed_f[, 1], bed_f[, 5], bed_f[, 2], bed_f[, 3], bed_f[, 6], bed_f[, 7], bed_f[, 8], bed_f[, 9],
      sep = ":"
    )
    all_hg38 <- inner_join(tmp_all, bed_f[, c("V4", "newID")], by = c("ID" = "V4"))
    all_hg38 <- all_hg38 |>
      dplyr::mutate(ID = newID) |>
      dplyr::select(!newID)
    write.table(
      all_hg38, glue("{dir_out}/DARTS/DNN/pred.{s}.hg38.txt"),
      quote = F, row.names = F, col.names = T, sep = "\t"
    )
  }
}

cmds[["rm coord files"]] <- glue(
  "rm -rf {dir_out}/DARTS/coord*"
)

for (s in c("SE", "RI", "A3SS", "A5SS")) {
  cmds[[glue("darts_bht_dnn {s}")]] <- glue(
    "Darts_BHT bayes_infer \\
    --rmats-count {rmats_count[s]} \\
    --od {dir_out}/DARTS/BHT_DNN/ \\
    --annot {rmats_anno[s]} \\
    --prior {dir_out}/DARTS/DNN/pred.{s}.hg38.txt \\
    -t {s} \\
    -c 0.05 \\
    --nthread {nproc} \\
    1>{dir_out}/DARTS/BHT_DNN/_log.bht_dnn.{s} 2>&1 && \\
    mv {dir_out}/DARTS/BHT_DNN/Darts_BHT.results.xlsx \\
    {dir_out}/DARTS/BHT_DNN/{s}.Darts_BHT.results.xlsx"
  )
}




