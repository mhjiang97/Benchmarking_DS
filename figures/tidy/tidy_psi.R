source("utils_tidy_psi.R")

dir_out <- "~/projects/AS/analysis/sim_1/"


##### rMATS #####
dir_rmats <- glue("{dir_out}/rMATS/")

list_rmats_psi <- psi_rmats(dir_rmats)
rmats_psi <- bind_rows(list_rmats_psi)

rmats <- readRDS(glue("{dir_rmats}/rmats.rds"))
rmats_psi$idx <- rmats$idx_truths[match(rmats_psi$myID, rmats$myID)]
saveRDS(rmats_psi, glue("{dir_rmats}/rmats_psi.rds"))


##### LeafCutter #####
dir_leafcutter <- glue("{dir_out}/LeafCutter/")

leafcutter_psi <- psi_leafcutter(glue("{dir_leafcutter}/cluster/leafcutter_perind.counts.gz"))
leafcutter <- readRDS(glue("{dir_leafcutter}/leafcutter.rds"))
leafcutter_psi$idx <- leafcutter$idx_truths[match(leafcutter_psi$myID, leafcutter$myID)]
saveRDS(leafcutter_psi, glue("{dir_leafcutter}/leafcutter_psi.rds"))


##### SplAdder #####
dir_spladder <- glue("{dir_out}/SplAdder/")

list_spladder_psi <- psi_spladder(dir_event = glue("{dir_spladder}/graphs/"))
spladder_psi <- bind_rows(list_spladder_psi)
nms <- str_replace_all(
  colnames(spladder_psi), pattern = ".SortedByCoord:psi", replacement = ""
)
colnames(spladder_psi) <- nms
spladder <- readRDS(glue("{dir_spladder}/spladder.rds"))
spladder_psi$idx <- spladder$idx_truths[match(spladder_psi$myID, spladder$myID)]
saveRDS(spladder_psi, glue("{dir_spladder}/spladder_psi.rds"))


##### SUPPA #####
dir_suppa <- glue("{dir_out}/SUPPA/")

list_suppa_psi <- psi_suppa(glue("{dir_suppa}/psi/"), s1 = "s1", s2 = "s2")
suppa_psi <- bind_rows(list_suppa_psi)
suppa <- readRDS(glue("{dir_suppa}/suppa.rds"))
suppa_psi$idx <- suppa$idx_truths[match(suppa_psi$myID, suppa$myID)]
saveRDS(suppa_psi, glue("{dir_suppa}/suppa_psi.rds"))


##### Whippet #####
dir_whippet <- glue("{dir_out}/Whippet/")

whippet_psi <- psi_whippet(glue("{dir_whippet}/quant/"), samples = paste0("sample", 1:8))
whippet <- readRDS(glue("{dir_whippet}/whippet.rds"))
whippet_psi$idx <- whippet$idx_truths_2[match(whippet_psi$myID, whippet$myID)]
saveRDS(whippet_psi, glue("{dir_whippet}/whippet_psi.rds"))


##### psichomics #####
dir_psichomics <- glue("{dir_out}/psichomics/")

psichomics_psi <- psi_psichomics(glue("{dir_psichomics}/psi.txt"))
psichomics <- readRDS(glue("{dir_psichomics}/psichomics.rds"))
psichomics_psi$idx <- psichomics$idx_truths[match(psichomics_psi$myID, psichomics$myID)]
saveRDS(psichomics_psi, glue("{dir_psichomics}/psichomics_psi.rds"))




