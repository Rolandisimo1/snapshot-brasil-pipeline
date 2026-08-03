## =============================================================================
## Module 6 — Species Interactions
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## SELF-CONTAINED: this script builds every table and figure in the rendered
## report from the raw combined-dataset files plus a small set of reference
## lookup tables that are themselves the output of one-time upstream work
## (species range mask, GJAM residual correlations from Module 4, and the
## array-level detection-history/covariate checkpoints from Module 5/6's
## array aggregation). Two classes of model are involved:
##
##   1. The 30 abundance-mediated interaction fits (paired N-mixture models,
##      NIMBLE/MCMC) are each expensive (up to 300,000 iterations across
##      three escalating retry rounds; the full batch takes many hours
##      sequentially, and the sandbox this pipeline was built in cannot run
##      more than one such fit at a time). By default this script LOADS the
##      precomputed coefficient table shipped alongside it rather than
##      refitting; set RUN_LIVE_PAIRS <- TRUE below to refit every pair from
##      scratch using fit_pair_array.R (shipped in mod6_code/), which this
##      script will source and call once per pair.
##   2. The GJAM residual correlation matrix (all 44 species) is Module 4's
##      output, not refit here; it is loaded as a reference lookup.
##
## INPUTS (place in DATA_DIR, default "data/", alongside this script):
##   all_pairs_array_level_coefs.csv   -- full posterior summary (mean/sd/lcl/ucl/Rhat)
##                                        for every parameter of all 30 array-level pairs
##   convergence_diagnostics_summary.csv -- per-pair max Rhat and convergence flag
##   gjam_residual_correlation.csv     -- 44x44 species residual correlation matrix
##                                        (Module 4 GJAM output, used as-is)
##   species_common_names.csv          -- scientific -> common name lookup
##   pair_topology.csv                 -- the 20 original pairs' detection fractions
##                                        and topology type, for reference
## Note: this module fits on Module 5's array-level covariate table, itself the same
## 11 covariates as the camera-level file subset down to the reduced 6-covariate
## ARRAY_SET at fit time -- not a separately reduced file.
## Set RUN_LIVE_PAIRS <- TRUE (and provide mod6_code/fit_pair_array.R plus the
## pair_inputs/*.rds files it expects) to refit from scratch instead.
## =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})

DATA_DIR <- "data"
FIGS_DIR <- "figs"
dir.create(FIGS_DIR, showWarnings = FALSE)

RUN_LIVE_PAIRS <- FALSE   # TRUE re-fits all 30 pairs via NIMBLE/MCMC (many hours); FALSE loads precomputed

## =============================================================================
## STEP 1: LOAD RESULTS (precomputed by default; refit if RUN_LIVE_PAIRS)
## =============================================================================
if (RUN_LIVE_PAIRS) {
  source(file.path("mod6_code", "fit_pair_array.R"))  # expects mod6_array/pair_inputs/*.rds
  stop("Live refit path invoked -- see fit_pair_array.R and the module's methods note ",
       "for the full escalating-retry MCMC procedure; this path is not run automatically ",
       "because a single pair can take hours and pairs must be fit strictly sequentially.")
} else {
  all_pairs <- fread(file.path(DATA_DIR, "all_pairs_array_level_coefs.csv"))
  conv <- fread(file.path(DATA_DIR, "convergence_diagnostics_summary.csv"))
  cat("Loaded precomputed array-level fits:", uniqueN(all_pairs$pair_idx, na.rm=TRUE) * 3, "topology groups,",
      "converged:", sum(conv$converged), "of", nrow(conv), "pairs\n")
}

common_names <- fread(file.path(DATA_DIR, "species_common_names.csv"))
name_map <- setNames(common_names$common, common_names$species)
to_common <- function(sci) ifelse(sci %in% names(name_map), name_map[sci], sci)

## =============================================================================
## STEP 2: THE INTERACTION COEFFICIENT (gamma0) FOR ALL 30 PAIRS
## =============================================================================
## gamma0 is the driver's log-abundance effect on the responder's log-abundance
## in the paired N-mixture model -- the actual interaction estimate (gamma1,
## gamma2, gamma3 are temperature-interaction terms tested separately; see the
## module text for why they are not the headline quantity).
gamma0 <- all_pairs[param == "gamma0"]
gamma0 <- merge(gamma0, conv[, .(driver, responder, source, converged)],
                 by = c("driver", "responder", "source"), all.x = TRUE)
gamma0[, sig := (lcl > 0) | (ucl < 0)]
gamma0[, driver_c := to_common(driver)]
gamma0[, responder_c := to_common(responder)]
cat("Significant pairs (95% CI excludes 0):", sum(gamma0$sig), "of", nrow(gamma0), "\n")

sig_pairs <- gamma0[sig == TRUE, .(driver_c, responder_c, source, mean, lcl, ucl)]
fwrite(sig_pairs, "gamma0_significant_pairs.csv")

## =============================================================================
## STEP 3: FOREST PLOT -- ALL 30 PAIRS, RANKED
## =============================================================================
gamma0[, label := paste0(driver_c, " \u2192 ", responder_c)]
gamma0_plot <- gamma0[order(mean)]
gamma0_plot[, y := .I]
gamma0_plot[, col := fifelse(!sig, "grey80", fifelse(mean > 0, "#2166ac", "#b2182b"))]

p_forest <- ggplot(gamma0_plot, aes(y = reorder(label, mean))) +
  geom_segment(aes(x = lcl, xend = ucl, yend = reorder(label, mean), color = col), linewidth = 1) +
  geom_point(aes(x = mean, color = col), size = 2) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_color_identity() +
  labs(x = "Interaction coefficient (\u03b30): driver's log-abundance effect on responder's log-abundance",
       y = NULL, title = "Abundance-mediated interaction strength, all 30 array-level pairs") +
  theme_minimal(base_size = 9)
ggsave(file.path(FIGS_DIR, "gamma0_forest_all_pairs.png"), p_forest, width = 9, height = 11, dpi = 150)
cat("Saved figs/gamma0_forest_all_pairs.png\n")

## =============================================================================
## STEP 4: THE ABUNDANCE-MEDIATED INTERACTION NETWORK
## =============================================================================
suppressMessages(library(igraph))
edges_abu <- gamma0[, .(from = driver_c, to = responder_c, sig, sign = sign(mean))]
g_abu <- graph_from_data_frame(edges_abu, directed = TRUE)
E(g_abu)$color <- ifelse(!edges_abu$sig, "grey85", ifelse(edges_abu$sign > 0, "#2166ac", "#b2182b"))
E(g_abu)$width <- ifelse(edges_abu$sig, 2.5, 0.6)
E(g_abu)$lty   <- ifelse(edges_abu$sig, 1, 2)
set.seed(42)
layout_abu <- layout_with_fr(g_abu)

png(file.path(FIGS_DIR, "interaction_network_abundance.png"), width = 1600, height = 1300, res = 140)
par(mar = c(4, 1, 3, 1))
plot(g_abu, layout = layout_abu, vertex.size = 14, vertex.color = "#4575b4",
     vertex.label.cex = 0.75, vertex.label.color = "black", edge.arrow.size = 0.5,
     main = "Abundance-mediated species interaction network \u2014 array level, 30 driver\u2192responder pairs")
legend("bottomleft", legend = c("Significant, positive", "Significant, negative", "Not significant"),
       col = c("#2166ac", "#b2182b", "grey85"), lty = c(1,1,2), lwd = c(2.5,2.5,0.6), bty = "n", cex = 0.8)
dev.off()
cat("Saved figs/interaction_network_abundance.png\n")

## =============================================================================
## STEP 5: THE GJAM RESIDUAL CO-OCCURRENCE NETWORK
## =============================================================================
gjam_corr <- fread(file.path(DATA_DIR, "gjam_residual_correlation.csv"))
sp_gjam <- gjam_corr$species
mat_gjam <- as.matrix(gjam_corr[, -"species"])
rownames(mat_gjam) <- sp_gjam

THRESH <- 0.4
edges_gjam <- data.table()
for (i in seq_along(sp_gjam)) for (j in seq_len(i-1)) {
  v <- mat_gjam[i, j]
  if (!is.na(v) && abs(v) > THRESH) {
    edges_gjam <- rbind(edges_gjam, data.table(from = to_common(sp_gjam[i]), to = to_common(sp_gjam[j]), r = v))
  }
}
cat("GJAM edges above |r|>", THRESH, ":", nrow(edges_gjam), "of", choose(length(sp_gjam), 2), "possible pairs\n")

g_gjam <- graph_from_data_frame(edges_gjam, directed = FALSE)
E(g_gjam)$width <- abs(edges_gjam$r) * 4
E(g_gjam)$color <- "#2166ac"
set.seed(11)
layout_gjam <- layout_with_fr(g_gjam)

png(file.path(FIGS_DIR, "gjam_network_full.png"), width = 1400, height = 1200, res = 140)
par(mar = c(1, 1, 3, 1))
plot(g_gjam, layout = layout_gjam, vertex.size = 16, vertex.color = "#66c2a5",
     vertex.label.cex = 0.8, vertex.label.color = "black",
     main = paste0("GJAM residual co-occurrence network (|r| > ", THRESH, ", habitat effects removed)"))
dev.off()
cat("Saved figs/gjam_network_full.png\n")

cat("\n=== Module 6 script complete: all figures + tables built from precomputed fits. ===\n")
cat("See module6_species_interactions.qmd for the narrated analysis and citations.\n")
