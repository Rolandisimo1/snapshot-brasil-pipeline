## =============================================================================
## Module 5 — Prediction & Mapping
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## SELF-CONTAINED: this script turns each species' fitted occupancy model
## (from Module 3/4) into a national prediction across a 17,766-cell grid,
## screens which species are worth mapping (AUC), flags where the model is
## extrapolating (MESS), applies the IUCN range mask, and turns the
## continuous occupancy surface into a three-outcome present/absent/
## not-assessed map at both the camera level (11 covariates) and the array
## level (6 covariates).
##
## Required input files (in DATA_DIR, default "data" alongside this script):
##   national_grid_covariates_11cov.csv  -- 17,766-cell national grid, the same 11
##                                           camera-level covariates as final_covariates.csv,
##                                           extracted from the identical Earth Engine
##                                           sources/definitions, plus MESS/extrapolating
##   all_species_pgocc_coefs.csv         -- 44-species camera-level PGOcc coefficients
##                                           (fit in Module 4, reused here as-is)
##   array_covariates.csv                -- 60-array covariate table (same 11 covariates as
##                                           the camera-level file; the array-level model below
##                                           subsets it down to the reduced 6-covariate
##                                           ARRAY_SET -- forest, savanna, native vegetation,
##                                           temperature, precipitation, distance to road --
##                                           via ARRAY_COVARS, since the array's coarser spatial
##                                           footprint doesn't support all 11 reliably)
##   ysum_array.csv / K_array.csv         -- 60 x 44 array-level detection-count matrices
##                                           (occasion-collapsed; reference outputs from the
##                                           array-level detection-history assembler)
##   array_ids.csv                       -- the 60 array IDs, in ysum/K row order
##   species_list.json                   -- the 44 modeled species/genera
##   grid_inrange_44sp.csv                -- per-cell IUCN range mask (44 species x 17,766 grid
##                                           cells), output of the Module 2.2 range-map exercise
##                                           applied to grid coordinates instead of camera sites
##   final_covariates.csv                 -- camera-site covariates, used to recover the
##                                           training-scale standardization constants AND (with
##                                           the two files below) to compute AUC/thresholds
##   final_deployments.csv, final_detections_mammals.csv -- raw deployment/detection tables,
##                                           needed to rebuild site-level presence/absence for
##                                           the AUC and threshold computation
##   species_range_mask.csv                -- site x species IUCN range mask (camera-level;
##                                           1 = candidate, 0 = structural zero), used to
##                                           restrict the AUC/threshold fit to in-range sites
##
## The array-level detection matrices and grid range mask are reference outputs from
## earlier modules' own assemblers and cannot be rebuilt from the three main combined-
## dataset files alone -- they are shipped alongside this script as fixed lookup tables.
##
## Requires: data.table, unmarked, jsonlite, ggplot2.

suppressMessages({
  library(data.table)
  library(unmarked)
  library(jsonlite)
})

DATA_DIR <- "data"
FIGS_DIR <- "figs"
dir.create(FIGS_DIR, showWarnings = FALSE)

CAMERA_COVARS <- c("forest_100m","savanna_100m","pasture_100m","cropland_100m",
                   "native_veg_1000m","temp_mean_C","precip_annual_mm","in_pa",
                   "ghsl_built_5000m","dist_road_m","dist_water_m")
ARRAY_COVARS  <- c("forest_100m","savanna_100m","native_veg_1000m",
                   "temp_mean_C","precip_annual_mm","dist_road_m")

## =============================================================================
## STEP 1: LOAD THE GRID AND CAMERA-LEVEL COEFFICIENTS
## =============================================================================
grid <- fread(file.path(DATA_DIR, "national_grid_covariates_11cov.csv"))
cat("Grid cells:", nrow(grid), "\n")
cat("Grid cells extrapolating (11-covariate MESS < 0):", sum(grid$extrapolating, na.rm=TRUE),
    sprintf("(%.1f%% of the grid)\n", 100*mean(grid$extrapolating, na.rm=TRUE)))

pgocc <- fread(file.path(DATA_DIR, "all_species_pgocc_coefs.csv"))
species <- fromJSON(file.path(DATA_DIR, "species_list.json"))
fc <- fread(file.path(DATA_DIR, "final_covariates.csv"))
fc_complete <- fc[complete.cases(fc[, ..CAMERA_COVARS])]

## Recover the training-scale standardization constants (mean/sd) exactly as
## used when the PGOcc models were fit, so the grid is standardized on the
## SAME scale as the coefficients, not the grid's own mean/sd.
znames <- c("forest_z","savanna_z","pasture_z","cropland_z","nveg_z","temp_z",
            "precip_z","in_pa","built_z","distroad_z","distwater_z")
site_scaler <- list()
for (i in seq_along(CAMERA_COVARS)) {
  zc <- znames[i]; raw <- CAMERA_COVARS[i]
  if (zc == "in_pa") next
  site_scaler[[zc]] <- list(raw = raw, mean = mean(fc_complete[[raw]]), sd = sd(fc_complete[[raw]]))
}

## =============================================================================
## STEP 2: PROJECT CAMERA-LEVEL PSI ACROSS THE GRID
## =============================================================================
sigmoid <- function(x) 1 / (1 + exp(-x))

pgocc_wide <- dcast(pgocc, species ~ param, value.var = "mean")
setDF(pgocc_wide); rownames(pgocc_wide) <- pgocc_wide$species

for (zc in znames) {
  if (zc == "in_pa") { grid[[paste0("z_", zc)]] <- grid$in_pa; next }
  s <- site_scaler[[zc]]
  grid[[paste0("z_", zc)]] <- (grid[[s$raw]] - s$mean) / s$sd
}
Xgrid <- as.matrix(grid[, paste0("z_", znames), with = FALSE])

psi_camera <- matrix(NA_real_, nrow = nrow(grid), ncol = length(species),
                      dimnames = list(grid$sid, species))
for (sci in species) {
  b0 <- pgocc_wide[sci, "(Intercept)"]
  betas <- as.numeric(pgocc_wide[sci, znames])
  lp <- b0 + Xgrid %*% betas
  psi_camera[, sci] <- sigmoid(lp)
}
cat("Camera-level psi computed for", length(species), "species x", nrow(grid), "grid cells\n")

## =============================================================================
## STEP 3: REFIT THE ARRAY-LEVEL SINGLE-SPECIES OCCUPANCY MODEL (6-COVARIATE SET)
## =============================================================================
## The array-level model is refit from scratch here (unlike the camera-level
## model, which reuses Module 4's saved coefficients) because it was never
## saved from an earlier module -- Module 3/4 only fit camera-level single-
## species occupancy.
arr_cov <- fread(file.path(DATA_DIR, "array_covariates.csv"))
array_ids <- fread(file.path(DATA_DIR, "array_ids.csv"), header = FALSE)$V1
arr_cov <- arr_cov[match(array_ids, array_id)]
arr_scaler <- list()
Xz_arr <- data.table(array_id = array_ids)
for (c in ARRAY_COVARS) {
  m <- mean(arr_cov[[c]]); s <- sd(arr_cov[[c]])
  arr_scaler[[c]] <- list(mean = m, sd = s)
  Xz_arr[[paste0(c, "_z")]] <- (arr_cov[[c]] - m) / s
}
Xz_arr_df <- as.data.frame(Xz_arr[, -1]); rownames(Xz_arr_df) <- array_ids

ysum_arr <- as.matrix(fread(file.path(DATA_DIR, "ysum_array.csv"), header = FALSE))
K_arr <- as.matrix(fread(file.path(DATA_DIR, "K_array.csv"), header = FALSE))

occ_form <- as.formula(paste("~1 ~", paste(paste0(ARRAY_COVARS, "_z"), collapse = "+")))
arr_coefs <- list(); arr_fitsum <- list(); arr_psi_pred <- list()
for (i in seq_along(species)) {
  sci <- species[i]
  n_occ_max <- max(K_arr[, i])
  y <- matrix(NA_real_, nrow = 60, ncol = n_occ_max)
  for (r in 1:60) {
    k <- K_arr[r, i]; ys <- ysum_arr[r, i]
    if (k > 0) { y[r, 1:k] <- 0; if (ys > 0) y[r, 1:ys] <- 1 }
  }
  keep <- colSums(!is.na(y)) > 0; y <- y[, keep, drop = FALSE]
  naive <- sum(rowSums(y, na.rm = TRUE) > 0)
  if (naive < 3) { arr_fitsum[[sci]] <- data.frame(species=sci, converged=FALSE, note="too_few"); next }
  umf <- tryCatch(unmarkedFrameOccu(y = y, siteCovs = Xz_arr_df), error = function(e) NULL)
  m <- tryCatch(occu(occ_form, umf), error = function(e) NULL)
  if (is.null(m) || is.na(logLik(m))) { arr_fitsum[[sci]] <- data.frame(species=sci, converged=FALSE, note="failed"); next }
  se <- SE(m, type = "state")
  unstable <- any(se > 5, na.rm = TRUE)
  arr_coefs[[sci]] <- coef(m, type = "state")
  arr_fitsum[[sci]] <- data.frame(species=sci, converged=TRUE, unstable_se=unstable, note="ok")
  ps <- predict(m, type = "state")$Predicted
  arr_psi_pred[[sci]] <- data.frame(species=sci, array_id=rownames(Xz_arr_df), psi=ps,
                                     detected = as.integer(rowSums(y, na.rm=TRUE) > 0))
}
arr_fitsum_df <- rbindlist(arr_fitsum, fill = TRUE)
projectable_arr <- arr_fitsum_df[converged == TRUE & unstable_se == FALSE, species]
cat("Array-level: ", sum(arr_fitsum_df$converged), "of", length(species), "converged;",
    length(projectable_arr), "stable + projectable\n")

## Project array-level psi across the grid for the projectable species
for (c in ARRAY_COVARS) grid[[paste0("za_", c)]] <- (grid[[c]] - arr_scaler[[c]]$mean) / arr_scaler[[c]]$sd
Xgrid_arr <- as.matrix(grid[, paste0("za_", ARRAY_COVARS), with = FALSE])
psi_array <- matrix(NA_real_, nrow = nrow(grid), ncol = length(projectable_arr),
                     dimnames = list(grid$sid, projectable_arr))
for (sci in projectable_arr) {
  betas <- arr_coefs[[sci]]
  b0 <- betas["psi(Int)"]
  bx <- betas[paste0("psi(", ARRAY_COVARS, "_z)")]
  lp <- b0 + Xgrid_arr %*% bx
  psi_array[, sci] <- sigmoid(lp)
}

## =============================================================================
## STEP 4: AUC SCREEN, THRESHOLD, RANGE MASK, AND THE THREE-OUTCOME STATUS MAP
## =============================================================================

## --- 4.1: rebuild site-level detection histories to compute training AUC ---
dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
det <- fread(file.path(DATA_DIR, "final_detections_mammals.csv"))
det[, sci := fifelse(!is.na(sci_mdd), sci_mdd, scientific_name)]
det[, date := as.IDate(date)]
dep[, start_date := as.IDate(start_date)]
dep[, end_date := as.IDate(end_date)]
site_window <- dep[, .(start = min(start_date), end = max(end_date)), by = site_id]

site_ids_new <- fc_complete[site_id %in% site_window$site_id, site_id]
sw <- site_window[match(site_ids_new, site_id)]

det_m <- merge(det[!is.na(date)], site_window, by = "site_id")
det_m <- det_m[date >= start & date <= end]

rmask <- fread(file.path(DATA_DIR, "species_range_mask.csv"))
setnames(rmask, 1, "site_id")
rmask <- rmask[match(site_ids_new, site_id)]

## --- 4.2: per-species empirical AUC + Youden (MaxSens+Spec) threshold ---
roc_auc <- function(y, p) {
  # Mann-Whitney U formulation: AUC = (sum of ranks of the positive class - n1*(n1+1)/2) / (n1*n0).
  # y and p must stay in the SAME (original) order throughout -- do not permute one without the other.
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(p)
  sum_r1 <- sum(ranks[y == 1])
  (sum_r1 - n1 * (n1 + 1) / 2) / (n1 * n0)
}
youden_thr <- function(y, p) {
  if (length(unique(y)) < 2) return(NA_real_)
  thrs <- sort(unique(p))
  best_j <- -Inf; best_t <- NA_real_
  for (t in thrs) {
    pred <- p >= t
    sens <- sum(pred & y == 1) / sum(y == 1)
    spec <- sum(!pred & y == 0) / sum(y == 0)
    j <- sens + spec - 1
    if (j > best_j) { best_j <- j; best_t <- t }
  }
  best_t
}

fc_sites <- fc_complete[match(site_ids_new, site_id)]
Xsite_cols <- lapply(znames, function(zc) {
  if (zc == "in_pa") return(fc_sites$in_pa)
  s <- site_scaler[[zc]]
  (fc_sites[[s$raw]] - s$mean) / s$sd
})
Xsite <- do.call(cbind, Xsite_cols)
colnames(Xsite) <- znames

auc_results <- data.table(species = species, auc_camera = NA_real_, thr_camera = NA_real_)
for (i in seq_along(species)) {
  target_sci <- species[i]  # named distinctly from det_m's own "sci" column to avoid
                             # a data.table scoping bug where `sci == sci` would silently
                             # compare the column to itself (always TRUE) instead of to
                             # this loop variable
  present_sites <- unique(det_m[sci == target_sci]$site_id)
  y <- as.integer(site_ids_new %in% present_sites)
  in_range <- if (target_sci %in% names(rmask)) rmask[[target_sci]] != 0 else rep(TRUE, length(site_ids_new))
  in_range[is.na(in_range)] <- TRUE
  b0 <- pgocc_wide[target_sci, "(Intercept)"]; betas <- as.numeric(pgocc_wide[target_sci, znames])
  psi_site <- sigmoid(b0 + Xsite %*% betas)
  ym <- y[in_range]; pm <- psi_site[in_range]
  auc_results[i, auc_camera := roc_auc(ym, pm)]
  auc_results[i, thr_camera := youden_thr(ym, pm)]
}
cat("Median camera-level AUC:", round(median(auc_results$auc_camera, na.rm=TRUE), 3), "\n")

## --- 4.3: range mask + MESS -> three-outcome status per grid cell ---
grid_inrange <- fread(file.path(DATA_DIR, "grid_inrange_44sp.csv"))
## grid_inrange only has a row for cells with a resolvable range-mask value (fewer than
## the full grid, e.g. cells outside every species' assessed extent) -- reindex it to
## the grid's own sid order so a species column lines up cell-for-cell, treating any
## grid cell absent from grid_inrange as NOT in range (conservative: unassessed rather
## than silently assumed present).
grid_inrange_aligned <- grid_inrange[match(grid$sid, sid)]

status_camera <- matrix(0L, nrow = nrow(grid), ncol = length(species),
                         dimnames = list(grid$sid, species))
for (i in seq_along(species)) {
  sci <- species[i]
  thr <- auc_results$thr_camera[i]
  in_range_grid <- if (sci %in% names(grid_inrange_aligned)) as.logical(grid_inrange_aligned[[sci]]) else rep(TRUE, nrow(grid))
  in_range_grid[is.na(in_range_grid)] <- FALSE
  extrap <- grid$extrapolating; extrap[is.na(extrap)] <- TRUE
  psi_i <- psi_camera[, sci]
  status_camera[, i] <- ifelse(!in_range_grid | extrap | is.na(thr), 0L, ifelse(psi_i >= thr, 2L, 1L))
}
pct_not_assessed <- mean(colMeans(status_camera == 0))
cat("Average share of grid cells 'not assessed' across all", length(species), "species:",
    sprintf("%.1f%%\n", 100*pct_not_assessed))

fwrite(auc_results, "all_species_auc_camera.csv")
fwrite(as.data.table(status_camera, keep.rownames = "sid"), "grid_status_camera.csv")

cat("\n=== Module 5 script complete: camera-level psi + AUC/threshold/status for all",
    length(species), "species,", length(projectable_arr), "array-level projections computed. ===\n")
cat("See module5_prediction_mapping.qmd for the narrated analysis and all diagnostic figures.\n")
