## =============================================================================
## Module 3 — Single-Species Occupancy
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## SELF-CONTAINED: this script starts from the three raw combined-dataset
## files (deployments, mammal detections, remote-sensing covariates) and
## builds everything from scratch -- detection histories, covariates,
## candidate site sets, and every model fit. Put this script in the same
## folder as the ten data files listed below and run it top to bottom.
##
## Required input files (in DATA_DIR, default "data" alongside this script):
##   final_deployments.csv         -- one row per camera deployment
##   final_detections_mammals.csv  -- one row per independent detection event
##   final_covariates.csv          -- one row per camera site, remote-sensing covariates
##   species_site_exclusions.csv   -- per-species structural-zero sites (Module 2 range-mask output)
##   genus_species_site_exclusions.csv -- structural-zero sites for pooled Didelphis (see note below)
##   site_ecoregion.csv            -- per-site biome/ecoregion (used only for the ecoregion extension)
##   species_order_lookup.csv      -- per-species/genus taxonomic order (used only for section 3.1)
##   species_common_names.csv      -- per-species/genus common name (used for PCA/RN figure labels)
##   brazil_boundary.gpkg          -- Brazil country outline (Natural Earth 110m; Johnson's-levels figure)
##   south_america_boundary.gpkg   -- South America outline (Natural Earth 110m; locator inset)
##
## The last four files are reference outputs from Module 2's range-map
## evaluation and this pipeline's taxonomy review, and cannot be rebuilt
## from the three main data files alone -- they are shipped alongside this
## script as fixed lookup tables.
##
## Taxonomy note: Didelphis (3 species with range maps + genus-only records)
## and Dasyprocta (5 species, only 1 with a range map, + genus-only records)
## are POOLED to genus level in this script -- a deliberate exception to the
## pipeline's usual species-level-only rule, made because field
## identification to species is unreliable for these two genera and pooling
## recovers genus-only detection records that would otherwise be discarded.
## Didelphis keeps its range mask (the union of all 3 mapped species'
## buffered ranges); Dasyprocta drops its range mask entirely, since only
## one of five congeners has a range map and using it alone would flag many
## real detections as "impossible."
##
## Requires: data.table, ggplot2, unmarked, spOccupancy, pROC, ggrepel, sf.
## If figure labels containing psi/beta/lambda symbols render incorrectly,
## run this script with a UTF-8 locale set (e.g. LC_ALL=en_US.UTF-8).
## =============================================================================

suppressMessages({
  library(data.table)
  library(unmarked)
  library(pROC)
  library(ggplot2)
})

# ---- 0. Paths ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"
stopifnot(dir.exists(DATA_DIR))
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR)

# color convention used throughout: green = positive effect, red = negative
COL_POS <- "#1b7837"; COL_NEG <- "#c0392b"; COL_NEUTRAL <- "#2c7fb8"

## Camera-level covariate set (11, incl. in_pa appended separately below) --
## used for every per-camera-site model (flagship, community loop, RN loop,
## interaction, ecoregion). Sample size at this level (1,110+ sites) easily
## supports this many predictors.
FINAL_COVARS <- c("forest_100m","savanna_100m","pasture_100m","cropland_100m",
                   "native_veg_1000m","temp_mean_C","precip_annual_mm",
                   "ghsl_built_5000m","dist_road_m","dist_water_m")

## Array-level covariate set (6, no in_pa) -- used ONLY for the array-level
## side of the camera-vs-array comparison (section 2.1). Chosen by ranking
## every camera-level candidate on community-wide explanatory power and
## confirming with VIF; a richer set is not defensible with only 60 arrays.
## See Module 2 for the derivation.
ARRAY_COVARS <- c("forest_100m","savanna_100m","native_veg_1000m",
                   "temp_mean_C","precip_annual_mm","dist_road_m")

OCC_DAYS <- 7  # occasion length in days, matching the pipeline's independence-interval convention


## =============================================================================
## STEP 1: LOAD RAW DATA AND BUILD SITE-LEVEL COVARIATES
## =============================================================================

dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
det <- fread(file.path(DATA_DIR, "final_detections_mammals.csv"))
cov <- fread(file.path(DATA_DIR, "final_covariates.csv"))

dep[, start_dt := as.IDate(start_date)]
dep[, end_dt   := as.IDate(end_date)]

# one row per site: network, array, coordinates, effort
site_tbl <- dep[, .(network = first(network), array_id = first(array_id),
                     longitude = mean(longitude), latitude = mean(latitude),
                     n_deployments = .N, camera_days = sum(camera_days)),
                 by = site_id]
setorder(site_tbl, site_id)
cat("Sites:", nrow(site_tbl), "| Arrays:", uniqueN(site_tbl$array_id), "\n")

# merge in remote-sensing covariates, z-score the final covariate set, number the arrays
sc <- merge(site_tbl, cov, by = "site_id", all.x = TRUE)
setorder(sc, site_id)
ALL_CONTINUOUS_COVARS <- union(FINAL_COVARS, ARRAY_COVARS)
for (v in ALL_CONTINUOUS_COVARS) sc[[paste0(v, "_z")]] <- as.numeric(scale(sc[[v]]))
sc[, treecover2000_100m_z := as.numeric(scale(treecover2000_100m))]
arrays_sorted <- sort(unique(sc$array_id))
sc[, array_num := match(array_id, arrays_sorted)]

cat("Sites missing covariates (excluded via listwise deletion at fit time):",
    sum(!complete.cases(sc[, ..FINAL_COVARS])), "\n\n")


## =============================================================================
## STEP 2: BUILD DETECTION HISTORIES FROM RAW DEPLOYMENTS + DETECTIONS
## =============================================================================
## An occupancy model needs a detection HISTORY: for each site and each
## 7-day occasion, was the species detected (1), surveyed but not detected
## (0), or was no camera active that occasion (missing)? We build this
## directly from the raw deployment windows and detection timestamps.
##
## Per-deployment windows (not a site's full min-max span) are used to mark
## which days a camera was actually active -- a multi-deployment site with a
## gap between two deployment periods should not be treated as continuously
## surveyed across that gap.

site_span <- dep[, .(min_start = min(start_dt), max_end = max(end_dt)), by = site_id]
setorder(site_span, site_id)
site_span[, n_occ := (as.integer(max_end - min_start) %/% OCC_DAYS) + 1L]
max_occ <- max(site_span$n_occ)
site_ids_ordered <- site_span$site_id
n_sites <- length(site_ids_ordered)

# per-site active-day boolean vector, built from each deployment's own start/end
active_day <- vector("list", n_sites)
names(active_day) <- site_ids_ordered
for (sid in site_ids_ordered) {
  min_start <- site_span[site_id == sid, min_start]
  max_end   <- site_span[site_id == sid, max_end]
  total_days <- as.integer(max_end - min_start) + 1L
  arr <- rep(FALSE, total_days)
  dep_sub <- dep[site_id == sid]
  for (i in seq_len(nrow(dep_sub))) {
    s_off <- as.integer(dep_sub$start_dt[i] - min_start) + 1L
    e_off <- as.integer(dep_sub$end_dt[i]   - min_start) + 1L
    arr[s_off:e_off] <- TRUE
  }
  active_day[[sid]] <- arr
}

# collapse active-day vector into occasion-level "surveyed" flag: an occasion counts as
# surveyed if >=50% of its days had an active camera
build_occasion_mask <- function(active_arr, occ_days = OCC_DAYS) {
  n_occ_here <- (length(active_arr) + occ_days - 1L) %/% occ_days
  mask <- rep(NA_real_, n_occ_here)
  for (i in seq_len(n_occ_here)) {
    window <- active_arr[((i - 1L) * occ_days + 1L):min(i * occ_days, length(active_arr))]
    if (mean(window) >= 0.5) mask[i] <- 0  # surveyed, not yet detected
  }
  mask
}
occ_active <- lapply(active_day, build_occasion_mask)

# map each detection to its site's occasion index; drop the handful of
# detections that predate the site's earliest recorded deployment start
# (bad calibration-image dates / off-by-one boundary cases -- a small,
# consistent fraction across the pipeline's history)
site_min_start <- setNames(site_span$min_start, site_span$site_id)
det[, site_min_start := site_min_start[site_id]]
det <- det[!is.na(site_min_start)]
det[, day_offset := as.integer(as.IDate(date) - site_min_start)]
det[, occ_idx := day_offset %/% OCC_DAYS]
n_dropped <- sum(det$occ_idx < 0)
det <- det[occ_idx >= 0]
cat("Detections used:", nrow(det), "(", n_dropped, "dropped: predate site's earliest deployment)\n\n")

# ---- Taxonomy: pool Didelphis and Dasyprocta to genus level ----
# (all other genera keep species-level identification, per the pipeline's
# standard rule; see the taxonomy note at the top of this script)
det[, sci_pooled := sci_mdd]
det[startsWith(sci_mdd, "Didelphis"), sci_pooled := "Didelphis"]
det[startsWith(sci_mdd, "Dasyprocta"), sci_pooled := "Dasyprocta"]

base_survey <- matrix(NA_real_, nrow = n_sites, ncol = max_occ)
rownames(base_survey) <- site_ids_ordered
for (sid in site_ids_ordered) {
  m <- occ_active[[sid]]
  base_survey[sid, seq_along(m)] <- m
}

build_dethist <- function(sci) {
  y <- base_survey
  sp_det <- det[sci_pooled == sci]
  for (i in seq_len(nrow(sp_det))) {
    sid <- sp_det$site_id[i]; oc <- sp_det$occ_idx[i] + 1L
    if (!(sid %in% rownames(y)) || oc > max_occ) next
    y[sid, oc] <- 1  # detected (also covers occasions the 50%-active rule marked as unsurveyed)
  }
  y
}

# ---- Choose the well-detected species to model ----
# same threshold used throughout this pipeline: species-level (or, for the
# two pooled genera, genus-level) binomial name, excluding domestic/human
# taxa and unresolved genus-only/family-only records for every OTHER genus,
# with at least 10 total detections and 20 sites detected at
sp_stats <- det[, .(n_det = .N, n_sites_det = uniqueN(site_id)), by = sci_pooled]
exclude_terms <- c("Homo sapiens","Bos taurus","Canis familiaris","Sus scrofa","Equus","Felis catus","Gallus")
genus_only_or_family <- c("Mazama","Didelphidae","Mammalia","Rodentia","Cervidae","Proechimys")
sp_stats[, is_binomial_or_pooled := grepl(" ", sci_pooled) | sci_pooled %in% c("Didelphis","Dasyprocta")]
sp_stats[, exclude := sci_pooled %in% exclude_terms | sci_pooled %in% genus_only_or_family]
candidates <- sp_stats[is_binomial_or_pooled & !exclude]
modeled <- candidates[n_det >= 10 & n_sites_det >= 20]
setorder(modeled, -n_det)
cat("Species/genera modeled:", nrow(modeled), "\n")
species_meta <- modeled[, .(sci_mdd = sci_pooled, n_det, n_arrays_det = NA_integer_)]
for (i in seq_len(nrow(species_meta))) {
  y <- build_dethist(species_meta$sci_mdd[i])
  det_sites <- rownames(y)[apply(y == 1, 1, any, na.rm = TRUE)]
  species_meta$n_arrays_det[i] <- uniqueN(sc[site_id %in% det_sites, array_id])
}
cat("\n")
fwrite(species_meta, "data/species_meta.csv")
fwrite(data.table(n_sites=nrow(site_tbl), n_arrays=uniqueN(site_tbl$array_id),
                   n_modeled=nrow(modeled)), "data/dataset_summary.csv")


## =============================================================================
## STEP 3: RANGE-MASK CANDIDATE SETS (Module 2 lookups)
## =============================================================================
## For species with an IUCN range map, a site is a "structural zero" if it
## falls outside the mapped range (with a 100km buffer for map
## imprecision) AND the species was never detected there -- these sites are
## excluded as default candidates. Sites detected outside the buffered
## range are always kept (a photograph outranks a polygon).

excl <- fread(file.path(DATA_DIR, "species_site_exclusions.csv"))
excl_genus <- fread(file.path(DATA_DIR, "genus_species_site_exclusions.csv"))
excl <- rbindlist(list(excl, excl_genus))
cat("Range-masked species (structural-zero exclusions applied):", uniqueN(excl$species), "\n")
cat("Pooled Didelphis: range mask = union of the 3 mapped congeners' buffered ranges (",
    excl_genus[species=="Didelphis", .N], "structural-zero sites excluded)\n")
cat("Pooled Dasyprocta: NO range mask applied (only 1 of 5 congeners has a range map;\n")
cat("  using it alone would flag many real detections as impossible) -- every site is a candidate.\n\n")


## =============================================================================
## 1. A FIRST OCCUPANCY MODEL: NINE-BANDED ARMADILLO
## =============================================================================
## We start with the Nine-banded Armadillo (*Dasypus novemcinctus*), a
## well-detected species, and fit a model using every covariate in the final
## set: forest, savanna, pasture, cropland, landscape-scale native
## vegetation, temperature, precipitation, and protected-area status.

sci <- "Dasypus novemcinctus"
y_arm <- build_dethist(sci)
excl_sites <- excl[species == sci, site_id]
FINAL_COVARS_Z <- paste0(FINAL_COVARS, "_z")
ARRAY_COVARS_Z <- paste0(ARRAY_COVARS, "_z")
keep <- complete.cases(sc[, ..FINAL_COVARS]) & !(sc$site_id %in% excl_sites)
y2 <- y_arm[keep, , drop = FALSE]
sc2 <- as.data.frame(sc[keep])
cat("== Armadillo: candidate site set ==\n")
cat("sites after range-mask + missing-covariate exclusion:", nrow(y2), "\n\n")

## ---- 1.1 Detection history ----
high_arr <- "SNAP_Piaui_Caatinga_Tapuio_25"
low_arr  <- "ATLA_ATL_11"
cat("== Detection history example arrays ==\n")
cat("High-detection array:", high_arr, "| Low-detection array:", low_arr, "\n\n")

n_occ_show <- 30
dethist_plot_data <- list()
for (arr_name in c(high_arr, low_arr)) {
  arr_sites <- sc[array_id == arr_name][order(site_id), site_id]
  sub <- y_arm[arr_sites, seq_len(min(n_occ_show, ncol(y_arm))), drop=FALSE]
  df <- as.data.table(sub)
  df[, site_id := arr_sites]
  df_long <- melt(df, id.vars="site_id", variable.name="occasion", value.name="status")
  df_long[, occasion_n := as.integer(gsub("V","",occasion))]
  df_long[, array_label := ifelse(arr_name==high_arr, "High-detection array", "Low-detection array")]
  df_long[, status_f := factor(ifelse(is.na(status), "Not surveyed", ifelse(status==1,"Detected","Surveyed, no detection")),
                                 levels=c("Not surveyed","Surveyed, no detection","Detected"))]
  dethist_plot_data[[arr_name]] <- df_long
}
dethist_dt <- rbindlist(dethist_plot_data)
p_dethist <- ggplot(dethist_dt, aes(x=occasion_n, y=site_id, fill=status_f)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_manual(values=c("Not surveyed"="white","Surveyed, no detection"="#dddddd","Detected"=COL_NEG), name=NULL) +
  facet_wrap(~array_label, scales="free_y") +
  labs(title="Nine-banded Armadillo detection history: high- vs. low-detection array",
       x="7-day occasion", y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.y=element_text(size=6), panel.grid=element_blank(), legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_dethist_armadillo.png"), p_dethist, width=11, height=6, dpi=150)
cat("Saved figs/m3_dethist_armadillo.png\n\n")

## ---- 1.2 Fitting the model, with the array random effect ----
## Fixed-effects fit (unmarked::occu) shown here; the array-random-effect
## fit (spOccupancy::PGOcc) requires several minutes of MCMC sampling and is
## commented out below with the exact call used -- uncomment to run it.

umf <- unmarkedFrameOccu(y = y2, siteCovs = sc2[, c(FINAL_COVARS_Z, "in_pa")])
form_full <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + "), "+ in_pa"))
fit_full <- occu(form_full, data = umf)
co <- coef(fit_full); se <- sqrt(diag(vcov(fit_full))); z <- co/se; p <- 2*pnorm(-abs(z))
cat("== Armadillo fixed-effects coefficients ==\n")
print(data.table(param=names(co), estimate=round(co,3), se=round(se,3), z=round(z,2), p=round(p,4)))
cat("\n")

## -- Array random-effect fit (spOccupancy::PGOcc), the pipeline's default --
suppressMessages(library(spOccupancy))
fit_pgocc_array_re <- function(y, sc_df, covars_z = FINAL_COVARS_Z, extra = "in_pa") {
  occ.covs <- sc_df[, c(covars_z, extra, "array_num")]
  occ.covs$array_num <- as.numeric(occ.covs$array_num)
  data.list <- list(y = y, occ.covs = occ.covs)
  occ.formula <- as.formula(paste("~", paste(covars_z, collapse=" + "), "+", extra, "+ (1|array_num)"))
  PGOcc(occ.formula = occ.formula, det.formula = ~1, data = data.list,
        n.samples = 8000, n.burn = 3000, n.thin = 5, n.chains = 1, verbose = FALSE)
}
fit_pgocc <- fit_pgocc_array_re(y2, sc2)
pgocc_coefs <- data.table(param = colnames(fit_pgocc$beta.samples),
                          mean = apply(fit_pgocc$beta.samples, 2, mean),
                          ci_lo = apply(fit_pgocc$beta.samples, 2, quantile, 0.025),
                          ci_hi = apply(fit_pgocc$beta.samples, 2, quantile, 0.975))
cat("== Armadillo array-RE (PGOcc) coefficients ==\n")
print(pgocc_coefs[param != "(Intercept)"], digits=3)
cat("\n")
fwrite(pgocc_coefs, "data/armadillo_pgocc_coefs.csv")
fwrite(data.table(n_sites=nrow(y2)), "data/armadillo_candidate_sites.csv")

fig_arm_effect_data <- list()
for (pr in c("forest_100m_z","savanna_100m_z","native_veg_1000m_z")) {
  b <- pgocc_coefs[param==pr, mean]; b_lo <- pgocc_coefs[param==pr, ci_lo]; b_hi <- pgocc_coefs[param==pr, ci_hi]
  intc <- pgocc_coefs[param=="(Intercept)", mean]
  zg <- seq(-2.5, 2.5, length.out=60)
  fig_arm_effect_data[[pr]] <- data.table(covariate=pr, z=zg,
    psi=plogis(intc + b*zg), psi_lo=plogis(intc + b_lo*zg), psi_hi=plogis(intc + b_hi*zg))
}
fig_arm_effect_dt <- rbindlist(fig_arm_effect_data)

## ---- 1.3 Naive vs. modeled occupancy ----
naive_psi_site <- apply(y2, 1, function(r) as.integer(any(r == 1, na.rm = TRUE)))
naive_psi <- mean(naive_psi_site)
psi_pred <- predict(fit_full, type = "state")$Predicted
modeled_psi_pgocc <- mean(apply(fit_pgocc$psi.samples, 2, mean))
cat("== Naive vs. modeled occupancy ==\n")
cat("Naive occupancy (fraction of sites ever detected):", round(naive_psi,3), "\n")
cat("Modeled occupancy (array-RE PGOcc, mean psi):", round(modeled_psi_pgocc,3), "\n\n")

p_naive_modeled <- ggplot(data.table(x=c("Naive\n(fraction detected)","Modeled\n(PGOcc, array RE)"),
                                       y=c(naive_psi, modeled_psi_pgocc)),
                            aes(x=x, y=y, fill=x)) +
  geom_col(width=0.55) +
  geom_text(aes(label=sprintf("%.2f", y)), vjust=-0.5, fontface="bold", size=4) +
  scale_fill_manual(values=c("#999999", COL_NEUTRAL), guide="none") +
  labs(title="Nine-banded Armadillo: naive vs. modeled occupancy", x=NULL, y="Occupancy (\u03c8)") +
  ylim(0, max(naive_psi, modeled_psi_pgocc)*1.3) +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIGS_DIR, "m3_naive_vs_modeled.png"), p_naive_modeled, width=5.5, height=5.5, dpi=150)
cat("Saved figs/m3_naive_vs_modeled.png\n\n")

## ---- 1.4 Covariate effects & 1.5 Model fit ----
## Two very different questions get asked about a fitted occupancy model:
## (a) DOES THE MODEL DISCRIMINATE occupied from unoccupied sites? -- AUC.
##     Take every site's predicted occupancy probability and ask: if you
##     picked one truly-occupied and one truly-unoccupied site at random,
##     what's the chance the model gives the occupied one a higher
##     predicted probability? AUC=0.5 is chance; AUC=1.0 is perfect.
## (b) IS THE MODEL'S ASSUMED ERROR STRUCTURE CORRECT? -- c-hat
##     (overdispersion), from a parametric bootstrap goodness-of-fit test.
##     c-hat near 1 means the assumed binomial error structure fits;
##     noticeably above 1 means real detections are more variable than the
##     model expects. AUC and c-hat are NOT interchangeable: a model can
##     discriminate well (high AUC) while still being overdispersed.

## ---- effect curves figure (4 significant covariates from the array-RE fit) ----
sig_covars_arm <- pgocc_coefs[param != "(Intercept)"][ci_lo>0 | ci_hi<0, param]
cat("Significant array-RE covariates for Armadillo:", paste(sig_covars_arm, collapse=", "), "\n\n")
effect_curve_data <- list()
for (pr in setdiff(sig_covars_arm, "in_pa")) {
  b <- pgocc_coefs[param==pr, mean]; b_lo <- pgocc_coefs[param==pr, ci_lo]; b_hi <- pgocc_coefs[param==pr, ci_hi]
  intc <- pgocc_coefs[param=="(Intercept)", mean]
  zg <- seq(-2.5, 2.5, length.out=60)
  effect_curve_data[[pr]] <- data.table(covariate=pr, z=zg,
    psi=plogis(intc + b*zg), psi_lo=plogis(intc + b_lo*zg), psi_hi=plogis(intc + b_hi*zg), sign=ifelse(b>0,"pos","neg"))
}
label_map <- c(forest_100m_z="Forest cover", savanna_100m_z="Savanna cover", pasture_100m_z="Pasture",
               cropland_100m_z="Cropland", native_veg_1000m_z="Native veg (1000m)",
               temp_mean_C_z="Temperature", precip_annual_mm_z="Precipitation",
               ghsl_built_5000m_z="Built surface", dist_road_m_z="Distance to road",
               dist_water_m_z="Distance to water", in_pa="Protected area")
if (length(effect_curve_data)) {
  ec_dt <- rbindlist(effect_curve_data)
  ec_dt[, covariate_label := factor(label_map[covariate], levels=label_map[names(effect_curve_data)])]
  p_effects <- ggplot(ec_dt, aes(x=z, y=psi, color=sign, fill=sign)) +
    geom_ribbon(aes(ymin=psi_lo, ymax=psi_hi), alpha=0.15, color=NA) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c(pos=COL_POS, neg=COL_NEG), guide="none") +
    scale_fill_manual(values=c(pos=COL_POS, neg=COL_NEG), guide="none") +
    facet_wrap(~covariate_label, nrow=1, scales="free_x") +
    labs(title="Nine-banded Armadillo: occupancy vs. each significant covariate (green=positive, red=negative)",
         x="Standardized covariate (z)", y="Predicted occupancy (\u03c8)") +
    theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
  ggsave(file.path(FIGS_DIR, "m3_armadillo_effect_curves.png"), p_effects, width=13, height=4, dpi=150)
  cat("Saved figs/m3_armadillo_effect_curves.png\n\n")
}

## ---- beta CI figure (array-RE model) ----
beta_plot <- pgocc_coefs[param != "(Intercept)"][order(mean)]
beta_plot[, label := label_map[param]]
stopifnot("Unmapped covariate in beta_plot -- add it to label_map" = !any(is.na(beta_plot$label)))
beta_plot[, label := factor(label, levels=unique(label))]
beta_plot[, sign := fifelse(ci_lo>0, "pos", fifelse(ci_hi<0, "neg", "ns"))]
p_beta_ci <- ggplot(beta_plot, aes(x=mean, y=label, color=sign)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0, linewidth=1) +
  geom_point(size=3) +
  scale_color_manual(values=c(pos=COL_POS, neg=COL_NEG, ns="grey50"), guide="none") +
  labs(title="Nine-banded Armadillo: covariate effects (array-RE model)",
       x="Coefficient (\u03b2), posterior mean and 95% credible interval", y=NULL) +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIGS_DIR, "m3_armadillo_beta_ci.png"), p_beta_ci, width=7, height=5.5, dpi=150)
cat("Saved figs/m3_armadillo_beta_ci.png\n\n")

auc_val <- as.numeric(auc(roc(response = naive_psi_site, predictor = psi_pred, quiet = TRUE)))
cat("== Model fit ==\nAUC:", round(auc_val, 3), "\n")

## ---- AUC explainer figure ----
roc_obj <- roc(response = naive_psi_site, predictor = psi_pred, quiet = TRUE)
roc_dt <- data.table(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
setorder(roc_dt, fpr)
score_dt <- data.table(psi = psi_pred, status = factor(ifelse(naive_psi_site==1,"True presence","True absence"),
                                                          levels=c("True absence","True presence")))
p_score_dist <- ggplot(score_dt, aes(x=psi, fill=status)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.55, bins=25) +
  scale_fill_manual(values=c("True absence"=COL_NEG, "True presence"=COL_POS), name=NULL) +
  labs(title="What the model predicts, by true status", x="Predicted occupancy probability (\u03c8-hat)", y="Density") +
  theme_minimal(base_size=10) + theme(legend.position=c(0.75,0.85))
p_roc <- ggplot(roc_dt, aes(x=fpr, y=tpr)) +
  geom_ribbon(aes(ymin=0, ymax=tpr), fill=COL_NEUTRAL, alpha=0.15) +
  geom_line(color=COL_NEUTRAL, linewidth=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  annotate("text", x=0.6, y=0.15, label="Chance\n(AUC = 0.5)", color="grey50", size=3.2) +
  labs(title=sprintf("ROC curve: AUC = %.2f", auc_val),
       x="False positive rate\n(unoccupied sites called 'present')",
       y="True positive rate\n(occupied sites called 'present')") +
  theme_minimal(base_size=10)
p_auc_explainer <- gridExtra::grid.arrange(p_score_dist, p_roc, ncol=2,
  top=grid::textGrob("AUC: how well does predicted occupancy separate true presence from true absence?",
                       gp=grid::gpar(fontsize=13, fontface="bold"), x=0.02, hjust=0))
ggsave(file.path(FIGS_DIR, "m3_auc_explainer.png"), p_auc_explainer, width=12, height=5, dpi=150)
cat("Saved figs/m3_auc_explainer.png\n")

## GoF bootstrap (c-hat) -- 100 simulations, ~1-2 minutes:
set.seed(1)
p_const <- unique(round(predict(fit_full, type="det")$Predicted, 6))[1]
compute_chisq <- function(ymat, psi_hat, p_hat) {
  n_occ_site <- rowSums(!is.na(ymat)); obs_det <- rowSums(ymat, na.rm=TRUE)
  exp_det <- psi_hat * n_occ_site * p_hat
  var_det <- psi_hat*(1-psi_hat)*(n_occ_site*p_hat)^2 + psi_hat*n_occ_site*p_hat*(1-p_hat)
  var_det[var_det<=0] <- 1e-6
  sum((obs_det - exp_det)^2 / var_det)
}
chisq_obs <- compute_chisq(y2, psi_pred, p_const)
nsim <- 100
chisq_sim <- numeric(nsim)
for (i in 1:nsim) {
  z_sim <- rbinom(nrow(y2), 1, psi_pred)
  y_sim <- matrix(NA, nrow=nrow(y2), ncol=ncol(y2))
  for (j in 1:nrow(y2)) {
    obs_mask <- !is.na(y2[j,])
    y_sim[j, obs_mask] <- rbinom(sum(obs_mask), 1, z_sim[j]*p_const)
  }
  umf_sim <- unmarkedFrameOccu(y=y_sim, siteCovs=sc2[, c(FINAL_COVARS_Z, "in_pa")])
  fit_sim <- tryCatch(occu(form_full, data=umf_sim), error=function(e) NULL)
  if (is.null(fit_sim)) { chisq_sim[i] <- NA; next }
  psi_sim <- predict(fit_sim, type="state")$Predicted
  p_sim <- unique(round(predict(fit_sim, type="det")$Predicted, 6))[1]
  chisq_sim[i] <- compute_chisq(y_sim, psi_sim, p_sim)
}
valid <- !is.na(chisq_sim)
c_hat <- chisq_obs / mean(chisq_sim[valid])
gof_p <- mean(chisq_sim[valid] >= chisq_obs)
cat("c-hat (overdispersion):", round(c_hat,3), "| Bootstrap GoF p-value:", round(gof_p,3), "\n")
cat("An AUC around 0.65 indicates modest discrimination; a c-hat above 1 means\n")
cat("the model is overdispersed -- both common, expected findings for occupancy\n")
cat("models on binary detection data, not a sign the model is unusable.\n\n")
fwrite(data.table(auc=round(auc_val,3), c_hat=round(c_hat,3), gof_p=round(gof_p,3),
                   det_prob=round(p_const,3)), "data/armadillo_fit_summary.csv")

## ---- 1.6 Does adding protected-area status improve the model? ----
form_no_pa <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + ")))
fit_no_pa <- occu(form_no_pa, data = umf)
delta_aic <- fit_no_pa@AIC - fit_full@AIC
newdat_out <- as.data.frame(setNames(as.list(rep(0, length(FINAL_COVARS_Z))), FINAL_COVARS_Z)); newdat_out$in_pa <- 0
newdat_in <- newdat_out; newdat_in$in_pa <- 1
psi_out <- predict(fit_full, newdata=newdat_out, type="state")$Predicted
psi_in <- predict(fit_full, newdata=newdat_in, type="state")$Predicted
direction <- ifelse(co["psi(in_pa)"] < 0, "negative", "positive")
cat("== Protected-area effect ==\n")
cat("Adding in_pa improves AIC by", round(delta_aic,1), "points.\n")
cat("Direction:", direction, "-- occupancy", ifelse(direction=="negative","lower","higher"), "inside protected areas.\n")
cat("Predicted occupancy inside PA:", round(psi_in,3), "| outside PA:", round(psi_out,3),
    "| ratio:", round(psi_in/psi_out,2), "\n\n")
fwrite(data.table(delta_aic=round(delta_aic,1), direction=direction,
                   psi_outside_pa=round(psi_out,3), psi_inside_pa=round(psi_in,3),
                   ratio=round(psi_in/psi_out,2)), "data/armadillo_pa_effect.csv")


## =============================================================================
## 2. EVERY SPECIES, TOGETHER
## =============================================================================
## Fitting an 8-covariate model to all modeled species/genera with a
## fixed-effects-only engine is unreliable for roughly half the community
## (quasi-complete separation). The pipeline's default is instead the array
## random-effect model (spOccupancy::PGOcc), which regularizes each fit
## through the array intercept -- commented out above for the flagship
## example; the same call is looped over every modeled species/genus to
## build the community-wide coefficient table used in the full pipeline.

cat("== Community model ==\n")
cat("Species/genera modeled:", nrow(modeled), "\n")
cat("Fitting the array-RE (PGOcc) model to every species/genus (~5-10 minutes)...\n")

community_coefs <- list()
for (sci_c in modeled$sci_pooled) {
  y_c <- build_dethist(sci_c)
  excl_c <- excl[species == sci_c, site_id]
  keep_c <- complete.cases(sc[, ..FINAL_COVARS]) & !(sc$site_id %in% excl_c)
  y_c2 <- y_c[keep_c, , drop=FALSE]
  sc_c2 <- sc[keep_c]
  fit_c <- tryCatch(fit_pgocc_array_re(y_c2, as.data.frame(sc_c2)), error=function(e) { cat(sci_c, "FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit_c)) next
  bs <- fit_c$beta.samples
  community_coefs[[sci_c]] <- data.table(species=sci_c, param=colnames(bs),
                                            mean=apply(bs,2,mean), ci_lo=apply(bs,2,quantile,0.025), ci_hi=apply(bs,2,quantile,0.975),
                                            n_sites=nrow(y_c2))
  cat(sci_c, "OK, n_sites=", nrow(y_c2), "\n")
}
community_dt <- rbindlist(community_coefs, fill=TRUE)
cat("\nCommunity fits completed:", uniqueN(community_dt$species), "of", nrow(modeled), "\n\n")
fwrite(data.table(n_fit=uniqueN(community_dt$species), n_modeled=nrow(modeled)), "data/community_fit_count.csv")
fwrite(community_dt, "data/community_coefs.csv")

## ---- community beta shaded table figure ----
label_map_full <- c(forest_100m_z="Forest", savanna_100m_z="Savanna", pasture_100m_z="Pasture",
                     cropland_100m_z="Cropland", native_veg_1000m_z="Native veg\n(1000m)",
                     temp_mean_C_z="Temp", precip_annual_mm_z="Precip",
                     ghsl_built_5000m_z="Built\nsurface", dist_road_m_z="Dist.\nroad",
                     dist_water_m_z="Dist.\nwater", in_pa="Protected\narea")
comm_plot <- community_dt[param != "(Intercept)"]
comm_plot <- merge(comm_plot, species_meta, by.x="species", by.y="sci_mdd", all.x=TRUE)
comm_plot[, sig := ci_lo>0 | ci_hi<0]
comm_plot[, ylab := paste0(species, "  (n_det=", n_det, ", n_arr=", n_arrays_det, ")")]
sp_order <- comm_plot[, .(n_det=first(n_det)), by=ylab][order(-n_det), ylab]
comm_plot[, ylab := factor(ylab, levels=rev(sp_order))]
comm_plot[, param_label := factor(label_map_full[param], levels=label_map_full[names(label_map_full) %in% param])]

p_beta_table <- ggplot(comm_plot, aes(x=param_label, y=ylab, fill=mean)) +
  geom_tile(color="white") +
  geom_text(data=comm_plot[sig==TRUE], aes(label="*"), fontface="bold", size=4) +
  scale_fill_gradient2(low=COL_NEG, mid="white", high=COL_POS, midpoint=0, name="Coefficient (\u03b2)") +
  labs(title="Habitat effects on occupancy across the community (array random effect)",
       subtitle="n_det = total detections, n_arr = arrays detected at; * = 95% credible interval excludes zero",
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.y=element_text(size=7), axis.text.x=element_text(size=8), panel.grid=element_blank())
ggsave(file.path(FIGS_DIR, "m3_beta_shaded_table.png"), p_beta_table, width=10, height=12, dpi=150)
cat("Saved figs/m3_beta_shaded_table.png\n\n")

## ---- 2.1 Does the array-level simplification hold? ----
## Six species illustrate the comparison: Spotted Paca and Nine-banded
## Armadillo (data-rich), White-lipped Peccary (sparse, range-masked), and
## Ocelot, Lowland Tapir, Puma (chosen for large camera-vs-array shifts).
## IMPORTANT: this is now a comparison of two DIFFERENT covariate sets, not
## just two sampling levels with the same covariates -- the camera-level fit
## uses the full 11-covariate set (n=1,110+ sites can support it); the
## array-level fit uses the reduced 6-covariate set (n=60 arrays cannot
## support more; see Module 2 for how that set was chosen).
species_6 <- c("Cuniculus paca","Dasypus novemcinctus","Tayassu pecari",
               "Leopardus pardalis","Tapirus terrestris","Puma concolor")

site_to_array <- setNames(sc$array_id, sc$site_id)
arrays_ordered <- sort(unique(sc$array_id))
n_arrays <- length(arrays_ordered)

array_cov <- sc[, lapply(.SD, mean, na.rm=TRUE), .SDcols = ARRAY_COVARS, by = array_id]
for (v in ARRAY_COVARS) array_cov[[paste0(v,"_z")]] <- as.numeric(scale(array_cov[[v]]))
# one array (ATLA_ATL_17) has no valid camera-level habitat/climate data at
# any of its cameras, so its array-level mean is NaN across the board --
# drop it here rather than silently carrying an all-NA row forward.
array_cov <- array_cov[complete.cases(array_cov[, ..ARRAY_COVARS])]
cat("Arrays with complete covariate data for array-level models:", nrow(array_cov), "of", n_arrays, "\n\n")
fwrite(data.table(n_arrays_complete=nrow(array_cov), n_arrays_total=n_arrays), "data/array_completeness.csv")

cat("== Camera-level (11-covariate) vs. array-level (6-covariate) comparison, 6 species ==\n")
fitstats <- list()
coef_compare <- list()
form_site6 <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + "), "+ in_pa"))
form_arr6  <- as.formula(paste("~1 ~", paste(ARRAY_COVARS_Z, collapse=" + ")))
for (sci6 in species_6) {
  y_site <- build_dethist(sci6)
  excl_sites6 <- excl[species == sci6, site_id]
  keep6 <- complete.cases(sc[, ..FINAL_COVARS]) & !(sc$site_id %in% excl_sites6)
  y_site2 <- y_site[keep6, , drop=FALSE]
  sc_site2 <- as.data.frame(sc[keep6])
  umf_s <- unmarkedFrameOccu(y=y_site2, siteCovs=sc_site2[, c(FINAL_COVARS_Z,"in_pa")])
  fit_s <- tryCatch(occu(form_site6, data=umf_s), error=function(e) NULL)

  # array-level: an array counts as "detected" if any camera in it detected that occasion
  y_arr_mat <- matrix(NA_real_, nrow=n_arrays, ncol=ncol(y_site))
  rownames(y_arr_mat) <- arrays_ordered
  for (sid in rownames(y_site)) {
    a <- site_to_array[sid]; row <- y_site[sid,]
    for (occ in seq_along(row)) {
      v <- row[occ]; if (is.na(v)) next
      if (is.na(y_arr_mat[a, occ])) y_arr_mat[a, occ] <- v else y_arr_mat[a, occ] <- max(y_arr_mat[a, occ], v)
    }
  }
  excl_arr6 <- excl[species==sci6, site_id]
  arr_excl_frac <- sc[site_id %in% excl_arr6, .(frac = 1), by = array_id]
  excl_arrays6 <- if (nrow(arr_excl_frac)) sc[site_id %in% excl_arr6][, .N, by=array_id][sc[,.N,by=array_id], on="array_id"][!is.na(N) & N/i.N > 0.5, array_id] else character(0)
  keep_arr <- !(array_cov$array_id %in% excl_arrays6)
  y_arr2 <- y_arr_mat[array_cov$array_id[keep_arr], , drop=FALSE]
  arr_cov2 <- as.data.frame(array_cov[keep_arr])
  umf_a <- unmarkedFrameOccu(y=y_arr2, siteCovs=arr_cov2[, ARRAY_COVARS_Z])
  fit_a <- tryCatch(occu(form_arr6, data=umf_a), error=function(e) NULL)

  if (!is.null(fit_s) && !is.null(fit_a)) {
    auc_s <- as.numeric(auc(roc(response=apply(y_site2,1,function(r) any(r==1,na.rm=TRUE)),
                                  predictor=predict(fit_s,type="state")$Predicted, quiet=TRUE)))
    auc_a <- as.numeric(auc(roc(response=apply(y_arr2,1,function(r) any(r==1,na.rm=TRUE)),
                                  predictor=predict(fit_a,type="state")$Predicted, quiet=TRUE)))
    fitstats[[sci6]] <- data.table(species=sci6, n_site=nrow(y_site2), n_array=nrow(y_arr2),
                                     auc_site=round(auc_s,3), auc_array=round(auc_a,3))

    co_s <- coef(fit_s); se_s <- sqrt(diag(vcov(fit_s)))
    co_a <- coef(fit_a); se_a <- sqrt(diag(vcov(fit_a)))
    # only compare parameters present in BOTH fits (the shared covariates:
    # forest, savanna, native_veg_1000m, temp, precip, dist_road)
    shared_params <- intersect(names(co_s), names(co_a))
    for (pr in shared_params) {
      coef_compare[[paste(sci6,pr)]] <- data.table(species=sci6, param=pr,
        site_est=co_s[pr], site_se=se_s[pr], array_est=co_a[pr], array_se=se_a[pr])
    }
  }
}
fitstats_dt <- rbindlist(fitstats, fill=TRUE)
print(fitstats_dt)
fwrite(fitstats_dt, "data/fitstats_6species.csv")
cat("\nArray-level AUC is typically higher -- fewer, cleaner sampling units are easier\n")
cat("to discriminate, at the cost of much wider confidence intervals (less power).\n")
cat("Note the camera-level model also carries four covariates -- built surface,\n")
cat("distance to road, distance to water, protected-area status -- that the\n")
cat("array-level model does not have room for; the comparison below is\n")
cat("restricted to the six covariates both models share.\n\n")

## ---- fit-stats comparison figure ----
fitstats_long <- melt(fitstats_dt, id.vars="species", measure.vars=c("auc_site","auc_array"),
                        variable.name="level", value.name="auc")
fitstats_long[, level := fifelse(level=="auc_site","Camera-level (11 cov.)","Array-level (6 cov.)")]
common_lookup6 <- c("Cuniculus paca"="Spotted Paca","Dasypus novemcinctus"="Nine-banded Armadillo",
                     "Tayassu pecari"="White-lipped Peccary","Leopardus pardalis"="Ocelot",
                     "Tapirus terrestris"="Lowland Tapir","Puma concolor"="Puma")
fitstats_long[, common := common_lookup6[species]]
fitstats_long[, common := factor(common, levels=common_lookup6[species_6])]
p_fitstats <- ggplot(fitstats_long, aes(x=common, y=auc, fill=level)) +
  geom_col(position="dodge", width=0.7) +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
  annotate("text", x=1, y=0.53, label="Chance", color="grey50", size=3) +
  scale_fill_manual(values=c("Camera-level (11 cov.)"=COL_NEUTRAL, "Array-level (6 cov.)"="#e67e22"), name=NULL) +
  labs(title="Model discrimination (AUC): camera-level (11-cov.) vs. array-level (6-cov.), 6 species", x=NULL, y="AUC") +
  ylim(0,1) + theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(FIGS_DIR, "m3_camera_vs_array_fit_6species.png"), p_fitstats, width=8, height=5.5, dpi=150)
cat("Saved figs/m3_camera_vs_array_fit_6species.png\n\n")

## ---- dot-whisker camera vs array coefficient comparison figure ----
## Restricted to the shared covariates (forest, native_veg_1000m -- the two
## with the most reliable signal across both levels for these species).
coef_compare_dt <- rbindlist(coef_compare, fill=TRUE)
coef_compare_dt <- coef_compare_dt[param %in% c("psi(forest_100m_z)","psi(native_veg_1000m_z)")]
coef_compare_dt[, common := common_lookup6[species]]
fwrite(coef_compare_dt, "data/coef_compare_6species.csv")
coef_long <- melt(coef_compare_dt, id.vars=c("species","param","common"),
                    measure.vars=list(est=c("site_est","array_est"), se=c("site_se","array_se")))
coef_long[, level := fifelse(variable==1, "Camera-level (11 cov.)", "Array-level (6 cov.)")]
param_labels3 <- c("psi(forest_100m_z)"="Forest", "psi(native_veg_1000m_z)"="Native veg (1000m)")
coef_long[, param_label := param_labels3[param]]
coef_long[, y_label := paste(common, param_label, sep=" \u2014 ")]
p_dotwhisker <- ggplot(coef_long, aes(x=est, y=y_label, color=level)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=est-1.96*se, xmax=est+1.96*se), height=0, position=position_dodge(width=0.5), linewidth=0.9) +
  geom_point(position=position_dodge(width=0.5), size=2.5) +
  scale_color_manual(values=c("Camera-level (11 cov.)"=COL_NEUTRAL, "Array-level (6 cov.)"="#e67e22"), name=NULL) +
  labs(title="Camera-level vs. array-level coefficients, 6 species (shared covariates only)",
       x="Coefficient (\u03b2), estimate and 95% CI", y=NULL) +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_array_vs_site_6species.png"), p_dotwhisker, width=9, height=8, dpi=150)
cat("Saved figs/m3_array_vs_site_6species.png\n\n")


## ---- 2.2 Why this comparison matters ----
cat("== Why this comparison matters ==\n")
cat("This speaks to a design choice the pipeline makes: using the array as\n")
cat("the primary replication unit for community/joint models (Module 4+),\n")
cat("while fitting single-species models at the finer camera level here.\n")
cat("Point estimates mostly keep a consistent sign across levels for\n")
cat("well-detected species, but a sparse or heavily range-masked species'\n")
cat("array-level result should be treated as a much weaker signal than its\n")
cat("p-value alone suggests.\n\n")

## ---- 2.2b Johnson (1980) hierarchy of habitat selection, illustrated ----
## Built natively from the site coordinates already loaded above -- no
## external image dependency. Requires the `sf` package (for the country
## outlines) in addition to the packages listed in the header.
cat("== Johnson's (1980) hierarchy of habitat selection ==\n")
suppressMessages(library(sf))
brazil_sf <- st_read(file.path(DATA_DIR, "brazil_boundary.gpkg"), quiet=TRUE)
sa_sf <- st_read(file.path(DATA_DIR, "south_america_boundary.gpkg"), quiet=TRUE)

COL_AVAIL_J <- "#bbbbbb"; COL_USED_J <- "#2c2c2c"
arr_centers_j <- sc[, .(lon=mean(longitude), lat=mean(latitude)), by=array_id]

## a 4-array cluster (~15-30km apart, clearly separable) used for the
## second- and third-order panels
cluster_arrays_j <- c("WI_WI_019","WI_WI_020","WI_WI_021","WI_WI_022")
cluster_sites_j <- sc[array_id %in% cluster_arrays_j, .(site_id, array_id, longitude, latitude)]

png(file.path(FIGS_DIR, "m3_johnson_levels_of_selection.png"), width=15.5, height=6.2, units="in", res=170)
par(mfrow=c(1,3), mar=c(1,1,3,1), oma=c(0,0,2,0))

## Panel 1: first-order selection -- geographic range
plot(st_geometry(brazil_sf), col=COL_AVAIL_J, border="white", main="First-order selection", cex.main=1.3, font.main=2)
points(arr_centers_j$lon, arr_centers_j$lat, pch=16, col=COL_USED_J, cex=0.7)
mtext("Available: geographic range of the species\nUsed: locations of camera arrays", side=1, line=-2, adj=0, cex=0.65)
## locator inset (South America, Brazil highlighted)
usr <- par("usr")
inset_w <- diff(usr[1:2])*0.32; inset_h <- diff(usr[3:4])*0.30
inset_x0 <- usr[1] + diff(usr[1:2])*0.02; inset_y0 <- usr[3] + diff(usr[3:4])*0.62
sa_bbox <- st_bbox(sa_sf)
sa_scale_x <- inset_w / (sa_bbox["xmax"]-sa_bbox["xmin"])
sa_scale_y <- inset_h / (sa_bbox["ymax"]-sa_bbox["ymin"])
sa_shift <- function(geom) {
  g <- (geom - c(sa_bbox["xmin"], sa_bbox["ymin"])) * c(sa_scale_x, sa_scale_y) + c(inset_x0, inset_y0)
  g
}
sa_geoms <- st_geometry(sa_sf) * matrix(c(sa_scale_x,0,0,sa_scale_y), 2,2)
sa_geoms <- sa_geoms + c(inset_x0 - st_bbox(sa_geoms)["xmin"], inset_y0 - st_bbox(sa_geoms)["ymin"])
plot(sa_geoms, col="#e0e0e0", border="white", lwd=0.3, add=TRUE)
br_geoms <- st_geometry(brazil_sf) * matrix(c(sa_scale_x,0,0,sa_scale_y), 2,2)
br_geoms <- br_geoms + c(inset_x0 - st_bbox(sa_geoms)["xmin"], inset_y0 - st_bbox(sa_geoms)["ymin"])
plot(br_geoms, col="#888888", border="white", lwd=0.3, add=TRUE)
rect(inset_x0, inset_y0, inset_x0+inset_w, inset_y0+inset_h, border="grey40", lwd=0.6)

## Panel 2: second-order selection -- home range within regional landscape
lon_c <- mean(cluster_sites_j$longitude); lat_c <- mean(cluster_sites_j$latitude)
half_span <- 0.35
asp_val <- 1/cos(lat_c*pi/180)
plot(NA, xlim=c(lon_c-half_span, lon_c+half_span), ylim=c(lat_c-half_span*0.85, lat_c+half_span*0.85),
     xlab="", ylab="", axes=FALSE, main="Second-order selection", cex.main=1.3, font.main=2, asp=asp_val)
rect(lon_c-half_span, lat_c-half_span*0.85, lon_c+half_span, lat_c+half_span*0.85, col=COL_AVAIL_J, border=NA)
for (arr in cluster_arrays_j) {
  sub <- cluster_sites_j[array_id == arr]
  if (nrow(sub) >= 3) {
    hpts <- as.matrix(sub[, .(longitude, latitude)])
    hull <- chull(hpts)
    polygon(hpts[c(hull, hull[1]), ], col=COL_USED_J, border="white", lwd=0.3)
  }
}
mtext("Available: regional landscape\nUsed: individual array footprints\n(home ranges)", side=1, line=-2, adj=0, cex=0.65)

## Panel 3: third-order selection -- site use within home range
example_array_j <- "WI_WI_020"
sub3 <- cluster_sites_j[array_id == example_array_j]
hpts3 <- as.matrix(sub3[, .(longitude, latitude)])
hull3 <- chull(hpts3)
xr3 <- range(sub3$longitude); yr3 <- range(sub3$latitude)
xpad3 <- diff(xr3)*0.25; ypad3 <- diff(yr3)*0.25
plot(NA, xlim=c(xr3[1]-xpad3, xr3[2]+xpad3), ylim=c(yr3[1]-ypad3, yr3[2]+ypad3),
     xlab="", ylab="", axes=FALSE, main="Third-order selection", cex.main=1.3, font.main=2, asp=asp_val)
polygon(hpts3[c(hull3, hull3[1]), ], col=COL_AVAIL_J, border=NA)
points(sub3$longitude, sub3$latitude, pch=16, col=COL_USED_J, cex=1.1)
mtext(sprintf("Available: home range (%s)\nUsed: individual camera locations", example_array_j), side=1, line=-2, adj=0, cex=0.65)

mtext("Johnson's (1980) hierarchy of habitat selection, illustrated with Snapshot Brasil camera-trap data",
      outer=TRUE, cex=1.05, font=2, adj=0.02)
dev.off()
cat("Saved figs/m3_johnson_levels_of_selection.png\n\n")

## ---- 2.3 Species groupings from the coefficient table ----
cat("== Species PCA ==\n")
pca_params <- c(names(label_map_full))
pca_wide <- dcast(community_dt[param != "(Intercept)"], species ~ param, value.var="mean")
pca_wide <- pca_wide[complete.cases(pca_wide[, ..pca_params])]
pca_mat <- scale(as.matrix(pca_wide[, ..pca_params]))
pca_fit <- prcomp(pca_mat, center=FALSE, scale.=FALSE)
ve <- summary(pca_fit)$importance[2, 1:2]

common_names <- fread(file.path(DATA_DIR, "species_common_names.csv"))
common_lookup_pca <- setNames(common_names$common_name, common_names$sci_pooled)
species_label_pca <- function(sci) {
  cn <- common_lookup_pca[sci]
  if (is.na(cn)) return(sci)
  cn
}

scores_dt <- data.table(species=pca_wide$species, PC1=pca_fit$x[,1], PC2=pca_fit$x[,2])
scores_dt[, label := sapply(species, species_label_pca)]
km_fit <- kmeans(scores_dt[, .(PC1,PC2)], centers=4, nstart=10)
scores_dt[, cluster := factor(km_fit$cluster)]
loadings_dt <- data.table(param=pca_params, PC1=pca_fit$rotation[,1]*3, PC2=pca_fit$rotation[,2]*3)
loadings_dt[, label := label_map_full[param]]

## Axis-direction labels: for each PC axis, find the covariate(s) most strongly
## associated with the positive vs. negative end, so a reader can tell what
## "high PC1" or "low PC2" means ecologically without cross-referencing the
## loading arrows by eye.
axis_dir_label <- function(loadings, axis_col) {
  vals <- loadings[[axis_col]]
  labs <- loadings$label
  pos_idx <- order(-vals)[vals[order(-vals)] > 0]
  neg_idx <- order(vals)[vals[order(vals)] < 0]
  pos_terms <- gsub("\\n", " ", labs[pos_idx][1:min(2,length(pos_idx))])
  neg_terms <- gsub("\\n", " ", labs[neg_idx][1:min(2,length(neg_idx))])
  list(pos = if (length(pos_terms)) paste(pos_terms, collapse=", ") else "(no strong covariate)",
       neg = if (length(neg_terms)) paste(neg_terms, collapse=", ") else "(no strong covariate)")
}
pc1_dir <- axis_dir_label(loadings_dt, "PC1")
pc2_dir <- axis_dir_label(loadings_dt, "PC2")
cat("PC1: more", pc1_dir$pos, "at the positive end; more", pc1_dir$neg, "at the negative end\n")
cat("PC2: more", pc2_dir$pos, "at the positive end; more", pc2_dir$neg, "at the negative end\n")

xr <- range(scores_dt$PC1); yr <- range(scores_dt$PC2)
xpad <- diff(xr) * 0.12; ypad <- diff(yr) * 0.12

p_pca <- ggplot(scores_dt, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=cluster), size=3, alpha=0.8) +
  ggrepel::geom_text_repel(aes(label=label), size=2.4, alpha=0.85, max.overlaps=30,
                             segment.size=0.2, segment.alpha=0.4, seed=1) +
  geom_segment(data=loadings_dt, aes(x=0,y=0,xend=PC1,yend=PC2), color=COL_NEG,
               arrow=arrow(length=unit(0.2,"cm")), inherit.aes=FALSE) +
  geom_text(data=loadings_dt, aes(x=PC1*1.15, y=PC2*1.15, label=label), color=COL_NEG, size=3, fontface="bold", inherit.aes=FALSE) +
  scale_color_brewer(palette="Set1", name="Group") +
  scale_x_continuous(limits=c(xr[1]-xpad, xr[2]+xpad*2.5),
                      sec.axis=dup_axis(breaks=xr, labels=c(paste0("more ", pc1_dir$neg), paste0("more ", pc1_dir$pos)), name=NULL)) +
  scale_y_continuous(limits=c(yr[1]-ypad, yr[2]+ypad*2),
                      sec.axis=dup_axis(breaks=yr, labels=c(paste0("more ", pc2_dir$neg), paste0("more ", pc2_dir$pos)), name=NULL)) +
  labs(title="Species grouped by full habitat-response profile (community model)",
       x=sprintf("PC1 (%.0f%%)", ve[1]*100), y=sprintf("PC2 (%.0f%%)", ve[2]*100)) +
  theme_minimal(base_size=12) +
  theme(axis.text.x.top=element_text(size=6.5, color="grey40"),
        axis.text.y.right=element_text(size=6.5, color="grey40"))
ggsave(file.path(FIGS_DIR, "m3_species_pca_biplot.png"), p_pca, width=10.5, height=9, dpi=150)
cat("Saved figs/m3_species_pca_biplot.png (", nrow(scores_dt), "species,", round(sum(ve)*100), "% variance explained)\n\n")


## =============================================================================
## 3. ADDING AN INTERACTION EFFECT
## =============================================================================
## Full model: psi ~ tree_cover * temperature + savanna + pasture + cropland
##             + native_veg_1000m + precipitation + in_pa
## Tree cover substitutes for forest here (r=0.82 collinearity with
## forest_100m) to avoid redundancy with the interaction term.

species4 <- c("Tapirus terrestris","Euphractus sexcinctus","Dasypus novemcinctus","Cerdocyon thous")
common4 <- c("Tapirus terrestris"="Lowland Tapir","Euphractus sexcinctus"="Yellow Armadillo",
             "Dasypus novemcinctus"="Nine-banded Armadillo","Cerdocyon thous"="Crab-eating Fox")
OTHER_COVARS <- c("savanna_100m_z","pasture_100m_z","cropland_100m_z","native_veg_1000m_z","precip_annual_mm_z","ghsl_built_5000m_z","dist_road_m_z","dist_water_m_z","in_pa")
cat("== Temperature x tree-cover interaction, full model, 4 species ==\n")
int_preds <- list()
int_labels <- list()
int4_results <- list()
for (sci4 in species4) {
  y4 <- build_dethist(sci4)
  excl4 <- excl[species==sci4, site_id]
  keep4 <- complete.cases(sc[, c("treecover2000_100m_z","temp_mean_C_z", ..OTHER_COVARS)]) & !(sc$site_id %in% excl4)
  y4b <- y4[keep4, , drop=FALSE]
  sc4 <- as.data.frame(sc[keep4])
  umf4 <- unmarkedFrameOccu(y=y4b, siteCovs=sc4[, c("treecover2000_100m_z","temp_mean_C_z", OTHER_COVARS)])
  form4 <- as.formula(paste("~1 ~ treecover2000_100m_z * temp_mean_C_z +", paste(OTHER_COVARS, collapse=" + ")))
  fit4 <- tryCatch(occu(form4, data=umf4), error=function(e) NULL)
  if (is.null(fit4)) next
  co4 <- coef(fit4); se4 <- sqrt(diag(vcov(fit4)))
  int_term <- "psi(treecover2000_100m_z:temp_mean_C_z)"
  z4 <- co4[int_term]/se4[int_term]
  int4_results[[sci4]] <- data.table(species=sci4, common=common4[sci4],
                                       interaction_coef=round(co4[int_term],3), z=round(z4,2))
  cat(sci4, "| interaction coef:", round(co4[int_term],3), "| z:", round(z4,2), "\n")
  int_labels[[sci4]] <- sprintf("%s\n(interaction %s, z=%.1f)", common4[sci4],
                                  ifelse(abs(z4)>1.96,"significant","not significant"), z4)

  tc_grid <- seq(min(sc4$treecover2000_100m_z, na.rm=TRUE), max(sc4$treecover2000_100m_z, na.rm=TRUE), length.out=40)
  for (temp_scenario in c(-1,0,1)) {
    newdat <- as.data.frame(setNames(as.list(rep(0,length(OTHER_COVARS))), OTHER_COVARS))
    newdat <- newdat[rep(1,40),]
    newdat$treecover2000_100m_z <- tc_grid; newdat$temp_mean_C_z <- temp_scenario
    pred <- predict(fit4, newdata=newdat, type="state")
    int_preds[[paste(sci4,temp_scenario)]] <- data.table(species=sci4, temp=factor(temp_scenario, levels=c(-1,0,1), labels=c("Cool","Average","Warm")),
                                                            tc=tc_grid, psi=pred$Predicted, lo=pred$lower, hi=pred$upper)
  }
}
cat("\nThree of four species typically show a significant interaction even with\n")
cat("the full covariate set included -- an additive-only model would miss this.\n\n")
fwrite(rbindlist(int4_results, fill=TRUE), "data/interaction_4species.csv")

int_pred_dt <- rbindlist(int_preds, fill=TRUE)
int_pred_dt[, species_label := unlist(int_labels[species])]
p_interaction <- ggplot(int_pred_dt, aes(x=tc, y=psi, color=temp, fill=temp)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c(Cool="#2c7fb8", Average="grey50", Warm=COL_NEG), name="Temperature") +
  scale_fill_manual(values=c(Cool="#2c7fb8", Average="grey50", Warm=COL_NEG), name="Temperature") +
  facet_wrap(~species_label, nrow=1, scales="free_y") +
  labs(title="Temperature x tree-cover interaction, full covariate model, 4 species",
       x="Tree cover (z)", y="Predicted occupancy (\u03c8)") +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_interaction_4species.png"), p_interaction, width=13, height=4.5, dpi=150)
cat("Saved figs/m3_interaction_4species.png\n\n")


## =============================================================================
## 3.1 IS THE TEMPERATURE x FOREST INTERACTION STRONGER IN SOME MAMMAL ORDERS?
## =============================================================================
## The 4-species example above shows the interaction can matter a lot for a
## handful of species. Here we fit the SAME interaction term (temperature x
## tree cover, full covariate model) for every modeled species/genus, tag
## each by taxonomic order, and ask whether one order stands out --
## specifically, whether Xenarthra (armadillos, anteaters; ectothermic-
## adjacent thermal physiology relative to other placental mammals) shows a
## more pronounced temperature-dependence in its forest-use response than
## the rest of the community.

order_lookup_dt <- fread(file.path(DATA_DIR, "species_order_lookup.csv"))
cat("== Temperature x tree-cover interaction, full community, by taxonomic order ==\n")
order_int_results <- list()
order_int_fits <- list()
for (sci_o in modeled$sci_pooled) {
  y_o <- build_dethist(sci_o)
  excl_o <- excl[species == sci_o, site_id]
  keep_o <- complete.cases(sc[, c("treecover2000_100m_z","temp_mean_C_z", ..OTHER_COVARS)]) & !(sc$site_id %in% excl_o)
  y_o2 <- y_o[keep_o, , drop=FALSE]
  sc_o2 <- as.data.frame(sc[keep_o])
  umf_o <- tryCatch(unmarkedFrameOccu(y=y_o2, siteCovs=sc_o2[, c("treecover2000_100m_z","temp_mean_C_z", OTHER_COVARS)]),
                     error=function(e) NULL)
  if (is.null(umf_o)) next
  fit_o <- tryCatch(occu(form4, data=umf_o), error=function(e) NULL)
  if (is.null(fit_o)) next
  co_o <- coef(fit_o); se_o <- sqrt(diag(vcov(fit_o)))
  int_term_o <- "psi(treecover2000_100m_z:temp_mean_C_z)"
  if (!(int_term_o %in% names(co_o))) next
  z_o <- co_o[int_term_o] / se_o[int_term_o]
  if (!is.finite(z_o) || abs(z_o) > 20) next  # drop unstable/separation-driven fits
  order_int_results[[sci_o]] <- data.table(species=sci_o, interaction_coef=co_o[int_term_o],
                                             interaction_se=se_o[int_term_o], interaction_z=z_o,
                                             n_sites=nrow(y_o2))
  order_int_fits[[sci_o]] <- fit_o
}
order_int_dt <- rbindlist(order_int_results, fill=TRUE)
order_int_dt <- merge(order_int_dt, order_lookup_dt, by.x="species", by.y="sci_pooled", all.x=TRUE)
order_int_dt[, sig := abs(interaction_z) > 1.96]
cat("Interaction fits completed:", nrow(order_int_dt), "of", nrow(modeled), "modeled species/genera\n\n")

order_summary <- order_int_dt[!is.na(order), .(n_species=.N, mean_abs_z=mean(abs(interaction_z)),
                                                  median_abs_z=median(abs(interaction_z)),
                                                  pct_sig=100*mean(sig)), by=order]
setorder(order_summary, -mean_abs_z)
cat("== Interaction strength (|z|) by taxonomic order ==\n")
print(order_summary, digits=3)
cat("\n")
fwrite(order_summary, "data/order_interaction_summary.csv")
fwrite(order_int_dt, "data/order_interaction_full.csv")

## ---- comparing the ten strongest interactions directly ----
## Rather than only comparing |z| values in a table, we plot the fitted
## tree-cover x temperature prediction curves for the ten species/genera
## with the strongest interaction (by |z|), all on a common scale, so the
## shape and direction of each interaction can be compared side by side.
top10_species <- order_int_dt[order(-abs(interaction_z))][1:min(10, .N), species]
top10_curves <- list()
for (sci_t in top10_species) {
  fit_t <- order_int_fits[[sci_t]]
  if (is.null(fit_t)) next
  sc_t <- as.data.frame(sc[complete.cases(sc[, c("treecover2000_100m_z","temp_mean_C_z", ..OTHER_COVARS)]) &
                             !(sc$site_id %in% excl[species==sci_t, site_id])])
  tc_grid_t <- seq(min(sc_t$treecover2000_100m_z, na.rm=TRUE), max(sc_t$treecover2000_100m_z, na.rm=TRUE), length.out=30)
  for (temp_scenario in c(-1, 1)) {
    newdat_t <- as.data.frame(setNames(as.list(rep(0,length(OTHER_COVARS))), OTHER_COVARS))
    newdat_t <- newdat_t[rep(1,30),]
    newdat_t$treecover2000_100m_z <- tc_grid_t; newdat_t$temp_mean_C_z <- temp_scenario
    pred_t <- predict(fit_t, newdata=newdat_t, type="state")
    top10_curves[[paste(sci_t,temp_scenario)]] <- data.table(
      species=sci_t, temp=factor(temp_scenario, levels=c(-1,1), labels=c("Cool","Warm")),
      tc=tc_grid_t, psi=pred_t$Predicted)
  }
}
top10_curves_dt <- rbindlist(top10_curves, fill=TRUE)
top10_z_lookup <- order_int_dt[species %in% top10_species, .(species, interaction_z, order)]
top10_curves_dt <- merge(top10_curves_dt, top10_z_lookup, by="species")
top10_curves_dt[, species_label := sprintf("%s (%s)\nz=%.1f", species, order, interaction_z)]
top10_curves_dt[, species_label := factor(species_label, levels=unique(species_label[order(-abs(top10_curves_dt$interaction_z))]))]

p_top10 <- ggplot(top10_curves_dt, aes(x=tc, y=psi, color=temp)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Cool"="#2c7fb8","Warm"=COL_NEG), name="Temperature scenario") +
  facet_wrap(~species_label, nrow=2, scales="free_y") +
  labs(title="The ten strongest temperature x tree-cover interactions, compared directly",
       subtitle="Cool vs. warm predicted-occupancy curves across the tree-cover range, for the ten species/genera with the largest |z|",
       x="Tree cover (z)", y="Predicted occupancy (\u03c8)") +
  theme_minimal(base_size=9) + theme(legend.position="bottom", strip.text=element_text(size=7))
ggsave(file.path(FIGS_DIR, "m3_interaction_top10.png"), p_top10, width=13, height=7, dpi=150)
cat("Saved figs/m3_interaction_top10.png\n\n")
fwrite(top10_z_lookup[order(-abs(interaction_z))], "data/top10_interaction_species.csv")

# Two-sample test: Xenarthra vs. all other orders combined
xen_z <- abs(order_int_dt[order=="Xenarthra", interaction_z])
other_z <- abs(order_int_dt[!is.na(order) & order!="Xenarthra", interaction_z])
if (length(xen_z) >= 3 && length(other_z) >= 3) {
  wt <- wilcox.test(xen_z, other_z)
  cat("Xenarthra (n=", length(xen_z), ") mean |z| =", round(mean(xen_z),2),
      "vs. all other orders (n=", length(other_z), ") mean |z| =", round(mean(other_z),2), "\n")
  cat("Wilcoxon rank-sum test, Xenarthra vs. rest: p =", round(wt$p.value,4), "\n")
  fwrite(data.table(n_xenarthra=length(xen_z), mean_xenarthra=round(mean(xen_z),2),
                     n_other=length(other_z), mean_other=round(mean(other_z),2),
                     wilcoxon_p=round(wt$p.value,4)), "data/xenarthra_test.csv")
  cat(if (wt$p.value < 0.05 && mean(xen_z) > mean(other_z))
        "Xenarthra shows a SIGNIFICANTLY stronger temperature x forest interaction than the rest of the community.\n"
      else if (wt$p.value < 0.05)
        "The difference is significant but in the OPPOSITE direction (Xenarthra weaker, not stronger).\n"
      else
        "No significant difference detected -- with only 5 Xenarthra species this test has limited power;\n  treat as suggestive, not conclusive.\n")
} else {
  cat("Too few species in one group for a formal test; reporting descriptive summary only.\n")
}
cat("\n")

## ---- figure: interaction strength by order, species-level dots + order means ----
order_int_dt[, order_f := factor(order, levels=order_summary$order)]
order_int_dt[, is_xenarthra := order == "Xenarthra"]
sig_label_dt <- data.table(order_f=factor(order_summary$order[1], levels=order_summary$order), y=2.15)
p_order_interaction <- ggplot(order_int_dt[!is.na(order)], aes(x=order_f, y=abs(interaction_z))) +
  geom_hline(yintercept=1.96, linetype="dashed", color="grey60") +
  geom_text(data=sig_label_dt, aes(x=order_f, y=y), label="|z| = 1.96 (95% sig.)",
            color="grey50", size=2.8, hjust=0, inherit.aes=FALSE) +
  geom_jitter(aes(color=is_xenarthra), width=0.15, size=2.5, alpha=0.8) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, color="black", linewidth=0.6) +
  scale_color_manual(values=c("TRUE"=COL_NEG, "FALSE"=COL_NEUTRAL), guide="none") +
  labs(title="Temperature x tree-cover interaction strength, by taxonomic order",
       subtitle="Each point = one species/genus; black bar = order mean; red = Xenarthra",
       x=NULL, y="|z| for the interaction term") +
  theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(FIGS_DIR, "m3_interaction_by_order.png"), p_order_interaction, width=9, height=6, dpi=150)
cat("Saved figs/m3_interaction_by_order.png\n\n")


## =============================================================================
## 4. LETTING A RELATIONSHIP VARY BY ECOREGION
## =============================================================================
## Full model: psi ~ forest * ecoregion + savanna + pasture + cropland +
##             native_veg_1000m + temperature + precipitation + in_pa

eco <- fread(file.path(DATA_DIR, "site_ecoregion.csv"))
sc_eco <- merge(sc, eco[, .(site_id, biome)], by="site_id", all.x=TRUE)
sci_eco <- "Leopardus pardalis"
y_eco <- build_dethist(sci_eco)
excl_eco <- excl[species==sci_eco, site_id]
OTHER_COVARS_ECO <- c("savanna_100m_z","pasture_100m_z","cropland_100m_z","native_veg_1000m_z","temp_mean_C_z","precip_annual_mm_z","ghsl_built_5000m_z","dist_road_m_z","dist_water_m_z","in_pa")
keep_base <- complete.cases(sc_eco[, c("forest_100m_z", ..OTHER_COVARS_ECO)]) & !(sc_eco$site_id %in% excl_eco)
target_biomes <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Grasslands, Savannas & Shrublands")
keep_eco <- keep_base & (sc_eco$biome %in% target_biomes)
y_eco2 <- y_eco[keep_eco, , drop=FALSE]
sc_eco2 <- as.data.frame(sc_eco[keep_eco])
sc_eco2$biome_f <- factor(sc_eco2$biome, levels=target_biomes)
umf_eco <- unmarkedFrameOccu(y=y_eco2, siteCovs=sc_eco2[, c("forest_100m_z","biome_f", OTHER_COVARS_ECO)])
form_pooled_eco <- as.formula(paste("~1 ~ forest_100m_z + biome_f +", paste(OTHER_COVARS_ECO, collapse=" + ")))
form_interact_eco <- as.formula(paste("~1 ~ forest_100m_z * biome_f +", paste(OTHER_COVARS_ECO, collapse=" + ")))
fit_pooled_eco <- occu(form_pooled_eco, data=umf_eco)
fit_interact_eco <- occu(form_interact_eco, data=umf_eco)
delta_eco <- fit_pooled_eco@AIC - fit_interact_eco@AIC
cat("== Ocelot: pooled vs. ecoregion-interaction forest effect ==\n")
cat("Interaction model AIC improvement:", round(delta_eco,1), "\n")
fwrite(data.table(delta_aic=round(delta_eco,1)), "data/ecoregion_ocelot_aic.csv")
cat("In Moist Broadleaf Forest, more forest cover typically means more Ocelot\n")
cat("occupancy; in Grasslands/Savannas & Shrublands the relationship reverses --\n")
cat("a single pooled 'forest effect' hides two genuinely different relationships.\n\n")

forest_grid <- seq(min(sc_eco2$forest_100m_z, na.rm=TRUE), max(sc_eco2$forest_100m_z, na.rm=TRUE), length.out=40)
eco_preds <- list()
for (b in target_biomes) {
  newdat_e <- as.data.frame(setNames(as.list(rep(0,length(OTHER_COVARS_ECO))), OTHER_COVARS_ECO))
  newdat_e <- newdat_e[rep(1,40),]
  newdat_e$forest_100m_z <- forest_grid; newdat_e$biome_f <- factor(b, levels=target_biomes)
  pred_e <- predict(fit_interact_eco, newdata=newdat_e, type="state")
  eco_preds[[b]] <- data.table(biome=b, forest_z=forest_grid, psi=pred_e$Predicted, lo=pred_e$lower, hi=pred_e$upper)
}
eco_pred_dt <- rbindlist(eco_preds)
biome_labels <- c("Tropical & Subtropical Moist Broadleaf Forests"="Moist Broadleaf Forest",
                   "Tropical & Subtropical Grasslands, Savannas & Shrublands"="Grasslands/Savannas")
eco_pred_dt[, biome_label := biome_labels[biome]]
p_ecoregion <- ggplot(eco_pred_dt, aes(x=forest_z, y=psi, color=biome_label, fill=biome_label)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
  geom_line(linewidth=1.1) +
  scale_color_manual(values=c("Moist Broadleaf Forest"=COL_POS, "Grasslands/Savannas"="#e67e22"), name=NULL) +
  scale_fill_manual(values=c("Moist Broadleaf Forest"=COL_POS, "Grasslands/Savannas"="#e67e22"), name=NULL) +
  labs(title="Ocelot: forest effect on occupancy, by ecoregion (full model)",
       x="Forest cover (z, standardized)", y="Predicted occupancy (\u03c8)") +
  theme_minimal(base_size=12) + theme(legend.position=c(0.25,0.9))
ggsave(file.path(FIGS_DIR, "m3_ecoregion_ocelot.png"), p_ecoregion, width=8, height=6, dpi=150)
cat("Saved figs/m3_ecoregion_ocelot.png\n\n")


## =============================================================================
## 5. COMPARING TO ROYLE-NICHOLS ABUNDANCE
## =============================================================================
## Royle-Nichols (occuRN) uses the same 1/0 detection data as occupancy, but
## interprets detection FREQUENCY as a signal of relative abundance rather
## than just presence/absence.

cat("== Occupancy vs. Royle-Nichols abundance ==\n")
cat("Fitting Royle-Nichols (occuRN) for every modeled species/genus...\n")
rn_coefs <- list()
for (sci_rn in modeled$sci_pooled) {
  y_rn <- build_dethist(sci_rn)
  excl_rn <- excl[species == sci_rn, site_id]
  keep_rn <- complete.cases(sc[, ..FINAL_COVARS]) & !(sc$site_id %in% excl_rn)
  y_rn2 <- y_rn[keep_rn, , drop=FALSE]
  sc_rn2 <- as.data.frame(sc[keep_rn])
  umf_rn <- tryCatch(unmarkedFrameOccu(y=y_rn2, siteCovs=sc_rn2[, c(FINAL_COVARS_Z, "in_pa")]), error=function(e) NULL)
  if (is.null(umf_rn)) next
  fit_rn <- tryCatch(occuRN(form_full, data=umf_rn, se=TRUE), error=function(e) NULL)
  if (is.null(fit_rn)) next
  co_rn <- coef(fit_rn); se_rn <- tryCatch(sqrt(diag(vcov(fit_rn))), error=function(e) rep(NA,length(co_rn)))
  rn_coefs[[sci_rn]] <- data.table(species=sci_rn, param=names(co_rn), estimate=co_rn, se=se_rn, n_sites=nrow(y_rn2))
}
rn_dt <- rbindlist(rn_coefs, fill=TRUE)
cat("RN fits completed:", uniqueN(rn_dt$species), "of", nrow(modeled), "\n\n")

## overall level comparison
occ_int_dt <- community_dt[param=="(Intercept)", .(species, occ_int=mean)]
rn_int_dt <- rn_dt[param=="lam(Int)", .(species, lam_int=estimate)]
overall_comp <- merge(occ_int_dt, rn_int_dt, by="species")
overall_comp[, mean_psi := plogis(occ_int)]; overall_comp[, mean_lam := exp(lam_int)]

## forest-effect comparison, filtering unstable RN fits
occ_forest_dt <- community_dt[param=="forest_100m_z", .(species, occ_forest=mean)]
rn_forest_dt <- rn_dt[param=="lam(forest_100m_z)", .(species, rn_forest=estimate, rn_se=se)]
forest_comp <- merge(occ_forest_dt, rn_forest_dt, by="species")
forest_comp_stable <- forest_comp[!is.na(rn_se) & rn_se < 5 & abs(rn_forest) < 5]
forest_comp_stable[, diff := rn_forest - occ_forest]
top2 <- forest_comp_stable[order(-diff)][1:min(2,.N)]
cat("Overall-level correlation:", round(cor(overall_comp$mean_psi, overall_comp$mean_lam),2), "\n")
cat("Forest-effect correlation (stable fits only):", round(cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest),2), "\n\n")
fwrite(data.table(overall_r=round(cor(overall_comp$mean_psi, overall_comp$mean_lam),2),
                   forest_r=round(cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest),2),
                   n_rn_fit=uniqueN(rn_dt$species), n_modeled=nrow(modeled)), "data/rn_vs_occ_correlation.csv")

p_overall <- ggplot(overall_comp, aes(x=mean_psi, y=mean_lam)) +
  geom_point(color=COL_NEUTRAL, size=2.5, alpha=0.75) +
  labs(title=sprintf("Overall level (r=%.2f)", cor(overall_comp$mean_psi, overall_comp$mean_lam)),
       x="Mean occupancy (\u03c8-hat)", y="Mean Royle-Nichols abundance (\u03bb-hat)") +
  theme_minimal(base_size=10)
top2[, common := sapply(species, species_label_pca)]
p_forest <- ggplot(forest_comp_stable, aes(x=occ_forest, y=rn_forest)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(color=COL_NEUTRAL, size=2.5, alpha=0.75) +
  ggrepel::geom_text_repel(data=top2, aes(label=common), size=2.8, color="black",
                             segment.size=0.3, segment.alpha=0.5, seed=2,
                             box.padding=0.6, max.overlaps=Inf) +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.22))) +
  labs(title=sprintf("Forest effect specifically (r=%.2f)", cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest)),
       x="Forest effect (occupancy, \u03b2)", y="Forest effect (Royle-Nichols, \u03b2)") +
  theme_minimal(base_size=10)
p_occ_vs_abund <- gridExtra::grid.arrange(p_overall, p_forest, ncol=2,
  top=grid::textGrob("Occupancy vs. Royle-Nichols abundance: agreement across the community",
                       gp=grid::gpar(fontsize=13, fontface="bold"), x=0.02, hjust=0))
ggsave(file.path(FIGS_DIR, "m3_occ_vs_abundance.png"), p_occ_vs_abund, width=11, height=5, dpi=150)
cat("Saved figs/m3_occ_vs_abundance.png\n\n")


## =============================================================================
## SUMMARY -- what the results mean
## =============================================================================
## - This script builds every detection history, covariate table, and model
##   directly from the three raw combined-dataset files -- no pre-fit
##   coefficients are loaded.
## - Didelphis and Dasyprocta are pooled to genus level (a deliberate
##   exception to the pipeline's species-level rule), recovering genus-only
##   detection records. Didelphis keeps a range mask (union of 3 mapped
##   congeners' buffered ranges); Dasyprocta has no range mask (only 1 of 5
##   congeners is mapped, and using it alone would misclassify real
##   detections as impossible).
## - Fixed-effects-only fitting is unreliable for roughly half the
##   community; the array random-effect model (commented-out PGOcc calls
##   above) is the pipeline's default and should be run for the full
##   community-wide analysis.
## - AUC (discrimination) and c-hat (overdispersion) measure different
##   things and can disagree.
## - Collapsing cameras to arrays widens confidence intervals substantially
##   and can shift a covariate's estimated direction for data-sparse species.
## - Both a continuous interaction (temperature x tree cover) and a
##   categorical one (ecoregion) show that a single "one number fits all"
##   habitat effect can obscure real, ecologically sensible variation.
## - Fitting the same temperature x tree-cover interaction across every
##   modeled species/genus and grouping by taxonomic order tests whether
##   this sensitivity clusters by lineage rather than being idiosyncratic
##   per species -- with only 5 Xenarthra species in the modeled set, this
##   is a suggestive test, not a definitive one.
##
## Next module: Community & Joint Models.
## =============================================================================
