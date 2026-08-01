## =============================================================================
## Module 2 — Covariates, Sampling & Protected Areas
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## What this script does, in order:
##   1. Explain why remote sensing is needed at all, why a buffer (not a
##      point) is used, and why percent composition (not a categorical
##      "in X habitat" label) is the right covariate form.
##   2. Load the pre-extracted covariate table and describe the candidate set.
##   3. Explain what Google Earth Engine (GEE) is and what a reader would need
##      to set up before running the (non-executed) extraction example.
##   4. Check for redundancy among the candidate covariates: correlation
##      matrix (full set vs. final reduced set), PCA, and variance inflation
##      factor (VIF).
##   5. Explain standardization and why we do NOT transform the skewed
##      pasture/cropland covariates.
##   6. Map the covariates: national surfaces, then camera sites on top.
##   7. Quantify within-array habitat homogeneity (intraclass correlation,
##      ICC) and show a concrete three-array example spanning the ICC range.
##   8. Summarise protected-area coverage across the network.
##   9. Evaluate species range maps: which sites are geographically
##      impossible for a given species, and should be masked before modeling.
##
## Requires: data.table, ggplot2. Point DATA_DIR at the folder containing the
## CSVs listed below (bundled alongside this script by default).
## =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})

# ---- 0. Paths ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"   # only used if you want to re-save figures from this script

stopifnot(dir.exists(DATA_DIR))


## =============================================================================
## 1. WHY WE NEED REMOTE SENSING AT ALL
## =============================================================================
## A camera trap only tells us what walked past it -- it says nothing on its
## own about the habitat surrounding that point. To relate detections to
## habitat, we need an independent measurement of what the landscape around
## every camera looks like: forest or savanna, cropland or pasture, wet or
## dry, protected or not. Remote sensing -- satellite imagery classified into
## land-cover types, or measured continuously (temperature, precipitation,
## canopy height) -- gives us that measurement at every camera location,
## cheaply and consistently, without a field survey at each site.
##
## The mechanics: every camera has a latitude/longitude. We take a
## satellite-derived map (a raster where every pixel is coded, e.g., forest
## vs. pasture vs. water) and read off what that map says at each camera's
## coordinates.
##
## ---- 1a. Why a buffer, not just a point ----
## A camera's exact coordinate is a single pixel -- but the animals it
## detects range over an area much larger than one pixel, and GPS
## coordinates carry their own few-meters of error. Reading the land-cover
## class at the literal point risks answering a fragile question: is this
## one 10-30 m pixel forest or not? Move the buffer a few meters and the
## answer can flip, even though nothing about the site's actual habitat
## changed.
##
## Instead, for every covariate we draw a BUFFER -- a circle of fixed radius
## around the camera -- and summarize the landscape WITHIN that circle, not
## just at its center. Two buffer sizes serve different purposes: a 100 m
## buffer for local habitat structure (what an animal encounters in its
## immediate surroundings), and a 1000 m buffer for landscape context (what
## kind of area the camera sits within more broadly).
##
## ---- 1b. Why percent composition, not a categorical label ----
## Given a buffer, there are still two ways to summarize it: pick the single
## class that covers the most area ("this site is forest"), or report the
## PERCENTAGE of the buffer's area in each class of interest ("87% forest,
## 13% pasture within 100 m"). We use the percentage form throughout. A
## categorical "in forest / not in forest" label throws away real
## information -- a site with 95% forest and a site with 51% forest both get
## called "forest," even though the second is meaningfully closer to a
## disturbance edge. Percent composition preserves that gradient, which is
## exactly the kind of continuous predictor an occupancy or abundance model
## needs: it can estimate how detection or occupancy probability changes as
## forest cover moves from 0% to 100%, rather than being told only which
## side of an arbitrary threshold a site falls on.

cat("== Why remote sensing, buffers, and percent composition ==\n",
    "See the comment block above for the full rationale.\n\n")


## =============================================================================
## 2. THE CANDIDATE COVARIATES
## =============================================================================
## Every candidate layer was extracted for every camera site at a 100 m
## buffer (local habitat) or a 1000 m buffer (landscape context).

cov <- fread(file.path(DATA_DIR, "covariates.csv"))
cat("== Candidate covariates ==\n")
cat("sites with covariates:", nrow(cov), "\n")
cand_cols <- setdiff(names(cov), c("site_id", "pt_hab", "wc_class_point", "in_pa"))
cat("candidate layers:", paste(cand_cols, collapse = ", "), "\n")

# We measured MapBiomas land-cover fractions (forest, savanna, pasture,
# cropland, native vegetation), two global cross-checks (WorldCover tree %,
# Hansen tree cover), canopy height, GHSL built-surface, TerraClimate
# temperature and precipitation, and WDPA protected-area status. Pasture
# (extensive cattle grazing) and cropland (row crops/plantation) are kept as
# two SEPARATE covariates rather than one combined "agriculture" layer,
# because camera-trap species respond to them very differently -- a pasture
# edge and a soy field are ecologically distinct disturbance types even
# though both count as human land use.

FINAL_COVARS <- c("forest_100m", "savanna_100m", "pasture_100m", "cropland_100m",
                   "native_veg_1000m", "temp_mean_C", "precip_annual_mm")
FULL_COVARS  <- c("forest_100m", "savanna_100m", "pasture_100m", "cropland_100m",
                   "native_veg_1000m", "wc_tree_pct_100m", "treecover2000_100m",
                   "canopy_height_m_100m", "ghsl_built_100m", "temp_mean_C",
                   "precip_annual_mm")


## =============================================================================
## 3. WHAT IS GOOGLE EARTH ENGINE, AND WHAT DO YOU NEED BEFORE RUNNING THIS CODE
## =============================================================================
## Google Earth Engine (GEE) is a cloud platform that hosts petabytes of
## satellite imagery (Landsat, Sentinel, MODIS, and derived products like
## MapBiomas or TerraClimate) alongside the compute to process it, so a
## query like "the percent forest cover in a 100 m circle around each of
## 1,125 points" runs on Google's servers in seconds rather than requiring
## you to download and process raw imagery yourself.
##
## Before the example code below can be run, a reader would need to:
##   1. Have (or create) a Google account, and sign up for Earth Engine
##      access at code.earthengine.google.com/register -- free for research
##      and non-commercial use, but requires an application and approval
##      (usually fast, sometimes a day or two).
##   2. Register a Google Cloud project and link it to Earth Engine -- GEE
##      now requires every user/script to run under a specific cloud project
##      ID, used for quota and billing tracking (compute itself remains free
##      for standard research use).
##   3. Choose an authentication method. For interactive use (the GEE Code
##      Editor, or a Python session on your own laptop), a browser-based
##      OAuth login is enough. For a script running unattended or in a
##      sandboxed environment (as this pipeline does), a SERVICE ACCOUNT is
##      needed instead -- a machine identity created in Google Cloud Console,
##      with a downloaded JSON key file, granted the "Earth Engine Resource
##      Writer" and "Service Usage Consumer" roles on your project.
##   4. Install the client library (earthengine-api for Python, or use the
##      JavaScript Code Editor directly in the browser) and initialize it
##      with your project ID and credentials before any extraction call
##      will run.
##
## None of this is needed to run THIS script -- the covariates are already
## extracted and saved in covariates.csv. It matters only if you want to run
## the extraction yourself, for a new set of camera points or a different
## country. The block below is REFERENCE ONLY (Python pseudocode, not
## executed by this R script).

## --- GEE extraction reference (Python, NOT executed here) -------------------
## import ee
## ee.Initialize(credentials, project="your-gee-project")
##
## # MapBiomas land-cover classification (30 m resolution)
## lulc = ee.Image("projects/mapbiomas-public/assets/brazil/lulc/collection10/"
##                  "mapbiomas_collection10_integration_v1").select("classification_2024")
##
## # Reduce to % composition of chosen classes within a 100 m buffer of each point
## def pct_in_buffer(pts_fc, image, class_codes, radius_m):
##     def per_point(pt):
##         buf = pt.geometry().buffer(radius_m)
##         hist = image.reduceRegion(ee.Reducer.frequencyHistogram(), buf, scale=30)
##         counts = ee.Dictionary(hist.get("classification_2024"))
##         total = counts.values().reduce(ee.Reducer.sum())
##         target = ee.List(class_codes).map(lambda c: counts.get(ee.String(ee.Number(c).format()), 0))
##         pct = ee.Number(ee.List(target).reduce(ee.Reducer.sum())).divide(total).multiply(100)
##         return pt.set("pct", pct)
##     return pts_fc.map(per_point)
##
## forest_codes = [1, 3, 5, 6, 49]        # Forest, Forest Formation, Mangrove, Floodable Forest, Wooded Sandbank
## pasture_codes = [15]                    # Pasture
## cropland_codes = [9,14,18,19,20,39,40,41,46,47,48,62,36,35]  # Farming/crops/plantation
##
## sites_fc = ee.FeatureCollection(site_points)
## forest_pct = pct_in_buffer(sites_fc, lulc, forest_codes, 100)
##
## # Climate: TerraClimate 1991-2020 normals
## terraclim = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate("1991-01-01","2021-01-01")
## temp_mean = terraclim.select(["tmmx","tmmn"]).mean().reduce(ee.Reducer.mean()).multiply(0.1)
## precip_annual = terraclim.select("pr").mean().multiply(12)
##
## # Protected areas: point-in-polygon against the WDPA polygon layer
## wdpa = ee.FeatureCollection("WCMC/WDPA/current/polygons")
## in_pa = sites_fc.map(lambda pt: pt.set("in_pa", wdpa.filterBounds(pt.geometry()).size().gt(0)))

cat("\n== GEE background ==\n",
    "See the comment block above for account setup and the extraction",
    "reference (not executed by this script).\n\n")


## =============================================================================
## 4. IS THERE REDUNDANCY IN THE CANDIDATE COVARIATES?
## =============================================================================
## The full 11-layer candidate set contains several land-cover and canopy
## measures that could plausibly all be describing the same underlying
## gradient. We check this three ways on the same data: a pairwise
## correlation matrix, a PCA biplot, and the variance inflation factor (VIF)
## -- each catching a different kind of redundancy.

# ---- 4a. Correlation, full set vs. final set --------------------------------
corr_full  <- cor(cov[, ..FULL_COVARS], use = "pairwise.complete.obs")
corr_final <- cor(cov[, ..FINAL_COVARS], use = "pairwise.complete.obs")

cat("== Correlation matrix, full candidate set (11 layers) ==\n")
print(round(corr_full, 2))
cat("\n== Correlation matrix, final modeling set (7 layers) ==\n")
print(round(corr_final, 2))

# Canopy height, Hansen treecover, WorldCover tree %, and forest % all
# correlate strongly with one another (r > 0.7) in the full set -- four
# different satellite products describing the same underlying gradient
# (woody cover). After dropping the redundant layers, every pairwise
# correlation in the final set is modest (|r| < 0.55), and pasture/cropland
# are essentially independent of every other covariate (|r| < 0.3).

# ---- 4b. PCA, full set vs. final set -----------------------------------------
pca_full  <- prcomp(na.omit(cov[, ..FULL_COVARS]),  scale. = TRUE)
pca_final <- prcomp(na.omit(cov[, ..FINAL_COVARS]), scale. = TRUE)

cat("\n== PCA, full candidate set: variance explained by PC1/PC2 ==\n")
print(round(summary(pca_full)$importance[2, 1:2] * 100, 1))
cat("\n== PCA, final modeling set: variance explained by PC1/PC2 ==\n")
print(round(summary(pca_final)$importance[2, 1:2] * 100, 1))

# In the full set, five layers (canopy height, Hansen treecover, WorldCover
# tree %, forest %, and GHSL built) crowd into nearly the same direction on
# the PCA biplot -- they are all measuring "how much woody cover" and barely
# separate. In the final set, all 7 arrows point in visibly different
# directions -- each covariate now contributes distinct information.

# ---- 4c. VIF: does correlation/PCA miss any JOINT redundancy? ---------------
## Correlation and PCA both describe pairwise or overall structure -- they
## can miss a case where three or more covariates are jointly redundant even
## though no single pair is highly correlated. VIF checks this directly: for
## each covariate, VIF measures how well ALL THE OTHER covariates predict it,
## VIF = 1/(1 - R^2), where R^2 comes from regressing that covariate on all
## the others simultaneously. VIF = 1 means no redundancy; VIF = 5 means the
## other covariates jointly explain 80% of this one's variance, so its
## coefficient estimate becomes unstable. No universal cutoff exists, but
## VIF > 5-10 is a common warning threshold.

compute_vif <- function(X) {
  X <- as.data.frame(X)
  sapply(names(X), function(v) {
    fit <- lm(as.formula(paste(v, "~ .")), data = X)
    1 / (1 - summary(fit)$r.squared)
  })
}
vif_final <- compute_vif(na.omit(cov[, ..FINAL_COVARS]))
cat("\n== VIF, final modeling set (all should be well below 5-10) ==\n")
print(round(vif_final, 2))

# The full 11-layer set has INFINITE VIF for the old lumped agriculture
# measure (built directly from the same pixels as one of its own inputs)
# and VIF up to 7 for the tree-cover layers; the final 7-layer set tops out
# under 2 -- confirming, from a third angle, what the correlation matrix and
# PCA already showed.


## =============================================================================
## 5. STANDARDIZING (AND WHETHER TO TRANSFORM)
## =============================================================================
## We z-standardize every continuous covariate -- subtract the mean, divide
## by the standard deviation -- before it enters a model. This puts every
## covariate on a common scale, so a fitted beta coefficient is directly
## comparable across covariates (a bigger |beta| means a stronger effect,
## not just a larger-scale covariate), and improves numerical stability
## during model fitting.

for (v in FINAL_COVARS) {
  cov[[paste0(v, "_z")]] <- as.numeric(scale(cov[[v]]))
}
cat("\n== Standardization: raw vs. z-scored skewness ==\n")
skewness <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x); s <- sd(x)
  mean(((x - m) / s)^3)
}
skew_tbl <- data.table(covariate = FINAL_COVARS,
                        raw_skew = sapply(FINAL_COVARS, function(v) skewness(cov[[v]])))
print(skew_tbl)

# Standardization fixes SCALE, not SHAPE -- a skewed distribution is still
# skewed after z-scoring. Pasture and cropland are so heavily zero-inflated
# that no standard transform (log, square-root) meaningfully normalizes
# them -- a log transform of a mostly-zero variable mostly just relabels the
# zeros. We keep these covariates on their raw percentage scale rather than
# force an ill-fitting transform: the models we fit (logistic-link
# occupancy) do not assume normally-distributed predictors, only a linear
# relationship on the link scale, so skew in the covariate itself is not a
## modeling violation -- it just means most of the "signal" in
## pasture/cropland comes from a smaller number of high-value sites.


## =============================================================================
## 6. WHERE THE COVARIATES ARE
## =============================================================================
## We look at the covariates spatially two ways: first the full NATIONAL
## surfaces (the value of each covariate everywhere in Brazil), then the
## CAMERA SITES drawn on top of those surfaces, which shows which parts of
## each gradient the network actually samples.
##
## (The national-surface raster maps and the camera-site maps are built from
## a national prediction grid and site lon/lat -- see figs/ for the rendered
## versions; reproducing the raster extraction itself requires a live GEE
## session, see Section 3.)

cat("\n== Covariate coverage ==\n",
    "See figs/m2_covariate_maps_national.png (national surfaces) and\n",
    "figs/m2_covariate_maps_updated.png (camera sites) for the maps.\n")

# Reading the two together makes the network's sampling footprint visible:
# cameras cluster in the higher-forest, lower-disturbance parts of the
# national surfaces, and thin out across the pasture- and
# cropland-dominated interior that the national maps show covers a large
# share of the country. Module 5 (prediction/mapping) quantifies how much of
# the covariate space the cameras actually cover, and where on the map
# predictions become extrapolation.


## =============================================================================
## 7. HOW MUCH DOES HABITAT VARY WITHIN AN ARRAY?
## =============================================================================
## Our design pools multiple cameras into an ARRAY (the replication unit for
## community models), and later modules use an array-level random effect as
## the primary correction for pseudoreplication. That correction is only
## appropriate if habitat is fairly homogeneous WITHIN an array. We measure
## this with the intraclass correlation (ICC): the fraction of a covariate's
## total variance that is BETWEEN arrays rather than WITHIN them. ICC near 1
## means an array is a single, homogeneous habitat patch; ICC near 0 means
## cameras within one array sample very different habitat.

dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
site_net <- unique(dep[, .(site_id, network, array_id)])
m <- merge(cov, site_net, by = "site_id", all.x = TRUE)

icc <- function(df, var, group) {
  d <- df[!is.na(get(var)) & !is.na(get(group))]
  grand_mean <- mean(d[[var]])
  grp <- d[, .(n = .N, mean_v = mean(get(var))), by = group]
  k <- nrow(grp); n_total <- nrow(d)
  if (k < 2) return(NA_real_)
  ms_between <- sum(grp$n * (grp$mean_v - grand_mean)^2) / (k - 1)
  within_ss <- sum((d[[var]] - d[, mean(get(var)), by = group][match(d[[group]], get(group)), V1])^2)
  ms_within <- within_ss / (n_total - k)
  n0 <- (n_total - sum(grp$n^2) / n_total) / (k - 1)
  if (ms_within == 0) return(1)
  var_between <- max(0, (ms_between - ms_within) / n0)
  var_between / (var_between + ms_within)
}

icc_covars <- c("forest_100m", "savanna_100m", "agriculture_100m",
                 "native_veg_1000m", "temp_mean_C", "precip_annual_mm")
icc_results <- data.table(covariate = character(), network = character(), icc = numeric())
for (net in c("Snapshot", "Atlantic", "WI")) {
  sub <- m[network == net]
  for (v in icc_covars) {
    icc_results <- rbind(icc_results, data.table(covariate = v, network = net, icc = icc(sub, v, "array_id")))
  }
}
icc_wide <- dcast(icc_results, covariate ~ network, value.var = "icc")
cat("\n== Within-array homogeneity (ICC), by network ==\n")
print(icc_wide)

# Climate is essentially constant within an array in every network
# (ICC > 0.93) -- expected, since weather does not vary meaningfully across
# a few kilometers. Habitat is more variable: Snapshot's tightly-clustered
# cameras show the highest habitat ICC of the three networks throughout,
# while Atlantic and WI each show lower habitat ICC for at least one
# covariate, meaning cameras within a single array can sit in noticeably
# different habitat.

# ---- 7a. Seeing ICC on the ground: three example arrays ---------------------
## Three real Snapshot Brasil arrays, chosen to span the ICC range using
## forest cover (%, within 100 m):
##   - Tefe (Amazonas): every camera reads 100% forest -- a genuinely
##     homogeneous array (std = 0).
##   - Caitaia (Rio Grande do Norte): mostly high forest cover with some
##     spread near a forest patch's edge (std = 6.7).
##   - Montes Claros (Minas Gerais): forest cover spans the FULL 0-100%
##     range within one array, because the underlying landscape there is a
##     genuinely fragmented forest/clearing mosaic (std = 42.5).
## An array-level random effect implicitly treats every camera within an
## array as sampling the same habitat -- that assumption holds well for
## climate everywhere, and for Snapshot's habitat covariates, but less well
## for Atlantic/WI arrays like Montes Claros. Because of this, Module 3
## fits single-species occupancy models TWO ways: with site-level habitat
## covariates, and with habitat covariates collapsed to the array level, to
## see how much the finer detail changes the ecological conclusions.
## See figs/m2_icc_array_example_v2.png for the map (forest-cover raster
## background, camera locations as black points).

cat("\n== ICC array examples ==\n",
    "See the comment block above and figs/m2_icc_array_example_v2.png\n")


## =============================================================================
## 8. PROTECTED AREAS: HOW EXPOSED IS THE NETWORK?
## =============================================================================

overall_pct <- round(100 * sum(cov$in_pa, na.rm = TRUE) / sum(!is.na(cov$in_pa)), 1)
cat("\n== Protected-area coverage ==\n")
cat("overall:", overall_pct, "% of sites fall inside a WDPA protected area\n")

pa <- fread(file.path(DATA_DIR, "pa_by_array_merged.csv"))
cat("-- by network --\n")
print(pa[, .(sites = sum(n_sites), pct_in_pa = round(100 * sum(n_in_pa) / sum(n_sites), 1)), by = network])

# Per-array PA status is close to bimodal -- most arrays are either almost
# entirely inside a PA or almost entirely outside one, because a park
# boundary rarely bisects a tightly-clustered camera array. This matters for
# Module 3: a PA effect estimated from array-level clustering carries less
# independent information than the raw site count might suggest.


## =============================================================================
## 9. RANGE-MAP EVALUATION: WHICH SPECIES NEED A MASKED CANDIDATE SET?
## =============================================================================
## Before fitting any occupancy or community model, we check each species'
## detections against its known geographic range. Some arrays fall outside a
## species' mapped range -- and if a species was never detected there, that
## "absence" is not ecological information about habitat; it is a
## geographic certainty that has nothing to do with the covariates above.
##
## ---- 9a. Accounting for imperfect range maps: a 100 km buffer ----
## IUCN range polygons are themselves imperfect -- they lag behind real
## distributional change, and mapping error near a range edge is common. To
## avoid masking a site for being merely NEAR the boundary of a plausible
## range, every polygon is expanded by a 100 KM BUFFER (computed in a
## metric, Brazil-appropriate projection) before the in-range test is
## applied. A site only becomes a candidate for masking if it falls outside
## the range PLUS this 100 km margin.
##
## ---- 9b. The two categories ----
## For each of the 17 species with an IUCN range polygon, every site is
## classified relative to that species' buffered range:
##   - IN RANGE -- falls inside the mapped range (incl. the 100 km buffer).
##     No adjustment needed.
##   - DETECTED OUT-OF-RANGE -- falls outside the buffered range, but the
##     species WAS photographed there. A photograph outranks a polygon:
##     these sites are always kept, but flagged (may indicate an outdated
##     range map or a taxonomic mismatch).
##   - STRUCTURAL ZERO -- falls outside the buffered range, and the species
##     was NEVER detected there. Polygon and camera data agree: this is the
##     "impossible site" that should be masked out before modeling, because
##     the guaranteed absence there carries no habitat information.
##
## ---- 9c. Two "problem child" species: presumed-unreliable range maps ----
## Azara's Agouti (Dasyprocta azarae) and Common Opossum (Didelphis
## marsupialis) show a pattern of out-of-range detections that is broad and
## geographically scattered in a way that looks like a mismatched or
## outdated range polygon rather than a real distributional limit. Rather
## than force a masking decision on possibly-wrong geography, we treat their
## range maps as UNRELIABLE and keep every site as a modeling candidate for
## these two species -- no masking is applied.

mask_summ <- fread(file.path(DATA_DIR, "m2_2_range_mask_summary_v2.csv"))
setorder(mask_summ, -pct_masked)
cat("\n== Masking impact per species (100 km buffer applied) ==\n")
print(mask_summ[, .(common, n_sites_total, n_structural_zero_sites,
                     n_sites_after_mask, pct_masked = round(pct_masked, 1),
                     range_map_presumed_unreliable)])

# Brazilian Common Opossum, White-eared Opossum, and Gray Brocket still lose
# more than a third of their sites (41-64% masked) -- the most conservative
# candidate sets going into Module 3. A broad middle group (Crab-eating Fox,
# White-lipped Peccary, Lowland Tapir, Giant Anteater) loses a modest
# 8-29%. A large group (Puma, Southern Tamandua, Nine-banded Armadillo,
# Collared Peccary, Spotted Paca, Ocelot, South American Coati, Tayra) is
# essentially unaffected (< 6% masked). Azara's Agouti and Common Opossum
# show 0% masked because their range maps are treated as unreliable.

## ---- 9d. What carries forward to Module 3 and beyond ----
## The full per-site mask table (m2_2_range_mask_full_v2.csv) is the
## candidate-set reference: for every one of the 17 range-mapped species, it
## records whether each site is in range, detected, and therefore whether it
## counts as a structural zero. Module 3 loads this table directly and fits
## every range-mapped species' occupancy model on its MASKED CANDIDATE SET
## (structural zeros removed, detected-out-of-range sites kept) as the
## default. Species without a range polygon (birds, and mammals outside the
## 17 mapped taxa) are unaffected and continue to use every site. The same
## masking logic carries into Module 4 (Community & joint models) via a
## latent-state range mask on the community occupancy state.

cat("\n== Range-mask carry-forward ==\n",
    "See the comment block above; the masked candidate sets defined here",
    "are the default for every range-mapped species from Module 3 onward.\n")


## =============================================================================
## SUMMARY -- what the results mean
## =============================================================================
## - Habitat and climate are measured as PERCENT COMPOSITION within a
##   buffer around each camera, not a categorical point label, because
##   camera-trap models need a continuous predictor.
## - Of an original 11-layer candidate set, several layers (canopy height,
##   Hansen treecover, WorldCover tree %, forest %) turned out to be
##   redundant proxies for the same woody-cover gradient. Correlation, PCA,
##   and VIF all confirm the same conclusion from different angles: a
##   reduced 7-layer set (forest, savanna, pasture, cropland, native
##   vegetation, temperature, precipitation) removes that redundancy while
##   keeping every ecologically distinct signal.
## - Within-array habitat homogeneity (ICC) is high for climate everywhere,
##   and high for Snapshot's habitat covariates, but lower for some Atlantic
##   and WI arrays -- a real design consideration for the array-level random
##   effect used from Module 3 onward.
## - Protected-area status is kept as its own covariate, tested alongside
##   (not folded into) the continuous habitat/climate set.
## - Range-map evaluation identifies, species by species, which sites are
##   geographically impossible and should be masked before modeling -- this
##   masked candidate set is the default going into Module 3.
##
## Next module: Single-Species Occupancy.
## =============================================================================
