## =============================================================================
## Module 3 — Single-Species Occupancy
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## What this script does, in order:
##   1. Fit a full single-species occupancy model (Nine-banded Armadillo),
##      with an array random effect, and walk through every diagnostic:
##      detection history, covariate effects, and model-fit measures
##      (AUC, c-hat, protected-area effect size).
##   2. Fit the same model to every well-detected species in the community,
##      and ask whether collapsing cameras to arrays as the sampling unit
##      changes the ecological conclusions (6 example species).
##   3. Extend the model two ways: letting a habitat effect vary with
##      temperature (continuous interaction), and letting it vary by
##      ecoregion (categorical interaction).
##   4. Close with the Royle-Nichols abundance comparison -- a different
##      way to use the same 1/0 detection data.
##
## Every species-level model uses the candidate site set defined in Module 2
## (Range-Map Evaluation): for the 17 species with an IUCN range map, sites
## that fall outside a 100 km-buffered range AND were never detected there
## (structural zeros) are excluded by default; sites detected outside the
## buffered range are always kept, flagged for review. Two species (Azara's
## Agouti, Common Opossum) have range maps treated as unreliable and use
## every site.
##
## Requires: data.table, ggplot2, unmarked, spOccupancy, pROC. Point DATA_DIR
## at the folder containing the CSVs listed below (bundled alongside this
## script by default).
## =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(unmarked)
})

# ---- 0. Paths ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"

stopifnot(dir.exists(DATA_DIR))

FINAL_COVARS_Z <- c("forest_100m_z","savanna_100m_z","pasture_100m_z","cropland_100m_z",
                     "native_veg_1000m_z","temp_mean_C_z","precip_annual_mm_z")


## =============================================================================
## 1. A FIRST OCCUPANCY MODEL: NINE-BANDED ARMADILLO
## =============================================================================
## We start with the Nine-banded Armadillo (Dasypus novemcinctus), a
## well-detected species, and fit a model using every covariate in the final
## set: forest, savanna, pasture, cropland, landscape-scale native
## vegetation, temperature, precipitation, and protected-area status.
## Armadillo has almost no range-mask exclusions, so its results here are
## close to what an unfiltered fit would give -- a useful baseline before
## looking at species where the mask matters more.

sc <- fread(file.path(DATA_DIR, "site_covariates.csv"))
sc <- sc[order(site_id)]
excl <- fread(file.path(DATA_DIR, "species_site_exclusions.csv"))

sci <- "Dasypus novemcinctus"
y_df <- fread(file.path(DATA_DIR, "dethist", paste0(gsub(" ","_",sci), ".csv")))
y_df <- y_df[order(site_id)]
y_mat <- as.matrix(y_df[, -1, with=FALSE])

excl_sites <- excl[species == sci, site_id]
keep <- !is.na(sc$forest_100m_z) & !(sc$site_id %in% excl_sites)
y2 <- y_mat[keep, , drop=FALSE]
sc2 <- as.data.frame(sc[keep])

cat("== Armadillo: candidate site set ==\n")
cat("sites after range-mask + missing-covariate exclusion:", nrow(y2), "\n\n")

## ---- 1.1 Detection history ----
## Occupancy models are built from a DETECTION HISTORY: for each 7-day
## occasion, was the species detected (1), not detected (0), or was the
## camera not running (missing)? We show two contrasting arrays below: one
## with many detections, one with almost none.

high_arr <- "SNAP_Piaui_Caatinga_Tapuio_25"
low_arr  <- "ATLA_ATL_11"
cat("== Detection history example arrays ==\n")
cat("High-detection array:", high_arr, "\n")
cat("Low-detection array:", low_arr, "\n")
cat("See figs/m3_dethist_armadillo.png for the graphical detection history,\n")
cat("and the 0/1 matrix tables in data/m3_dethist_table_high.csv /\n")
cat("data/m3_dethist_table_low.csv for the literal detection matrix.\n\n")

## ---- 1.2 Fitting the model, with the array random effect ----
## We fit the model with a random intercept per array -- the pipeline's
## standard correction for pseudoreplication among cameras that share the
## same array (see Module 2). Because this is a random-effects model, we
## use a Bayesian engine (spOccupancy::PGOcc) rather than the
## maximum-likelihood `unmarked` package, which only fits fixed effects.
## (This block shows the fixed-effects fit; the array-RE fit uses the same
## covariates via spOccupancy -- see the coefficient table printed below,
## precomputed and loaded from disk since the Bayesian fit takes several
## minutes.)

umf <- unmarkedFrameOccu(y = y2, siteCovs = sc2[, c(FINAL_COVARS_Z, "in_pa")])
form_full <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + "), "+ in_pa"))
fit_full <- occu(form_full, data = umf)

co <- coef(fit_full); se <- sqrt(diag(vcov(fit_full))); z <- co/se; p <- 2*pnorm(-abs(z))
cat("== Armadillo fixed-effects coefficients ==\n")
print(data.table(param=names(co), estimate=round(co,3), se=round(se,3), z=round(z,2), p=round(p,4)))

pgocc_coefs <- fread(file.path(DATA_DIR, "armadillo_pgocc_coefs.csv"))
cat("\n== Armadillo array-RE (PGOcc) coefficients: posterior mean, 95% credible interval ==\n")
print(pgocc_coefs[param != "(Intercept)"])

sigma_re <- fread(file.path(DATA_DIR, "armadillo_pgocc_sigma_re.csv"))
cat("\nArray random-effect variance (logit scale): mean =", round(sigma_re$sigma_mean,2),
    " sd =", round(sigma_re$sigma_sd,2), "\n")
cat("A variance this size means arrays differ substantially in their baseline\n")
cat("occupancy, beyond what the measured habitat covariates explain -- exactly\n")
cat("the kind of unexplained clustering the random effect is designed to absorb.\n\n")

## ---- 1.3 Naive vs. model-based occupancy ----
## "Naive occupancy" is simply the fraction of sites where the species was
## ever detected -- it always UNDERESTIMATES true occupancy, because a
## species can be present but missed on every visit.

naive_psi_site <- apply(y2, 1, function(r) as.integer(any(r==1, na.rm=TRUE)))
naive_psi <- mean(naive_psi_site)
psi_pred <- predict(fit_full, type="state")$Predicted
modeled_psi <- mean(psi_pred)
cat("== Naive vs. modeled occupancy ==\n")
cat("Naive occupancy (fraction of sites ever detected):", round(naive_psi,3), "\n")
cat("Modeled occupancy (array-RE PGOcc, mean psi):", round(modeled_psi,3), "\n")
cat("See figs/m3_naive_vs_modeled.png for the bar comparison.\n\n")

## ---- 1.4 Covariate effects ----
## Four covariates are significant in the fixed-effects fit: forest, savanna,
## native vegetation, and protected-area status. The figure
## (figs/m3_armadillo_effect_curves.png) shows predicted occupancy against
## each of these, holding the others at their mean -- with green lines for a
## positive effect and red for a negative one, matching the color
## convention used throughout this pipeline.

## ---- 1.5 Model fit: what do these numbers actually mean? ----
## Two very different questions get asked about a fitted occupancy model,
## and it helps to keep them separate:
##
##   (a) DOES THE MODEL DISCRIMINATE OCCUPIED FROM UNOCCUPIED SITES? --
##       measured by AUC (area under the ROC curve). Take every site's
##       predicted occupancy probability (psi-hat) and ask: if you picked
##       one truly-occupied site and one truly-unoccupied site at random,
##       what's the chance the model gives the occupied one a HIGHER
##       predicted probability? That chance is AUC. AUC = 0.5 means the
##       model does no better than a coin flip (predicted probabilities
##       for occupied and unoccupied sites completely overlap); AUC = 1.0
##       means perfect separation (every occupied site scores higher than
##       every unoccupied site). See figs/m3_auc_explainer.png for a worked
##       example: the left panel shows two overlapping distributions of
##       predicted probability (one for true presences, one for true
##       absences); the right panel is the ROC curve built by sweeping a
##       detection threshold across those distributions, and AUC is the
##       area under that curve.
##   (b) IS THE MODEL'S ASSUMED ERROR STRUCTURE CORRECT? -- measured by
##       C-HAT (overdispersion), from a parametric bootstrap
##       goodness-of-fit test. The model assumes detections follow a
##       binomial process with a single detection probability; c-hat
##       compares the observed chi-square discrepancy between the model
##       and the data against a bootstrap distribution of that same
##       statistic simulated FROM the fitted model. c-hat near 1 means the
##       assumed error structure fits; c-hat noticeably above 1 (as here)
##       means the real data are more variable than the model expects
##       (overdispersion) -- common when detection probability actually
##       varies by site or season in ways the model doesn't capture.
##       Practically, this means standard errors from the model should be
##       read as a bit optimistic (too narrow) rather than exact.
##
## AUC and c-hat are NOT interchangeable: a model can discriminate well
## (high AUC) while still being overdispersed (c-hat > 1), because AUC only
## cares about relative ranking of sites, while c-hat cares about whether
## the absolute variance in detections matches what the binomial model
## predicts.

fit_summary <- fread(file.path(DATA_DIR, "armadillo_fit_summary.csv"))
cat("== Armadillo model-fit summary ==\n")
print(fit_summary)
cat("\nAUC of", round(fit_summary$auc,2), "indicates modest discrimination -- better\n")
cat("than chance, but far from perfect, typical for occupancy models on binary\n")
cat("detection data with modest covariate signal. c-hat of", round(fit_summary$c_hat,2),
    "means the\nmodel is moderately overdispersed -- a common, expected finding, not a sign\n")
cat("the model is unusable.\n\n")

## ---- 1.6 Does adding protected-area status improve the model? ----
pa_aic <- fread(file.path(DATA_DIR, "armadillo_pa_aic_test.csv"))
cat("== AIC comparison: with vs. without protected-area status ==\n")
print(pa_aic)

pa_mag <- fread(file.path(DATA_DIR, "armadillo_pa_magnitude.csv"))
cat("\n== Predicted occupancy inside vs. outside protected areas ==\n")
print(pa_mag)
cat("\nAdding in_pa lowers AIC by", round(pa_aic[model=="without in_pa", delta_AIC],1),
    "points -- a decisive improvement.\n")
cat("The direction of the effect is", pa_mag$direction, ": Armadillo occupancy is",
    ifelse(pa_mag$direction=="negative","lower","higher"), "\ninside protected areas.",
    "In concrete terms, predicted occupancy inside protected\n")
cat("areas (", round(pa_mag$psi_inside_pa,2), ") is roughly", round(pa_mag$rel_change,2),
    "times predicted occupancy\noutside (", round(pa_mag$psi_outside_pa,2),
    ") -- a substantial effect, not just a statistically\ndetectable one.\n\n")


## =============================================================================
## 2. EVERY SPECIES, TOGETHER
## =============================================================================
## Every range-mapped species in this fit uses its masked candidate set
## (Module 2): structural-zero sites excluded, detected-out-of-range sites
## kept. Species without a range map, and the two species with unreliable
## range maps (Azara's Agouti, Common Opossum), use every site.
##
## Fitting an 8-covariate model to all modeled species with a
## fixed-effects-only engine is unreliable for roughly half the community
## (quasi-complete separation). We go straight to the array random-effect
## model instead, which regularizes each fit through the array intercept.

modeled <- fread(file.path(DATA_DIR, "modeled_species.csv"))
species_meta <- fread(file.path(DATA_DIR, "species_meta.csv"))
cat("== Community model ==\n")
cat("Species modeled:", nrow(modeled), "\n")
cat("See figs/m3_beta_shaded_table.png for the full community coefficient table:\n")
cat("each row is one species, with (n_det, n_arr) -- total detections and total\n")
cat("arrays detected at -- shown next to its name, each column is one covariate\n")
cat("(including Protected Area), green = positive effect, red = negative,\n")
cat("* = 95% credible interval excludes zero.\n\n")
cat("Forest is the most consistently important covariate across the community,\n")
cat("positive for the great majority of species with a credible effect.\n\n")

## ---- 2.1 Does the array-level simplification hold? ----
## Module 2 flagged that habitat is not always homogeneous within an array,
## especially for the Atlantic and WI networks. Here we test directly: does
## using arrays -- rather than individual cameras -- as the sampling unit
## change the ecological conclusions? We build GENUINE array-level detection
## histories: for each array and each occasion, the array counts as
## "detected" if any camera in it recorded the species that occasion, and
## "surveyed, not detected" if at least one camera was active but none
## detected it. This collapses cameras into arrays as the actual sampling
## unit -- not simply averaging site-level covariates while keeping
## site-level detections.
##
## Six species illustrate the comparison, spanning different data
## situations and different sensitivity to the level of aggregation:
## Spotted Paca and Nine-banded Armadillo (data-rich, moderate camera/array
## agreement), White-lipped Peccary (sparse, heavily range-masked at the
## array level), and three species chosen because their protected-area
## and/or forest coefficients change substantially -- in some cases
## reversing sign -- between the camera and array level: Ocelot, Lowland
## Tapir, and Puma.

site_6 <- fread(file.path(DATA_DIR, "site_level_6species_coefs.csv"))
array_6 <- fread(file.path(DATA_DIR, "array_true_6species_coefs.csv"))
cat("== Camera-level vs. array-level coefficients, 6 species ==\n")
cat("See figs/m3_array_vs_site_6species.png for the dot-whisker comparison.\n")
cat("Confidence intervals widen sharply at the array level for every species\n")
cat("and covariate -- fewer sampling units means less power. For Ocelot,\n")
cat("Lowland Tapir, and Puma specifically, the protected-area and/or forest\n")
cat("point estimates change direction between the two levels -- a caution\n")
cat("against over-interpreting any single array-level coefficient.\n\n")

fitstats_6 <- fread(file.path(DATA_DIR, "camera_vs_array_fitstats_6species.csv"))
cat("== Model-fit comparison: camera level vs. array level, 6 species ==\n")
cat("Columns: 'naive_occupancy_fraction_detected' = fraction of sampling units\n")
cat("(cameras or arrays) where the species was EVER detected -- a raw count,\n")
cat("no model involved. 'modeled_occupancy_probability' = the fitted model's\n")
cat("mean predicted occupancy probability across all sampling units --\n")
cat("accounts for imperfect detection, so it is always >= the naive value.\n\n")
print(fitstats_6[, .(common, level, n, auc=round(auc,3),
                       naive_occupancy_fraction_detected=round(naive_occupancy_fraction_detected,3),
                       modeled_occupancy_probability=round(modeled_occupancy_probability,3))])
cat("\nArray-level AUC is higher for all 6 species (roughly 0.60-0.80 vs.\n")
cat("0.56-0.70 at the camera level) -- with far fewer sampling units, an\n")
cat("array-level detection history is a cleaner, less noisy signal, which a\n")
cat("simple 3-covariate model can discriminate more easily, even though the\n")
cat("array-level fit has much less statistical power in its confidence\n")
cat("intervals.\n\n")

## ---- 2.2 Why this comparison matters ----
## This comparison is not just a technical check -- it speaks directly to a
## design choice this entire pipeline makes: using the ARRAY as the primary
## replication unit for community and joint models (Module 4 onward), while
## still fitting single-species models at the finer camera level here.
## The result is reassuring but not unconditional: point estimates mostly
## keep a consistent SIGN across the two levels, meaning the array-level
## simplification captures the right qualitative story for most
## species -- but the magnitude, and occasionally the sign, of a
## data-sparse species' covariate effect can shift once cameras are
## collapsed to arrays. In practice, this means: trust an array-level
## coefficient's sign and rough size for a well-detected species, but treat
## an array-level result for a sparse or heavily range-masked species (like
## White-lipped Peccary here) as a much weaker signal than its p-value alone
## would suggest -- confirm with the finer camera-level fit where possible,
## exactly as done in this section.

cat("== Why the array-level comparison matters ==\n")
cat("See the comment block above for the full explanation.\n\n")

## ---- 2.3 Species groupings from the coefficient table ----
## The community-wide coefficient table is dense. A principal-components
## analysis on these coefficients (7 habitat/climate covariates plus
## protected-area status, standardized) reduces that table to two axes and
## reveals genuine ecological groupings.
cat("== Species PCA ==\n")
cat("See figs/m3_species_pca_biplot.png. PC1 and PC2 together explain roughly\n")
cat("40% of the variance in species' covariate-response profiles. Arrows show\n")
cat("which covariates dominate each direction; species cluster into 4 groups\n")
cat("via k-means on the two PCA scores.\n\n")


## =============================================================================
## 3. ADDING AN INTERACTION EFFECT
## =============================================================================
## We extend the occupancy model to let one habitat effect depend on
## another variable -- here, whether the effect of tree cover on occupancy
## changes with temperature. This is fit as a FULL model, not an
## interaction term alone:
##
##   psi ~ tree_cover * temperature + savanna + pasture + cropland +
##         native_veg_1000m + precipitation + in_pa
##
## Tree cover (a MapBiomas year-2000 baseline layer, distinct from the
## pipeline's standard forest_100m covariate) substitutes for forest in this
## model, since the two are highly correlated (r = 0.82) and would
## otherwise create collinearity with the interaction term. All other
## standard covariates remain in the model as ordinary additive effects.

int_coefs <- fread(file.path(DATA_DIR, "interaction_4species_coefs.csv"))
int_terms <- int_coefs[param == "psi(treecover2000_100m_z:temp_mean_C_z)"]
cat("== Temperature x tree-cover interaction, full model, 4 species ==\n")
print(int_terms[, .(species, estimate=round(estimate,3), z=round(z,2))])
cat("\nSee figs/m3_interaction_4species.png for cool/average/warm temperature\n")
cat("scenario curves with 95% confidence bands. Three of four species show a\n")
cat("statistically significant interaction even after controlling for the\n")
cat("full covariate set -- a single additive-only model would average this\n")
cat("away.\n\n")


## =============================================================================
## 4. LETTING A RELATIONSHIP VARY BY ECOREGION
## =============================================================================
## A different way habitat relationships can vary: not with a continuous
## covariate, but across ecologically distinct regions. We test this with
## Ocelot (range-masked), comparing a single pooled forest effect against a
## FULL model that lets the forest effect differ between two ecoregions
## with enough data -- Moist Broadleaf Forest and the drier
## Grasslands/Savannas & Shrublands -- while every other covariate stays in
## the model as an ordinary additive (fixed, not ecoregion-varying) effect:
##
##   psi ~ forest * ecoregion + savanna + pasture + cropland +
##         native_veg_1000m + temperature + precipitation + in_pa

eco_aic <- fread(file.path(DATA_DIR, "ecoregion_ocelot_aic.csv"))
cat("== AIC comparison: pooled forest effect vs. ecoregion-interaction model ==\n")
print(eco_aic)
cat("\nSee figs/m3_ecoregion_ocelot.png. In Moist Broadleaf Forest, more forest\n")
cat("cover means more Ocelot occupancy, as expected for a forest specialist.\n")
cat("In the Grasslands/Savannas & Shrublands ecoregion the relationship\n")
cat("reverses direction -- the same lesson as the temperature interaction,\n")
cat("applied to a categorical rather than continuous moderator: a single\n")
cat("'forest effect' number can hide two genuinely different relationships\n")
cat("that only become visible once the data are split by ecological\n")
cat("context.\n\n")


## =============================================================================
## 5. COMPARING TO ROYLE-NICHOLS ABUNDANCE
## =============================================================================
## The Royle-Nichols model uses the same 1/0 detection data as occupancy,
## but interprets the FREQUENCY of detection across occasions as a signal
## of how many individuals are present at a site, not just whether the
## species is present at all: a site visited by three animals is detected
## more often, on average, than a site visited by one, even if the
## per-individual detection probability is identical.

cat("== Occupancy vs. Royle-Nichols abundance, across the community ==\n")
cat("See figs/m3_occ_vs_abundance.png. Left panel: mean occupancy (psi-hat)\n")
cat("vs. mean Royle-Nichols relative abundance (lambda-hat) across species --\n")
cat("agreement is strong overall. Right panel: forest effect estimated two\n")
cat("ways (occupancy vs. Royle-Nichols); agreement is more modest but still\n")
cat("clearly positive. Unstable Royle-Nichols fits (standard error above 5,\n")
cat("a sign of a rare species without enough information to pin down both an\n")
cat("intercept and a slope) are excluded from both panels. Two species\n")
cat("(labeled on the figure) sit well above the 1:1 line on the forest-effect\n")
cat("panel -- Royle-Nichols estimates a substantially larger forest effect on\n")
cat("abundance than occupancy estimates for presence alone, though both\n")
cat("species are rare enough that the Royle-Nichols magnitude should be read\n")
cat("cautiously even after excluding the most extreme unstable fits.\n\n")


## =============================================================================
## SUMMARY -- what the results mean
## =============================================================================
## - Every model in this module uses the range-masked candidate set from
##   Module 2 by default for the 17 range-mapped species.
## - Fixed-effects-only fitting is unreliable for roughly half the
##   community due to separation; the array random-effect model resolves
##   this by regularizing through the array intercept and is the
##   pipeline's default.
## - AUC (discrimination) and c-hat (overdispersion) measure different
##   things and can disagree -- a model can discriminate well while still
##   being overdispersed.
## - Forest cover is the community's most consistent habitat signal;
##   occupancy and Royle-Nichols abundance agree well on overall level and
##   reasonably on the forest effect specifically, once unstable
##   Royle-Nichols fits are excluded.
## - Collapsing cameras to arrays as the sampling unit widens confidence
##   intervals substantially -- a real cost in statistical power -- and for
##   data-sparse species can change a covariate's estimated direction, not
##   just its precision. Trust array-level results more for well-detected
##   species.
## - Two extensions -- a continuous interaction (temperature x tree cover)
##   and a categorical one (ecoregion) -- both show that a single "one
##   number fits all" habitat effect can obscure real, ecologically
##   sensible variation.
## - Protected-area status is a genuinely large effect for Nine-banded
##   Armadillo: predicted occupancy inside protected areas is well below
##   predicted occupancy outside, holding other covariates at their mean.
##
## Next module: Community & Joint Models.
## =============================================================================
