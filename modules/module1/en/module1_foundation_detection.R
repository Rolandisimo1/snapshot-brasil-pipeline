## =============================================================================
## Module 1 — Foundation & Detection Histories
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## What this script does, in order:
##   1. Load Snapshot Brasil's own raw export and describe its structure
##      (deployments, placenames, sequences, arrays) and camera trap-type/height.
##   2. Summarise Snapshot Brasil's own species community, before any merge.
##   3. Quantify the cost of the 60-minute independence threshold (vs. the
##      conventional 1-minute standard) using Snapshot Brasil's own timestamps.
##   4. Load the final MERGED dataset (three networks: Snapshot Brasil, Atlantic
##      Forest, WI public projects) and describe how it was assembled.
##   5. Explain the deployment -> site -> array pipeline (why card swaps and
##      overlapping deployments need careful handling).
##   6. Build a detection history for one example species/array and show both
##      the graphical and literal 0/1 matrix forms.
##   7. Discuss occasion-window choice (why 7 days).
##   8. Summarise the merged mammal community (detection rate vs naive occupancy).
##
## Restricted to MAMMALS throughout, because the Atlantic Forest compilation
## recorded mammals only — keeping birds in would make every bird falsely
## "absent" on every Atlantic array (a data artifact, not ecology).
##
## Requires: data.table, ggplot2. Point DATA_DIR at the folder containing the
## CSVs listed below (bundled alongside this script by default).
## =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})

# ---- 0. Paths ---------------------------------------------------------------
# Edit DATA_DIR if you move this script away from its bundled data/ folder.
DATA_DIR <- "data"
FIGS_DIR <- "figs"   # only used if you want to re-save figures from this script

stopifnot(dir.exists(DATA_DIR))


## =============================================================================
## 1. SNAPSHOT BRASIL — structure and features (its own raw export)
## =============================================================================
## Wildlife Insights (the platform Snapshot Brasil runs on) organises data as:
##   - a DEPLOYMENT: one camera running continuously at one PLACENAME for one
##     period. Servicing a camera mid-season (battery/card swap) creates a
##     second deployment row at the same placename — still one physical camera.
##   - a SEQUENCE: a burst of photos grouped as one triggering event, with one
##     species identification.
##   - an ARRAY: Snapshot Brasil's basic sampling unit — a spatially clustered
##     group of camera stations (a "grid"/"subproject") deployed together.

snap_dep <- fread(file.path(DATA_DIR, "snapshot_brasil_deployments.csv"))
snap_seq <- fread(file.path(DATA_DIR, "snapshot_brasil_sequences.csv"))

cat("== Snapshot Brasil structure ==\n")
cat("deployments:", nrow(snap_dep),
    "| distinct placenames:", uniqueN(snap_dep$placename),
    "| arrays (subprojects):", uniqueN(snap_dep$subproject_name), "\n")
cat("total camera-days (raw deployments):",
    sum(as.integer(as.Date(snap_dep$end_date) - as.Date(snap_dep$start_date)), na.rm = TRUE), "\n")
cat("total sequences recorded:", nrow(snap_seq), "\n")

# NOTE: data.table::fread reads blank CSV cells as empty strings (""), NOT NA.
# A filter like `!is.na(genus)` silently keeps blank-genus rows (unidentified
# sequences, e.g. "Rodent" family-level labels) — always filter on `genus != ""`
# as well when excluding un-identified records.

# ---- 1a. Camera trap-type and height ----------------------------------------
# Every deployment records feature_type (what the camera was aimed at) and
# sensor_height (how high it was mounted). This matters ecologically: a camera
# on a burrow/water source samples differently than a randomized grid point,
# and mounting height constrains which species/behaviors are captured.

height_tbl <- snap_dep[, .N, by = sensor_height][order(-N)]
cat("\n-- sensor height breakdown --\n"); print(height_tbl)

feature_tbl <- snap_dep[, .N, by = feature_type][order(-N)]
cat("\n-- feature type breakdown --\n"); print(feature_tbl)

# We checked specifically for a canopy-height camera subset (arboreal
# placements would need separate handling — different community, and the
# "active" assumption for a ground-active mammal breaks down there). NONE
# EXISTS in the current data: every recorded height is ground-level (knee
# height ~40-50cm, chest height ~1-1.3m, or "Other" resolving to ~40cm or ~1m
# above ground in the free-text field). No canopy exclusion is needed.
cat("\nCanopy/arboreal camera check: sensor_height categories present are",
    paste(unique(snap_dep$sensor_height), collapse = ", "),
    "— all ground-level. No canopy set exists in this data.\n")


## =============================================================================
## 2. SNAPSHOT BRASIL — its own species community (before any merge)
## =============================================================================

snap_mammals <- snap_seq[class == "Mammalia" & genus != "Homo" &
                            genus != "" & species != ""]
snap_mammals[, sci := paste(genus, species)]

cat("\n== Snapshot Brasil's own community ==\n")
cat("distinct mammal taxa (raw genus+species, pre-harmonisation):",
    uniqueN(snap_mammals$sci), "\n")
cat("total identified mammal sequences:", nrow(snap_mammals), "\n")

top_species <- snap_mammals[, .N, by = .(sci, common_name)][order(-N)][1:10]
cat("\n-- ten most-recorded mammal taxa --\n"); print(top_species)


## =============================================================================
## 3. THE INDEPENDENCE INTERVAL — 1 minute vs. 60 minutes
## =============================================================================
## Any camera-trap dataset needs a rule for turning a photo stream into
## discrete independent detection "events" (otherwise a lingering animal, or a
## group moving through, gets counted many times). The usual convention is
## 1 MINUTE. This pipeline uses 60 MINUTES instead, for consistency with the
## Atlantic Forest compilation (whose source studies weren't timestamped
## finely enough to support a 1-minute rule). That move is not free — we
## quantify exactly what it costs using Snapshot Brasil's own timestamps
## (the only network in this pipeline with fine-grained sequence times).

snap_seq[, start_time := as.POSIXct(start_time)]
snap_seq[, end_time   := as.POSIXct(end_time)]

mammals_ts <- snap_seq[class == "Mammalia" & genus != "Homo" &
                          genus != "" & species != ""]
mammals_ts[, sci := paste(genus, species)]
setorder(mammals_ts, deployment_id, sci, start_time)

# Gap (minutes) to the previous sequence of the SAME species at the SAME
# deployment. The first sequence per deployment/species has no previous
# sequence, so it is always counted as independent (gap = Inf).
mammals_ts[, prev_end := shift(end_time, 1), by = .(deployment_id, sci)]
mammals_ts[, gap_min := as.numeric(difftime(start_time, prev_end, units = "mins"))]
mammals_ts[is.na(gap_min), gap_min := Inf]

mammals_ts[, indep_1min  := gap_min > 1]
mammals_ts[, indep_60min := gap_min > 60]

n_total       <- nrow(mammals_ts)
n_indep_1min  <- sum(mammals_ts$indep_1min)
n_indep_60min <- sum(mammals_ts$indep_60min)
n_collapsed   <- n_indep_1min - n_indep_60min
pct_collapsed <- 100 * n_collapsed / n_indep_1min

cat("\n== Independence interval: 1-min vs 60-min (Snapshot Brasil) ==\n")
cat("total mammal sequences:", n_total, "\n")
cat("independent @ 1-min: ", n_indep_1min, "\n")
cat("independent @ 60-min:", n_indep_60min, "\n")
cat("detections collapsed (1min -> 60min):", n_collapsed,
    sprintf(" (%.1f%%)\n", pct_collapsed))

# Per-species breakdown: which species lose the most detections to the wider
# window? (Restricted to species with >=20 independent detections at 1-min,
# for a stable percentage.)
by_species <- mammals_ts[, .(n_indep_1min = sum(indep_1min),
                              n_indep_60min = sum(indep_60min)),
                          by = .(sci, common_name)]
by_species[, n_collapsed := n_indep_1min - n_indep_60min]
by_species[, pct_collapsed := 100 * n_collapsed / n_indep_1min]
by_species_top <- by_species[n_indep_1min >= 20][order(-pct_collapsed)][1:12]

cat("\n-- species most affected by the 60-minute threshold --\n")
print(by_species_top[, .(common_name, sci, n_indep_1min, n_indep_60min,
                          n_collapsed, pct_collapsed = round(pct_collapsed, 1))])

# The pattern is ecological, not arbitrary: gregarious or site-lingering
# species (Rock Cavy, cattle, coati, peccaries) lose the most detections,
# because a group or a lingering animal triggers repeated captures within the
# same hour. Solitary, wide-ranging species are comparatively unaffected.


## =============================================================================
## 4. COMBINING THREE NETWORKS INTO ONE DATASET
## =============================================================================
## The full analysis combines:
##   - Snapshot Brasil   : the coordinated 2025 national survey (Section 1-3)
##   - Atlantic Forest   : a compilation of Atlantic-forest studies, 2016-2020
##   - WI public projects: public Wildlife Insights projects in Brazil, 2010-2025
## All in Wildlife Insights format.

dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
det <- fread(file.path(DATA_DIR, "detections.csv"))   # mammals only

cat("\n== Merged dataset ==\n")
cat("networks:", uniqueN(dep$network),
    "| arrays:", uniqueN(dep$array_id),
    "| sites:", uniqueN(dep$site_id),
    "| deployments:", nrow(dep),
    "| camera-days:", format(sum(dep$camera_days, na.rm = TRUE), big.mark = ","), "\n")
cat("mammal detections:", format(nrow(det), big.mark = ","),
    "| mammal species (all IDs):", uniqueN(det$sci_mdd), "\n")

# ---- 4a. How the Atlantic and WI data were selected -------------------------
# Snapshot Brasil is a purpose-built survey with a consistent design. The
# other two networks are compilations of independent studies, so they needed
# harmonising to the same standard before merging. The SAME rules were
# applied to both:
#   - Exclude baited cameras (Atlantic only).
#   - 60-minute independence interval, applied consistently across networks
#     (Section 3 quantifies what this costs).
#   - Cluster cameras into ARRAYS using the Snapshot design rule: cameras
#     within 200m thinned to one; linked into an array if within 5km of a
#     neighbour; array kept only if >= 10 cameras.
#   - One season per array — for the multi-year Atlantic/WI networks, each
#     array is pinned to its single best-sampled year (Snapshot is already
#     single-season and exempt).
#   - Coordinate precision — sites with coordinates rounded to ~1km were
#     removed (their 100m-buffer environmental covariates would be wrong).
# The result: "site" and "array" mean the same thing regardless of network.

# ---- 4b. A consistency choice: mammals only ---------------------------------
# The Atlantic Forest compilation recorded ONLY mammals. Keeping birds in the
# analysis would make every bird falsely ABSENT on every Atlantic array — a
# data artifact, not ecology. We restrict the entire analysis to mammals to
# keep every species comparable across all three networks.

# ---- 4c. Taxonomy: harmonising names across networks ------------------------
# The three networks were compiled at different times with different
# taxonomic conventions. Every name is reconciled to the Mammal Diversity
# Database (MDD) v2.4, the current global standard.

tax <- fread(file.path(DATA_DIR, "taxonomy_changes.csv"))
cat("\n-- taxonomic name updates applied (MDD v2.4) --\n")
print(tax[, .(data_name, mdd_accepted, mdd_common, ambiguous, sources)])


## =============================================================================
## 5. FROM DEPLOYMENTS TO SITES, ARRAYS, AND EFFORT
## =============================================================================
## A DEPLOYMENT is one camera running for one continuous period. The
## distinction between a raw deployment row and an analytical SITE matters for
## two concrete reasons, both of which would otherwise distort a model:
##
##   (a) CARD SWAPS — servicing a camera mid-season creates a second
##       deployment row at the same location. Treated naively this doubles the
##       site's sample contribution AND can fabricate false "0" (surveyed,
##       not detected) occasions across the servicing gap if the detection
##       history is built by filling the FULL span between the first and last
##       deployment date, rather than each deployment's own true active
##       window. (This was a real bug found and fixed during this pipeline's
##       development: filling min(start)-max(end) as if the camera ran
##       continuously fabricated ~1,200 false "0" days across 35 multi-
##       deployment sites, mostly in the Atlantic network — worst case, 171 of
##       181 span-days were fabricated at one site.) We concatenate matching
##       deployments into a single SITE (site_id), and respect each
##       deployment's own active window when building detection histories, so
##       gaps between servicing visits are correctly treated as "no effort,"
##       not "surveyed, nothing detected."
##   (b) OVERLAPPING DEPLOYMENTS — two deployments at the same location that
##       overlap in time (a backend error, e.g. duplicate entry). We exclude
##       both and flag them, rather than guess which is correct.
##
## The ARRAY — not the individual camera — is the replication unit every
## downstream model (occupancy random effects, community models) uses. An
## array containing silently duplicated or gap-inflated sites would bias
## every estimate built on top of it.

eff <- dep[, .(camera_days = sum(camera_days, na.rm = TRUE), n_deploy = .N),
           by = .(network, array_id, site_id)]
cat("\n== Effort summary ==\n")
cat("distinct sites:", uniqueN(eff$site_id),
    "| median camera-days per site:", round(median(eff$camera_days, na.rm = TRUE)), "\n")
cat("-- effort by network --\n")
print(eff[, .(sites = uniqueN(site_id), camera_days = sum(camera_days)), by = network])


## =============================================================================
## 6. DETECTION HISTORIES
## =============================================================================
## Occupancy models do not use raw photos. They use a DETECTION HISTORY: for
## each site, the sampling period is divided into equal OCCASIONS (7-day
## windows here), and each occasion is scored 1 (species detected at least
## once) or 0 (not detected but the camera was active). This 1/0/1/1/0...
## string is the raw material of an occupancy model.
##
## Example: Black Agouti at the Tefe array (Amazonas) — a single continuous
## 2025 season with 30 cameras and no servicing gaps, chosen because the
## picture is straightforward to read (no atypical mid-season data gaps).

dh_tbl <- fread(file.path(DATA_DIR, "m1_dethist_table.csv"), colClasses = "character")
cat("\n== Detection history example: Black Agouti, Tefe array ==\n")
cat("Literal 0/1 matrix (first 15 of 20 occasions, all 30 cameras;",
    "blank = camera not yet/no longer active):\n")
print(dh_tbl)

# See figs/m1_dethist_new_example.png for the graphical version of this same
# matrix (grey = active/no detection, red = detected, white = not active).
# Twenty of thirty cameras detected the species at least once; one camera
# shows an unusually long run of consecutive detections — a real behavioral
# signal (a resident individual using that spot repeatedly), not an artifact.


## =============================================================================
## 7. CHOOSING THE OCCASION WINDOW
## =============================================================================
## The occasion length is a MODELLING CHOICE, not a property of the data. It
## trades off two things:
##   - Too short (e.g. 1 day): most occasions are 0 simply because the animal
##     wasn't there THAT DAY, even if it uses the site. Detection probability
##     per occasion is very low and occupancy is underestimated.
##   - Too long (e.g. 30 days): almost every occasion becomes 1 for common
##     species, so the model loses the information that distinguishes sites,
##     and short deployments yield too few occasions to estimate detection.
##
## Different species can warrant different windows (a wide-ranging, low-
## density carnivore visiting every few weeks needs longer than a resident
## rodent photographed nightly), but a single common window keeps every
## species on the same footing in a multi-species pipeline.
##
## WE ADOPT A 7-DAY OCCASION WINDOW as the pipeline default: short enough
## that common species still vary across occasions, long enough that the
## short deployments common in a multi-network dataset still yield several
## occasions.

cat("\n== Occasion window ==\nPipeline default: 7-day occasions ",
    "(see comment block above for the reasoning).\n")


## =============================================================================
## 8. WHAT THE MERGED COMMUNITY LOOKS LIKE
## =============================================================================
## A first look at the mammal community across all three merged networks.
## For every species we compute two simple summaries:
##   - Detection rate    : detections per 100 camera-nights (raw signal strength)
##   - Naive occupancy   : % of sites where the species was detected at least
##                          once (uncorrected for imperfect detection)

modeled <- fread(file.path(DATA_DIR, "m1_modeled_species.csv"))
cat("\n== Merged community ==\n")
cat("modeled species (>=30 detections, >=8 sites):", nrow(modeled), "\n")
cat("(see figs/m1_xmas_mammals.png for the full detection-rate vs.",
    "naive-occupancy figure)\n")

# Naive occupancy always UNDERESTIMATES true occupancy — a species can be
# present but missed on every visit. Correcting this is the job of Module 3
# (single-species occupancy models).


## =============================================================================
## SUMMARY — what the results mean
## =============================================================================
## - Snapshot Brasil's own structure (deployments, placenames, sequences,
##   arrays) and camera trap-type/height are all standard ground-level
##   placements, with no canopy set to exclude.
## - The 60-minute independence threshold, adopted for consistency with the
##   Atlantic network, has a real, measurable cost concentrated in gregarious
##   and site-lingering species — a trade worth knowing about, not a hidden
##   assumption.
## - We have a single, harmonised, MAMMAL-ONLY dataset from three networks,
##   spanning a wide environmental gradient across Brazil.
## - The analytical units (sites, arrays, detection histories) mean the same
##   thing across all three sources — getting the deployment-to-site step
##   right (respecting real camera-active windows, not the full span between
##   servicing visits) matters because every downstream model treats the
##   array as its replication unit.
## - A 7-day occasion window is the pipeline default, balancing information
##   against the short deployments in the merged data.
##
## Next module: Covariates & Sampling Summary — habitat, climate, protected-
## area extraction and collinearity screening.
## =============================================================================
