## =============================================================================
## Module 7 — Temporal Interactions: Diel Overlap & Time-Since-Event (PAMM)
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## SELF-CONTAINED: this script rebuilds every table and figure in the rendered
## report from precomputed summary/curve tables, without re-running the
## computationally heavy activity kernels (activity::fitact, 1000 bootstraps
## each) or the PAMM hazard fits (mgcv piecewise-exponential GAMs) -- exactly
## matching the report's own stated convention ("shown as reference code, not
## re-run"). Put this script in the same folder as the data files listed below
## and run it top to bottom.
##
## Required input files (in DATA_DIR, default "data" alongside this script):
##   temporal_feasibility.csv        -- per-species timestamped-detection counts
##                                       and diel-overlap-bar pass/fail
##   hazard_feasibility_verdict.csv  -- per-driver best well-distributed responder
##                                       pairing and hazard-model feasibility verdict
##   diel_overlap.csv                -- Dhat4 overlap coefficients + bootstrap
##                                       difference tests for the focal species pairs
##   diel_density_curves.csv         -- activity density curves (0-24h) for the
##                                       dog/human/felid diel-overlap figures
##   allspecies_density.csv          -- activity density curves for all 44 modelled
##                                       species (community activity-stack figures)
##   activity_group_key.csv          -- species membership of the 4 PCA habitat groups
##   allspecies_meta.csv             -- species -> order + group + n, camera-ready
##                                       metadata for the overall activity stack
##   biome_meta.csv                  -- species -> biome + order + n, for the
##                                       per-biome activity stack
##   jaguar_activity_curve.csv       -- jaguar's own activity density (reference
##                                       curve overlaid on the peccary 2x2 figure)
##   peccary_jaguar_test.csv         -- Collared Peccary diel overlap test,
##                                       jaguar-present vs jaguar-absent sites
##   peccary_temp_test.csv           -- same test, hot vs cool sites (heat-hypothesis
##                                       control)
##   peccary_temp_curves.csv         -- Collared Peccary activity density, hot vs
##                                       cool sites
##   peccary_2x2_curves.csv          -- Collared Peccary activity density, jaguar x
##                                       temperature 2x2 design
##   peccary_2x2_meta.csv            -- per-cell n/locations/nocturnal-fraction for
##                                       the 2x2 design
##   human_nocturnality.csv          -- per-species nocturnal-fraction shift, low
##                                       vs high human-pressure sites
##   seasonal_shift.csv              -- per-species nocturnal-fraction shift, wet
##                                       vs dry season
##   lunar_effect.csv                -- per-species activity shift, dark vs bright
##                                       moon nights
##   pamm_results.csv                -- pooled + per-subproject PAMM hazard
##                                       coefficients (tsh_coef, HR/week, p) for
##                                       every focal driver-responder pair
##   influence_check.csv             -- per-camera influence check (drop top-2
##                                       cameras, refit) for the two significant
##                                       subproject-level PAMM signals
##   coati_piracicaba_percamera.csv  -- per-camera detection-rate breakdown for
##                                       the Piracicaba coati signal
##
## To actually REFIT the activity kernels or PAMM hazard models from raw
## detections (hours, not seconds -- shown as reference only), set RUN_LIVE
## below to TRUE and see mod7_code/*.R for the exact fitting scripts; this
## requires the `activity`, `overlap`, and `mgcv` R packages plus the raw
## per-pair paired-event-data (PED) tables, which are not shipped here since
## the fits are not intended to be re-run routinely.
##
## Requires: data.table, ggplot2, ggrepel, patchwork.

RUN_LIVE <- FALSE   # TRUE re-fits everything from raw data (hours) via mod7_code/*.R;
                    # FALSE (default) loads the precomputed tables above and rebuilds
                    # every figure/table from them, matching the report's own convention.

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---- 0. Paths ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"
stopifnot(dir.exists(DATA_DIR))
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR)

# green = positive/high, red = negative/low; a consistent handful of hues used
# throughout every figure in this script.
COL_POS   <- "#1a9850"
COL_NEG   <- "#d73027"
COL_NEUTRAL <- "#4575b4"
COL_DOG   <- "#6a3d9a"
COL_HUMAN <- "#e08214"

read_d <- function(f) fread(file.path(DATA_DIR, f))


## =============================================================================
## STEP 1: FEASIBILITY -- can this analysis honestly be run at all?
## =============================================================================
## Before fitting anything: count the qualifying events per species, and if a
## model would rest on a handful of events at one or two cameras, do not fit
## it. Reporting *why* an analysis was not run is itself a legitimate result.

feas <- read_d("temporal_feasibility.csv")
vd   <- read_d("hazard_feasibility_verdict.csv")
cat(sprintf("timestamped detections across the focal species table: %s\n", format(sum(feas$n_timestamped), big.mark=",")))
cat(sprintf("species clearing the diel-overlap bar (>=50 events): %d\n\n", sum(feas$n_timestamped >= 50)))
cat("Hazard-model feasibility by driver:\n")
print(vd[, .(driver, total_timestamped, best_responder, best_events, best_locations, hazard_feasible)])

# ---- feasibility figure: diel-overlap bar (all species) + hazard bar (3 drivers) ----
feas_top <- feas[order(-n_timestamped)][1:10]
feas_top$species <- factor(feas_top$species, levels=rev(feas_top$species))
pA <- ggplot(feas_top, aes(x=species, y=n_timestamped)) +
  geom_col(fill=COL_NEUTRAL) + coord_flip() +
  geom_hline(yintercept=50, color=COL_NEG, linetype="dashed") +
  geom_text(aes(label=n_timestamped), hjust=-0.15, size=3) +
  labs(x=NULL, y="timestamped detections (network-wide)",
       title="A. Diel-overlap feasibility") +
  theme_minimal(base_size=10) + expand_limits(y=max(feas_top$n_timestamped)*1.12)

vd$fill_col <- ifelse(vd$hazard_feasible, COL_POS, COL_NEG)
pB <- ggplot(vd, aes(x=driver, y=best_events, fill=fill_col)) +
  geom_col() + scale_fill_identity() +
  geom_text(aes(label=sprintf("%d ev / %d loc", best_events, best_locations)), vjust=-0.4, size=3) +
  labs(x=NULL, y="best well-distributed driver->native paired events",
       title="B. Hazard feasibility (best responder per driver)") +
  theme_minimal(base_size=10) + expand_limits(y=max(vd$best_events)*1.15)

p_feas <- pA + pB
ggsave(file.path(FIGS_DIR, "temporal_feasibility.png"), p_feas, width=13, height=6, dpi=150)
cat("\nSaved figs/temporal_feasibility.png\n")


## =============================================================================
## STEP 2: DIEL ACTIVITY OVERLAP
## =============================================================================
## Compares WHEN (solar time, 0-24h) each species is active, using Ridout &
## Linkie's Dhat4 overlap coefficient (0 = no overlap, 1 = identical timing).
## The density curves here are precomputed (activity::fitact, 1000 bootstraps,
## reference code in mod7_code/fit_diel.R); this step only re-plots them.

ov  <- read_d("diel_overlap.csv")
dc  <- read_d("diel_density_curves.csv")

dog_ov <- ov[sp1 == "Domestic Dog"][order(Dhat4)]
cat("Diel overlap (Dhat4) between the dog and each native (lower = less overlap):\n")
print(dog_ov[, .(native=sp2, n_native=n2, Dhat4, p_differ)])

# ---- dog-vs-native overlap: three natives spanning the overlap range ----
sel_native <- c("Spotted Paca", "Collared Peccary", "Gray Brocket")
dc_long <- melt(dc[, c("hour","Domestic Dog", sel_native), with=FALSE], id.vars="hour",
                 variable.name="species", value.name="density")
p_dog <- ggplot(dc_long, aes(x=hour, y=density, color=species)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=setNames(c(COL_DOG, "#238b45","#fdae61","#3182bd"),
                                      c("Domestic Dog", sel_native))) +
  labs(x="hour (solar time)", y="activity density",
       title="Diel activity overlap: Domestic Dog vs three natives") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR, "activity_overlap_dog.png"), p_dog, width=9, height=5.5, dpi=150)
cat("Saved figs/activity_overlap_dog.png\n")

# ---- Ocelot x Puma ----
dc_felid <- melt(dc[, c("hour","Puma","Ocelot"), with=FALSE], id.vars="hour",
                  variable.name="species", value.name="density")
p_felid <- ggplot(dc_felid, aes(x=hour, y=density, color=species)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Puma"="#e6550d","Ocelot"="#31a354")) +
  labs(x="hour (solar time)", y="activity density",
       title="Ocelot x Puma diel overlap (Dhat4=0.83, p<0.001)") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR, "activity_overlap_felids.png"), p_felid, width=8, height=5.5, dpi=150)
cat("Saved figs/activity_overlap_felids.png\n")


## =============================================================================
## STEP 3: THE WHOLE COMMUNITY'S DAY -- stacked activity, event-weighted
## =============================================================================
## Height = summed, event-weighted activity across species; colour = taxonomic
## order, shade = species within order. Built from precomputed per-species
## activity densities (mod7_data/allspecies_density.csv), weighted by each
## species' own network-wide detection count (n) so common species dominate
## the visible stack the way they dominate the camera record.

allsp_dens <- read_d("allspecies_density.csv")
allsp_meta <- read_d("allspecies_meta.csv")
group_key  <- read_d("activity_group_key.csv")

order_pal <- c(Artiodactyla="#3182bd", Carnivora="#e6550d", Rodentia="#31a354",
               Didelphimorphia="#756bb1", Cingulata="#636363", Pilosa="#a1d99b",
               Perissodactyla="#fd8d3c", Primates="#fdae6b", Lagomorpha="#bdbdbd")

build_stack <- function(dens_wide, meta, weight_col="n", title="") {
  # dens_wide: hour x species density matrix; meta: species -> order/group/n
  long <- melt(dens_wide, id.vars="hour", variable.name="species", value.name="density")
  long[, species := as.character(species)]
  long <- merge(long, meta[, .(species=sci, order, w=get(weight_col))], by="species")
  long[, weighted := density * w]
  agg <- long[, .(activity=sum(weighted)), by=.(hour, order)]
  agg[, order := factor(order, levels=names(order_pal))]
  ggplot(agg, aes(x=hour, y=activity, fill=order)) +
    geom_area(position="stack", alpha=0.85) +
    scale_fill_manual(values=order_pal, drop=FALSE) +
    labs(x="hour (solar time)", y="event-weighted activity (summed)", fill="Order", title=title) +
    theme_minimal(base_size=11)
}

p_overall <- build_stack(allsp_dens, allsp_meta, title="Total daily activity, all 44 modelled species")
ggsave(file.path(FIGS_DIR, "activity_stack_overall.png"), p_overall, width=10, height=6, dpi=150)
cat("Saved figs/activity_stack_overall.png\n")

# ---- by the 4 Module 3 PCA habitat groups ----
p_by_group_list <- list()
for (g in sort(unique(group_key$group))) {
  gk <- group_key[group == g]
  sp_in_group <- gk$species
  cols_present <- intersect(c("hour", sp_in_group), names(allsp_dens))
  if (length(cols_present) < 2) next
  sub_meta <- allsp_meta[sci %in% sp_in_group]
  lbl <- unique(gk$group_label)
  p_by_group_list[[as.character(g)]] <- build_stack(allsp_dens[, ..cols_present], sub_meta, title=lbl)
}
p_group <- wrap_plots(p_by_group_list, ncol=2) + plot_layout(guides="collect")
ggsave(file.path(FIGS_DIR, "activity_stack_bygroup.png"), p_group, width=13, height=9, dpi=150)
cat("Saved figs/activity_stack_bygroup.png\n")

# ---- key table: species membership of the four groups ----
group_summary <- group_key[, .N, by=.(group, group_label)][order(group)]
setnames(group_summary, "N", "species_n")
cat("\nSpecies membership of the four habitat community groups:\n")
print(group_summary)

# ---- by WWF biome ----
biome_meta <- read_d("biome_meta.csv")
biomes_to_show <- setdiff(unique(biome_meta$biome), "Mangrove")  # too few events, omitted per report
p_by_biome_list <- list()
for (bm_name in biomes_to_show) {
  bmeta <- biome_meta[biome == bm_name]
  sp_in_biome <- bmeta$sci
  cols_present <- intersect(c("hour", unique(bmeta$common)), names(allsp_dens))
  # allsp_density columns are scientific names, not common -- match on sci instead
  cols_present <- intersect(c("hour", sp_in_biome), names(allsp_dens))
  if (length(cols_present) < 2) next
  p_by_biome_list[[bm_name]] <- build_stack(allsp_dens[, ..cols_present], bmeta, title=bm_name)
}
p_biome <- wrap_plots(p_by_biome_list, ncol=2) + plot_layout(guides="collect")
ggsave(file.path(FIGS_DIR, "activity_stack_bybiome.png"), p_biome, width=13, height=9, dpi=150)
cat("Saved figs/activity_stack_bybiome.png\n")


## =============================================================================
## STEP 4: COLLARED PECCARY -- JAGUAR-PRESENT VS JAGUAR-ABSENT SITES
## =============================================================================
## The within-species contrast is the honest test for a behavioural response
## (vs an intrinsic difference in niche): does the SAME species act differently
## depending on whether its main predator is around?

jt <- read_d("peccary_jaguar_test.csv")
cat(sprintf("Peccary diel overlap by jaguar presence: Dhat4=%.2f (p=%.4f); nocturnal present %.0f%% vs absent %.0f%%\n",
            jt$Dhat4, jt$p_differ, 100*jt$noct_present, 100*jt$noct_absent))
tt2 <- read_d("peccary_temp_test.csv")
cat(sprintf("Same test, by temperature instead: Dhat4=%.2f (p=%.4f); nocturnal hot %.0f%% vs cool %.0f%%\n",
            tt2$Dhat4, tt2$p_differ, 100*tt2$noct_hot, 100*tt2$noct_cool))

# jaguar-present/absent + hot/cool contrast, side by side
pjt_curves <- read_d("peccary_2x2_curves.csv")  # reuses the jaguar dimension from the 2x2 design
jag_curve  <- read_d("jaguar_activity_curve.csv")
ptemp_curves <- read_d("peccary_temp_curves.csv")

d1 <- melt(pjt_curves[, .(hour, `Jaguar present`=jaguar_hot, `Jaguar absent`=nojaguar_hot)],
           id.vars="hour", variable.name="condition", value.name="density")
p_jag <- ggplot(d1, aes(x=hour, y=density, color=condition)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Jaguar present"="#e6550d", "Jaguar absent"="#31a354")) +
  labs(x="hour (solar time)", y="activity density", color=NULL,
       title="Collared Peccary: jaguar present vs absent") +
  theme_minimal(base_size=10)

d2 <- melt(ptemp_curves[, .(hour, Hot=hot, Cool=cool)], id.vars="hour",
           variable.name="condition", value.name="density")
p_temp <- ggplot(d2, aes(x=hour, y=density, color=condition)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Hot"="#e6550d", "Cool"="#3182bd")) +
  labs(x="hour (solar time)", y="activity density", color=NULL,
       title="Collared Peccary: hot vs cool sites (temperature control)") +
  theme_minimal(base_size=10)

p_jag_temp <- p_jag + p_temp
ggsave(file.path(FIGS_DIR, "peccary_jaguar_diel_contrast.png"), p_jag_temp, width=13, height=5.5, dpi=150)
cat("Saved figs/peccary_jaguar_diel_contrast.png\n")

# ---- 2x2 design: jaguar (yes/no) x temperature (hot/cool) ----
meta2x2 <- read_d("peccary_2x2_meta.csv")
cat("\n2x2 design cell sizes (n / locations / nocturnal fraction):\n")
print(meta2x2)

# The "jaguar present + cool" cell (44 peccary detections, 2 locations) is NOT
# treated as a real condition: those 2 "cool" sites sit at 23.2-23.4C, right on
# the hot/cool temperature split, and jaguars in this dataset are otherwise
# almost entirely absent from cool sites (only 2 of 244 jaguar detections
# network-wide). Rather than plot a density curve built from that thin,
# borderline sample -- which produces a spurious sharp spike, an artifact of
# kernel smoothing on ~44 points -- this panel is deliberately left blank with
# an explanatory annotation, matching the report's own stated design.
d2x2 <- melt(pjt_curves[, .(hour, nojaguar_cool, nojaguar_hot, jaguar_hot)],
             id.vars=c("hour"), variable.name="cell", value.name="density")
d2x2[, jaguar_status := ifelse(grepl("^nojaguar", cell), "Jaguar absent", "Jaguar present")]
d2x2[, temp := ifelse(grepl("hot$", cell), "Hot", "Cool")]
d2x2[, jaguar_status := factor(jaguar_status, levels=c("Jaguar absent","Jaguar present"))]
d2x2[, temp := factor(temp, levels=c("Cool","Hot"))]

# jag_curve's own "jaguar" column (a density value) must not share a name with the
# facet variable -- facet_grid() resolves panel membership per-layer, and a same-named
# column in a different layer's data (however unrelated) gets swept in as if it were
# the facet key, silently exploding into one facet per unique density value.
# The jaguar reference curve is only drawn in the 3 real panels; the blank
# "Jaguar present x Cool" panel gets no curves at all, matching the original design.
jag_curve_ref <- copy(jag_curve)
setnames(jag_curve_ref, "jaguar", "density_ref")
panel_combos <- data.table(
  jaguar_status = factor(c("Jaguar absent","Jaguar absent","Jaguar present"), levels=levels(d2x2$jaguar_status)),
  temp          = factor(c("Cool","Hot","Hot"), levels=levels(d2x2$temp))
)
jag_curve_full <- rbindlist(lapply(seq_len(nrow(panel_combos)), function(i) {
  cbind(jag_curve_ref[, .(hour, density_ref)],
        jaguar_status = panel_combos$jaguar_status[i], temp = panel_combos$temp[i])
}))

blank_note <- data.table(
  jaguar_status = factor("Jaguar present", levels=levels(d2x2$jaguar_status)),
  temp = factor("Cool", levels=levels(d2x2$temp)),
  x = 12, y = 0.25,
  label = "jaguars essentially absent from cool sites\n\nonly 2 of 244 jaguar detections network-wide come from cool sites;\nthe 2 'cool' peccary sites sit at 23.2-23.4C,\nright on the hot/cool split -- not a real condition"
)

p_2x2 <- ggplot(d2x2, aes(x=hour, y=density)) +
  geom_line(color=COL_NEUTRAL, linewidth=1) +
  geom_line(data=jag_curve_full, aes(x=hour, y=density_ref), color="black", linetype="dotted", linewidth=0.8) +
  geom_text(data=blank_note, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            size=3, color="grey40", lineheight=0.9) +
  facet_grid(jaguar_status ~ temp, drop=FALSE) +
  labs(x="hour (solar time)", y="Collared Peccary activity density",
       title="Collared Peccary activity: jaguar presence x temperature",
       subtitle="dotted line = jaguar's own activity curve (hot-forest sites, n=242)") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "peccary_jaguar_temp_2x2.png"), p_2x2, width=9, height=8, dpi=150)
cat("Saved figs/peccary_jaguar_temp_2x2.png\n")


## =============================================================================
## STEP 5: THREE MORE AXES THE MERGED DATA CAN CARRY
## =============================================================================
## (1) Human-induced nocturnality (Gaynor-style): does a native shift nocturnal
##     under high human pressure?
## (2) Seasonal activity shift: wet vs dry-season diel curves per species.
## (3) Lunar effect: do strongly-nocturnal species avoid bright moonlit nights?

# ---- (1) human-induced nocturnality ----
hn <- read_d("human_nocturnality.csv")
hn[, sig := p_differ < 0.05]
hn[, common := factor(common, levels=common[order(dnoct)])]
p_hn <- ggplot(hn, aes(x=common, y=dnoct, color=sig)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_segment(aes(xend=common, y=0, yend=dnoct)) +
  geom_point(size=3) +
  scale_color_manual(values=c("TRUE"=COL_NEG, "FALSE"="grey60"), guide="none") +
  coord_flip() +
  labs(x=NULL, y="change in nocturnal fraction (high - low human pressure)",
       title="Human-induced nocturnality: 12 species with sufficient data") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "human_nocturnality.png"), p_hn, width=9, height=6, dpi=150)
n_toward_night <- sum(hn$dnoct > 0); n_away <- sum(hn$dnoct < 0)
cat(sprintf("Human-induced nocturnality: %d of %d species shift toward night, %d away\n", n_toward_night, nrow(hn), n_away))
cat("Saved figs/human_nocturnality.png\n")

# ---- (2) seasonal shift ----
ss <- read_d("seasonal_shift.csv")
ss[, dnoct := noct_wet - noct_dry]
ss[, sig := p_differ < 0.05]
top6 <- ss[sig == TRUE][order(-abs(dnoct))][1:min(6,.N)]
cat(sprintf("\nSeasonal shift: %d of %d species show a significant wet/dry difference\n", sum(ss$sig), nrow(ss)))
print(top6[, .(common, n_wet, n_dry, Dhat4, dnoct)])
p_seasonal <- ggplot(top6, aes(x=reorder(common, dnoct), y=dnoct)) +
  geom_col(fill=COL_NEUTRAL) + coord_flip() +
  labs(x=NULL, y="nocturnal-fraction shift (wet - dry)",
       title="Seasonal activity shift: 6 largest significant changes") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "seasonal_shift.png"), p_seasonal, width=9, height=5, dpi=150)
cat("Saved figs/seasonal_shift.png\n")

# ---- (3) lunar effect ----
# illum_shift = each species' own mean detection-night illumination minus the
# community-night baseline (mean_illum_null); negative = lunar-phobic (fewer
# detections on bright nights), positive = lunar-philic. Significance is the
# `p` column (bootstrap test of illum_shift against 0).
le <- read_d("lunar_effect.csv")
sig_le <- le[p < 0.05]
cat(sprintf("\nLunar effect: %d of %d strongly-nocturnal species show a significant moon response\n", nrow(sig_le), nrow(le)))
print(sig_le[, .(common, n_night, n_bright, n_dark, illum_shift, p)])
le[, common := factor(common, levels=common[order(illum_shift)])]
p_lunar <- ggplot(le, aes(x=common, y=illum_shift, color=p<0.05)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_segment(aes(xend=common, y=0, yend=illum_shift)) +
  geom_point(size=3) +
  scale_color_manual(values=c("TRUE"=COL_NEG, "FALSE"="grey60"), guide="none") +
  coord_flip() +
  labs(x=NULL, y="moon-illumination shift vs community-night baseline (negative = lunar-phobic)",
       title="Lunar illumination effect: 16 strongly-nocturnal species") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "lunar_effect.png"), p_lunar, width=9, height=7, dpi=150)
cat("Saved figs/lunar_effect.png\n")


## =============================================================================
## STEP 6: TIME-SINCE-EVENT HAZARD MODELS (PAMM)
## =============================================================================
## Does a native's short-term detection rate change in the days after a human
## or dog passes the same camera? Piecewise-additive mixed model (binomial
## GAM): y ~ tsh + s(hour,cc) + subproject + s(camera,re). Coefficients here
## are loaded from the precomputed fits (mod7_code/fit_pamm.R, reference only).

res <- read_d("pamm_results.csv")

# ---- human: pooled, network-wide ----
pooled_h <- res[scope == "pooled" & grepl("Human", pair)]
cat("Pooled time-since-human hazard (no native shows a significant network-wide effect):\n")
print(pooled_h[, .(pair, events, HR_per_week, p)])
cat(sprintf("significant pooled human effects (p<0.05): %d of %d\n\n", sum(pooled_h$p < 0.05), nrow(pooled_h)))

pooled_h[, sig := p < 0.05]
p_pamm_human <- ggplot(pooled_h, aes(x=reorder(pair, HR_per_week), y=HR_per_week)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey50") +
  geom_point(size=3, color=COL_NEUTRAL) + coord_flip() +
  labs(x=NULL, y="hazard ratio per week (1 = no effect)",
       title="Time-since-human hazard: four focal natives, network-wide") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "pamm_human.png"), p_pamm_human, width=9, height=5, dpi=150)
cat("Saved figs/pamm_human.png\n")

# ---- Coati at Piracicaba: a local signal worth chasing, and its influence check ----
inf <- read_d("influence_check.csv")
coati_inf <- inf[grepl("Coati", signal)]
cat("\nPer-camera influence check, Piracicaba coati signal:\n")
print(coati_inf[, .(signal, full_p, drop2_p, robust)])

coati_pc <- read_d("coati_piracicaba_percamera.csv")
coati_pc_long <- melt(coati_pc, id.vars=c("camera","dets","bins"),
                       measure.vars=c("rate_early","rate_late"),
                       variable.name="period", value.name="rate")
coati_pc_long[, period := ifelse(period=="rate_early", "first days after human", "later")]
p_coati_a <- ggplot(pooled_h[grepl("Coati", pair)] , aes(x=1, y=HR_per_week)) +
  geom_point(size=4, color=COL_NEUTRAL) +
  annotate("text", x=1, y=pooled_h[grepl("Coati",pair)]$HR_per_week, label="pooled", vjust=-1.5) +
  labs(title="Coati detection rate\nvs days since human", x=NULL, y="HR per week") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_blank())
p_coati_b <- ggplot(coati_pc_long, aes(x=reorder(camera, -dets), y=rate, fill=period)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c("first days after human"=COL_HUMAN, "later"="grey60")) +
  labs(x=NULL, y="detections per 1,000 bins", fill=NULL,
       title="Per-camera breakdown: signal concentration") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=30, hjust=1))
p_coati <- p_coati_a + p_coati_b + plot_layout(widths=c(1,2))
ggsave(file.path(FIGS_DIR, "pamm_coati.png"), p_coati, width=11, height=5.5, dpi=150)
cat("Saved figs/pamm_coati.png\n")

# ---- Dogs and Spotted Paca ----
dog_res <- res[grepl("Dog", pair)]
cat("\nTime-since-dog model for Spotted Paca (pooled and best-sampled subproject):\n")
print(dog_res[, .(scope, events, HR_per_week, p)])

dog_pooled <- dog_res[scope=="pooled"]
dog_sub    <- dog_res[scope!="pooled"]
p_dog_pamm <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", color="grey50") +
  geom_point(data=dog_pooled, aes(x="pooled", y=HR_per_week), size=4, color=COL_DOG) +
  geom_point(data=dog_sub, aes(x="best subproject\n(dashed line)", y=HR_per_week), size=4, color=COL_DOG, shape=17) +
  labs(x=NULL, y="hazard ratio per week (1 = no effect)",
       title="Spotted Paca detection rate vs days since a dog") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR, "pamm_dog_paca.png"), p_dog_pamm, width=8, height=5.5, dpi=150)
cat("Saved figs/pamm_dog_paca.png\n")

cat("\n=== Module 7 script complete: feasibility screen, diel overlap, activity stacks,\n")
cat("    peccary/jaguar contrast, 3 extended axes, and PAMM results all rebuilt. ===\n")
cat("See module7_temporal_pamm.qmd for the full narrated analysis and citations.\n")
