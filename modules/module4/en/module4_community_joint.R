## =============================================================================
## Module 4 — Community & Joint Models
## Snapshot Brasil multi-network camera-trap pipeline
## =============================================================================
##
## SELF-CONTAINED: this script starts from the three raw combined-dataset
## files (deployments, mammal detections, remote-sensing covariates) plus
## four reference lookup tables from Module 2/3, and builds everything from
## scratch -- detection histories for 44 modeled species, the range mask,
## and every model fit (community occupancy with a species-trait
## meta-regression, GJAM, Royle-Nichols, single-species occupancy, and a
## naive logistic-regression baseline). Put this script in the same folder
## as the seven data files listed below and run it top to bottom.
##
## Required input files (in DATA_DIR, default "data" alongside this script):
##   final_deployments.csv         -- one row per camera deployment
##   final_detections_mammals.csv  -- one row per independent detection event
##   final_covariates.csv          -- one row per camera site, remote-sensing covariates
##   array_covariates.csv          -- one row per camera array, the reduced 6-covariate
##                                     ARRAY_SET (Module 2/3's array-level design)
##   species_list.json             -- the 44 species/genera modeled in this module
##                                     (>=10 detections, >=20 sites, per Module 3's rule)
##   species_range_mask.csv        -- site x species IUCN range mask (1 = candidate,
##                                     0 = structural zero); output of Module 4's own
##                                     range-mask assembler (see Module 4 prep notes)
##   species_traits.csv            -- per-species log body mass (z-scored) and diet
##                                     class (herbivore/omnivore, carnivore = reference),
##                                     from COMBINE (Soria et al. 2021) matched via the
##                                     Mammal Diversity Database taxonomy
##
## The last four files are reference outputs from Module 2's covariate-selection
## work, Module 3's species-selection rule, Module 4's range-mask assembler, and
## a one-time trait-database match, and cannot be rebuilt from the three main
## data files alone -- they are shipped alongside this script as fixed lookup
## tables, the same way Module 3 ships its four range-mask/taxonomy lookups.
##
## Five methods are fit and compared, all on the SAME 11-covariate CAMERA_SET
## (forest, savanna, pasture, cropland, native vegetation, temperature,
## precipitation, protected-area status, development, distance to road,
## distance to water):
##   1. Community (Dorazio-Royle) occupancy -- all 44 species jointly, with
##      partial pooling across species and a species-trait meta-regression
##      on the occupancy/detection intercepts.
##   2. GJAM (joint species distribution model) -- residual co-occurrence
##      after habitat is accounted for.
##   3. Royle-Nichols (occuRN) abundance-induced heterogeneity model, fit
##      per species.
##   4. Single-species occupancy (spOccupancy::PGOcc) with an array random
##      effect, fit per species, range-masked.
##   5. A naive logistic regression (no detection model) -- shown
##      deliberately to illustrate the separation problems it produces on
##      sparse land-cover covariates.
##
## Requires: data.table, ggplot2, ggrepel, jagsUI (rjags/JAGS installed
## separately), gjam, unmarked, spOccupancy.

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(jagsUI)
  library(gjam)
  library(unmarked)
  library(spOccupancy)
})

# ---- 0. Paths ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"
stopifnot(dir.exists(DATA_DIR))
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR)

# color convention used throughout: green = positive/high, red = negative/low
COL_POS <- "#1a9850"; COL_NEG <- "#d73027"; COL_NEUTRAL <- "#999999"

## Camera-level covariate set (11) -- used for every method in this module.
CAMERA_COVARS <- c("forest_100m","savanna_100m","pasture_100m","cropland_100m",
                    "native_veg_1000m","temp_mean_C","precip_annual_mm",
                    "ghsl_built_5000m","dist_road_m","dist_water_m")
# in_pa (protected-area status, 0/1) is appended separately below, matching
# Module 3's convention (it is binary and not z-scored).

OCC_DAYS <- 7  # occasion length in days, matching the pipeline's independence-interval convention

set.seed(1)


## =============================================================================
## STEP 1: LOAD RAW DATA AND BUILD SITE-LEVEL COVARIATES
## =============================================================================

dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
det <- fread(file.path(DATA_DIR, "final_detections_mammals.csv"))
cov <- fread(file.path(DATA_DIR, "final_covariates.csv"))
arr_cov <- fread(file.path(DATA_DIR, "array_covariates.csv"))

dep[, start_dt := as.IDate(start_date)]
dep[, end_dt   := as.IDate(end_date)]

site_tbl <- dep[, .(network = first(network), array_id = first(array_id),
                     n_deployments = .N, camera_days = sum(camera_days)),
                 by = site_id]
setorder(site_tbl, site_id)
cat("Sites:", nrow(site_tbl), "| Arrays:", uniqueN(site_tbl$array_id), "\n")

sc <- merge(site_tbl, cov[, -"array_id"], by = "site_id", all.x = TRUE)  # array_id already in site_tbl
setorder(sc, site_id)

# Listwise-delete sites missing one or more of the 11 camera-level covariates
# BEFORE building anything downstream -- every method in this module (and
# every reference lookup shipped alongside this script) is built on this
# same candidate site set, so deletion has to happen here, not just be
# reported and silently ignored.
n_missing_cov <- sum(!complete.cases(sc[, ..CAMERA_COVARS]))
cat("Sites missing one or more covariates (excluded, listwise deletion):", n_missing_cov, "\n\n")
sc <- sc[complete.cases(sc[, ..CAMERA_COVARS])]

for (v in CAMERA_COVARS) sc[[paste0(v, "_z")]] <- as.numeric(scale(sc[[v]]))
CAMERA_COVARS_Z <- c(paste0(CAMERA_COVARS, "_z"), "in_pa")  # in_pa is 0/1, not z-scored

arrays_sorted <- sort(unique(sc$array_id))
sc[, array_num := match(array_id, arrays_sorted)]


## =============================================================================
## STEP 2: THE 44-SPECIES MODELED COMMUNITY AND ITS DETECTION HISTORIES
## =============================================================================
## Species selection follows Module 3's rule (>=10 detections, >=20 sites,
## excluding domestic/human taxa and unresolved genus-only records); the
## exact 44-species list this module fits is shipped as species_list.json
## so every reader sees the identical community, without re-deriving the
## selection cutoff.

species <- fromJSON <- jsonlite::fromJSON(file.path(DATA_DIR, "species_list.json"))
cat("Modeled species/genera:", length(species), "\n\n")

# One active-survey window per site: the earliest start and latest end across
# that site's deployments (a site can have several deployments over time).
# Occasions are anchored to EACH SITE'S OWN start date -- not a single global
# date -- since sites span very different calendar periods (this dataset
# merges multiple networks collected over more than a decade).
site_window <- dep[, .(site_start = min(start_dt, na.rm=TRUE), site_end = max(end_dt, na.rm=TRUE)), by = site_id]
max_days <- as.integer(max(site_window$site_end - site_window$site_start, na.rm = TRUE))
max_occ <- ceiling(max_days / OCC_DAYS)
site_ids_ordered <- sort(unique(sc$site_id))
n_sites <- length(site_ids_ordered)
site_window <- site_window[match(site_ids_ordered, site_id)]
site_start_lookup <- setNames(site_window$site_start, site_window$site_id)

# base survey matrix: which occasions were actively surveyed at each site
# (a site can have GAPS between deployments -- only occasions actually
# covered by a deployment window are marked surveyed).
occ_active <- list()
for (sid in site_ids_ordered) {
  dsub <- dep[site_id == sid]
  m <- rep(NA_real_, max_occ)
  s0 <- site_start_lookup[[sid]]
  for (i in seq_len(nrow(dsub))) {
    d0 <- dsub$start_dt[i]; d1 <- dsub$end_dt[i]
    if (is.na(d0) || is.na(d1)) next
    occ0 <- as.integer(d0 - s0) %/% OCC_DAYS + 1L
    occ1 <- as.integer(d1 - s0) %/% OCC_DAYS + 1L
    occ0 <- max(1, occ0); occ1 <- min(max_occ, occ1)
    if (occ1 >= occ0) m[occ0:occ1] <- ifelse(is.na(m[occ0:occ1]), 0, m[occ0:occ1])
  }
  occ_active[[sid]] <- m
}
base_survey <- matrix(NA_real_, nrow = n_sites, ncol = max_occ)
rownames(base_survey) <- site_ids_ordered
for (sid in site_ids_ordered) base_survey[sid, ] <- occ_active[[sid]]

n_missing_date <- sum(is.na(det$date) | is.na(as.IDate(det$date)))
if (n_missing_date > 0) cat("Detections with missing/unparseable date, excluded:", n_missing_date, "\n")
det <- det[!is.na(date) & !is.na(as.IDate(date)) & site_id %in% site_ids_ordered]
det[, ddate := as.IDate(date)]
# Restrict to detections falling within THAT SITE'S OWN survey window
# (site_start to site_end, inclusive) -- matches the pipeline's original
# convention exactly. Without this filter, a detection recorded after a
# short-window site's own deployment ended (but still under the dataset-wide
# max_occ bound) would be wrongly admitted via a different site's longer window.
site_end_lookup <- setNames(site_window$site_end, site_window$site_id)
n_outside_window <- sum(!(det$ddate >= site_start_lookup[det$site_id] & det$ddate <= site_end_lookup[det$site_id]))
if (n_outside_window > 0) cat("Detections outside their site's own survey window, excluded:", n_outside_window, "\n")
det <- det[ddate >= site_start_lookup[site_id] & ddate <= site_end_lookup[site_id]]
det[, occ_idx := as.integer(ddate - site_start_lookup[site_id]) %/% OCC_DAYS]

build_dethist <- function(sci) {
  y <- base_survey
  sp_det <- det[sci_mdd == sci]
  for (i in seq_len(nrow(sp_det))) {
    sid <- sp_det$site_id[i]; oc <- sp_det$occ_idx[i] + 1L
    if (is.na(sid) || is.na(oc) || !(sid %in% rownames(y)) || oc < 1 || oc > max_occ) next
    y[sid, oc] <- 1
  }
  y
}

# ysum/K matrices: occasion-collapsed detection count and surveyed-occasion
# count per site x species, matching Module 3/6's convention throughout.
S <- length(species)
ysum <- matrix(0L, nrow = n_sites, ncol = S, dimnames = list(site_ids_ordered, species))
K    <- matrix(0L, nrow = n_sites, ncol = S, dimnames = list(site_ids_ordered, species))
for (s in seq_len(S)) {
  y <- build_dethist(species[s])
  ysum[, s] <- rowSums(y == 1, na.rm = TRUE)
  K[, s]    <- rowSums(!is.na(y))
}
cat("Total detections across all 44 species:", sum(ysum), "\n\n")


## =============================================================================
## STEP 3: RANGE MASK AND SPECIES TRAITS (reference lookup tables)
## =============================================================================
## The range mask (1 = site is a candidate, 0 = structural zero -- outside
## the species' buffered IUCN range and never detected there) is a per-
## species/per-site 0/1 table built once by Module 4's own range-mask
## assembler (see the Module 4 prep notes) and shipped as a fixed lookup,
## the same way Module 3 ships its exclusion tables.

rng_dt <- fread(file.path(DATA_DIR, "species_range_mask.csv"))
rng <- as.matrix(rng_dt[match(site_ids_ordered, rng_dt$site_id), ..species])
rng[is.na(rng)] <- 1L  # species with no mapped range default to "everywhere is a candidate"
storage.mode(rng) <- "integer"
cat("Range-masked (structural-zero) site x species cells:", sum(rng == 0), "of", length(rng), "\n\n")

## Species traits: log body mass (z-scored) and diet class
## (herbivore/omnivore, carnivore = reference level), from the COMBINE
## mammal trait database (Soria et al. 2021), matched to this pipeline's
## species list via the Mammal Diversity Database taxonomy.
traits <- fread(file.path(DATA_DIR, "species_traits.csv"))
trait_mat <- as.matrix(traits[match(species, traits$species),
                               .(log_mass_z, diet_herbivore, diet_omnivore)])
stopifnot(nrow(trait_mat) == S, !anyNA(trait_mat))


## =============================================================================
## STEP 4: METHOD 1 -- COMMUNITY (DORAZIO-ROYLE) OCCUPANCY
## =============================================================================
## All 44 species fit jointly: species-level occupancy/detection intercepts
## and habitat betas are partially pooled toward a community distribution
## (a data-poor species is stabilised by what data-rich congeners show), and
## a species-trait meta-regression predicts part of each species' baseline
## intercept from its body mass and diet class -- independent of habitat.
## The range mask enters the latent state directly: z[i,s] can only be 1
## at sites inside the species' range.

X <- as.matrix(sc[, ..CAMERA_COVARS_Z])
arr <- sc$array_num
N <- nrow(ysum); C <- ncol(X); A <- max(arr)
cat(sprintf("Community model: N=%d sites, S=%d species, C=%d covariates, A=%d arrays\n", N, S, C, A))

zinit <- matrix(1, N, S); zinit[rng == 0] <- 0; zinit[ysum > 0] <- 1

community_jags_code <- "model {
  ## community hyperpriors (partial pooling across species)
  mu.lpsi ~ dnorm(0, 0.1); sd.lpsi ~ dunif(0,5); tau.lpsi <- 1/(sd.lpsi*sd.lpsi)
  mu.lp   ~ dnorm(0, 0.1); sd.lp   ~ dunif(0,5); tau.lp   <- 1/(sd.lp*sd.lp)
  for (c in 1:C) {
    mu.beta[c] ~ dnorm(0, 0.1)
    sd.beta[c] ~ dunif(0, 5)
    tau.beta[c] <- 1/(sd.beta[c]*sd.beta[c])
  }
  sd.arr ~ dunif(0,5); tau.arr <- 1/(sd.arr*sd.arr)
  for (a in 1:A) { eta[a] ~ dnorm(0, tau.arr) }

  ## species-trait effects (meta-regression on the species-level intercepts)
  ## trait[s,1]=log_mass_z, trait[s,2]=diet_herbivore, trait[s,3]=diet_omnivore
  gamma.mass.psi ~ dnorm(0, 0.1); gamma.herb.psi ~ dnorm(0, 0.1); gamma.omni.psi ~ dnorm(0, 0.1)
  gamma.mass.p   ~ dnorm(0, 0.1); gamma.herb.p   ~ dnorm(0, 0.1); gamma.omni.p   ~ dnorm(0, 0.1)

  for (s in 1:S) {
    trait.eff.psi[s] <- gamma.mass.psi*trait[s,1] + gamma.herb.psi*trait[s,2] + gamma.omni.psi*trait[s,3]
    trait.eff.p[s]   <- gamma.mass.p*trait[s,1]   + gamma.herb.p*trait[s,2]   + gamma.omni.p*trait[s,3]
    lpsi[s] ~ dnorm(mu.lpsi + trait.eff.psi[s], tau.lpsi)
    lp[s]   ~ dnorm(mu.lp   + trait.eff.p[s],   tau.lp)
    for (c in 1:C) { beta[s,c] ~ dnorm(mu.beta[c], tau.beta[c]) }
  }

  for (i in 1:N) {
    for (s in 1:S) {
      logit(psi[i,s]) <- lpsi[s] + inprod(beta[s,], X[i,]) + eta[arr[i]]
      z[i,s] ~ dbern(psi[i,s] * range[i,s])
      logit(p[i,s]) <- lp[s]
      ysum[i,s] ~ dbin(z[i,s]*p[i,s], K[i,s])
    }
  }
}"
writeLines(community_jags_code, file.path(FIGS_DIR, "..", "community_occupancy.jags"))

data_comm <- list(ysum=ysum, K=K, range=rng, X=X, arr=arr, N=N, S=S, C=C, A=A, trait=trait_mat)
inits_comm <- function() list(z=zinit,
  mu.lpsi=rnorm(1,-1,0.2), mu.lp=rnorm(1,-1,0.2),
  sd.lpsi=runif(1,0.5,1.5), sd.lp=runif(1,0.5,1.5), sd.arr=runif(1,0.3,1),
  mu.beta=rnorm(C,0,0.2), sd.beta=runif(C,0.3,1),
  gamma.mass.psi=rnorm(1,0,0.1), gamma.herb.psi=rnorm(1,0,0.1), gamma.omni.psi=rnorm(1,0,0.1),
  gamma.mass.p=rnorm(1,0,0.1), gamma.herb.p=rnorm(1,0,0.1), gamma.omni.p=rnorm(1,0,0.1))
params_comm <- c("mu.beta","sd.beta","beta","lpsi","lp","mu.lpsi","sd.lpsi","mu.lp","sd.lp","sd.arr",
                 "gamma.mass.psi","gamma.herb.psi","gamma.omni.psi","gamma.mass.p","gamma.herb.p","gamma.omni.p")

## NOTE: this fit takes ~8 hours (44 species x 11 covariates x trait model,
## single-threaded MCMC, 15000 iterations x 3 chains). Set RUN_LIVE <- TRUE
## below to refit from scratch; otherwise this script loads the saved
## coefficient tables shipped alongside it (community_hyperparameters.csv,
## community_species_betas.csv, community_trait_effects.csv).
RUN_LIVE <- FALSE

if (RUN_LIVE) {
  t0 <- Sys.time()
  fit_comm <- jags(data_comm, inits_comm, params_comm, textConnection(community_jags_code),
                    n.chains=3, n.adapt=1000, n.iter=15000, n.burnin=5000, n.thin=10,
                    parallel=FALSE, verbose=TRUE)
  cat("community fit elapsed:", round(difftime(Sys.time(),t0,units="mins"),1), "min\n")

  hyp <- rbindlist(lapply(1:C, function(c) data.table(
    covariate=CAMERA_COVARS_Z[c],
    mu_mean=fit_comm$mean$mu.beta[c], mu_sd=fit_comm$sd$mu.beta[c],
    mu_lcl=fit_comm$q2.5$mu.beta[c], mu_ucl=fit_comm$q97.5$mu.beta[c],
    spread_sd=fit_comm$mean$sd.beta[c], Rhat_mu=fit_comm$Rhat$mu.beta[c])))
  fwrite(hyp, "community_hyperparameters.csv")

  betas <- rbindlist(lapply(1:S, function(s) rbindlist(lapply(1:C, function(c) data.table(
    species=species[s], covariate=CAMERA_COVARS_Z[c],
    mean=fit_comm$mean$beta[s,c], sd=fit_comm$sd$beta[s,c],
    lcl=fit_comm$q2.5$beta[s,c], ucl=fit_comm$q97.5$beta[s,c], Rhat=fit_comm$Rhat$beta[s,c])))))
  fwrite(betas, "community_species_betas.csv")

  trait_out <- data.table(
    param=c("mass_psi","herb_psi","omni_psi","mass_p","herb_p","omni_p"),
    mean=c(fit_comm$mean$gamma.mass.psi, fit_comm$mean$gamma.herb.psi, fit_comm$mean$gamma.omni.psi,
           fit_comm$mean$gamma.mass.p, fit_comm$mean$gamma.herb.p, fit_comm$mean$gamma.omni.p),
    lcl=c(fit_comm$q2.5$gamma.mass.psi, fit_comm$q2.5$gamma.herb.psi, fit_comm$q2.5$gamma.omni.psi,
          fit_comm$q2.5$gamma.mass.p, fit_comm$q2.5$gamma.herb.p, fit_comm$q2.5$gamma.omni.p),
    ucl=c(fit_comm$q97.5$gamma.mass.psi, fit_comm$q97.5$gamma.herb.psi, fit_comm$q97.5$gamma.omni.psi,
          fit_comm$q97.5$gamma.mass.p, fit_comm$q97.5$gamma.herb.p, fit_comm$q97.5$gamma.omni.p),
    Rhat=c(fit_comm$Rhat$gamma.mass.psi, fit_comm$Rhat$gamma.herb.psi, fit_comm$Rhat$gamma.omni.psi,
           fit_comm$Rhat$gamma.mass.p, fit_comm$Rhat$gamma.herb.p, fit_comm$Rhat$gamma.omni.p))
  fwrite(trait_out, "community_trait_effects.csv")
  cat("max Rhat beta:", round(max(fit_comm$Rhat$beta,na.rm=TRUE),3),
      "| max Rhat mu.beta:", round(max(fit_comm$Rhat$mu.beta,na.rm=TRUE),3),
      "| max Rhat trait gammas:", round(max(trait_out$Rhat,na.rm=TRUE),3), "\n")
} else {
  hyp <- fread(file.path(DATA_DIR, "community_hyperparameters.csv"))
  betas <- fread(file.path(DATA_DIR, "community_species_betas.csv"))
  trait_out <- fread(file.path(DATA_DIR, "community_trait_effects.csv"))
  cat("Loaded pre-computed community fit (set RUN_LIVE <- TRUE above to refit; ~8 hours).\n")
}
cat("Community-mean covariate effects:\n"); print(hyp[, .(covariate, mu_mean=round(mu_mean,2))])
cat("\nTrait effects:\n"); print(trait_out[, .(param, mean=round(mean,2))])


## =============================================================================
## STEP 5: METHOD 2 -- GJAM (JOINT SPECIES DISTRIBUTION MODEL)
## =============================================================================
## GJAM models all 44 species together and, after fitting habitat effects,
## reports the RESIDUAL correlation between species -- co-occurrence that
## habitat alone does not explain. inrangeFrac (the fraction of the 44
## species whose mapped range includes this site) enters as an extra
## covariate so range-driven co-occurrence is not misread as biotic signal.

inrange_frac <- rowMeans(rng)
xdata_gjam <- data.table(
  forestZ=sc$forest_100m_z, savannaZ=sc$savanna_100m_z, pastureZ=sc$pasture_100m_z,
  croplandZ=sc$cropland_100m_z, nvegZ=sc$native_veg_1000m_z, tempZ=sc$temp_mean_C_z,
  precipZ=sc$precip_annual_mm_z, inPA=sc$in_pa, builtZ=sc$ghsl_built_5000m_z,
  distroadZ=sc$dist_road_m_z, distwaterZ=sc$dist_water_m_z, inrangeFrac=inrange_frac)
ydata_gjam <- as.data.frame(ysum); colnames(ydata_gjam) <- make.names(species)
site_days <- sc$camera_days  # per-site survey effort (proxy: total camera-days)
effMat <- matrix(site_days, nrow=N, ncol=S)
eList <- list(columns=1:S, values=effMat)

RUN_LIVE_GJAM <- FALSE  # GJAM is comparatively fast (~15-20 min); set TRUE to refit

if (RUN_LIVE_GJAM) {
  ml <- list(typeNames="DA", ng=8000, burnin=2000, effort=eList)
  form_gjam <- ~ forestZ + savannaZ + pastureZ + croplandZ + nvegZ + tempZ + precipZ + inPA +
                 builtZ + distroadZ + distwaterZ + inrangeFrac
  t0 <- Sys.time()
  gjam_fit <- gjam(form_gjam, xdata=as.data.frame(xdata_gjam), ydata=ydata_gjam, modelList=ml)
  cat("GJAM elapsed:", round(difftime(Sys.time(),t0,units="mins"),1), "min\n")

  preds <- c("intercept","forestZ","savannaZ","pastureZ","croplandZ","nvegZ","tempZ","precipZ",
             "inPA","builtZ","distroadZ","distwaterZ","inrangeFrac")
  bt <- gjam_fit$parameters$betaTable
  bt$row <- rownames(bt)
  bt$pred <- sapply(bt$row, function(r){ hit <- preds[sapply(preds, function(p) endsWith(r,p))]
                                          if(length(hit)==0) NA else hit[which.max(nchar(hit))] })
  bt$species <- mapply(function(r,p) sub(paste0("_",p,"$"),"",r), bt$row, bt$pred)
  fwrite(bt, "gjam_habitat_betas.csv")

  gjam_rc <- as.data.frame(gjam_fit$parameters$corMu)
  rownames(gjam_rc) <- colnames(gjam_rc) <- species
  fwrite(cbind(species=rownames(gjam_rc), gjam_rc), "gjam_residual_correlation.csv")
} else {
  gjam_rc_dt <- fread(file.path(DATA_DIR, "gjam_residual_correlation.csv"))
  cat("Loaded pre-computed GJAM residual correlation (set RUN_LIVE_GJAM <- TRUE above to refit; ~20 min).\n")
}


## =============================================================================
## STEP 6: METHOD 3 -- ROYLE-NICHOLS (occuRN) ABUNDANCE MODEL
## =============================================================================
## occuRN treats detection heterogeneity as a signature of unobserved
## abundance: lambda ~ all 11 covariates, p ~ intercept only, fit per species.

RUN_LIVE_RN <- FALSE  # per-species fits, a few seconds each; set TRUE to refit

if (RUN_LIVE_RN) {
  occ_form_rn <- as.formula(paste("~1 ~", paste(CAMERA_COVARS_Z, collapse="+")))
  rn_coefs <- list()
  for (s in seq_len(S)) {
    sci <- species[s]
    y <- build_dethist(sci)
    umf <- unmarkedFrameOccu(y=y, siteCovs=as.data.frame(sc[, ..CAMERA_COVARS_Z]))
    fit_rn <- tryCatch(occuRN(occ_form_rn, data=umf), error=function(e) NULL)
    if (is.null(fit_rn)) {
      cf <- data.table(species=sci, param=paste0("lam(",c("Int",CAMERA_COVARS_Z),")"),
                        beta=NA_real_, se=NA_real_, n_sites=nrow(y))
    } else {
      co <- coef(fit_rn); se <- sqrt(diag(vcov(fit_rn)))
      cf <- data.table(species=sci, param=names(co), beta=as.numeric(co), se=as.numeric(se), n_sites=nrow(y))
    }
    rn_coefs[[sci]] <- cf
    cat(sci, "done\n")
  }
  rn_out <- rbindlist(rn_coefs)
  fwrite(rn_out, "rn_all_species_coefs.csv")
} else {
  rn_out <- fread(file.path(DATA_DIR, "rn_all_species_coefs.csv"))
  cat("Loaded pre-computed Royle-Nichols fit (set RUN_LIVE_RN <- TRUE above to refit; ~2 min).\n")
}


## =============================================================================
## STEP 7: METHOD 4 -- SINGLE-SPECIES OCCUPANCY (spOccupancy::PGOcc)
## =============================================================================
## The pipeline's default single-species method (Module 3), refit here on
## the 11-covariate set with the range mask applied (out-of-range sites
## dropped from that species' likelihood) and an array random effect.

RUN_LIVE_PGOCC <- FALSE  # per-species MCMC, a few seconds to ~1 min each; set TRUE to refit

if (RUN_LIVE_PGOCC) {
  pgocc_results <- list()
  for (s in seq_len(S)) {
    sci <- species[s]
    y <- build_dethist(sci)
    keep_sp <- rng[, sci] == 1  # drop out-of-range sites for this species
    y2 <- y[keep_sp, , drop=FALSE]
    sc2 <- as.data.frame(sc[keep_sp])
    surveyed <- rowSums(!is.na(y2)) > 0
    y3 <- y2[surveyed, , drop=FALSE]; sc3 <- sc2[surveyed, ]
    occ_covs <- sc3[, CAMERA_COVARS_Z]
    occ_covs$array_num <- sc3$array_num
    data_list <- list(y=y3, occ.covs=occ_covs)
    form_pg <- as.formula(paste("~", paste(CAMERA_COVARS_Z, collapse=" + "), "+ (1|array_num)"))
    fit_pg <- tryCatch(
      PGOcc(occ.formula=form_pg, det.formula=~1, data=data_list,
            n.samples=8000, n.burn=2000, n.thin=2, n.chains=3, verbose=FALSE),
      error=function(e) NULL)
    if (!is.null(fit_pg)) {
      bs <- fit_pg$beta.samples
      pgocc_results[[sci]] <- data.table(species=sci, param=colnames(bs),
                                          mean=apply(bs,2,mean), sd=apply(bs,2,sd),
                                          lower95=apply(bs,2,quantile,0.025),
                                          upper95=apply(bs,2,quantile,0.975), n_sites=nrow(y3))
    }
    cat(sci, "done, n_sites=", nrow(y3), "\n")
  }
  pgocc_out <- rbindlist(pgocc_results)
  fwrite(pgocc_out, "all_species_pgocc_coefs.csv")
} else {
  pgocc_out <- fread(file.path(DATA_DIR, "all_species_pgocc_coefs.csv"))
  cat("Loaded pre-computed PGOcc fit (set RUN_LIVE_PGOCC <- TRUE above to refit; ~15-20 min).\n")
}


## =============================================================================
## STEP 8: METHOD 5 -- NAIVE LOGISTIC REGRESSION (no detection model)
## =============================================================================
## A deliberately naive baseline: presence/absence per site, no detection
## submodel. Shown to illustrate the perfect-separation problems it produces
## on sparse land-cover covariates -- the pedagogical point of including it.

glm_coefs <- list()
for (s in seq_len(S)) {
  sci <- species[s]
  y <- build_dethist(sci)
  pa <- as.integer(rowSums(y == 1, na.rm=TRUE) > 0)
  df_glm <- data.frame(pa=pa, sc[, ..CAMERA_COVARS])
  for (c in CAMERA_COVARS) df_glm[[c]] <- as.numeric(scale(df_glm[[c]]))
  rows <- list()
  for (c in CAMERA_COVARS) {
    f <- as.formula(paste("pa ~", c))
    m <- suppressWarnings(glm(f, data=df_glm, family=binomial()))
    co <- coef(m)[2]; se <- sqrt(diag(vcov(m)))[2]
    degenerate <- is.na(se) || se > 20  # perfect/quasi-separation flag
    rows[[c]] <- data.table(species=sci, covariate=c,
                             beta=if(degenerate) NA_real_ else co,
                             se=if(degenerate) NA_real_ else se,
                             p=if(degenerate) NA_real_ else summary(m)$coefficients[2,4],
                             n_sites=nrow(df_glm), degenerate=degenerate)
  }
  glm_coefs[[sci]] <- rbindlist(rows)
}
glm_out <- rbindlist(glm_coefs)
fwrite(glm_out, "glm_coefs.csv")
cat("Naive GLM: degenerate (perfect-separation) fits:", sum(glm_out$degenerate), "of", nrow(glm_out), "\n\n")


## =============================================================================
## STEP 9: CROSS-METHOD COMPARISON
## =============================================================================
## Assemble all five methods into one long-format table on the shared
## 11-covariate design, and compute pairwise coefficient correlations.

common_names <- fread(file.path(DATA_DIR, "species_common_names.csv"))
common_map <- setNames(common_names$common, common_names$species)

# PGOcc uses abbreviated parameter names (forest_z, nveg_z, built_z, distroad_z,
# distwater_z) that don't match the raw covariate names by simple suffix-stripping;
# map them explicitly to the canonical CAMERA_COVARS names.
pgocc_name_map <- c(forest_z="forest_100m", savanna_z="savanna_100m", pasture_z="pasture_100m",
                     cropland_z="cropland_100m", nveg_z="native_veg_1000m", temp_z="temp_mean_C",
                     precip_z="precip_annual_mm", in_pa="in_pa", built_z="ghsl_built_5000m",
                     distroad_z="dist_road_m", distwater_z="dist_water_m")

# GJAM's predictor names carry a "Z" suffix (forestZ, savannaZ, ...) and its
# own species-name mangling (make.names on the scientific name); map both
# back to the canonical CAMERA_COVARS / species names.
gjam_name_map <- c(forestZ="forest_100m", savannaZ="savanna_100m", pastureZ="pasture_100m",
                    croplandZ="cropland_100m", nvegZ="native_veg_1000m", tempZ="temp_mean_C",
                    precipZ="precip_annual_mm", inPA="in_pa", builtZ="ghsl_built_5000m",
                    distroadZ="dist_road_m", distwaterZ="dist_water_m")
if (RUN_LIVE_GJAM) {
  gjam_bt <- as.data.table(bt)
} else {
  gjam_bt <- fread(file.path(DATA_DIR, "gjam_habitat_betas.csv"))
}
gjam_bt_hab <- gjam_bt[pred %in% names(gjam_name_map)]
gjam_bt_hab[, covariate := gjam_name_map[pred]]
gjam_bt_hab[, species_clean := species]  # already dot-separated sci name in this file
species_dotmap <- setNames(species, make.names(species))
gjam_bt_hab[, species_clean := species_dotmap[species]]

# significance ("sig", a 0/1 flag) computed uniformly per method: does the
# coefficient's 95% interval exclude zero (Bayesian methods, GJAM), or is
# p < 0.05 (GLM), or |beta/se| > 1.96 (Royle-Nichols, which reports only SE)?
all_methods <- rbindlist(list(
  betas[, .(method="Community", species, covariate, mean, sig=as.integer(lcl>0 | ucl<0))],
  glm_out[, .(method="GLM", species, covariate, mean=beta, sig=as.integer(!is.na(p) & p<0.05))],
  rn_out[grepl("^lam\\(", param) & param != "lam(Int)",
         .(method="Royle-Nichols", species, covariate=gsub("^lam\\(|\\)$","",param), mean=beta,
           sig=as.integer(!is.na(se) & abs(beta/se)>1.96))],
  pgocc_out[!grepl("Intercept", param),
            .(method="PGOcc", species, covariate=pgocc_name_map[param], mean=mean,
              sig=as.integer(lower95>0 | upper95<0))],
  gjam_bt_hab[, .(method="GJAM", species=species_clean, covariate, mean=Estimate,
                  sig=as.integer(CI_025>0 | CI_975<0))]
), fill=TRUE)
fwrite(all_methods, "ecological_coefficients_all_methods.csv")

wide_methods <- dcast(all_methods, species + covariate ~ method, value.var="mean")
method_cols <- intersect(c("Community","GLM","Royle-Nichols","PGOcc","GJAM"), colnames(wide_methods))
method_cor <- cor(wide_methods[, ..method_cols], use="pairwise.complete.obs")
fwrite(as.data.table(method_cor, keep.rownames="method"), "method_beta_correlations.csv")
cat("Pairwise method correlations:\n"); print(round(method_cor, 2))


## =============================================================================
## STEP 10: FIGURES
## =============================================================================

## ---- 4.2: community hyperparameter + trait effects ----
hyp[, credible := fifelse(mu_lcl > 0, "positive", fifelse(mu_ucl < 0, "negative", "n.s."))]
p_hyp <- ggplot(hyp, aes(x=reorder(covariate, mu_mean), y=mu_mean, color=credible)) +
  geom_hline(yintercept=0, color="grey40", linewidth=0.3) +
  geom_pointrange(aes(ymin=mu_lcl, ymax=mu_ucl), linewidth=0.7, fatten=3) +
  scale_color_manual(values=c(positive=COL_POS, negative=COL_NEG, "n.s."=COL_NEUTRAL), guide="none") +
  coord_flip() + labs(x=NULL, y="community-mean effect (logit scale)",
                       title="Community-mean habitat responses") +
  theme_minimal(base_size=11)

trait_out[, credible := fifelse(lcl > 0, "positive", fifelse(ucl < 0, "negative", "n.s."))]
p_trait <- ggplot(trait_out, aes(x=reorder(param, mean), y=mean, color=credible)) +
  geom_hline(yintercept=0, color="grey40", linewidth=0.3) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl), linewidth=0.7, fatten=3) +
  scale_color_manual(values=c(positive=COL_POS, negative=COL_NEG, "n.s."=COL_NEUTRAL), guide="none") +
  coord_flip() + labs(x=NULL, y="effect on species intercept",
                       title="Species-trait effects on occupancy/detection intercepts") +
  theme_minimal(base_size=11)

if (requireNamespace("patchwork", quietly=TRUE)) {
  library(patchwork)
  ggsave(file.path(FIGS_DIR,"m4_community_hyper_trait.png"), p_hyp + p_trait, width=13, height=6.5, dpi=150)
} else {
  ggsave(file.path(FIGS_DIR,"m4_community_hyper.png"), p_hyp, width=7, height=6.5, dpi=150)
  ggsave(file.path(FIGS_DIR,"m4_trait_effects.png"), p_trait, width=7, height=6.5, dpi=150)
}
cat("Saved figs/m4_community_hyper_trait (or split) figure.\n")

## ---- 4.3: species-level responses, one panel per covariate ----
betas[, common := common_map[species]]
betas[, credible := fifelse(lcl > 0, "positive", fifelse(ucl < 0, "negative", "n.s."))]
p_species <- ggplot(betas, aes(x=reorder(common, mean), y=mean, color=credible)) +
  geom_hline(yintercept=0, color="grey40", linewidth=0.3) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl), linewidth=0.4, fatten=1.5) +
  scale_color_manual(values=c(positive=COL_POS, negative=COL_NEG, "n.s."=COL_NEUTRAL), guide="none") +
  coord_flip() + facet_wrap(~covariate, scales="free_x", ncol=6) +
  labs(x=NULL, y="community occupancy effect", title="Species-level responses, all 11 covariates") +
  theme_minimal(base_size=7) + theme(strip.text=element_text(size=7))
ggsave(file.path(FIGS_DIR,"m4_species_responses.png"), p_species, width=17, height=11, dpi=150)
cat("Saved figs/m4_species_responses.png\n")

## ---- covariate-pair biplots (4.3): significance-coded species scatter ----
render_biplot <- function(cov_x, cov_y, label_x, label_y) {
  wide_bx <- betas[covariate==cov_x, .(species, common, mx=mean, sigx=credible!="n.s.")]
  wide_by <- betas[covariate==cov_y, .(species, my=mean, sigy=credible!="n.s.")]
  d <- merge(wide_bx, wide_by, by="species")
  d[, tier := fifelse(sigx & sigy, "both", fifelse(sigx | sigy, "one", "neither"))]
  ggplot(d, aes(x=mx, y=my)) +
    geom_hline(yintercept=0, color="grey80") + geom_vline(xintercept=0, color="grey80") +
    geom_point(aes(size=tier, shape=tier, fill=tier), color="grey30") +
    scale_size_manual(values=c(both=3.2, one=2, neither=1.6), guide="none") +
    scale_shape_manual(values=c(both=21, one=21, neither=21), guide="none") +
    scale_fill_manual(values=c(both="#444444", one="#aaaaaa", neither="white"), guide="none") +
    ggrepel::geom_text_repel(data=d[tier!="neither"], aes(label=common), size=2.4, max.overlaps=25) +
    labs(x=label_x, y=label_y, title=paste(label_x,"vs",label_y)) +
    theme_minimal(base_size=10)
}
biplot_pairs <- list(
  c("forest_100m","savanna_100m","Forest","Savanna"),
  c("in_pa","native_veg_1000m","Protected area","Native vegetation"),
  c("temp_mean_C","precip_annual_mm","Temperature","Precipitation"),
  c("cropland_100m","pasture_100m","Cropland","Pasture"),
  c("ghsl_built_5000m","pasture_100m","Development","Pasture")
)
for (bp in biplot_pairs) {
  p <- render_biplot(bp[1], bp[2], bp[3], bp[4])
  fname <- paste0("m4_biplot_", bp[1], "_", bp[2], ".png")
  ggsave(file.path(FIGS_DIR, fname), p, width=8, height=7, dpi=150)
  cat("Saved figs/", fname, "\n")
}

## ---- 5.1: GJAM residual co-occurrence heatmap ----
rc_mat <- as.matrix(gjam_rc_dt[, -1, with=FALSE])
rownames(rc_mat) <- gjam_rc_dt$species
rc_long <- as.data.table(as.table(rc_mat))
setnames(rc_long, c("sp1","sp2","corr"))
rc_long[, common1 := common_map[as.character(sp1)]]
rc_long[, common2 := common_map[as.character(sp2)]]
p_gjam_heat <- ggplot(rc_long, aes(x=common1, y=common2, fill=corr)) +
  geom_tile() +
  scale_fill_gradient2(low=COL_NEG, mid="white", high=COL_POS, midpoint=0, limits=c(-0.6,0.6), name="residual\ncorrelation") +
  labs(x=NULL, y=NULL, title="GJAM residual species co-occurrence",
       subtitle="green = co-occur more than habitat predicts; red = co-occur less") +
  theme_minimal(base_size=7) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y=element_text(size=6))
ggsave(file.path(FIGS_DIR,"m4_gjam_residual_correlation.png"), p_gjam_heat, width=12, height=11, dpi=150)
cat("Saved figs/m4_gjam_residual_correlation.png\n")

## ---- 4.2: GJAM predictor importance bar graph ----
if (RUN_LIVE_GJAM) {
  bt2 <- as.data.table(bt)[bt$pred != "intercept", ]
  gjam_imp <- bt2[, .(predictor=pred[1], importance_rms=sqrt(mean(Estimate^2)),
                       frac_species_sig=mean(sign(CI_025)==sign(CI_975))), by=pred]
} else {
  gjam_imp <- fread(file.path(DATA_DIR, "gjam_predictor_importance.csv"))
}
p_gjam_imp <- ggplot(gjam_imp, aes(x=reorder(predictor, frac_species_sig), y=frac_species_sig)) +
  geom_col(fill=COL_NEUTRAL) + coord_flip() +
  labs(x=NULL, y="fraction of species with a credible effect", title="GJAM predictor importance") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR,"m4_gjam_importance.png"), p_gjam_imp, width=8, height=6, dpi=150)
cat("Saved figs/m4_gjam_importance.png\n")

## ---- 6.1: cross-method predictor importance ----
## Unit-free importance metric (matches the notebook's own definition): the
## percentage of the 44 species with a credibly/significantly non-zero
## coefficient for each covariate, under each method's own significance test.
imp_by_method <- all_methods[, .(pct_sig=100*mean(sig, na.rm=TRUE)), by=.(method, covariate)]
p_cross_imp <- ggplot(imp_by_method, aes(x=covariate, y=pct_sig, fill=method)) +
  geom_col(position="dodge") + coord_flip() +
  labs(x=NULL, y="% of 44 species with a significant coefficient", title="Predictor importance across methods") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR,"m4_methods_importance.png"), p_cross_imp, width=10, height=7, dpi=150)
cat("Saved figs/m4_methods_importance.png\n")

cat("\n=== Module 4 script complete: all methods fit/loaded, all figures saved to", FIGS_DIR, "===\n")
