## =============================================================================
## Módulo 5 — Predição e Mapeamento
## Pipeline multi-rede de armadilhas fotográficas do Snapshot Brasil
## =============================================================================
##
## AUTOCONTIDO: este script transforma o modelo de ocupação ajustado de cada espécie
## (do Módulo 3/4) em uma predição nacional em uma grade de 17,766 células,
## faz a triagem de quais espécies vale a pena mapear (AUC), sinaliza onde o modelo está
## extrapolando (MESS), aplica a máscara de distribuição (IUCN) e transforma a
## superfície contínua de ocupação em um mapa de três resultados presente/ausente/
## não avaliado, tanto no nível de câmera (11 covariáveis) quanto no nível de array
## (6 covariáveis).
##
## Arquivos de entrada necessários (em DATA_DIR, padrão "data" junto a este script):
## national_grid_covariates_11cov.csv  -- grade nacional de 17,766 células, as mesmas 11
## covariáveis de nível de câmera do final_covariates.csv,
## extraídas das mesmas fontes/definições do Earth Engine,
## além de MESS/extrapolação
## all_species_pgocc_coefs.csv         -- coeficientes PGOcc de nível de câmera de 44 espécies
## (ajustados no Módulo 4, reaproveitados aqui como estão)
## array_covariates.csv                -- tabela de covariáveis de 60 arrays (as mesmas 11 covariáveis
## do arquivo de nível de câmera; o modelo de nível de array abaixo
## as reduz ao conjunto reduzido de 6 covariáveis
## ARRAY_SET -- floresta, savana, vegetação nativa,
## temperatura, precipitação, distância até estrada --
## via ARRAY_COVARS, já que a resolução espacial mais grosseira do array
## não sustenta todas as 11 de forma confiável)
## ysum_array.csv / K_array.csv         -- matrizes de contagem de detecção de nível de array 60 x 44
## (colapsadas por ocasião; saídas de referência do
## montador de históricos de detecção de nível de array)
## array_ids.csv                       -- os 60 IDs de array, na ordem das linhas de ysum/K
## species_list.json                   -- as 44 espécies/gêneros modelados
## grid_inrange_44sp.csv                -- máscara de distribuição (IUCN) por célula (44 espécies x 17,766
## células da grade), saída do exercício de mapa de distribuição do Módulo 2.2
## aplicado às coordenadas da grade em vez dos sítios de câmera
## final_covariates.csv                 -- covariáveis dos sítios de câmera, usadas para recuperar
## as constantes de padronização em escala de treinamento E (com
## os dois arquivos abaixo) para calcular AUC/limiares
## final_deployments.csv, final_detections_mammals.csv -- tabelas brutas de implantação/detecção,
## necessárias para reconstruir a presença/ausência em nível de sítio para
## o cálculo de AUC e limiar
## species_range_mask.csv                -- máscara de distribuição (IUCN) sítio x espécie (nível de câmera;
## 1 = candidato, 0 = zero estrutural), usada para
## restringir o ajuste de AUC/limiar aos sítios dentro da distribuição
##
## As matrizes de detecção de nível de array e a máscara de distribuição da grade são saídas de referência dos
## próprios montadores de módulos anteriores e não podem ser reconstruídas apenas a partir dos três
## arquivos principais do conjunto de dados combinado -- elas são fornecidas junto com este script como tabelas de consulta fixas.
##
## Requer: data.table, unmarked, jsonlite, ggplot2.

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
## PASSO 1: CARREGAR A GRADE E OS COEFICIENTES DE NÍVEL DE CÂMERA
## =============================================================================
grid <- fread(file.path(DATA_DIR, "national_grid_covariates_11cov.csv"))
cat("Grid cells:", nrow(grid), "\n")
cat("Grid cells extrapolating (11-covariate MESS < 0):", sum(grid$extrapolating, na.rm=TRUE),
    sprintf("(%.1f%% of the grid)\n", 100*mean(grid$extrapolating, na.rm=TRUE)))

pgocc <- fread(file.path(DATA_DIR, "all_species_pgocc_coefs.csv"))
species <- fromJSON(file.path(DATA_DIR, "species_list.json"))
fc <- fread(file.path(DATA_DIR, "final_covariates.csv"))
fc_complete <- fc[complete.cases(fc[, ..CAMERA_COVARS])]

## Recuperar as constantes de padronização de escala do treinamento (média/dp) exatamente como
## usadas quando os modelos PGOcc foram ajustados, para que a grade seja padronizada na
## MESMA escala dos coeficientes, e não na média/dp da própria grade.
znames <- c("forest_z","savanna_z","pasture_z","cropland_z","nveg_z","temp_z",
            "precip_z","in_pa","built_z","distroad_z","distwater_z")
site_scaler <- list()
for (i in seq_along(CAMERA_COVARS)) {
  zc <- znames[i]; raw <- CAMERA_COVARS[i]
  if (zc == "in_pa") next
  site_scaler[[zc]] <- list(raw = raw, mean = mean(fc_complete[[raw]]), sd = sd(fc_complete[[raw]]))
}

## =============================================================================
## ETAPA 2: PROJETAR PSI NO NÍVEL DE CÂMERA EM TODA A GRADE
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
## PASSO 3: REAJUSTAR O MODELO DE OCUPAÇÃO DE ESPÉCIE ÚNICA EM NÍVEL DE ARRAY (CONJUNTO DE 6 COVARIÁVEIS)
## =============================================================================
## O modelo em nível de array é reajustado do zero aqui (diferente do modelo
## em nível de câmera, que reutiliza os coeficientes salvos do Module 4) porque ele nunca foi
## salvo em um módulo anterior -- Module 3/4 ajustou apenas ocupação de espécie
## única em nível de câmera.
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

## Projetar o psi de nível de array pela grade para a espécie projetável
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
## PASSO 4: TRIAGEM DE AUC, LIMIAR, MÁSCARA DE DISTRIBUIÇÃO (IUCN) E O MAPA DE STATUS DE TRÊS RESULTADOS
## =============================================================================

## --- 4.1: reconstruir os históricos de detecção em nível de site para calcular a AUC de treinamento ---
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

## --- 4.2: AUC empírico por espécie + limiar de Youden (MaxSens+Spec) ---
roc_auc <- function(y, p) {
  # Formulação de Mann-Whitney U: AUC = (soma dos ranks da classe positiva - n1*(n1+1)/2) / (n1*n0).
  # y e p devem permanecer na MESMA ordem (original) durante todo o processo -- não permute um sem o outro.
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
  target_sci <- species[i]  # # nomeado de forma distinta da própria coluna "sci" de det_m para evitar
                             # um bug de escopo do data.table em que `sci == sci` compararia silenciosamente
                             # a coluna com ela mesma (sempre TRUE) em vez de com
                             # esta variável do loop
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

## --- 4.3: máscara de distribuição (IUCN) + MESS -> status de três resultados por célula da grade ---
grid_inrange <- fread(file.path(DATA_DIR, "grid_inrange_44sp.csv"))
## grid_inrange só tem uma linha para células com um valor de máscara de distribuição (IUCN) resolvível (menos do que
## a grade completa, por exemplo, células fora da extensão avaliada de todas as espécies) -- reindexe-o para
## a própria ordem de sid da grade, para que uma coluna de espécie se alinhe célula a célula, tratando qualquer
## célula da grade ausente de grid_inrange como NÃO dentro da distribuição (conservador: não avaliada em vez
## de silenciosamente assumida como presente).
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
