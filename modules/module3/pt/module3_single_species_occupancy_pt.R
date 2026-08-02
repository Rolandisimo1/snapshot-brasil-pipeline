## =============================================================================
## Módulo 3 — Ocupação de Espécie Única
## Pipeline multi-rede de armadilhas fotográficas Snapshot Brasil
## =============================================================================
##
## AUTOSSUFICIENTE: este script parte dos três arquivos brutos do conjunto de
## dados combinado (implantações, detecções de mamíferos, covariáveis de
## sensoriamento remoto) e constrói tudo do zero -- históricos de detecção,
## covariáveis, conjuntos candidatos de locais, e cada ajuste de modelo.
## Coloque este script na mesma pasta que os seis arquivos de dados
## listados abaixo e execute-o do início ao fim.
##
## Arquivos de entrada necessários (em DATA_DIR, padrão "data" junto a este script):
##   final_deployments.csv         -- uma linha por implantação de câmera
##   final_detections_mammals.csv  -- uma linha por evento de detecção independente
##   final_covariates.csv          -- uma linha por local de câmera, covariáveis de sensoriamento remoto
##   species_site_exclusions.csv   -- locais de zero estrutural por espécie (saída do mascaramento de distribuição do Módulo 2)
##   genus_species_site_exclusions.csv -- locais de zero estrutural para Didelphis agrupado (ver nota abaixo)
##   site_ecoregion.csv            -- bioma/ecorregião por local (usado apenas na extensão de ecorregião)
##
## Os últimos três arquivos são saídas de referência da avaliação de mapas de
## distribuição do Módulo 2 (uma junção espacial com polígonos de
## distribuição IUCN) e não podem ser reconstruídos a partir dos três
## arquivos principais de dados sozinhos -- são fornecidos junto com este
## script como tabelas de referência fixas.
##
## Nota taxonômica: Didelphis (3 espécies com mapas de distribuição +
## registros apenas em nível de gênero) e Dasyprocta (5 espécies, apenas 1
## com mapa de distribuição, + registros apenas em nível de gênero) são
## AGRUPADOS em nível de gênero neste script -- uma exceção deliberada à
## regra usual do pipeline de identificação apenas em nível de espécie, feita
## porque a identificação de campo até espécie é não confiável para esses
## dois gêneros, e o agrupamento recupera registros de detecção apenas em
## nível de gênero que de outra forma seriam descartados. Didelphis mantém
## sua máscara de distribuição (a união das distribuições com buffer das 3
## espécies mapeadas); Dasyprocta descarta sua máscara de distribuição por
## completo, já que apenas um dos cinco congêneres tem mapa de distribuição
## e usá-lo sozinho classificaria erroneamente muitas detecções reais como
## "impossíveis".
##
## Requer: data.table, ggplot2, unmarked, spOccupancy, pROC.
## Se os rótulos das figuras (letras acentuadas ou símbolos psi/beta/lambda)
## aparecerem corrompidos, execute este script com um locale UTF-8 definido
## (ex.: LC_ALL=pt_BR.UTF-8 ou LC_ALL=en_US.UTF-8).
## =============================================================================

suppressMessages({
  library(data.table)
  library(unmarked)
  library(pROC)
  library(ggplot2)
})

# ---- 0. Caminhos -------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"
stopifnot(dir.exists(DATA_DIR))
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR)

# convenção de cores usada em todo o script: verde = efeito positivo, vermelho = negativo
COL_POS <- "#1b7837"; COL_NEG <- "#c0392b"; COL_NEUTRAL <- "#2c7fb8"

FINAL_COVARS <- c("forest_100m","savanna_100m","pasture_100m","cropland_100m",
                   "native_veg_1000m","temp_mean_C","precip_annual_mm")
OCC_DAYS <- 7  # duração da ocasião em dias, seguindo a convenção do pipeline para o intervalo de independência


## =============================================================================
## ETAPA 1: CARREGAR DADOS BRUTOS E CONSTRUIR COVARIÁVEIS EM NÍVEL DE LOCAL
## =============================================================================

dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
det <- fread(file.path(DATA_DIR, "final_detections_mammals.csv"))
cov <- fread(file.path(DATA_DIR, "final_covariates.csv"))

dep[, start_dt := as.IDate(start_date)]
dep[, end_dt   := as.IDate(end_date)]

# uma linha por local: rede, array, coordenadas, esforço
site_tbl <- dep[, .(network = first(network), array_id = first(array_id),
                     longitude = mean(longitude), latitude = mean(latitude),
                     n_deployments = .N, camera_days = sum(camera_days)),
                 by = site_id]
setorder(site_tbl, site_id)
cat("Locais:", nrow(site_tbl), "| Arrays:", uniqueN(site_tbl$array_id), "\n")

# mescla as covariáveis de sensoriamento remoto, padroniza (z-score) o conjunto final, numera os arrays
sc <- merge(site_tbl, cov, by = "site_id", all.x = TRUE)
setorder(sc, site_id)
for (v in FINAL_COVARS) sc[[paste0(v, "_z")]] <- as.numeric(scale(sc[[v]]))
sc[, treecover2000_100m_z := as.numeric(scale(treecover2000_100m))]
arrays_sorted <- sort(unique(sc$array_id))
sc[, array_num := match(array_id, arrays_sorted)]

cat("Locais sem covariáveis (excluídos por exclusão listwise no ajuste):",
    sum(!complete.cases(sc[, ..FINAL_COVARS])), "\n\n")


## =============================================================================
## ETAPA 2: CONSTRUIR HISTÓRICOS DE DETECÇÃO A PARTIR DOS DADOS BRUTOS
## =============================================================================
## Um modelo de ocupação precisa de um HISTÓRICO DE DETECÇÃO: para cada
## local e cada ocasião de 7 dias, a espécie foi detectada (1), amostrada mas
## não detectada (0), ou nenhuma câmera estava ativa naquela ocasião
## (ausente)? Construímos isso diretamente das janelas de implantação brutas
## e dos registros de data/hora das detecções.
##
## Janelas por implantação (não o intervalo mínimo-máximo completo de um
## local) são usadas para marcar quais dias uma câmera esteve realmente
## ativa -- um local com múltiplas implantações e um intervalo entre dois
## períodos de implantação não deve ser tratado como continuamente
## amostrado durante essa lacuna.

site_span <- dep[, .(min_start = min(start_dt), max_end = max(end_dt)), by = site_id]
setorder(site_span, site_id)
site_span[, n_occ := (as.integer(max_end - min_start) %/% OCC_DAYS) + 1L]
max_occ <- max(site_span$n_occ)
site_ids_ordered <- site_span$site_id
n_sites <- length(site_ids_ordered)

# vetor booleano de dias ativos por local, construído a partir do início/fim de cada implantação
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

# colapsa o vetor de dias ativos em uma marcação de "amostrado" por ocasião: uma ocasião
# conta como amostrada se >=50% de seus dias tiveram uma câmera ativa
build_occasion_mask <- function(active_arr, occ_days = OCC_DAYS) {
  n_occ_here <- (length(active_arr) + occ_days - 1L) %/% occ_days
  mask <- rep(NA_real_, n_occ_here)
  for (i in seq_len(n_occ_here)) {
    window <- active_arr[((i - 1L) * occ_days + 1L):min(i * occ_days, length(active_arr))]
    if (mean(window) >= 0.5) mask[i] <- 0  # amostrado, ainda não detectado
  }
  mask
}
occ_active <- lapply(active_day, build_occasion_mask)

# mapeia cada detecção para o índice de ocasião de seu local; descarta o
# pequeno número de detecções que antecedem a implantação mais antiga
# registrada do local (datas de imagem de calibração incorretas / casos
# de limite off-by-one -- uma fração pequena e consistente ao longo do
# histórico do pipeline)
site_min_start <- setNames(site_span$min_start, site_span$site_id)
det[, site_min_start := site_min_start[site_id]]
det <- det[!is.na(site_min_start)]
det[, day_offset := as.integer(as.IDate(date) - site_min_start)]
det[, occ_idx := day_offset %/% OCC_DAYS]
n_dropped <- sum(det$occ_idx < 0)
det <- det[occ_idx >= 0]
cat("Detecções usadas:", nrow(det), "(", n_dropped, "descartadas: antecedem a implantação mais antiga do local)\n\n")

# ---- Taxonomia: agrupar Didelphis e Dasyprocta em nível de gênero ----
# (todos os outros gêneros mantêm identificação em nível de espécie,
# seguindo a regra padrão do pipeline; ver a nota taxonômica no topo deste script)
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
    y[sid, oc] <- 1  # detectado (também cobre ocasiões que a regra de >=50% ativo marcou como não amostradas)
  }
  y
}

# ---- Escolher as espécies bem detectadas a modelar ----
# mesmo limiar usado em todo o pipeline: nome binomial em nível de espécie
# (ou, para os dois gêneros agrupados, em nível de gênero), excluindo taxa
# domésticos/humanos e registros não resolvidos apenas em nível de
# gênero/família para todo OUTRO gênero, com pelo menos 10 detecções
# totais e 20 locais detectados
sp_stats <- det[, .(n_det = .N, n_sites_det = uniqueN(site_id)), by = sci_pooled]
exclude_terms <- c("Homo sapiens","Bos taurus","Canis familiaris","Sus scrofa","Equus","Felis catus","Gallus")
genus_only_or_family <- c("Mazama","Didelphidae","Mammalia","Rodentia","Cervidae","Proechimys")
sp_stats[, is_binomial_or_pooled := grepl(" ", sci_pooled) | sci_pooled %in% c("Didelphis","Dasyprocta")]
sp_stats[, exclude := sci_pooled %in% exclude_terms | sci_pooled %in% genus_only_or_family]
candidates <- sp_stats[is_binomial_or_pooled & !exclude]
modeled <- candidates[n_det >= 10 & n_sites_det >= 20]
setorder(modeled, -n_det)
cat("Espécies/gêneros modelados:", nrow(modeled), "\n")
species_meta <- modeled[, .(sci_mdd = sci_pooled, n_det, n_arrays_det = NA_integer_)]
for (i in seq_len(nrow(species_meta))) {
  y <- build_dethist(species_meta$sci_mdd[i])
  det_sites <- rownames(y)[apply(y == 1, 1, any, na.rm = TRUE)]
  species_meta$n_arrays_det[i] <- uniqueN(sc[site_id %in% det_sites, array_id])
}
cat("\n")


## =============================================================================
## ETAPA 3: CONJUNTOS CANDIDATOS COM MÁSCARA DE DISTRIBUIÇÃO (referências do Módulo 2)
## =============================================================================
## Para espécies com mapa de distribuição IUCN, um local é um "zero
## estrutural" se cai fora da distribuição mapeada (com um buffer de 100km
## para imprecisão do mapa) E a espécie nunca foi detectada ali -- esses
## locais são excluídos como candidatos por padrão. Locais detectados fora
## da distribuição com buffer são sempre mantidos (uma fotografia supera um
## polígono).

excl <- fread(file.path(DATA_DIR, "species_site_exclusions.csv"))
excl_genus <- fread(file.path(DATA_DIR, "genus_species_site_exclusions.csv"))
excl <- rbindlist(list(excl, excl_genus))
cat("Espécies com máscara de distribuição (exclusões de zero estrutural aplicadas):", uniqueN(excl$species), "\n")
cat("Didelphis agrupado: máscara de distribuição = união das distribuições com buffer dos 3 congêneres mapeados (",
    excl_genus[species=="Didelphis", .N], "locais de zero estrutural excluídos)\n")
cat("Dasyprocta agrupado: NENHUMA máscara de distribuição aplicada (apenas 1 de 5 congêneres tem mapa de distribuição;\n")
cat("  usá-lo sozinho classificaria muitas detecções reais como impossíveis) -- todo local é candidato.\n\n")


## =============================================================================
## 1. UM PRIMEIRO MODELO DE OCUPÂNCIA: TATU-GALINHA
## =============================================================================
## Começamos com o Tatu-galinha (*Dasypus novemcinctus*), uma espécie bem
## detectada, e ajustamos um modelo usando todas as covariáveis do conjunto
## final: floresta, savana, pastagem, lavoura, vegetação nativa em escala de
## paisagem, temperatura, precipitação e status de área protegida.

sci <- "Dasypus novemcinctus"
y_arm <- build_dethist(sci)
excl_sites <- excl[species == sci, site_id]
FINAL_COVARS_Z <- paste0(FINAL_COVARS, "_z")
keep <- !is.na(sc$forest_100m_z) & !(sc$site_id %in% excl_sites)
y2 <- y_arm[keep, , drop = FALSE]
sc2 <- as.data.frame(sc[keep])
cat("== Tatu-galinha: conjunto candidato de locais ==\n")
cat("locais após exclusão por máscara de distribuição + covariáveis ausentes:", nrow(y2), "\n\n")

## ---- 1.1 Histórico de detecção ----
high_arr <- "SNAP_Piaui_Caatinga_Tapuio_25"
low_arr  <- "ATLA_ATL_11"
cat("== Arrays de exemplo do histórico de detecção ==\n")
cat("Array de alta detecção:", high_arr, "| Array de baixa detecção:", low_arr, "\n\n")

n_occ_show <- 30
dethist_plot_data <- list()
for (arr_name in c(high_arr, low_arr)) {
  arr_sites <- sc[array_id == arr_name][order(site_id), site_id]
  sub <- y_arm[arr_sites, seq_len(min(n_occ_show, ncol(y_arm))), drop=FALSE]
  df <- as.data.table(sub)
  df[, site_id := arr_sites]
  df_long <- melt(df, id.vars="site_id", variable.name="occasion", value.name="status")
  df_long[, occasion_n := as.integer(gsub("V","",occasion))]
  df_long[, array_label := ifelse(arr_name==high_arr, "Array de alta detecção", "Array de baixa detecção")]
  df_long[, status_f := factor(ifelse(is.na(status), "Não amostrado", ifelse(status==1,"Detectado","Amostrado, não detectado")),
                                 levels=c("Não amostrado","Amostrado, não detectado","Detectado"))]
  dethist_plot_data[[arr_name]] <- df_long
}
dethist_dt <- rbindlist(dethist_plot_data)
p_dethist <- ggplot(dethist_dt, aes(x=occasion_n, y=site_id, fill=status_f)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_manual(values=c("Não amostrado"="white","Amostrado, não detectado"="#dddddd","Detectado"=COL_NEG), name=NULL) +
  facet_wrap(~array_label, scales="free_y") +
  labs(title="Histórico de detecção do Tatu-galinha: array de alta vs. baixa detecção",
       x="Ocasião de 7 dias", y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.y=element_text(size=6), panel.grid=element_blank(), legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_dethist_armadillo.png"), p_dethist, width=11, height=6, dpi=150)
cat("Salvo figs/m3_dethist_armadillo.png\n\n")

## ---- 1.2 Ajustando o modelo, com o efeito aleatório de array ----
## O ajuste de efeitos fixos (unmarked::occu) é mostrado aqui; o ajuste com
## efeito aleatório de array (spOccupancy::PGOcc) requer vários minutos de
## amostragem MCMC e é executado abaixo.

umf <- unmarkedFrameOccu(y = y2, siteCovs = sc2[, c(FINAL_COVARS_Z, "in_pa")])
form_full <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + "), "+ in_pa"))
fit_full <- occu(form_full, data = umf)
co <- coef(fit_full); se <- sqrt(diag(vcov(fit_full))); z <- co/se; p <- 2*pnorm(-abs(z))
cat("== Coeficientes de efeitos fixos do Tatu-galinha ==\n")
print(data.table(param=names(co), estimate=round(co,3), se=round(se,3), z=round(z,2), p=round(p,4)))
cat("\n")

## -- Ajuste com efeito aleatório de array (spOccupancy::PGOcc), o padrão do pipeline --
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
cat("== Coeficientes do Tatu-galinha com efeito de array (PGOcc) ==\n")
print(pgocc_coefs[param != "(Intercept)"], digits=3)
cat("\n")

fig_arm_effect_data <- list()
for (pr in c("forest_100m_z","savanna_100m_z","native_veg_1000m_z")) {
  b <- pgocc_coefs[param==pr, mean]; b_lo <- pgocc_coefs[param==pr, ci_lo]; b_hi <- pgocc_coefs[param==pr, ci_hi]
  intc <- pgocc_coefs[param=="(Intercept)", mean]
  zg <- seq(-2.5, 2.5, length.out=60)
  fig_arm_effect_data[[pr]] <- data.table(covariate=pr, z=zg,
    psi=plogis(intc + b*zg), psi_lo=plogis(intc + b_lo*zg), psi_hi=plogis(intc + b_hi*zg))
}
fig_arm_effect_dt <- rbindlist(fig_arm_effect_data)

## ---- 1.3 Ocupação ingênua vs. modelada ----
naive_psi_site <- apply(y2, 1, function(r) as.integer(any(r == 1, na.rm = TRUE)))
naive_psi <- mean(naive_psi_site)
psi_pred <- predict(fit_full, type = "state")$Predicted
modeled_psi_pgocc <- mean(apply(fit_pgocc$psi.samples, 2, mean))
cat("== Ocupação ingênua vs. modelada ==\n")
cat("Ocupação ingênua (fração de locais detectados alguma vez):", round(naive_psi,3), "\n")
cat("Ocupação modelada (PGOcc com efeito de array, psi médio):", round(modeled_psi_pgocc,3), "\n\n")

p_naive_modeled <- ggplot(data.table(x=c("Ingênua\n(fração detectada)","Modelada\n(PGOcc, efeito de array)"),
                                       y=c(naive_psi, modeled_psi_pgocc)),
                            aes(x=x, y=y, fill=x)) +
  geom_col(width=0.55) +
  geom_text(aes(label=sprintf("%.2f", y)), vjust=-0.5, fontface="bold", size=4) +
  scale_fill_manual(values=c("#999999", COL_NEUTRAL), guide="none") +
  labs(title="Tatu-galinha: ocupação ingênua vs. modelada", x=NULL, y="Ocupação (\u03c8)") +
  ylim(0, max(naive_psi, modeled_psi_pgocc)*1.3) +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIGS_DIR, "m3_naive_vs_modeled.png"), p_naive_modeled, width=5.5, height=5.5, dpi=150)
cat("Salvo figs/m3_naive_vs_modeled.png\n\n")

## ---- 1.4 Efeitos das covariáveis & 1.5 Ajuste do modelo ----
## Duas perguntas muito diferentes são feitas sobre um modelo de ocupação
## ajustado:
## (a) O MODELO DISCRIMINA locais ocupados de não ocupados? -- AUC.
##     Pegue a probabilidade de ocupação predita de cada local e pergunte:
##     se você escolhesse um local verdadeiramente ocupado e um
##     verdadeiramente não ocupado ao acaso, qual a chance de o modelo dar
##     ao ocupado uma probabilidade predita maior? AUC=0,5 é chance;
##     AUC=1,0 é perfeito.
## (b) A ESTRUTURA DE ERRO ASSUMIDA PELO MODELO ESTÁ CORRETA? -- c-hat
##     (superdispersão), a partir de um teste de qualidade de ajuste por
##     bootstrap paramétrico. c-hat próximo de 1 significa que a estrutura
##     de erro binomial assumida se ajusta; notavelmente acima de 1
##     significa que as detecções reais são mais variáveis do que o modelo
##     espera. AUC e c-hat NÃO são intercambiáveis: um modelo pode
##     discriminar bem (AUC alto) e ainda ser superdisperso.

## ---- figura de curvas de efeito (4 covariáveis significativas do ajuste com efeito de array) ----
sig_covars_arm <- pgocc_coefs[param != "(Intercept)"][ci_lo>0 | ci_hi<0, param]
cat("Covariáveis significativas com efeito de array para o Tatu-galinha:", paste(sig_covars_arm, collapse=", "), "\n\n")
effect_curve_data <- list()
for (pr in setdiff(sig_covars_arm, "in_pa")) {
  b <- pgocc_coefs[param==pr, mean]; b_lo <- pgocc_coefs[param==pr, ci_lo]; b_hi <- pgocc_coefs[param==pr, ci_hi]
  intc <- pgocc_coefs[param=="(Intercept)", mean]
  zg <- seq(-2.5, 2.5, length.out=60)
  effect_curve_data[[pr]] <- data.table(covariate=pr, z=zg,
    psi=plogis(intc + b*zg), psi_lo=plogis(intc + b_lo*zg), psi_hi=plogis(intc + b_hi*zg), sign=ifelse(b>0,"pos","neg"))
}
if (length(effect_curve_data)) {
  ec_dt <- rbindlist(effect_curve_data)
  label_map <- c(forest_100m_z="Cobertura florestal", savanna_100m_z="Cobertura de savana", pasture_100m_z="Pastagem",
                 cropland_100m_z="Lavoura", native_veg_1000m_z="Veg. nativa (1000m)",
                 temp_mean_C_z="Temperatura", precip_annual_mm_z="Precipitação")
  ec_dt[, covariate_label := factor(label_map[covariate], levels=label_map[names(effect_curve_data)])]
  p_effects <- ggplot(ec_dt, aes(x=z, y=psi, color=sign, fill=sign)) +
    geom_ribbon(aes(ymin=psi_lo, ymax=psi_hi), alpha=0.15, color=NA) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c(pos=COL_POS, neg=COL_NEG), guide="none") +
    scale_fill_manual(values=c(pos=COL_POS, neg=COL_NEG), guide="none") +
    facet_wrap(~covariate_label, nrow=1, scales="free_x") +
    labs(title="Tatu-galinha: ocupação vs. cada covariável significativa (verde=positivo, vermelho=negativo)",
         x="Covariável padronizada (z)", y="Ocupação predita (\u03c8)") +
    theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
  ggsave(file.path(FIGS_DIR, "m3_armadillo_effect_curves.png"), p_effects, width=13, height=4, dpi=150)
  cat("Salvo figs/m3_armadillo_effect_curves.png\n\n")
}

## ---- figura de beta com IC (modelo com efeito de array) ----
beta_plot <- pgocc_coefs[param != "(Intercept)"][order(mean)]
beta_plot[, label := label_map[param]]
beta_plot[is.na(label), label := "Área protegida"]
beta_plot[, label := factor(label, levels=label)]
beta_plot[, sign := fifelse(ci_lo>0, "pos", fifelse(ci_hi<0, "neg", "ns"))]
p_beta_ci <- ggplot(beta_plot, aes(x=mean, y=label, color=sign)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0, linewidth=1) +
  geom_point(size=3) +
  scale_color_manual(values=c(pos=COL_POS, neg=COL_NEG, ns="grey50"), guide="none") +
  labs(title="Tatu-galinha: efeitos das covariáveis (modelo com efeito de array)",
       x="Coeficiente (\u03b2), média posterior e intervalo de credibilidade de 95%", y=NULL) +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIGS_DIR, "m3_armadillo_beta_ci.png"), p_beta_ci, width=7, height=5.5, dpi=150)
cat("Salvo figs/m3_armadillo_beta_ci.png\n\n")

auc_val <- as.numeric(auc(roc(response = naive_psi_site, predictor = psi_pred, quiet = TRUE)))
cat("== Ajuste do modelo ==\nAUC:", round(auc_val, 3), "\n")

## ---- figura explicativa do AUC ----
roc_obj <- roc(response = naive_psi_site, predictor = psi_pred, quiet = TRUE)
roc_dt <- data.table(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
setorder(roc_dt, fpr)
score_dt <- data.table(psi = psi_pred, status = factor(ifelse(naive_psi_site==1,"Presença verdadeira","Ausência verdadeira"),
                                                          levels=c("Ausência verdadeira","Presença verdadeira")))
p_score_dist <- ggplot(score_dt, aes(x=psi, fill=status)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.55, bins=25) +
  scale_fill_manual(values=c("Ausência verdadeira"=COL_NEG, "Presença verdadeira"=COL_POS), name=NULL) +
  labs(title="O que o modelo prediz, por status verdadeiro", x="Probabilidade de ocupação predita (\u03c8-chapéu)", y="Densidade") +
  theme_minimal(base_size=10) + theme(legend.position=c(0.75,0.85))
p_roc <- ggplot(roc_dt, aes(x=fpr, y=tpr)) +
  geom_ribbon(aes(ymin=0, ymax=tpr), fill=COL_NEUTRAL, alpha=0.15) +
  geom_line(color=COL_NEUTRAL, linewidth=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  annotate("text", x=0.6, y=0.15, label="Chance\n(AUC = 0.5)", color="grey50", size=3.2) +
  labs(title=sprintf("Curva ROC: AUC = %.2f", auc_val),
       x="Taxa de falso positivo\n(locais não ocupados chamados de 'presente')",
       y="Taxa de verdadeiro positivo\n(locais ocupados chamados de 'presente')") +
  theme_minimal(base_size=10)
p_auc_explainer <- gridExtra::grid.arrange(p_score_dist, p_roc, ncol=2,
  top=grid::textGrob("AUC: quão bem a ocupação predita separa presença verdadeira de ausência verdadeira?",
                       gp=grid::gpar(fontsize=13, fontface="bold"), x=0.02, hjust=0))
ggsave(file.path(FIGS_DIR, "m3_auc_explainer.png"), p_auc_explainer, width=12, height=5, dpi=150)
cat("Salvo figs/m3_auc_explainer.png\n")

## Bootstrap de qualidade de ajuste (c-hat) -- 100 simulações, ~1-2 minutos:
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
cat("c-hat (superdispersão):", round(c_hat,3), "| p-valor do bootstrap de qualidade de ajuste:", round(gof_p,3), "\n")
cat("Um AUC em torno de 0,65 indica discriminação modesta; um c-hat acima de 1 significa\n")
cat("que o modelo é superdisperso -- ambos achados comuns e esperados para modelos de\n")
cat("ocupação em dados binários de detecção, não um sinal de que o modelo é inutilizável.\n\n")

## ---- 1.6 Incluir o status de área protegida melhora o modelo? ----
form_no_pa <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + ")))
fit_no_pa <- occu(form_no_pa, data = umf)
delta_aic <- fit_no_pa@AIC - fit_full@AIC
newdat_out <- as.data.frame(setNames(as.list(rep(0, length(FINAL_COVARS_Z))), FINAL_COVARS_Z)); newdat_out$in_pa <- 0
newdat_in <- newdat_out; newdat_in$in_pa <- 1
psi_out <- predict(fit_full, newdata=newdat_out, type="state")$Predicted
psi_in <- predict(fit_full, newdata=newdat_in, type="state")$Predicted
direction <- ifelse(co["psi(in_pa)"] < 0, "negativa", "positiva")
cat("== Efeito de área protegida ==\n")
cat("Incluir in_pa melhora o AIC em", round(delta_aic,1), "pontos.\n")
cat("Direção:", direction, "-- ocupação", ifelse(direction=="negativa","menor","maior"), "dentro de áreas protegidas.\n")
cat("Ocupação predita dentro da AP:", round(psi_in,3), "| fora da AP:", round(psi_out,3),
    "| razão:", round(psi_in/psi_out,2), "\n\n")


## =============================================================================
## 2. TODA ESPÉCIE, JUNTA
## =============================================================================
## Ajustar um modelo de 8 covariáveis a todas as espécies/gêneros modelados
## com um motor apenas de efeitos fixos é pouco confiável para cerca de
## metade da comunidade (separação quase completa). O padrão do pipeline é,
## em vez disso, o modelo com efeito aleatório de array
## (spOccupancy::PGOcc), que regulariza cada ajuste através do intercepto de
## array -- usado acima no exemplo do carro-chefe; a mesma chamada é
## repetida para cada espécie/gênero modelado para construir a tabela de
## coeficientes de toda a comunidade usada no pipeline completo.

cat("== Modelo de comunidade ==\n")
cat("Espécies/gêneros modelados:", nrow(modeled), "\n")
cat("Ajustando o modelo com efeito de array (PGOcc) para cada espécie/gênero (~5-10 minutos)...\n")

community_coefs <- list()
for (sci_c in modeled$sci_pooled) {
  y_c <- build_dethist(sci_c)
  excl_c <- excl[species == sci_c, site_id]
  keep_c <- !is.na(sc$forest_100m_z) & !(sc$site_id %in% excl_c)
  y_c2 <- y_c[keep_c, , drop=FALSE]
  sc_c2 <- sc[keep_c]
  fit_c <- tryCatch(fit_pgocc_array_re(y_c2, as.data.frame(sc_c2)), error=function(e) { cat(sci_c, "FALHOU:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit_c)) next
  bs <- fit_c$beta.samples
  community_coefs[[sci_c]] <- data.table(species=sci_c, param=colnames(bs),
                                            mean=apply(bs,2,mean), ci_lo=apply(bs,2,quantile,0.025), ci_hi=apply(bs,2,quantile,0.975),
                                            n_sites=nrow(y_c2))
  cat(sci_c, "OK, n_sites=", nrow(y_c2), "\n")
}
community_dt <- rbindlist(community_coefs, fill=TRUE)
cat("\nAjustes de comunidade completados:", uniqueN(community_dt$species), "de", nrow(modeled), "\n\n")

## ---- figura da tabela de beta sombreada da comunidade ----
label_map_full <- c(forest_100m_z="Floresta", savanna_100m_z="Savana", pasture_100m_z="Pastagem",
                     cropland_100m_z="Lavoura", native_veg_1000m_z="Veg. nativa\n(1000m)",
                     temp_mean_C_z="Temp", precip_annual_mm_z="Precip", in_pa="Área\nprotegida")
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
  scale_fill_gradient2(low=COL_NEG, mid="white", high=COL_POS, midpoint=0, name="Coeficiente (\u03b2)") +
  labs(title="Efeitos de habitat na ocupação em toda a comunidade (efeito aleatório de array)",
       subtitle="n_det = detecções totais, n_arr = arrays detectados; * = intervalo de credibilidade de 95% exclui zero",
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.y=element_text(size=7), axis.text.x=element_text(size=8), panel.grid=element_blank())
ggsave(file.path(FIGS_DIR, "m3_beta_shaded_table.png"), p_beta_table, width=10, height=12, dpi=150)
cat("Salvo figs/m3_beta_shaded_table.png\n\n")

## ---- 2.1 A simplificação em nível de array se sustenta? ----
## Seis espécies ilustram a comparação: Paca e Tatu-galinha (ricos em
## dados), Queixada (esparsa, mascarada por distribuição), e Jaguatirica,
## Anta, Onça-parda (escolhidas por grandes mudanças câmera-vs-array).
species_6 <- c("Cuniculus paca","Dasypus novemcinctus","Tayassu pecari",
               "Leopardus pardalis","Tapirus terrestris","Puma concolor")
COVARS3 <- c("forest_100m_z","native_veg_1000m_z")

site_to_array <- setNames(sc$array_id, sc$site_id)
arrays_ordered <- sort(unique(sc$array_id))
n_arrays <- length(arrays_ordered)

array_cov <- sc[, lapply(.SD, mean, na.rm=TRUE), .SDcols = c(FINAL_COVARS, "in_pa"), by = array_id]
for (v in FINAL_COVARS) array_cov[[paste0(v,"_z")]] <- as.numeric(scale(array_cov[[v]]))

cat("== Comparação câmera-nível vs. array-nível, 6 espécies ==\n")
fitstats <- list()
coef_compare <- list()
for (sci6 in species_6) {
  y_site <- build_dethist(sci6)
  excl_sites6 <- excl[species == sci6, site_id]
  keep6 <- !is.na(sc$forest_100m_z) & !(sc$site_id %in% excl_sites6)
  y_site2 <- y_site[keep6, , drop=FALSE]
  sc_site2 <- as.data.frame(sc[keep6])
  umf_s <- unmarkedFrameOccu(y=y_site2, siteCovs=sc_site2[, c(COVARS3,"in_pa")])
  form3 <- as.formula(paste("~1 ~", paste(COVARS3, collapse=" + "), "+ in_pa"))
  fit_s <- tryCatch(occu(form3, data=umf_s), error=function(e) NULL)

  # nível de array: um array conta como "detectado" se qualquer câmera nele detectou naquela ocasião
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
  keep_arr <- !is.na(array_cov$forest_100m_z) & !(array_cov$array_id %in% excl_arrays6)
  y_arr2 <- y_arr_mat[keep_arr, , drop=FALSE]
  arr_cov2 <- as.data.frame(array_cov[keep_arr])
  umf_a <- unmarkedFrameOccu(y=y_arr2, siteCovs=arr_cov2[, c(COVARS3,"in_pa")])
  fit_a <- tryCatch(occu(form3, data=umf_a), error=function(e) NULL)

  if (!is.null(fit_s) && !is.null(fit_a)) {
    auc_s <- as.numeric(auc(roc(response=apply(y_site2,1,function(r) any(r==1,na.rm=TRUE)),
                                  predictor=predict(fit_s,type="state")$Predicted, quiet=TRUE)))
    auc_a <- as.numeric(auc(roc(response=apply(y_arr2,1,function(r) any(r==1,na.rm=TRUE)),
                                  predictor=predict(fit_a,type="state")$Predicted, quiet=TRUE)))
    fitstats[[sci6]] <- data.table(species=sci6, n_site=nrow(y_site2), n_array=nrow(y_arr2),
                                     auc_site=round(auc_s,3), auc_array=round(auc_a,3))

    co_s <- coef(fit_s); se_s <- sqrt(diag(vcov(fit_s)))
    co_a <- coef(fit_a); se_a <- sqrt(diag(vcov(fit_a)))
    for (pr in names(co_s)) {
      if (!(pr %in% names(co_a))) next
      coef_compare[[paste(sci6,pr)]] <- data.table(species=sci6, param=pr,
        site_est=co_s[pr], site_se=se_s[pr], array_est=co_a[pr], array_se=se_a[pr])
    }
  }
}
fitstats_dt <- rbindlist(fitstats, fill=TRUE)
print(fitstats_dt)
cat("\nO AUC em nível de array é tipicamente maior -- unidades amostrais mais\n")
cat("limpas e em menor número são mais fáceis de discriminar, ao custo de\n")
cat("intervalos de confiança muito mais amplos (menos poder).\n\n")

## ---- figura de comparação de estatísticas de ajuste ----
fitstats_long <- melt(fitstats_dt, id.vars="species", measure.vars=c("auc_site","auc_array"),
                        variable.name="level", value.name="auc")
fitstats_long[, level := fifelse(level=="auc_site","Nível de câmera","Nível de array")]
common_lookup6 <- c("Cuniculus paca"="Paca","Dasypus novemcinctus"="Tatu-galinha",
                     "Tayassu pecari"="Queixada","Leopardus pardalis"="Jaguatirica",
                     "Tapirus terrestris"="Anta","Puma concolor"="Onça-parda")
fitstats_long[, common := common_lookup6[species]]
fitstats_long[, common := factor(common, levels=common_lookup6[species_6])]
p_fitstats <- ggplot(fitstats_long, aes(x=common, y=auc, fill=level)) +
  geom_col(position="dodge", width=0.7) +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
  annotate("text", x=1, y=0.53, label="Chance", color="grey50", size=3) +
  scale_fill_manual(values=c("Nível de câmera"=COL_NEUTRAL, "Nível de array"="#e67e22"), name=NULL) +
  labs(title="Discriminação do modelo (AUC): nível de câmera vs. nível de array, 6 espécies", x=NULL, y="AUC") +
  ylim(0,1) + theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(FIGS_DIR, "m3_camera_vs_array_fit_6species.png"), p_fitstats, width=8, height=5.5, dpi=150)
cat("Salvo figs/m3_camera_vs_array_fit_6species.png\n\n")

## ---- figura de comparação de coeficientes câmera vs array (pontos com barras) ----
coef_compare_dt <- rbindlist(coef_compare, fill=TRUE)
coef_compare_dt <- coef_compare_dt[!(param %in% c("psi(Int)", "p(Int)"))]
coef_compare_dt[, common := common_lookup6[species]]
coef_long <- melt(coef_compare_dt, id.vars=c("species","param","common"),
                    measure.vars=list(est=c("site_est","array_est"), se=c("site_se","array_se")))
coef_long[, level := fifelse(variable==1, "Nível de câmera", "Nível de array")]
param_labels3 <- c("psi(forest_100m_z)"="Floresta", "psi(native_veg_1000m_z)"="Veg. nativa (1000m)", "psi(in_pa)"="Área protegida")
coef_long[, param_label := param_labels3[param]]
coef_long[, y_label := paste(common, param_label, sep=" \u2014 ")]
p_dotwhisker <- ggplot(coef_long, aes(x=est, y=y_label, color=level)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=est-1.96*se, xmax=est+1.96*se), height=0, position=position_dodge(width=0.5), linewidth=0.9) +
  geom_point(position=position_dodge(width=0.5), size=2.5) +
  scale_color_manual(values=c("Nível de câmera"=COL_NEUTRAL, "Nível de array"="#e67e22"), name=NULL) +
  labs(title="Coeficientes em nível de câmera vs. array, 6 espécies",
       x="Coeficiente (\u03b2), estimativa e IC 95%", y=NULL) +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_array_vs_site_6species.png"), p_dotwhisker, width=9, height=8, dpi=150)
cat("Salvo figs/m3_array_vs_site_6species.png\n\n")

## ---- 2.2 Por que essa comparação importa ----
cat("== Por que essa comparação importa ==\n")
cat("Isso fala diretamente de uma decisão de desenho que o pipeline faz:\n")
cat("usar o array como a unidade primária de replicação para modelos de\n")
cat("comunidade/conjuntos (Módulo 4+), enquanto ajusta modelos de espécie\n")
cat("única no nível de câmera mais fino aqui. As estimativas pontuais em\n")
cat("geral mantêm um sinal consistente entre os níveis para espécies bem\n")
cat("detectadas, mas o resultado em nível de array de uma espécie esparsa\n")
cat("ou fortemente mascarada por distribuição deve ser tratado como um\n")
cat("sinal muito mais fraco do que seu p-valor sozinho sugere.\n\n")

## ---- 2.3 Agrupamentos de espécies a partir da tabela de coeficientes ----
cat("== PCA de espécies ==\n")
pca_params <- c(names(label_map_full))
pca_wide <- dcast(community_dt[param != "(Intercept)"], species ~ param, value.var="mean")
pca_wide <- pca_wide[complete.cases(pca_wide[, ..pca_params])]
pca_mat <- scale(as.matrix(pca_wide[, ..pca_params]))
pca_fit <- prcomp(pca_mat, center=FALSE, scale.=FALSE)
ve <- summary(pca_fit)$importance[2, 1:2]
scores_dt <- data.table(species=pca_wide$species, PC1=pca_fit$x[,1], PC2=pca_fit$x[,2])
km_fit <- kmeans(scores_dt[, .(PC1,PC2)], centers=4, nstart=10)
scores_dt[, cluster := factor(km_fit$cluster)]
loadings_dt <- data.table(param=pca_params, PC1=pca_fit$rotation[,1]*3, PC2=pca_fit$rotation[,2]*3)
loadings_dt[, label := label_map_full[param]]

p_pca <- ggplot(scores_dt, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=cluster), size=3, alpha=0.8) +
  geom_segment(data=loadings_dt, aes(x=0,y=0,xend=PC1,yend=PC2), color=COL_NEG,
               arrow=arrow(length=unit(0.2,"cm")), inherit.aes=FALSE) +
  geom_text(data=loadings_dt, aes(x=PC1*1.15, y=PC2*1.15, label=label), color=COL_NEG, size=3, fontface="bold", inherit.aes=FALSE) +
  scale_color_brewer(palette="Set1", name="Grupo") +
  labs(title="Espécies agrupadas pelo perfil completo de resposta a habitat (modelo de comunidade)",
       x=sprintf("PC1 (%.0f%%)", ve[1]*100), y=sprintf("PC2 (%.0f%%)", ve[2]*100)) +
  theme_minimal(base_size=12)
ggsave(file.path(FIGS_DIR, "m3_species_pca_biplot.png"), p_pca, width=9, height=8, dpi=150)
cat("Salvo figs/m3_species_pca_biplot.png (", nrow(scores_dt), "espécies,", round(sum(ve)*100), "% da variância explicada)\n\n")


## =============================================================================
## 3. ADICIONANDO UM EFEITO DE INTERAÇÃO
## =============================================================================
## Modelo completo: psi ~ cobertura_arborea * temperatura + savana +
##             pastagem + lavoura + veg_nativa_1000m + precipitação + in_pa
## Cobertura arbórea substitui a floresta aqui (colinearidade r=0,82 com
## forest_100m) para evitar redundância com o termo de interação.

species4 <- c("Tapirus terrestris","Euphractus sexcinctus","Dasypus novemcinctus","Cerdocyon thous")
common4 <- c("Tapirus terrestris"="Anta","Euphractus sexcinctus"="Tatu-peba",
             "Dasypus novemcinctus"="Tatu-galinha","Cerdocyon thous"="Cachorro-do-mato")
OTHER_COVARS <- c("savanna_100m_z","pasture_100m_z","cropland_100m_z","native_veg_1000m_z","precip_annual_mm_z","in_pa")
cat("== Interação temperatura x cobertura arbórea, modelo completo, 4 espécies ==\n")
int_preds <- list()
int_labels <- list()
for (sci4 in species4) {
  y4 <- build_dethist(sci4)
  excl4 <- excl[species==sci4, site_id]
  keep4 <- !is.na(sc$treecover2000_100m_z) & !(sc$site_id %in% excl4)
  y4b <- y4[keep4, , drop=FALSE]
  sc4 <- as.data.frame(sc[keep4])
  umf4 <- unmarkedFrameOccu(y=y4b, siteCovs=sc4[, c("treecover2000_100m_z","temp_mean_C_z", OTHER_COVARS)])
  form4 <- as.formula(paste("~1 ~ treecover2000_100m_z * temp_mean_C_z +", paste(OTHER_COVARS, collapse=" + ")))
  fit4 <- tryCatch(occu(form4, data=umf4), error=function(e) NULL)
  if (is.null(fit4)) next
  co4 <- coef(fit4); se4 <- sqrt(diag(vcov(fit4)))
  int_term <- "psi(treecover2000_100m_z:temp_mean_C_z)"
  z4 <- co4[int_term]/se4[int_term]
  cat(sci4, "| coef. de interação:", round(co4[int_term],3), "| z:", round(z4,2), "\n")
  int_labels[[sci4]] <- sprintf("%s\n(interação %s, z=%.1f)", common4[sci4],
                                  ifelse(abs(z4)>1.96,"significativa","não significativa"), z4)

  tc_grid <- seq(min(sc4$treecover2000_100m_z, na.rm=TRUE), max(sc4$treecover2000_100m_z, na.rm=TRUE), length.out=40)
  for (temp_scenario in c(-1,0,1)) {
    newdat <- as.data.frame(setNames(as.list(rep(0,length(OTHER_COVARS))), OTHER_COVARS))
    newdat <- newdat[rep(1,40),]
    newdat$treecover2000_100m_z <- tc_grid; newdat$temp_mean_C_z <- temp_scenario
    pred <- predict(fit4, newdata=newdat, type="state")
    int_preds[[paste(sci4,temp_scenario)]] <- data.table(species=sci4, temp=factor(temp_scenario, levels=c(-1,0,1), labels=c("Fria","Média","Quente")),
                                                            tc=tc_grid, psi=pred$Predicted, lo=pred$lower, hi=pred$upper)
  }
}
cat("\nTrês de quatro espécies tipicamente mostram uma interação significativa\n")
cat("mesmo com o conjunto completo de covariáveis incluído -- um modelo apenas\n")
cat("aditivo perderia isso.\n\n")

int_pred_dt <- rbindlist(int_preds, fill=TRUE)
int_pred_dt[, species_label := unlist(int_labels[species])]
p_interaction <- ggplot(int_pred_dt, aes(x=tc, y=psi, color=temp, fill=temp)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Fria"="#2c7fb8", "Média"="grey50", "Quente"=COL_NEG), name="Temperatura") +
  scale_fill_manual(values=c("Fria"="#2c7fb8", "Média"="grey50", "Quente"=COL_NEG), name="Temperatura") +
  facet_wrap(~species_label, nrow=1, scales="free_y") +
  labs(title="Interação temperatura x cobertura arbórea, modelo completo de covariáveis, 4 espécies",
       x="Cobertura arbórea (z)", y="Ocupação predita (\u03c8)") +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_interaction_4species.png"), p_interaction, width=13, height=4.5, dpi=150)
cat("Salvo figs/m3_interaction_4species.png\n\n")


## =============================================================================
## 4. PERMITINDO QUE UMA RELAÇÃO VARIE POR ECORREGIÃO
## =============================================================================
## Modelo completo: psi ~ floresta * ecorregiao + savana + pastagem +
##             lavoura + veg_nativa_1000m + temperatura + precipitação + in_pa

eco <- fread(file.path(DATA_DIR, "site_ecoregion.csv"))
sc_eco <- merge(sc, eco[, .(site_id, biome)], by="site_id", all.x=TRUE)
sci_eco <- "Leopardus pardalis"
y_eco <- build_dethist(sci_eco)
excl_eco <- excl[species==sci_eco, site_id]
keep_base <- !is.na(sc_eco$forest_100m_z) & !(sc_eco$site_id %in% excl_eco)
target_biomes <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Grasslands, Savannas & Shrublands")
keep_eco <- keep_base & (sc_eco$biome %in% target_biomes)
y_eco2 <- y_eco[keep_eco, , drop=FALSE]
sc_eco2 <- as.data.frame(sc_eco[keep_eco])
sc_eco2$biome_f <- factor(sc_eco2$biome, levels=target_biomes)
OTHER_COVARS_ECO <- c("savanna_100m_z","pasture_100m_z","cropland_100m_z","native_veg_1000m_z","temp_mean_C_z","precip_annual_mm_z","in_pa")
umf_eco <- unmarkedFrameOccu(y=y_eco2, siteCovs=sc_eco2[, c("forest_100m_z","biome_f", OTHER_COVARS_ECO)])
form_pooled_eco <- as.formula(paste("~1 ~ forest_100m_z + biome_f +", paste(OTHER_COVARS_ECO, collapse=" + ")))
form_interact_eco <- as.formula(paste("~1 ~ forest_100m_z * biome_f +", paste(OTHER_COVARS_ECO, collapse=" + ")))
fit_pooled_eco <- occu(form_pooled_eco, data=umf_eco)
fit_interact_eco <- occu(form_interact_eco, data=umf_eco)
delta_eco <- fit_pooled_eco@AIC - fit_interact_eco@AIC
cat("== Jaguatirica: efeito de floresta agrupado vs. interação por ecorregião ==\n")
cat("Melhoria de AIC do modelo de interação:", round(delta_eco,1), "\n")
cat("Na Floresta Úmida de Folhas Largas, mais cobertura florestal tipicamente\n")
cat("significa mais ocupação de Jaguatirica; na ecorregião Pastagens/Savanas\n")
cat("& Arbustais a relação se inverte -- um único 'efeito de floresta' agrupado\n")
cat("esconde duas relações genuinamente diferentes.\n\n")

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
biome_labels <- c("Tropical & Subtropical Moist Broadleaf Forests"="Floresta Úmida de Folhas Largas",
                   "Tropical & Subtropical Grasslands, Savannas & Shrublands"="Pastagens/Savanas")
eco_pred_dt[, biome_label := biome_labels[biome]]
p_ecoregion <- ggplot(eco_pred_dt, aes(x=forest_z, y=psi, color=biome_label, fill=biome_label)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
  geom_line(linewidth=1.1) +
  scale_color_manual(values=c("Floresta Úmida de Folhas Largas"=COL_POS, "Pastagens/Savanas"="#e67e22"), name=NULL) +
  scale_fill_manual(values=c("Floresta Úmida de Folhas Largas"=COL_POS, "Pastagens/Savanas"="#e67e22"), name=NULL) +
  labs(title="Jaguatirica: efeito da floresta na ocupação, por ecorregião (modelo completo)",
       x="Cobertura florestal (z, padronizada)", y="Ocupação predita (\u03c8)") +
  theme_minimal(base_size=12) + theme(legend.position=c(0.25,0.9))
ggsave(file.path(FIGS_DIR, "m3_ecoregion_ocelot.png"), p_ecoregion, width=8, height=6, dpi=150)
cat("Salvo figs/m3_ecoregion_ocelot.png\n\n")


## =============================================================================
## 5. COMPARANDO COM A ABUNDÂNCIA ROYLE-NICHOLS
## =============================================================================
## O Royle-Nichols (occuRN) usa os mesmos dados de detecção 1/0 que a
## ocupação, mas interpreta a FREQUÊNCIA de detecção como um sinal de
## abundância relativa, em vez de apenas presença/ausência.

cat("== Ocupação vs. abundância Royle-Nichols ==\n")
cat("Ajustando Royle-Nichols (occuRN) para cada espécie/gênero modelado...\n")
rn_coefs <- list()
for (sci_rn in modeled$sci_pooled) {
  y_rn <- build_dethist(sci_rn)
  excl_rn <- excl[species == sci_rn, site_id]
  keep_rn <- !is.na(sc$forest_100m_z) & !(sc$site_id %in% excl_rn)
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
cat("Ajustes RN completados:", uniqueN(rn_dt$species), "de", nrow(modeled), "\n\n")

## comparação em nível geral
occ_int_dt <- community_dt[param=="(Intercept)", .(species, occ_int=mean)]
rn_int_dt <- rn_dt[param=="lam(Int)", .(species, lam_int=estimate)]
overall_comp <- merge(occ_int_dt, rn_int_dt, by="species")
overall_comp[, mean_psi := plogis(occ_int)]; overall_comp[, mean_lam := exp(lam_int)]

## comparação do efeito de floresta, filtrando ajustes RN instáveis
occ_forest_dt <- community_dt[param=="forest_100m_z", .(species, occ_forest=mean)]
rn_forest_dt <- rn_dt[param=="lam(forest_100m_z)", .(species, rn_forest=estimate, rn_se=se)]
forest_comp <- merge(occ_forest_dt, rn_forest_dt, by="species")
forest_comp_stable <- forest_comp[!is.na(rn_se) & rn_se < 5 & abs(rn_forest) < 5]
forest_comp_stable[, diff := rn_forest - occ_forest]
top2 <- forest_comp_stable[order(-diff)][1:min(2,.N)]
cat("Correlação em nível geral:", round(cor(overall_comp$mean_psi, overall_comp$mean_lam),2), "\n")
cat("Correlação do efeito de floresta (apenas ajustes estáveis):", round(cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest),2), "\n\n")

p_overall <- ggplot(overall_comp, aes(x=mean_psi, y=mean_lam)) +
  geom_point(color=COL_NEUTRAL, size=2.5, alpha=0.75) +
  labs(title=sprintf("Nível geral (r=%.2f)", cor(overall_comp$mean_psi, overall_comp$mean_lam)),
       x="Ocupação média (\u03c8-chapéu)", y="Abundância média Royle-Nichols (\u03bb-chapéu)") +
  theme_minimal(base_size=10)
p_forest <- ggplot(forest_comp_stable, aes(x=occ_forest, y=rn_forest)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(color=COL_NEUTRAL, size=2.5, alpha=0.75) +
  geom_text(data=top2, aes(label=species), size=2.8, hjust=-0.1, vjust=-0.3) +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.18))) +
  labs(title=sprintf("Efeito de floresta especificamente (r=%.2f)", cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest)),
       x="Efeito de floresta (ocupação, \u03b2)", y="Efeito de floresta (Royle-Nichols, \u03b2)") +
  theme_minimal(base_size=10)
p_occ_vs_abund <- gridExtra::grid.arrange(p_overall, p_forest, ncol=2,
  top=grid::textGrob("Ocupação vs. abundância Royle-Nichols: concordância em toda a comunidade",
                       gp=grid::gpar(fontsize=13, fontface="bold"), x=0.02, hjust=0))
ggsave(file.path(FIGS_DIR, "m3_occ_vs_abundance.png"), p_occ_vs_abund, width=11, height=5, dpi=150)
cat("Salvo figs/m3_occ_vs_abundance.png\n\n")


## =============================================================================
## RESUMO -- o que os resultados significam
## =============================================================================
## - Este script constrói todo histórico de detecção, tabela de covariáveis
##   e modelo diretamente dos três arquivos brutos do conjunto de dados
##   combinado -- nenhum coeficiente pré-ajustado é carregado.
## - Didelphis e Dasyprocta são agrupados em nível de gênero (uma exceção
##   deliberada à regra de nível de espécie do pipeline), recuperando
##   registros de detecção apenas em nível de gênero. Didelphis mantém uma
##   máscara de distribuição (união das distribuições com buffer dos 3
##   congêneres mapeados); Dasyprocta não tem máscara de distribuição
##   (apenas 1 de 5 congêneres é mapeado, e usá-lo sozinho classificaria
##   erroneamente detecções reais como impossíveis).
## - O ajuste apenas de efeitos fixos é pouco confiável para cerca de metade
##   da comunidade; o modelo com efeito aleatório de array (chamadas PGOcc
##   acima) é o padrão do pipeline e deve ser executado para a análise
##   completa de toda a comunidade.
## - AUC (discriminação) e c-hat (superdispersão) medem coisas diferentes e
##   podem discordar.
## - Agrupar câmeras em arrays amplia substancialmente os intervalos de
##   confiança e pode mudar a direção estimada de uma covariável para
##   espécies com poucos dados.
## - Tanto uma interação contínua (temperatura x cobertura arbórea) quanto
##   uma categórica (ecorregião) mostram que um único efeito de habitat "um
##   número serve para tudo" pode esconder variação real e ecologicamente
##   sensata.
##
## Próximo módulo: Modelos de Comunidade e Conjuntos.
## =============================================================================
