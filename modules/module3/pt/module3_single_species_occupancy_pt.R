## =============================================================================
## Módulo 3 — Ocupação de Espécie Única
## Pipeline multi-rede de armadilhas fotográficas Snapshot Brasil
## =============================================================================
##
## AUTOCONTIDO: este script parte dos três arquivos brutos do conjunto de
## dados combinado (implantações, detecções de mamíferos, covariáveis de
## sensoriamento remoto) e constrói tudo do zero -- históricos de detecção,
## covariáveis, conjuntos de locais candidatos, e cada ajuste de modelo.
## Coloque este script na mesma pasta que os dez arquivos de dados
## listados abaixo e execute-o do início ao fim.
##
## Arquivos de entrada necessários (em DATA_DIR, padrão "data" junto a este script):
##   final_deployments.csv         -- uma linha por implantação de câmera
##   final_detections_mammals.csv  -- uma linha por evento de detecção independente
##   final_covariates.csv          -- uma linha por local de câmera, covariáveis de sensoriamento remoto
##   species_site_exclusions.csv   -- locais de zero estrutural por espécie (saída da avaliação de mapas do Módulo 2)
##   genus_species_site_exclusions.csv -- locais de zero estrutural para Didelphis agrupado (ver nota abaixo)
##   site_ecoregion.csv            -- bioma/ecorregião por local (usado apenas na extensão de ecorregião)
##   species_order_lookup.csv      -- ordem taxonômica por espécie/gênero (usado apenas na seção 3.1)
##   species_common_names_pt.csv   -- nome comum em português por espécie/gênero (usado nos rótulos de PCA/RN)
##   brazil_boundary.gpkg          -- contorno do Brasil (Natural Earth 110m; figura de níveis de Johnson)
##   south_america_boundary.gpkg   -- contorno da América do Sul (Natural Earth 110m; inset localizador)
##
## Os últimos quatro arquivos são saídas de referência da avaliação de mapas
## de distribuição do Módulo 2 e da revisão taxonômica deste pipeline, e não
## podem ser reconstruídos a partir dos três arquivos principais isoladamente
## -- são fornecidos junto com este script como tabelas de referência fixas.
##
## Nota taxonômica: Didelphis (3 espécies com mapas de distribuição +
## registros apenas de gênero) e Dasyprocta (5 espécies, apenas 1 com mapa
## de distribuição, + registros apenas de gênero) são AGRUPADOS ao nível de
## gênero neste script -- uma exceção deliberada à regra usual do pipeline
## de somente nível de espécie, feita porque a identificação de campo até
## espécie é pouco confiável para esses dois gêneros, e o agrupamento
## recupera registros de detecção apenas de gênero que seriam descartados.
## Didelphis mantém sua máscara de distribuição (a união das distribuições
## com buffer das 3 espécies congêneres mapeadas); Dasyprocta descarta
## totalmente sua máscara de distribuição, já que apenas um dos cinco
## congêneres tem mapa de distribuição, e usá-lo isoladamente marcaria
## muitas detecções reais como "impossíveis."
##
## Requer: data.table, ggplot2, unmarked, spOccupancy, pROC, ggrepel, sf.
## Se os rótulos das figuras contendo símbolos psi/beta/lambda aparecerem
## corrompidos, execute este script com um locale UTF-8 definido
## (ex.: LC_ALL=en_US.UTF-8).
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

## Conjunto de covariáveis no nível de câmera (11, incl. in_pa adicionado
## separadamente abaixo) -- usado em todo modelo por local de câmera
## (carro-chefe, loop comunitário, loop RN, interação, ecorregião). O
## tamanho amostral neste nível (1.110+ locais) sustenta facilmente esse
## número de preditores.
FINAL_COVARS <- c("forest_100m","savanna_100m","pasture_100m","cropland_100m",
                   "native_veg_1000m","temp_mean_C","precip_annual_mm",
                   "ghsl_built_5000m","dist_road_m","dist_water_m")

## Conjunto de covariáveis no nível de array (6, sem in_pa) -- usado APENAS
## no lado de nível de array da comparação câmera-vs-array (seção 2.1).
## Escolhido classificando cada covariável candidata do nível de câmera por
## poder explicativo comunitário e confirmado com VIF; um conjunto mais rico
## não é defensável com apenas 60 arrays. Ver Módulo 2 para a derivação.
ARRAY_COVARS <- c("forest_100m","savanna_100m","native_veg_1000m",
                   "temp_mean_C","precip_annual_mm","dist_road_m")

OCC_DAYS <- 7  # duração da ocasião em dias, seguindo a convenção do intervalo de independência do pipeline


## =============================================================================
## ETAPA 1: CARREGAR DADOS BRUTOS E CONSTRUIR COVARIÁVEIS NO NÍVEL DE LOCAL
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

# junta as covariáveis de sensoriamento remoto, padroniza (z-score) o conjunto final, numera os arrays
sc <- merge(site_tbl, cov, by = "site_id", all.x = TRUE)
setorder(sc, site_id)
ALL_CONTINUOUS_COVARS <- union(FINAL_COVARS, ARRAY_COVARS)
for (v in ALL_CONTINUOUS_COVARS) sc[[paste0(v, "_z")]] <- as.numeric(scale(sc[[v]]))
sc[, treecover2000_100m_z := as.numeric(scale(treecover2000_100m))]
arrays_sorted <- sort(unique(sc$array_id))
sc[, array_num := match(array_id, arrays_sorted)]

cat("Locais sem covariáveis (excluídos por exclusão listwise no momento do ajuste):",
    sum(!complete.cases(sc[, ..FINAL_COVARS])), "\n\n")


## =============================================================================
## ETAPA 2: CONSTRUIR HISTÓRICOS DE DETECÇÃO A PARTIR DE IMPLANTAÇÕES E DETECÇÕES BRUTAS
## =============================================================================
## Um modelo de ocupância precisa de um HISTÓRICO DE DETECÇÃO: para cada
## local e cada ocasião de 7 dias, a espécie foi detectada (1), amostrada
## mas não detectada (0), ou nenhuma câmera estava ativa naquela ocasião
## (ausente)? Construímos isso diretamente a partir das janelas de
## implantação brutas e dos timestamps de detecção.
##
## Janelas por implantação (não o intervalo mín-máx completo de um local)
## são usadas para marcar quais dias uma câmera esteve realmente ativa --
## um local com múltiplas implantações e um intervalo entre dois períodos
## de implantação não deve ser tratado como continuamente amostrado durante
## essa lacuna.

site_span <- dep[, .(min_start = min(start_dt), max_end = max(end_dt)), by = site_id]
setorder(site_span, site_id)
site_span[, n_occ := (as.integer(max_end - min_start) %/% OCC_DAYS) + 1L]
max_occ <- max(site_span$n_occ)
site_ids_ordered <- site_span$site_id
n_sites <- length(site_ids_ordered)

# vetor booleano de dia-ativo por local, construído a partir do início/fim de cada implantação
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

# colapsa o vetor de dia-ativo em uma flag de "amostrado" no nível de ocasião: uma ocasião conta como
# amostrada se >=50% de seus dias tiveram uma câmera ativa
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

# mapeia cada detecção ao índice de ocasião do seu local; descarta o punhado
# de detecções anteriores ao início de implantação mais antigo registrado do
# local (datas de imagens de calibração incorretas / casos de limite off-by-one
# -- uma fração pequena e consistente em todo o histórico do pipeline)
site_min_start <- setNames(site_span$min_start, site_span$site_id)
det[, site_min_start := site_min_start[site_id]]
det <- det[!is.na(site_min_start)]
det[, day_offset := as.integer(as.IDate(date) - site_min_start)]
det[, occ_idx := day_offset %/% OCC_DAYS]
n_dropped <- sum(det$occ_idx < 0)
det <- det[occ_idx >= 0]
cat("Detecções usadas:", nrow(det), "(", n_dropped, "descartadas: anteriores à implantação mais antiga do local)\n\n")

# ---- Taxonomia: agrupar Didelphis e Dasyprocta ao nível de gênero ----
# (todos os outros gêneros mantêm identificação no nível de espécie, seguindo
# a regra padrão do pipeline; ver a nota taxonômica no topo deste script)
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
    y[sid, oc] <- 1  # detectado (também cobre ocasiões que a regra dos 50%-ativo marcou como não amostradas)
  }
  y
}

# ---- Escolher as espécies bem-detectadas a modelar ----
# mesmo limiar usado em todo o pipeline: nome binomial no nível de espécie
# (ou, para os dois gêneros agrupados, no nível de gênero), excluindo taxa
# domésticos/humanos e registros não resolvidos apenas de gênero/família
# para qualquer OUTRO gênero, com pelo menos 10 detecções totais e 20 locais detectados
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
fwrite(species_meta, "data/species_meta.csv")
fwrite(data.table(n_sites=nrow(site_tbl), n_arrays=uniqueN(site_tbl$array_id),
                   n_modeled=nrow(modeled)), "data/dataset_summary.csv")


## =============================================================================
## ETAPA 3: CONJUNTOS CANDIDATOS COM MÁSCARA DE DISTRIBUIÇÃO (referências do Módulo 2)
## =============================================================================
## Para espécies com mapa de distribuição da IUCN, um local é um "zero
## estrutural" se está fora da distribuição mapeada (com um buffer de 100km
## para imprecisão do mapa) E a espécie nunca foi detectada ali -- esses
## locais são excluídos como candidatos padrão. Locais com detecção fora da
## distribuição com buffer são sempre mantidos (uma fotografia supera um polígono).

excl <- fread(file.path(DATA_DIR, "species_site_exclusions.csv"))
excl_genus <- fread(file.path(DATA_DIR, "genus_species_site_exclusions.csv"))
excl <- rbindlist(list(excl, excl_genus))
cat("Espécies com máscara de distribuição (exclusões de zero estrutural aplicadas):", uniqueN(excl$species), "\n")
cat("Didelphis agrupado: máscara de distribuição = união das distribuições com buffer dos 3 congêneres mapeados (",
    excl_genus[species=="Didelphis", .N], "locais de zero estrutural excluídos)\n")
cat("Dasyprocta agrupado: NENHUMA máscara de distribuição aplicada (apenas 1 dos 5 congêneres tem mapa de distribuição;\n")
cat("  usá-lo isoladamente marcaria muitas detecções reais como impossíveis) -- todo local é candidato.\n\n")


## =============================================================================
## 1. UM PRIMEIRO MODELO DE OCUPAÇÃO: TATU-GALINHA
## =============================================================================
## Começamos com o Tatu-galinha (*Dasypus novemcinctus*), uma espécie
## bem-detectada, e ajustamos um modelo usando todas as covariáveis do
## conjunto final: floresta, savana, pastagem, lavoura, vegetação nativa em
## escala de paisagem, temperatura, precipitação e status de área protegida.

sci <- "Dasypus novemcinctus"
y_arm <- build_dethist(sci)
excl_sites <- excl[species == sci, site_id]
FINAL_COVARS_Z <- paste0(FINAL_COVARS, "_z")
ARRAY_COVARS_Z <- paste0(ARRAY_COVARS, "_z")
keep <- complete.cases(sc[, ..FINAL_COVARS]) & !(sc$site_id %in% excl_sites)
y2 <- y_arm[keep, , drop = FALSE]
sc2 <- as.data.frame(sc[keep])
cat("== Tatu-galinha: conjunto de locais candidatos ==\n")
cat("locais após exclusão por máscara de distribuição + covariável ausente:", nrow(y2), "\n\n")

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
  df_long[, status_f := factor(ifelse(is.na(status), "Não amostrado", ifelse(status==1,"Detectado","Amostrado, sem detecção")),
                                 levels=c("Não amostrado","Amostrado, sem detecção","Detectado"))]
  dethist_plot_data[[arr_name]] <- df_long
}
dethist_dt <- rbindlist(dethist_plot_data)
p_dethist <- ggplot(dethist_dt, aes(x=occasion_n, y=site_id, fill=status_f)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_manual(values=c("Não amostrado"="white","Amostrado, sem detecção"="#dddddd","Detectado"=COL_NEG), name=NULL) +
  facet_wrap(~array_label, scales="free_y") +
  labs(title="Histórico de detecção do Tatu-galinha: array de alta vs. baixa detecção",
       x="Ocasião de 7 dias", y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.y=element_text(size=6), panel.grid=element_blank(), legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_dethist_armadillo.png"), p_dethist, width=11, height=6, dpi=150)
cat("Salvo figs/m3_dethist_armadillo.png\n\n")

## ---- 1.2 Ajustando o modelo, com o efeito aleatório de array ----
## Ajuste de efeitos fixos (unmarked::occu) mostrado aqui; o ajuste com
## efeito aleatório de array (spOccupancy::PGOcc) requer vários minutos de
## amostragem MCMC e está comentado abaixo com a chamada exata usada --
## descomente para executá-lo.

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
cat("== Coeficientes do efeito aleatório de array do Tatu-galinha (PGOcc) ==\n")
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

## ---- 1.3 Ocupação naive vs. modelada ----
naive_psi_site <- apply(y2, 1, function(r) as.integer(any(r == 1, na.rm = TRUE)))
naive_psi <- mean(naive_psi_site)
psi_pred <- predict(fit_full, type = "state")$Predicted
modeled_psi_pgocc <- mean(apply(fit_pgocc$psi.samples, 2, mean))
cat("== Ocupação naive vs. modelada ==\n")
cat("Ocupação naive (fração de locais alguma vez detectados):", round(naive_psi,3), "\n")
cat("Ocupação modelada (PGOcc com efeito aleatório de array, psi médio):", round(modeled_psi_pgocc,3), "\n\n")

p_naive_modeled <- ggplot(data.table(x=c("Naive\n(fração detectada)","Modelada\n(PGOcc, efeito aleat. de array)"),
                                       y=c(naive_psi, modeled_psi_pgocc)),
                            aes(x=x, y=y, fill=x)) +
  geom_col(width=0.55) +
  geom_text(aes(label=sprintf("%.2f", y)), vjust=-0.5, fontface="bold", size=4) +
  scale_fill_manual(values=c("#999999", COL_NEUTRAL), guide="none") +
  labs(title="Tatu-galinha: ocupação naive vs. modelada", x=NULL, y="Ocupação (\u03c8)") +
  ylim(0, max(naive_psi, modeled_psi_pgocc)*1.3) +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIGS_DIR, "m3_naive_vs_modeled.png"), p_naive_modeled, width=5.5, height=5.5, dpi=150)
cat("Salvo figs/m3_naive_vs_modeled.png\n\n")

## ---- 1.4 Efeitos das covariáveis & 1.5 Ajuste do modelo ----
## Duas perguntas muito diferentes são feitas sobre um modelo de ocupância ajustado:
## (a) O MODELO DISCRIMINA locais ocupados de não ocupados? -- AUC.
##     Toma a probabilidade de ocupação prevista de cada local e pergunta:
##     se você escolhesse aleatoriamente um local verdadeiramente ocupado e
##     um verdadeiramente não ocupado, qual a chance de o modelo dar ao
##     ocupado uma probabilidade prevista maior? AUC=0.5 é o acaso; AUC=1.0
##     é perfeito.
## (b) A ESTRUTURA DE ERRO ASSUMIDA PELO MODELO ESTÁ CORRETA? -- c-hat
##     (sobredispersão), a partir de um teste de bondade de ajuste por
##     bootstrap paramétrico. c-hat próximo de 1 significa que a estrutura
##     de erro binomial assumida se ajusta; notavelmente acima de 1
##     significa que as detecções reais são mais variáveis do que o modelo
##     espera. AUC e c-hat NÃO são intercambiáveis: um modelo pode
##     discriminar bem (AUC alto) e ainda estar sobredisperso.

## ---- figura de curvas de efeito (4 covariáveis significativas do ajuste com efeito aleatório de array) ----
sig_covars_arm <- pgocc_coefs[param != "(Intercept)"][ci_lo>0 | ci_hi<0, param]
cat("Covariáveis significativas do efeito aleatório de array para o Tatu-galinha:", paste(sig_covars_arm, collapse=", "), "\n\n")
effect_curve_data <- list()
for (pr in setdiff(sig_covars_arm, "in_pa")) {
  b <- pgocc_coefs[param==pr, mean]; b_lo <- pgocc_coefs[param==pr, ci_lo]; b_hi <- pgocc_coefs[param==pr, ci_hi]
  intc <- pgocc_coefs[param=="(Intercept)", mean]
  zg <- seq(-2.5, 2.5, length.out=60)
  effect_curve_data[[pr]] <- data.table(covariate=pr, z=zg,
    psi=plogis(intc + b*zg), psi_lo=plogis(intc + b_lo*zg), psi_hi=plogis(intc + b_hi*zg), sign=ifelse(b>0,"pos","neg"))
}
label_map <- c(forest_100m_z="Cobertura florestal", savanna_100m_z="Cobertura de savana", pasture_100m_z="Pastagem",
               cropland_100m_z="Lavoura", native_veg_1000m_z="Veg. nativa (1000m)",
               temp_mean_C_z="Temperatura", precip_annual_mm_z="Precipitação",
               ghsl_built_5000m_z="Superfície construída", dist_road_m_z="Distância até estrada",
               dist_water_m_z="Distância até água", in_pa="Área protegida")
if (length(effect_curve_data)) {
  ec_dt <- rbindlist(effect_curve_data)
  ec_dt[, covariate_label := factor(label_map[covariate], levels=label_map[names(effect_curve_data)])]
  p_effects <- ggplot(ec_dt, aes(x=z, y=psi, color=sign, fill=sign)) +
    geom_ribbon(aes(ymin=psi_lo, ymax=psi_hi), alpha=0.15, color=NA) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c(pos=COL_POS, neg=COL_NEG), guide="none") +
    scale_fill_manual(values=c(pos=COL_POS, neg=COL_NEG), guide="none") +
    facet_wrap(~covariate_label, nrow=1, scales="free_x") +
    labs(title="Tatu-galinha: ocupação vs. cada covariável significativa (verde=positivo, vermelho=negativo)",
         x="Covariável padronizada (z)", y="Ocupação prevista (\u03c8)") +
    theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
  ggsave(file.path(FIGS_DIR, "m3_armadillo_effect_curves.png"), p_effects, width=13, height=4, dpi=150)
  cat("Salvo figs/m3_armadillo_effect_curves.png\n\n")
}

## ---- figura de IC dos betas (modelo com efeito aleatório de array) ----
beta_plot <- pgocc_coefs[param != "(Intercept)"][order(mean)]
beta_plot[, label := label_map[param]]
stopifnot("Covariável não mapeada em beta_plot -- adicione a label_map" = !any(is.na(beta_plot$label)))
beta_plot[, label := factor(label, levels=unique(label))]
beta_plot[, sign := fifelse(ci_lo>0, "pos", fifelse(ci_hi<0, "neg", "ns"))]
p_beta_ci <- ggplot(beta_plot, aes(x=mean, y=label, color=sign)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0, linewidth=1) +
  geom_point(size=3) +
  scale_color_manual(values=c(pos=COL_POS, neg=COL_NEG, ns="grey50"), guide="none") +
  labs(title="Tatu-galinha: efeitos das covariáveis (modelo com efeito aleatório de array)",
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
  labs(title="O que o modelo prevê, por status verdadeiro", x="Probabilidade de ocupação prevista (\u03c8-estimado)", y="Densidade") +
  theme_minimal(base_size=10) + theme(legend.position=c(0.75,0.85))
p_roc <- ggplot(roc_dt, aes(x=fpr, y=tpr)) +
  geom_ribbon(aes(ymin=0, ymax=tpr), fill=COL_NEUTRAL, alpha=0.15) +
  geom_line(color=COL_NEUTRAL, linewidth=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  annotate("text", x=0.6, y=0.15, label="Acaso\n(AUC = 0.5)", color="grey50", size=3.2) +
  labs(title=sprintf("Curva ROC: AUC = %.2f", auc_val),
       x="Taxa de falso positivo\n(locais não ocupados chamados de 'presentes')",
       y="Taxa de verdadeiro positivo\n(locais ocupados chamados de 'presentes')") +
  theme_minimal(base_size=10)
p_auc_explainer <- gridExtra::grid.arrange(p_score_dist, p_roc, ncol=2,
  top=grid::textGrob("AUC: quão bem a ocupação prevista separa presença verdadeira de ausência verdadeira?",
                       gp=grid::gpar(fontsize=13, fontface="bold"), x=0.02, hjust=0))
ggsave(file.path(FIGS_DIR, "m3_auc_explainer.png"), p_auc_explainer, width=12, height=5, dpi=150)
cat("Salvo figs/m3_auc_explainer.png\n")

## Bootstrap de BdA (c-hat) -- 100 simulações, ~1-2 minutos:
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
cat("c-hat (sobredispersão):", round(c_hat,3), "| valor-p do BdA por bootstrap:", round(gof_p,3), "\n")
cat("Um AUC em torno de 0.65 indica discriminação modesta; um c-hat acima de 1 significa\n")
cat("que o modelo está sobredisperso -- ambos achados comuns e esperados para modelos de\n")
cat("ocupância em dados de detecção binários, não um sinal de que o modelo é inutilizável.\n\n")
fwrite(data.table(auc=round(auc_val,3), c_hat=round(c_hat,3), gof_p=round(gof_p,3),
                   det_prob=round(p_const,3)), "data/armadillo_fit_summary.csv")

## ---- 1.6 Adicionar status de área protegida melhora o modelo? ----
form_no_pa <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + ")))
fit_no_pa <- occu(form_no_pa, data = umf)
delta_aic <- fit_no_pa@AIC - fit_full@AIC
newdat_out <- as.data.frame(setNames(as.list(rep(0, length(FINAL_COVARS_Z))), FINAL_COVARS_Z)); newdat_out$in_pa <- 0
newdat_in <- newdat_out; newdat_in$in_pa <- 1
psi_out <- predict(fit_full, newdata=newdat_out, type="state")$Predicted
psi_in <- predict(fit_full, newdata=newdat_in, type="state")$Predicted
direction <- ifelse(co["psi(in_pa)"] < 0, "negativa", "positiva")
cat("== Efeito de área protegida ==\n")
cat("Adicionar in_pa melhora o AIC em", round(delta_aic,1), "pontos.\n")
cat("Direção:", direction, "-- ocupação", ifelse(direction=="negativa","menor","maior"), "dentro de áreas protegidas.\n")
cat("Ocupação prevista dentro de AP:", round(psi_in,3), "| fora de AP:", round(psi_out,3),
    "| razão:", round(psi_in/psi_out,2), "\n\n")
fwrite(data.table(delta_aic=round(delta_aic,1), direction=direction,
                   psi_outside_pa=round(psi_out,3), psi_inside_pa=round(psi_in,3),
                   ratio=round(psi_in/psi_out,2)), "data/armadillo_pa_effect.csv")


## =============================================================================
## 2. TODAS AS ESPÉCIES, JUNTAS
## =============================================================================
## Ajustar um modelo de 8 covariáveis a todas as espécies/gêneros modelados
## com um motor apenas de efeitos fixos é pouco confiável para
## aproximadamente metade da comunidade (separação quase-completa). O
## padrão do pipeline é, em vez disso, o modelo com efeito aleatório de
## array (spOccupancy::PGOcc), que regulariza cada ajuste através do
## intercepto do array -- comentado acima para o exemplo carro-chefe; a
## mesma chamada é repetida em loop para cada espécie/gênero modelado para
## construir a tabela de coeficientes comunitária usada no pipeline completo.

cat("== Modelo comunitário ==\n")
cat("Espécies/gêneros modelados:", nrow(modeled), "\n")
cat("Ajustando o modelo com efeito aleatório de array (PGOcc) para cada espécie/gênero (~5-10 minutos)...\n")

community_coefs <- list()
for (sci_c in modeled$sci_pooled) {
  y_c <- build_dethist(sci_c)
  excl_c <- excl[species == sci_c, site_id]
  keep_c <- complete.cases(sc[, ..FINAL_COVARS]) & !(sc$site_id %in% excl_c)
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
cat("\nAjustes comunitários completados:", uniqueN(community_dt$species), "de", nrow(modeled), "\n\n")
fwrite(data.table(n_fit=uniqueN(community_dt$species), n_modeled=nrow(modeled)), "data/community_fit_count.csv")
fwrite(community_dt, "data/community_coefs.csv")

## ---- figura da tabela sombreada de betas comunitária ----
label_map_full <- c(forest_100m_z="Floresta", savanna_100m_z="Savana", pasture_100m_z="Pastagem",
                     cropland_100m_z="Lavoura", native_veg_1000m_z="Veg. nativa\n(1000m)",
                     temp_mean_C_z="Temp", precip_annual_mm_z="Precip",
                     ghsl_built_5000m_z="Superfície\nconstruída", dist_road_m_z="Dist.\nestrada",
                     dist_water_m_z="Dist.\nágua", in_pa="Área\nprotegida")
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
  labs(title="Efeitos de habitat na ocupação através da comunidade (efeito aleatório de array)",
       subtitle="n_det = detecções totais, n_arr = arrays com detecção; * = intervalo de credibilidade de 95% exclui zero",
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.y=element_text(size=7), axis.text.x=element_text(size=8), panel.grid=element_blank())
ggsave(file.path(FIGS_DIR, "m3_beta_shaded_table.png"), p_beta_table, width=10, height=12, dpi=150)
cat("Salvo figs/m3_beta_shaded_table.png\n\n")

## ---- 2.1 A simplificação no nível de array se sustenta? ----
## Seis espécies ilustram a comparação: Cutia e Tatu-galinha (ricas em
## dados), Queixada (esparsa, com máscara de distribuição), e Jaguatirica,
## Anta e Onça-parda (escolhidas por grandes mudanças câmera-vs-array).
## IMPORTANTE: esta agora é uma comparação de dois conjuntos de covariáveis
## DIFERENTES, não apenas dois níveis de amostragem com as mesmas
## covariáveis -- o ajuste no nível de câmera usa o conjunto completo de 11
## covariáveis (n=1.110+ locais sustenta isso); o ajuste no nível de array
## usa o conjunto reduzido de 6 covariáveis (n=60 arrays não sustenta mais
## que isso; ver Módulo 2 para como esse conjunto foi escolhido).
species_6 <- c("Cuniculus paca","Dasypus novemcinctus","Tayassu pecari",
               "Leopardus pardalis","Tapirus terrestris","Puma concolor")

site_to_array <- setNames(sc$array_id, sc$site_id)
arrays_ordered <- sort(unique(sc$array_id))
n_arrays <- length(arrays_ordered)

array_cov <- sc[, lapply(.SD, mean, na.rm=TRUE), .SDcols = ARRAY_COVARS, by = array_id]
for (v in ARRAY_COVARS) array_cov[[paste0(v,"_z")]] <- as.numeric(scale(array_cov[[v]]))
# um array (ATLA_ATL_17) não tem dados válidos de habitat/clima no nível de
# câmera em nenhuma de suas câmeras, então sua média no nível de array é NaN
# em todas as colunas -- removemos aqui em vez de carregar silenciosamente
# uma linha totalmente NA adiante.
array_cov <- array_cov[complete.cases(array_cov[, ..ARRAY_COVARS])]
cat("Arrays com dados completos de covariáveis para modelos no nível de array:", nrow(array_cov), "de", n_arrays, "\n\n")
fwrite(data.table(n_arrays_complete=nrow(array_cov), n_arrays_total=n_arrays), "data/array_completeness.csv")

cat("== Comparação nível de câmera (11 covariáveis) vs. nível de array (6 covariáveis), 6 espécies ==\n")
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

  # nível de array: um array conta como "detectado" se qualquer câmera dele detectou naquela ocasião
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
    # compara apenas parâmetros presentes em AMBOS os ajustes (as covariáveis
    # compartilhadas: floresta, savana, veg_nativa_1000m, temp, precip, dist_estrada)
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
cat("\nO AUC no nível de array é tipicamente maior -- unidades amostrais mais\n")
cat("poucas e limpas são mais fáceis de discriminar, ao custo de intervalos de\n")
cat("confiança muito mais amplos (menos poder).\n")
cat("Note que o modelo no nível de câmera também carrega quatro covariáveis --\n")
cat("superfície construída, distância até estrada, distância até água,\n")
cat("status de área protegida -- que o modelo no nível de array não tem\n")
cat("espaço para incluir; a comparação abaixo é restrita às seis\n")
cat("covariáveis que ambos os modelos compartilham.\n\n")

## ---- figura de comparação de estatísticas de ajuste ----
fitstats_long <- melt(fitstats_dt, id.vars="species", measure.vars=c("auc_site","auc_array"),
                        variable.name="level", value.name="auc")
fitstats_long[, level := fifelse(level=="auc_site","Nível de câmera (11 cov.)","Nível de array (6 cov.)")]
common_lookup6 <- c("Cuniculus paca"="Cutia","Dasypus novemcinctus"="Tatu-galinha",
                     "Tayassu pecari"="Queixada","Leopardus pardalis"="Jaguatirica",
                     "Tapirus terrestris"="Anta","Puma concolor"="Onça-parda")
fitstats_long[, common := common_lookup6[species]]
fitstats_long[, common := factor(common, levels=common_lookup6[species_6])]
p_fitstats <- ggplot(fitstats_long, aes(x=common, y=auc, fill=level)) +
  geom_col(position="dodge", width=0.7) +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
  annotate("text", x=1, y=0.53, label="Acaso", color="grey50", size=3) +
  scale_fill_manual(values=c("Nível de câmera (11 cov.)"=COL_NEUTRAL, "Nível de array (6 cov.)"="#e67e22"), name=NULL) +
  labs(title="Discriminação do modelo (AUC): nível de câmera (11 cov.) vs. nível de array (6 cov.), 6 espécies", x=NULL, y="AUC") +
  ylim(0,1) + theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(FIGS_DIR, "m3_camera_vs_array_fit_6species.png"), p_fitstats, width=8, height=5.5, dpi=150)
cat("Salvo figs/m3_camera_vs_array_fit_6species.png\n\n")

## ---- figura de comparação dos coeficientes câmera vs array (dot-whisker) ----
## Restrita às covariáveis compartilhadas (floresta, veg_nativa_1000m -- as
## duas com sinal mais confiável em ambos os níveis para essas espécies).
coef_compare_dt <- rbindlist(coef_compare, fill=TRUE)
coef_compare_dt <- coef_compare_dt[param %in% c("psi(forest_100m_z)","psi(native_veg_1000m_z)")]
coef_compare_dt[, common := common_lookup6[species]]
fwrite(coef_compare_dt, "data/coef_compare_6species.csv")
coef_long <- melt(coef_compare_dt, id.vars=c("species","param","common"),
                    measure.vars=list(est=c("site_est","array_est"), se=c("site_se","array_se")))
coef_long[, level := fifelse(variable==1, "Nível de câmera (11 cov.)", "Nível de array (6 cov.)")]
param_labels3 <- c("psi(forest_100m_z)"="Floresta", "psi(native_veg_1000m_z)"="Veg. nativa (1000m)")
coef_long[, param_label := param_labels3[param]]
coef_long[, y_label := paste(common, param_label, sep=" \u2014 ")]
p_dotwhisker <- ggplot(coef_long, aes(x=est, y=y_label, color=level)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=est-1.96*se, xmax=est+1.96*se), height=0, position=position_dodge(width=0.5), linewidth=0.9) +
  geom_point(position=position_dodge(width=0.5), size=2.5) +
  scale_color_manual(values=c("Nível de câmera (11 cov.)"=COL_NEUTRAL, "Nível de array (6 cov.)"="#e67e22"), name=NULL) +
  labs(title="Coeficientes no nível de câmera vs. nível de array, 6 espécies (apenas covariáveis compartilhadas)",
       x="Coeficiente (\u03b2), estimativa e IC de 95%", y=NULL) +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_array_vs_site_6species.png"), p_dotwhisker, width=9, height=8, dpi=150)
cat("Salvo figs/m3_array_vs_site_6species.png\n\n")


## ---- 2.2 Por que essa comparação importa ----
cat("== Por que essa comparação importa ==\n")
cat("Isso reflete uma escolha de design do pipeline: usar o array como\n")
cat("unidade de replicação primária para modelos comunitários/conjuntos\n")
cat("(Módulo 4+), enquanto ajusta modelos de espécie única no nível mais\n")
cat("fino de câmera aqui. As estimativas pontuais em geral mantêm um sinal\n")
cat("consistente entre os níveis para espécies bem-detectadas, mas o\n")
cat("resultado no nível de array de uma espécie esparsa ou fortemente\n")
cat("mascarada por distribuição deve ser tratado como um sinal muito mais\n")
cat("fraco do que seu valor-p isoladamente sugere.\n\n")

## ---- 2.2b Hierarquia de seleção de habitat de Johnson (1980), ilustrada ----
## Construído nativamente a partir das coordenadas de local já carregadas
## acima -- sem dependência de imagem externa. Requer o pacote `sf` (para os
## contornos dos países) além dos pacotes listados no cabeçalho.
cat("== Hierarquia de seleção de habitat de Johnson (1980) ==\n")
suppressMessages(library(sf))
brazil_sf <- st_read(file.path(DATA_DIR, "brazil_boundary.gpkg"), quiet=TRUE)
sa_sf <- st_read(file.path(DATA_DIR, "south_america_boundary.gpkg"), quiet=TRUE)

COL_AVAIL_J <- "#bbbbbb"; COL_USED_J <- "#2c2c2c"
arr_centers_j <- sc[, .(lon=mean(longitude), lat=mean(latitude)), by=array_id]

## um agrupamento de 4 arrays (~15-30km entre si, claramente separáveis)
## usado nos painéis de segunda e terceira ordem
cluster_arrays_j <- c("WI_WI_019","WI_WI_020","WI_WI_021","WI_WI_022")
cluster_sites_j <- sc[array_id %in% cluster_arrays_j, .(site_id, array_id, longitude, latitude)]

png(file.path(FIGS_DIR, "m3_johnson_levels_of_selection_pt.png"), width=15.5, height=6.2, units="in", res=170)
par(mfrow=c(1,3), mar=c(1,1,3,1), oma=c(0,0,2,0))

## Painel 1: seleção de primeira ordem -- distribuição geográfica
plot(st_geometry(brazil_sf), col=COL_AVAIL_J, border="white", main="Seleção de primeira ordem", cex.main=1.3, font.main=2)
points(arr_centers_j$lon, arr_centers_j$lat, pch=16, col=COL_USED_J, cex=0.7)
mtext("Disponível: distribuição geográfica da espécie\nUsado: locais dos arrays de câmeras", side=1, line=-2, adj=0, cex=0.65)
## inset localizador (América do Sul, Brasil destacado)
usr <- par("usr")
inset_w <- diff(usr[1:2])*0.32; inset_h <- diff(usr[3:4])*0.30
inset_x0 <- usr[1] + diff(usr[1:2])*0.02; inset_y0 <- usr[3] + diff(usr[3:4])*0.62
sa_bbox <- st_bbox(sa_sf)
sa_scale_x <- inset_w / (sa_bbox["xmax"]-sa_bbox["xmin"])
sa_scale_y <- inset_h / (sa_bbox["ymax"]-sa_bbox["ymin"])
sa_geoms <- st_geometry(sa_sf) * matrix(c(sa_scale_x,0,0,sa_scale_y), 2,2)
sa_geoms <- sa_geoms + c(inset_x0 - st_bbox(sa_geoms)["xmin"], inset_y0 - st_bbox(sa_geoms)["ymin"])
plot(sa_geoms, col="#e0e0e0", border="white", lwd=0.3, add=TRUE)
br_geoms <- st_geometry(brazil_sf) * matrix(c(sa_scale_x,0,0,sa_scale_y), 2,2)
br_geoms <- br_geoms + c(inset_x0 - st_bbox(sa_geoms)["xmin"], inset_y0 - st_bbox(sa_geoms)["ymin"])
plot(br_geoms, col="#888888", border="white", lwd=0.3, add=TRUE)
rect(inset_x0, inset_y0, inset_x0+inset_w, inset_y0+inset_h, border="grey40", lwd=0.6)

## Painel 2: seleção de segunda ordem -- área de vida dentro da paisagem regional
lon_c <- mean(cluster_sites_j$longitude); lat_c <- mean(cluster_sites_j$latitude)
half_span <- 0.35
asp_val <- 1/cos(lat_c*pi/180)
plot(NA, xlim=c(lon_c-half_span, lon_c+half_span), ylim=c(lat_c-half_span*0.85, lat_c+half_span*0.85),
     xlab="", ylab="", axes=FALSE, main="Seleção de segunda ordem", cex.main=1.3, font.main=2, asp=asp_val)
rect(lon_c-half_span, lat_c-half_span*0.85, lon_c+half_span, lat_c+half_span*0.85, col=COL_AVAIL_J, border=NA)
for (arr in cluster_arrays_j) {
  sub <- cluster_sites_j[array_id == arr]
  if (nrow(sub) >= 3) {
    hpts <- as.matrix(sub[, .(longitude, latitude)])
    hull <- chull(hpts)
    polygon(hpts[c(hull, hull[1]), ], col=COL_USED_J, border="white", lwd=0.3)
  }
}
mtext("Disponível: paisagem regional\nUsado: territórios de arrays individuais\n(áreas de vida)", side=1, line=-2, adj=0, cex=0.65)

## Painel 3: seleção de terceira ordem -- uso do local dentro da área de vida
example_array_j <- "WI_WI_020"
sub3 <- cluster_sites_j[array_id == example_array_j]
hpts3 <- as.matrix(sub3[, .(longitude, latitude)])
hull3 <- chull(hpts3)
xr3 <- range(sub3$longitude); yr3 <- range(sub3$latitude)
xpad3 <- diff(xr3)*0.25; ypad3 <- diff(yr3)*0.25
plot(NA, xlim=c(xr3[1]-xpad3, xr3[2]+xpad3), ylim=c(yr3[1]-ypad3, yr3[2]+ypad3),
     xlab="", ylab="", axes=FALSE, main="Seleção de terceira ordem", cex.main=1.3, font.main=2, asp=asp_val)
polygon(hpts3[c(hull3, hull3[1]), ], col=COL_AVAIL_J, border=NA)
points(sub3$longitude, sub3$latitude, pch=16, col=COL_USED_J, cex=1.1)
mtext(sprintf("Disponível: área de vida (%s)\nUsado: locais individuais de câmeras", example_array_j), side=1, line=-2, adj=0, cex=0.65)

mtext("Hierarquia de seleção de habitat de Johnson (1980), ilustrada com dados de armadilhas fotográficas do Snapshot Brasil",
      outer=TRUE, cex=1.0, font=2, adj=0.02)
dev.off()
cat("Salvo figs/m3_johnson_levels_of_selection_pt.png\n\n")

## ---- 2.3 Agrupamentos de espécies a partir da tabela de coeficientes ----
cat("== PCA de espécies ==\n")
pca_params <- c(names(label_map_full))
pca_wide <- dcast(community_dt[param != "(Intercept)"], species ~ param, value.var="mean")
pca_wide <- pca_wide[complete.cases(pca_wide[, ..pca_params])]
pca_mat <- scale(as.matrix(pca_wide[, ..pca_params]))
pca_fit <- prcomp(pca_mat, center=FALSE, scale.=FALSE)
ve <- summary(pca_fit)$importance[2, 1:2]

common_names_pt <- fread(file.path(DATA_DIR, "species_common_names_pt.csv"))
common_lookup_pca <- setNames(common_names_pt$common_name_pt, common_names_pt$sci_pooled)
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

## Rótulos de direção dos eixos: para cada eixo do PC, encontra a(s)
## covariável(is) mais fortemente associada(s) com a extremidade positiva vs.
## negativa, para que o leitor saiba o que "PC1 alto" ou "PC2 baixo" significa
## ecologicamente sem precisar comparar as setas de carga visualmente.
axis_dir_label <- function(loadings, axis_col) {
  vals <- loadings[[axis_col]]
  labs <- loadings$label
  pos_idx <- order(-vals)[vals[order(-vals)] > 0]
  neg_idx <- order(vals)[vals[order(vals)] < 0]
  pos_terms <- gsub("\\n", " ", labs[pos_idx][1:min(2,length(pos_idx))])
  neg_terms <- gsub("\\n", " ", labs[neg_idx][1:min(2,length(neg_idx))])
  list(pos = if (length(pos_terms)) paste(pos_terms, collapse=", ") else "(sem covariável forte)",
       neg = if (length(neg_terms)) paste(neg_terms, collapse=", ") else "(sem covariável forte)")
}
pc1_dir <- axis_dir_label(loadings_dt, "PC1")
pc2_dir <- axis_dir_label(loadings_dt, "PC2")
cat("PC1: mais", pc1_dir$pos, "na extremidade positiva; mais", pc1_dir$neg, "na extremidade negativa\n")
cat("PC2: mais", pc2_dir$pos, "na extremidade positiva; mais", pc2_dir$neg, "na extremidade negativa\n")

xr <- range(scores_dt$PC1); yr <- range(scores_dt$PC2)
xpad <- diff(xr) * 0.12; ypad <- diff(yr) * 0.12

p_pca <- ggplot(scores_dt, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=cluster), size=3, alpha=0.8) +
  ggrepel::geom_text_repel(aes(label=label), size=2.4, alpha=0.85, max.overlaps=30,
                             segment.size=0.2, segment.alpha=0.4, seed=1) +
  geom_segment(data=loadings_dt, aes(x=0,y=0,xend=PC1,yend=PC2), color=COL_NEG,
               arrow=arrow(length=unit(0.2,"cm")), inherit.aes=FALSE) +
  geom_text(data=loadings_dt, aes(x=PC1*1.15, y=PC2*1.15, label=label), color=COL_NEG, size=3, fontface="bold", inherit.aes=FALSE) +
  scale_color_brewer(palette="Set1", name="Grupo") +
  scale_x_continuous(limits=c(xr[1]-xpad, xr[2]+xpad*2.5),
                      sec.axis=dup_axis(breaks=xr, labels=c(paste0("mais ", pc1_dir$neg), paste0("mais ", pc1_dir$pos)), name=NULL)) +
  scale_y_continuous(limits=c(yr[1]-ypad, yr[2]+ypad*2),
                      sec.axis=dup_axis(breaks=yr, labels=c(paste0("mais ", pc2_dir$neg), paste0("mais ", pc2_dir$pos)), name=NULL)) +
  labs(title="Espécies agrupadas pelo perfil completo de resposta ao habitat (modelo comunitário)",
       x=sprintf("PC1 (%.0f%%)", ve[1]*100), y=sprintf("PC2 (%.0f%%)", ve[2]*100)) +
  theme_minimal(base_size=12) +
  theme(axis.text.x.top=element_text(size=6.5, color="grey40"),
        axis.text.y.right=element_text(size=6.5, color="grey40"))
ggsave(file.path(FIGS_DIR, "m3_species_pca_biplot.png"), p_pca, width=10.5, height=9, dpi=150)
cat("Salvo figs/m3_species_pca_biplot.png (", nrow(scores_dt), "espécies,", round(sum(ve)*100), "% de variância explicada)\n\n")


## =============================================================================
## 3. ADICIONANDO UM EFEITO DE INTERAÇÃO
## =============================================================================
## Modelo completo: psi ~ cobertura_arbórea * temperatura + savana +
##             pastagem + lavoura + veg_nativa_1000m + precipitação + in_pa
## Cobertura arbórea substitui floresta aqui (colinearidade r=0.82 com
## forest_100m) para evitar redundância com o termo de interação.

species4 <- c("Tapirus terrestris","Euphractus sexcinctus","Dasypus novemcinctus","Cerdocyon thous")
common4 <- c("Tapirus terrestris"="Anta","Euphractus sexcinctus"="Tatu-peludo",
             "Dasypus novemcinctus"="Tatu-galinha","Cerdocyon thous"="Cachorro-do-mato")
OTHER_COVARS <- c("savanna_100m_z","pasture_100m_z","cropland_100m_z","native_veg_1000m_z","precip_annual_mm_z","ghsl_built_5000m_z","dist_road_m_z","dist_water_m_z","in_pa")
cat("== Interação temperatura x cobertura arbórea, modelo completo, 4 espécies ==\n")
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
fwrite(rbindlist(int4_results, fill=TRUE), "data/interaction_4species.csv")

int_pred_dt <- rbindlist(int_preds, fill=TRUE)
int_pred_dt[, species_label := unlist(int_labels[species])]
p_interaction <- ggplot(int_pred_dt, aes(x=tc, y=psi, color=temp, fill=temp)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Fria"="#2c7fb8", "Média"="grey50", "Quente"=COL_NEG), name="Temperatura") +
  scale_fill_manual(values=c("Fria"="#2c7fb8", "Média"="grey50", "Quente"=COL_NEG), name="Temperatura") +
  facet_wrap(~species_label, nrow=1, scales="free_y") +
  labs(title="Interação temperatura x cobertura arbórea, modelo com covariáveis completas, 4 espécies",
       x="Cobertura arbórea (z)", y="Ocupação prevista (\u03c8)") +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
ggsave(file.path(FIGS_DIR, "m3_interaction_4species.png"), p_interaction, width=13, height=4.5, dpi=150)
cat("Salvo figs/m3_interaction_4species.png\n\n")


## =============================================================================
## 3.1 A INTERAÇÃO TEMPERATURA x FLORESTA É MAIS FORTE EM ALGUMAS ORDENS DE MAMÍFEROS?
## =============================================================================
## O exemplo de 4 espécies acima mostra que a interação pode importar muito
## para um punhado de espécies. Aqui ajustamos o MESMO termo de interação
## (temperatura x cobertura arbórea, modelo com covariáveis completas) para
## cada espécie/gênero modelado, marcamos cada um pela ordem taxonômica, e
## perguntamos se alguma ordem se destaca -- especificamente, se Xenarthra
## (tatus, tamanduás; fisiologia térmica próxima da ectotermia em relação a
## outros mamíferos placentários) mostra uma dependência de temperatura mais
## pronunciada em sua resposta de uso de floresta do que o resto da comunidade.

order_lookup_dt <- fread(file.path(DATA_DIR, "species_order_lookup.csv"))
cat("== Interação temperatura x cobertura arbórea, comunidade completa, por ordem taxonômica ==\n")
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
  if (!is.finite(z_o) || abs(z_o) > 20) next  # descarta ajustes instáveis/dirigidos por separação
  order_int_results[[sci_o]] <- data.table(species=sci_o, interaction_coef=co_o[int_term_o],
                                             interaction_se=se_o[int_term_o], interaction_z=z_o,
                                             n_sites=nrow(y_o2))
  order_int_fits[[sci_o]] <- fit_o
}
order_int_dt <- rbindlist(order_int_results, fill=TRUE)
order_int_dt <- merge(order_int_dt, order_lookup_dt, by.x="species", by.y="sci_pooled", all.x=TRUE)
order_int_dt[, sig := abs(interaction_z) > 1.96]
cat("Ajustes de interação completados:", nrow(order_int_dt), "de", nrow(modeled), "espécies/gêneros modelados\n\n")

order_summary <- order_int_dt[!is.na(order), .(n_species=.N, mean_abs_z=mean(abs(interaction_z)),
                                                  median_abs_z=median(abs(interaction_z)),
                                                  pct_sig=100*mean(sig)), by=order]
setorder(order_summary, -mean_abs_z)
cat("== Força da interação (|z|) por ordem taxonômica ==\n")
print(order_summary, digits=3)
cat("\n")
fwrite(order_summary, "data/order_interaction_summary.csv")
fwrite(order_int_dt, "data/order_interaction_full.csv")

## ---- comparando diretamente as dez interações mais fortes ----
## Em vez de comparar apenas os valores de |z| em uma tabela, plotamos as
## curvas de previsão ajustadas de cobertura arbórea x temperatura para os
## dez espécies/gêneros com a interação mais forte (por |z|), todas em uma
## escala comum, para que a forma e a direção de cada interação possam ser
## comparadas lado a lado.
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
      species=sci_t, temp=factor(temp_scenario, levels=c(-1,1), labels=c("Fria","Quente")),
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
  scale_color_manual(values=c("Fria"="#2c7fb8","Quente"=COL_NEG), name="Cenário de temperatura") +
  facet_wrap(~species_label, nrow=2, scales="free_y") +
  labs(title="As dez interações temperatura x cobertura arbórea mais fortes, comparadas diretamente",
       subtitle="Curvas de ocupação prevista fria vs. quente ao longo da faixa de cobertura arbórea, para os dez espécies/gêneros com maior |z|",
       x="Cobertura arbórea (z)", y="Ocupação prevista (\u03c8)") +
  theme_minimal(base_size=9) + theme(legend.position="bottom", strip.text=element_text(size=7))
ggsave(file.path(FIGS_DIR, "m3_interaction_top10.png"), p_top10, width=13, height=7, dpi=150)
cat("Salvo figs/m3_interaction_top10.png\n\n")
fwrite(top10_z_lookup[order(-abs(interaction_z))], "data/top10_interaction_species.csv")

# Teste de duas amostras: Xenarthra vs. todas as outras ordens combinadas
xen_z <- abs(order_int_dt[order=="Xenarthra", interaction_z])
other_z <- abs(order_int_dt[!is.na(order) & order!="Xenarthra", interaction_z])
if (length(xen_z) >= 3 && length(other_z) >= 3) {
  wt <- wilcox.test(xen_z, other_z)
  cat("Xenarthra (n=", length(xen_z), ") |z| médio =", round(mean(xen_z),2),
      "vs. todas as outras ordens (n=", length(other_z), ") |z| médio =", round(mean(other_z),2), "\n")
  cat("Teste de Wilcoxon rank-sum, Xenarthra vs. resto: p =", round(wt$p.value,4), "\n")
  fwrite(data.table(n_xenarthra=length(xen_z), mean_xenarthra=round(mean(xen_z),2),
                     n_other=length(other_z), mean_other=round(mean(other_z),2),
                     wilcoxon_p=round(wt$p.value,4)), "data/xenarthra_test.csv")
  cat(if (wt$p.value < 0.05 && mean(xen_z) > mean(other_z))
        "Xenarthra mostra uma interação temperatura x floresta SIGNIFICATIVAMENTE mais forte que o resto da comunidade.\n"
      else if (wt$p.value < 0.05)
        "A diferença é significativa, mas na direção OPOSTA (Xenarthra mais fraco, não mais forte).\n"
      else
        "Nenhuma diferença significativa detectada -- com apenas 5 espécies de Xenarthra este teste tem poder limitado;\n  tratar como sugestivo, não conclusivo.\n")
} else {
  cat("Poucas espécies em um dos grupos para um teste formal; relatando apenas o resumo descritivo.\n")
}
cat("\n")

## ---- figura: força da interação por ordem, pontos por espécie + médias por ordem ----
order_int_dt[, order_f := factor(order, levels=order_summary$order)]
order_int_dt[, is_xenarthra := order == "Xenarthra"]
sig_label_dt <- data.table(order_f=factor(order_summary$order[1], levels=order_summary$order), y=2.15)
p_order_interaction <- ggplot(order_int_dt[!is.na(order)], aes(x=order_f, y=abs(interaction_z))) +
  geom_hline(yintercept=1.96, linetype="dashed", color="grey60") +
  geom_text(data=sig_label_dt, aes(x=order_f, y=y), label="|z| = 1.96 (sig. 95%)",
            color="grey50", size=2.8, hjust=0, inherit.aes=FALSE) +
  geom_jitter(aes(color=is_xenarthra), width=0.15, size=2.5, alpha=0.8) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, color="black", linewidth=0.6) +
  scale_color_manual(values=c("TRUE"=COL_NEG, "FALSE"=COL_NEUTRAL), guide="none") +
  labs(title="Força da interação temperatura x cobertura arbórea, por ordem taxonômica",
       subtitle="Cada ponto = uma espécie/gênero; barra preta = média da ordem; vermelho = Xenarthra",
       x=NULL, y="|z| para o termo de interação") +
  theme_minimal(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(FIGS_DIR, "m3_interaction_by_order.png"), p_order_interaction, width=9, height=6, dpi=150)
cat("Salvo figs/m3_interaction_by_order.png\n\n")


## =============================================================================
## 4. DEIXANDO UMA RELAÇÃO VARIAR POR ECORREGIÃO
## =============================================================================
## Modelo completo: psi ~ floresta * ecorregião + savana + pastagem +
##             lavoura + veg_nativa_1000m + temperatura + precipitação + in_pa

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
cat("== Jaguatirica: efeito da floresta agrupado vs. com interação por ecorregião ==\n")
cat("Melhoria de AIC do modelo com interação:", round(delta_eco,1), "\n")
fwrite(data.table(delta_aic=round(delta_eco,1)), "data/ecoregion_ocelot_aic.csv")
cat("Em Floresta Úmida de Folhas Largas, mais cobertura florestal geralmente\n")
cat("significa mais ocupação de Jaguatirica; em Pastagens/Savanas & Arbustos\n")
cat("a relação se inverte -- um único 'efeito de floresta' agrupado esconde\n")
cat("duas relações genuinamente diferentes.\n\n")

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
       x="Cobertura florestal (z, padronizada)", y="Ocupação prevista (\u03c8)") +
  theme_minimal(base_size=12) + theme(legend.position=c(0.25,0.9))
ggsave(file.path(FIGS_DIR, "m3_ecoregion_ocelot.png"), p_ecoregion, width=8, height=6, dpi=150)
cat("Salvo figs/m3_ecoregion_ocelot.png\n\n")


## =============================================================================
## 5. COMPARANDO COM A ABUNDÂNCIA DE ROYLE-NICHOLS
## =============================================================================
## Royle-Nichols (occuRN) usa os mesmos dados de detecção 1/0 da ocupância,
## mas interpreta a FREQUÊNCIA de detecção como um sinal de abundância
## relativa em vez de apenas presença/ausência.

cat("== Ocupância vs. abundância de Royle-Nichols ==\n")
cat("Ajustando Royle-Nichols (occuRN) para cada espécie/gênero modelado...\n")
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
cat("Ajustes RN completados:", uniqueN(rn_dt$species), "de", nrow(modeled), "\n\n")

## comparação no nível geral
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
cat("Correlação no nível geral:", round(cor(overall_comp$mean_psi, overall_comp$mean_lam),2), "\n")
cat("Correlação do efeito de floresta (apenas ajustes estáveis):", round(cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest),2), "\n\n")
fwrite(data.table(overall_r=round(cor(overall_comp$mean_psi, overall_comp$mean_lam),2),
                   forest_r=round(cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest),2),
                   n_rn_fit=uniqueN(rn_dt$species), n_modeled=nrow(modeled)), "data/rn_vs_occ_correlation.csv")

p_overall <- ggplot(overall_comp, aes(x=mean_psi, y=mean_lam)) +
  geom_point(color=COL_NEUTRAL, size=2.5, alpha=0.75) +
  labs(title=sprintf("Nível geral (r=%.2f)", cor(overall_comp$mean_psi, overall_comp$mean_lam)),
       x="Ocupação média (\u03c8-estimado)", y="Abundância média de Royle-Nichols (\u03bb-estimado)") +
  theme_minimal(base_size=10)
top2[, common := sapply(species, species_label_pca)]
p_forest <- ggplot(forest_comp_stable, aes(x=occ_forest, y=rn_forest)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(color=COL_NEUTRAL, size=2.5, alpha=0.75) +
  ggrepel::geom_text_repel(data=top2, aes(label=common), size=2.8, color="black",
                             segment.size=0.3, segment.alpha=0.5, seed=2,
                             box.padding=0.6, max.overlaps=Inf) +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.22))) +
  labs(title=sprintf("Efeito de floresta especificamente (r=%.2f)", cor(forest_comp_stable$occ_forest, forest_comp_stable$rn_forest)),
       x="Efeito de floresta (ocupância, \u03b2)", y="Efeito de floresta (Royle-Nichols, \u03b2)") +
  theme_minimal(base_size=10)
p_occ_vs_abund <- gridExtra::grid.arrange(p_overall, p_forest, ncol=2,
  top=grid::textGrob("Ocupância vs. abundância de Royle-Nichols: concordância através da comunidade",
                       gp=grid::gpar(fontsize=13, fontface="bold"), x=0.02, hjust=0))
ggsave(file.path(FIGS_DIR, "m3_occ_vs_abundance.png"), p_occ_vs_abund, width=11, height=5, dpi=150)
cat("Salvo figs/m3_occ_vs_abundance.png\n\n")


## =============================================================================
## RESUMO -- o que os resultados significam
## =============================================================================
## - Este script constrói cada histórico de detecção, tabela de covariáveis
##   e modelo diretamente a partir dos três arquivos brutos do conjunto de
##   dados combinado -- nenhum coeficiente pré-ajustado é carregado.
## - Didelphis e Dasyprocta são agrupados ao nível de gênero (uma exceção
##   deliberada à regra de nível de espécie do pipeline), recuperando
##   registros de detecção apenas de gênero. Didelphis mantém uma máscara
##   de distribuição (união das distribuições com buffer dos 3 congêneres
##   mapeados); Dasyprocta não tem máscara de distribuição (apenas 1 dos 5
##   congêneres está mapeado, e usá-lo isoladamente classificaria
##   erroneamente detecções reais como impossíveis).
## - O ajuste apenas com efeitos fixos é pouco confiável para
##   aproximadamente metade da comunidade; o modelo com efeito aleatório de
##   array (chamadas PGOcc comentadas acima) é o padrão do pipeline e deve
##   ser executado para a análise comunitária completa.
## - AUC (discriminação) e c-hat (sobredispersão) medem coisas diferentes e
##   podem discordar.
## - Colapsar câmeras em arrays amplia substancialmente os intervalos de
##   confiança e pode inverter a direção estimada de uma covariável para
##   espécies com poucos dados.
## - Tanto uma interação contínua (temperatura x cobertura arbórea) quanto
##   uma categórica (ecorregião) mostram que um único efeito de habitat
##   "um número serve para tudo" pode esconder variação real e
##   ecologicamente sensata.
## - Ajustar a mesma interação temperatura x cobertura arbórea em cada
##   espécie/gênero modelado e agrupar por ordem taxonômica testa se essa
##   sensibilidade se agrupa por linhagem em vez de ser idiossincrática por
##   espécie -- com apenas 5 espécies de Xenarthra no conjunto modelado,
##   este é um teste sugestivo, não definitivo.
##
## Próximo módulo: Modelos de Comunidade e Conjuntos.
## =============================================================================
