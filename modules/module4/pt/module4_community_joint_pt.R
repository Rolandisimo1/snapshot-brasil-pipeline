## =============================================================================
## Módulo 4 — Modelos de Comunidade e Conjuntos
## Pipeline de câmeras-trap multi-rede do Snapshot Brasil
## =============================================================================
##
## 1. AUTOCONTIDO: este script parte dos três arquivos brutos de dataset combinado
## 2. (implantações, detecções de mamíferos, covariáveis de sensoriamento remoto) mais
## 3. quatro tabelas de referência do Módulo 2/3, e constrói tudo do
## 4. zero -- históricos de detecção para 44 espécies modeladas, a máscara de distribuição,
## 5. e todos os ajustes de modelo (ocupação de comunidade com uma
## 6. meta-regressão de traço de espécie, GJAM, Royle-Nichols, ocupação de espécie única, e uma
## 7. linha de base ingênua de regressão logística). Coloque este script na mesma pasta
## 8. dos sete arquivos de dados listados abaixo e execute-o de cima a baixo.
##
## 9. Arquivos de entrada necessários (em DATA_DIR, padrão "data" ao lado deste script):
## 10. final_deployments.csv         -- uma linha por implantação de câmera
## 11. final_detections_mammals.csv  -- uma linha por evento de detecção independente
## 12. final_covariates.csv          -- uma linha por local de câmera, covariáveis de sensoriamento remoto
## 13. array_covariates.csv          -- uma linha por array de câmeras, o conjunto reduzido de 6 covariáveis
## 14. ARRAY_SET (design em nível de array do Módulo 2/3)
## 15. species_list.json             -- as 44 espécies/gêneros modelados neste módulo
## 16. (>=10 detecções, >=20 locais, conforme a regra do Módulo 3)
## 17. species_range_mask.csv        -- máscara de distribuição local x espécie da IUCN (1 = candidato,
## 18. 0 = zero estrutural); saída do próprio montador de
## 19. máscara de distribuição do Módulo 4 (ver notas de preparo do Módulo 4)
## 20. species_traits.csv            -- massa corporal logarítmica por espécie (padronizada em z) e classe
## 21. de dieta (herbívoro/onívoro, carnívoro = referência),
## 22. do COMBINE (Soria et al. 2021) associado via a
## 23. taxonomia da Mammal Diversity Database
##
## 24. Os últimos quatro arquivos são saídas de referência do trabalho de seleção de covariáveis
## 25. do Módulo 2, da regra de seleção de espécies do Módulo 3, do montador de máscara de distribuição do Módulo 4, e
## 26. de uma associação única com banco de dados de traços, e não podem ser reconstruídos a partir dos três
## 27. arquivos de dados principais isoladamente -- eles são fornecidos junto com este script como tabelas
## 28. de referência fixas, da mesma forma que o Módulo 3 fornece suas quatro tabelas de referência de máscara de distribuição/taxonomia.
##
## 29. Cinco métodos são ajustados e comparados, todos no MESMO CAMERA_SET de 11 covariáveis
## 30. (floresta, savana, pastagem, cultivo agrícola, vegetação nativa, temperatura,
## 31. precipitação, status de área protegida, desenvolvimento, distância à estrada,
## 32. distância à água):
## 33. 1. Ocupação de comunidade (Dorazio-Royle) -- todas as 44 espécies conjuntamente, com
## 34. pooling parcial entre espécies e uma meta-regressão de traço de espécie
## 35. nos interceptos de ocupação/detecção.
## 36. 2. GJAM (modelo conjunto de distribuição de espécies) -- co-ocorrência residual
## 37. após o habitat ser contabilizado.
## 38. 3. Modelo de heterogeneidade induzida por abundância Royle-Nichols (occuRN), ajustado
## 39. por espécie.
## 40. 4. Ocupação de espécie única (spOccupancy::PGOcc) com um efeito aleatório
## 41. de array, ajustado por espécie, com máscara de distribuição.
## 42. 5. Uma regressão logística ingênua (sem modelo de detecção) -- mostrada
## 43. deliberadamente para ilustrar os problemas de separação que ela produz em
## 44. covariáveis esparsas de cobertura do solo.
##
## 45. Requer: data.table, ggplot2, ggrepel, jagsUI (rjags/JAGS instalado
## 46. separadamente), gjam, unmarked, spOccupancy.

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(jagsUI)
  library(gjam)
  library(unmarked)
  library(spOccupancy)
})

# Linha traduzida: ---- 0. Caminhos ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"
stopifnot(dir.exists(DATA_DIR))
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR)

# Verde = positivo/alto, vermelho = negativo/baixo
COL_POS <- "#1a9850"; COL_NEG <- "#d73027"; COL_NEUTRAL <- "#999999"

## Conjunto de covariáveis no nível de câmera (11) -- usado para todos os métodos deste módulo.
CAMERA_COVARS <- c("forest_100m","savanna_100m","pasture_100m","cropland_100m",
                    "native_veg_1000m","temp_mean_C","precip_annual_mm",
                    "ghsl_built_5000m","dist_road_m","dist_water_m")
# in_pa (status de área protegida, 0/1) é anexado separadamente abaixo, seguindo
# a convenção do Módulo 3 (é binário e não é padronizado em z-score).

OCC_DAYS <- 7  # duração da ocasião em dias, seguindo a convenção de intervalo de independência do pipeline

set.seed(1)


## =============================================================================
## PASSO 1: CARREGAR OS DADOS BRUTOS E CONSTRUIR AS COVARIÁVEIS EM NÍVEL DE SÍTIO
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

sc <- merge(site_tbl, cov[, -"array_id"], by = "site_id", all.x = TRUE)  # "array_id" já presente em site_tbl
setorder(sc, site_id)

# 1. Excluir listwise os sites com uma ou mais das 11 covariáveis em nível de câmera ausentes
# 2. ANTES de construir qualquer coisa a jusante -- todo método neste módulo (e
# 3. toda consulta de referência distribuída junto com este script) é construído sobre este
# 4. mesmo conjunto de sites candidatos, então a exclusão precisa acontecer aqui, não apenas ser
# 5. relatada e silenciosamente ignorada.
n_missing_cov <- sum(!complete.cases(sc[, ..CAMERA_COVARS]))
cat("Sites missing one or more covariates (excluded, listwise deletion):", n_missing_cov, "\n\n")
sc <- sc[complete.cases(sc[, ..CAMERA_COVARS])]

for (v in CAMERA_COVARS) sc[[paste0(v, "_z")]] <- as.numeric(scale(sc[[v]]))
CAMERA_COVARS_Z <- c(paste0(CAMERA_COVARS, "_z"), "in_pa")  # Line 1: in_pa é 0/1, não padronizado (z-score)

arrays_sorted <- sort(unique(sc$array_id))
sc[, array_num := match(array_id, arrays_sorted)]


## =============================================================================
## PASSO 2: A COMUNIDADE MODELADA DE 44 ESPÉCIES E SEUS HISTÓRICOS DE DETECÇÃO
## =============================================================================
## A seleção de espécies segue a regra do Módulo 3 (>=10 detecções, >=20 sites,
## excluindo táxons domésticos/humanos e registros não resolvidos apenas em nível de gênero); a
## lista exata de 44 espécies que este módulo ajusta é distribuída como species_list.json
## para que todo leitor veja a mesma comunidade, sem precisar rederivar o
## critério de corte da seleção.

species <- fromJSON <- jsonlite::fromJSON(file.path(DATA_DIR, "species_list.json"))
cat("Modeled species/genera:", length(species), "\n\n")

# 1. Uma janela de amostragem ativa por site: a data de início mais antiga e a data de término mais recente entre
# 2. as implantações daquele site (um site pode ter várias implantações ao longo do tempo).
# 3. As ocasiões são ancoradas na PRÓPRIA DATA DE INÍCIO DE CADA SITE -- e não em uma única data
# 4. global -- já que os sites abrangem períodos de calendário muito diferentes (este conjunto de dados
# 5. reúne múltiplas redes coletadas ao longo de mais de uma década).
site_window <- dep[, .(site_start = min(start_dt, na.rm=TRUE), site_end = max(end_dt, na.rm=TRUE)), by = site_id]
max_days <- as.integer(max(site_window$site_end - site_window$site_start, na.rm = TRUE))
max_occ <- ceiling(max_days / OCC_DAYS)
site_ids_ordered <- sort(unique(sc$site_id))
n_sites <- length(site_ids_ordered)
site_window <- site_window[match(site_ids_ordered, site_id)]
site_start_lookup <- setNames(site_window$site_start, site_window$site_id)

# 1. matriz base do levantamento: quais ocasiões foram efetivamente amostradas em cada sítio
# 2. (um sítio pode ter LACUNAS entre as instalações -- apenas as ocasiões efetivamente
# 3. cobertas por uma janela de instalação são marcadas como amostradas).
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
# 1. Restringir às detecções que caem dentro da PRÓPRIA JANELA DE AMOSTRAGEM DO SITE
# 2. (site_start até site_end, inclusive) -- corresponde exatamente à convenção
# 3. original do pipeline. Sem esse filtro, uma detecção registrada após o
# 4. término da implantação de um site com janela curta (mas ainda dentro do limite
# 5. max_occ do conjunto de dados) seria admitida indevidamente através da janela mais longa de outro site.
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

# 1. Matrizes ysum/K: contagem de detecções colapsada por ocasião e ocasiões amostradas
# 2. contagem por site x espécie, seguindo a convenção do Módulo 3/6 ao longo de todo o texto.
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
## 1. ETAPA 3: MÁSCARA DE DISTRIBUIÇÃO E TRAÇOS DAS ESPÉCIES (tabelas de referência)
## =============================================================================
## 2. A máscara de distribuição (1 = local é candidato, 0 = zero estrutural -- fora
## 3. da distribuição da IUCN com buffer da espécie e nunca detectada ali) é uma tabela
## 4. 0/1 por espécie/por local construída uma única vez pelo próprio montador de máscara
## 5. de distribuição do Módulo 4 (ver as notas de preparo do Módulo 4) e distribuída como uma referência fixa,
## 6. da mesma forma que o Módulo 3 distribui suas tabelas de exclusão.

rng_dt <- fread(file.path(DATA_DIR, "species_range_mask.csv"))
rng <- as.matrix(rng_dt[match(site_ids_ordered, rng_dt$site_id), ..species])
rng[is.na(rng)] <- 1L  # Espécies sem distribuição mapeada assumem por padrão "todo lugar é um candidato"
storage.mode(rng) <- "integer"
cat("Range-masked (structural-zero) site x species cells:", sum(rng == 0), "of", length(rng), "\n\n")

## 1. Traços das espécies: log da massa corporal (padronizado por z-score) e classe de dieta
## 2. (herbívoro/onívoro, carnívoro = nível de referência), do banco de dados de traços de mamíferos COMBINE
## 3. (Soria et al. 2021), pareado com a lista de espécies deste pipeline
## 4. por meio da taxonomia do Mammal Diversity Database.
traits <- fread(file.path(DATA_DIR, "species_traits.csv"))
trait_mat <- as.matrix(traits[match(species, traits$species),
                               .(log_mass_z, diet_herbivore, diet_omnivore)])
stopifnot(nrow(trait_mat) == S, !anyNA(trait_mat))


## =============================================================================
## PASSO 4: MÉTODO 1 -- OCUPAÇÃO DE COMUNIDADE (DORAZIO-ROYLE)
## =============================================================================
## Todas as 44 espécies ajustadas conjuntamente: os interceptos de ocupação/detecção em nível de espécie
## e os betas de habitat são parcialmente agrupados em direção a uma distribuição de comunidade
## (uma espécie com poucos dados é estabilizada pelo que congêneres com muitos dados mostram), e
## uma meta-regressão de traços de espécie prevê parte do intercepto basal de cada espécie
## a partir de sua massa corporal e classe de dieta -- independentemente do habitat.
## A máscara de distribuição entra diretamente no estado latente: z[i,s] só pode ser 1
## em sítios dentro da distribuição da espécie.

X <- as.matrix(sc[, ..CAMERA_COVARS_Z])
arr <- sc$array_num
N <- nrow(ysum); C <- ncol(X); A <- max(arr)
cat(sprintf("Community model: N=%d sites, S=%d species, C=%d covariates, A=%d arrays\n", N, S, C, A))

zinit <- matrix(1, N, S); zinit[rng == 0] <- 0; zinit[ysum > 0] <- 1

community_jags_code <- "model {
  ## Hiperpriors da comunidade (agrupamento parcial entre espécies)
  mu.lpsi ~ dnorm(0, 0.1); sd.lpsi ~ dunif(0,5); tau.lpsi <- 1/(sd.lpsi*sd.lpsi)
  mu.lp   ~ dnorm(0, 0.1); sd.lp   ~ dunif(0,5); tau.lp   <- 1/(sd.lp*sd.lp)
  for (c in 1:C) {
    mu.beta[c] ~ dnorm(0, 0.1)
    sd.beta[c] ~ dunif(0, 5)
    tau.beta[c] <- 1/(sd.beta[c]*sd.beta[c])
  }
  sd.arr ~ dunif(0,5); tau.arr <- 1/(sd.arr*sd.arr)
  for (a in 1:A) { eta[a] ~ dnorm(0, tau.arr) }

  ## efeitos de traços de espécie (meta-regressão nos interceptos em nível de espécie)
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

## NOTA: este ajuste leva ~8 horas (44 espécies x 11 covariáveis x modelo de traço,
## MCMC single-threaded, 15000 iterações x 3 cadeias). Defina RUN_LIVE <- TRUE
## abaixo para reajustar do zero; caso contrário, este script carrega as
## tabelas de coeficientes salvas que o acompanham (community_hyperparameters.csv,
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
## 1. PASSO 5: MÉTODO 2 -- GJAM (MODELO CONJUNTO DE DISTRIBUIÇÃO DE ESPÉCIES)
## =============================================================================
## 2. O GJAM modela as 44 espécies em conjunto e, após ajustar os efeitos de hábitat,
## 3. reporta a correlação RESIDUAL entre espécies -- co-ocorrência que
## 4. o hábitat sozinho não explica. inrangeFrac (a fração das 44
## 5. espécies cuja distribuição mapeada inclui este local) entra como uma covariável
## 6. adicional para que a co-ocorrência causada pela distribuição não seja interpretada erroneamente como sinal biótico.

inrange_frac <- rowMeans(rng)
xdata_gjam <- data.table(
  forestZ=sc$forest_100m_z, savannaZ=sc$savanna_100m_z, pastureZ=sc$pasture_100m_z,
  croplandZ=sc$cropland_100m_z, nvegZ=sc$native_veg_1000m_z, tempZ=sc$temp_mean_C_z,
  precipZ=sc$precip_annual_mm_z, inPA=sc$in_pa, builtZ=sc$ghsl_built_5000m_z,
  distroadZ=sc$dist_road_m_z, distwaterZ=sc$dist_water_m_z, inrangeFrac=inrange_frac)
ydata_gjam <- as.data.frame(ysum); colnames(ydata_gjam) <- make.names(species)
site_days <- sc$camera_days  # Esforço amostral por sítio (proxy: total de dias-câmera)
effMat <- matrix(site_days, nrow=N, ncol=S)
eList <- list(columns=1:S, values=effMat)

RUN_LIVE_GJAM <- FALSE  # GJAM é comparativamente rápido (~15-20 min); defina como TRUE para reajustar

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
## 1. PASSO 6: MÉTODO 3 -- MODELO DE ABUNDÂNCIA ROYLE-NICHOLS (occuRN)
## =============================================================================
## 2. O occuRN trata a heterogeneidade de detecção como uma assinatura de abundância
## 3. não observada: lambda ~ todas as 11 covariáveis, p ~ apenas intercepto, ajustado por espécie.

RUN_LIVE_RN <- FALSE  # 4. ajustes por espécie, poucos segundos cada; definir como TRUE para reajustar

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
## 1. PASSO 7: MÉTODO 4 -- OCUPAÇÃO DE ESPÉCIE ÚNICA (spOccupancy::PGOcc)
## =============================================================================
## 2. O método padrão de espécie única do pipeline (Módulo 3), reajustado aqui
## 3. no conjunto de 11 covariáveis com a máscara de distribuição aplicada (sítios fora
## 4. da distribuição removidos da verossimilhança daquela espécie) e um efeito aleatório de array.

RUN_LIVE_PGOCC <- FALSE  # 5. MCMC por espécie, alguns segundos a ~1 min cada; defina TRUE para reajustar

if (RUN_LIVE_PGOCC) {
  pgocc_results <- list()
  for (s in seq_len(S)) {
    sci <- species[s]
    y <- build_dethist(sci)
    keep_sp <- rng[, sci] == 1  # Sites fora da área de distribuição desta espécie são removidos
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
## 1. PASSO 8: MÉTODO 5 -- REGRESSÃO LOGÍSTICA INGÊNUA (sem modelo de detecção)
## =============================================================================
## 2. Uma referência deliberadamente ingênua: presença/ausência por local, sem submodelo
## 3. de detecção. Apresentada para ilustrar os problemas de separação perfeita que produz
## 4. em covariáveis esparsas de uso do solo -- o ponto pedagógico de incluí-la.

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
    degenerate <- is.na(se) || se > 20  # Ocupação/quasi-separação perfeita
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
## PASSO 9: COMPARAÇÃO ENTRE MÉTODOS
## =============================================================================
## Reúna os cinco métodos em uma tabela única no formato longo, usando o
## design compartilhado de 11 covariáveis, e calcule as correlações par a par dos coeficientes.

common_names <- fread(file.path(DATA_DIR, "species_common_names.csv"))
common_map <- setNames(common_names$common, common_names$species)

# 1. O PGOcc usa nomes abreviados de parâmetros (forest_z, nveg_z, built_z, distroad_z,
# 2. distwater_z) que não correspondem aos nomes brutos das covariáveis por simples remoção de sufixo;
# 3. mapeie-os explicitamente para os nomes canônicos de CAMERA_COVARS.
pgocc_name_map <- c(forest_z="forest_100m", savanna_z="savanna_100m", pasture_z="pasture_100m",
                     cropland_z="cropland_100m", nveg_z="native_veg_1000m", temp_z="temp_mean_C",
                     precip_z="precip_annual_mm", in_pa="in_pa", built_z="ghsl_built_5000m",
                     distroad_z="dist_road_m", distwater_z="dist_water_m")

# 1. Os nomes de preditores do GJAM carregam um sufixo "Z" (forestZ, savannaZ, ...) e sua
# 2. própria transformação dos nomes de espécies (make.names sobre o nome científico); mapeie ambos
# 3. de volta para os nomes canônicos de CAMERA_COVARS / espécies.
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
gjam_bt_hab[, species_clean := species]  # Dot-separated no nome científico já presente neste arquivo
species_dotmap <- setNames(species, make.names(species))
gjam_bt_hab[, species_clean := species_dotmap[species]]

# 1. significância ("sig", uma flag 0/1) calculada de forma uniforme por método: o
# 2. intervalo de 95% do coeficiente exclui zero (métodos Bayesianos, GJAM), ou
# 3. p < 0,05 (GLM), ou |beta/se| > 1.96 (Royle-Nichols, que reporta apenas o SE)?
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
## ETAPA 10: FIGURAS
## =============================================================================

## ---- 4.2: hiperparâmetro da comunidade + efeitos de traço ----
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

## ---- 4.3: respostas em nível de espécie, um painel por covariável ----
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

## ---- gráficos de dispersão bivariados de pares de covariáveis (4.3): dispersão de espécies codificada por significância ----
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

## ---- 5.1: Heatmap de correlação residual de co-ocorrência do GJAM ----
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

## 4.2: Gráfico de barras de importância dos preditores do GJAM ----
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

## 1. ---- 6.1: importância dos preditores entre métodos ----
## 2. Métrica de importância sem unidade (segue a mesma definição usada no notebook): a
## 3. porcentagem das 44 espécies com um coeficiente credível/significativamente diferente de zero
## 4. para cada covariável, segundo o teste de significância próprio de cada método.
imp_by_method <- all_methods[, .(pct_sig=100*mean(sig, na.rm=TRUE)), by=.(method, covariate)]
p_cross_imp <- ggplot(imp_by_method, aes(x=covariate, y=pct_sig, fill=method)) +
  geom_col(position="dodge") + coord_flip() +
  labs(x=NULL, y="% of 44 species with a significant coefficient", title="Predictor importance across methods") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR,"m4_methods_importance.png"), p_cross_imp, width=10, height=7, dpi=150)
cat("Saved figs/m4_methods_importance.png\n")

cat("\n=== Module 4 script complete: all methods fit/loaded, all figures saved to", FIGS_DIR, "===\n")
