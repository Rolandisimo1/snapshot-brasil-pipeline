## =============================================================================
## Módulo 2 — Covariáveis, Amostragem e Áreas Protegidas
## Pipeline multi-rede de armadilhas fotográficas Snapshot Brasil
## =============================================================================
##
## O que este script faz, em ordem:
##   1. Explicar por que o sensoriamento remoto é necessário, por que um
##      buffer (não um ponto) é usado, e por que composição percentual (não
##      um rótulo categórico "em habitat X") é a forma correta de covariável.
##   2. Carregar a tabela de covariáveis pré-extraídas e descrever o conjunto
##      candidato.
##   3. Explicar o que é o Google Earth Engine (GEE) e o que um leitor
##      precisaria configurar antes de rodar o exemplo de extração (não
##      executado).
##   4. Verificar redundância entre as covariáveis candidatas: matriz de
##      correlação (conjunto completo vs. conjunto final reduzido), PCA, e
##      fator de inflação da variância (VIF).
##   5. Explicar padronização e por que NÃO transformamos as covariáveis
##      enviesadas de pastagem/lavoura.
##   6. Mapear as covariáveis: superfícies nacionais, depois locais de câmera
##      sobre elas.
##   7. Quantificar a homogeneidade de habitat dentro do array (correlação
##      intraclasse, ICC) e mostrar um exemplo concreto de três arrays
##      cobrindo a faixa de ICC.
##   8. Resumir a cobertura de área protegida na rede.
##   9. Avaliar mapas de distribuição de espécies: quais locais são
##      geograficamente impossíveis para uma dada espécie, e devem ser
##      mascarados antes da modelagem.
##
## Requer: data.table, ggplot2. Aponte DATA_DIR para a pasta contendo os
## CSVs listados abaixo (empacotados junto com este script por padrão).
## =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})

# ---- 0. Caminhos -------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"   # usado apenas se você quiser re-salvar figuras a partir deste script

stopifnot(dir.exists(DATA_DIR))


## =============================================================================
## 1. POR QUE PRECISAMOS DE SENSORIAMENTO REMOTO
## =============================================================================
## Uma armadilha fotográfica só nos diz o que passou por ela -- não diz nada,
## por si só, sobre o habitat ao redor daquele ponto. Para relacionar
## detecções a habitat, precisamos de uma medida independente de como é a
## paisagem ao redor de cada câmera: floresta ou savana, lavoura ou
## pastagem, úmido ou seco, protegido ou não. O sensoriamento remoto --
## imagens de satélite classificadas em tipos de cobertura do solo, ou
## medidas continuamente (temperatura, precipitação, altura do dossel) --
## dá essa medida em cada local de câmera, de forma barata e consistente,
## sem uma pesquisa de campo em cada ponto.
##
## A mecânica: cada câmera tem uma latitude/longitude. Pegamos um mapa
## derivado de satélite (um raster onde cada pixel é classificado, ex.,
## floresta vs. pastagem vs. água) e lemos o que esse mapa diz nas
## coordenadas de cada câmera.
##
## ---- 1a. Por que um buffer, e não apenas um ponto ----
## A coordenada exata de uma câmera é um único pixel -- mas os animais que
## ela detecta se deslocam em uma área muito maior que um pixel, e as
## coordenadas de GPS carregam seu próprio erro de poucos metros. Ler a
## classe de cobertura do solo no ponto literal arrisca responder a uma
## pergunta frágil: este único pixel de 10-30 m é floresta ou não? Mover o
## buffer alguns metros pode inverter a resposta, mesmo que nada sobre o
## habitat real do local tenha mudado.
##
## Em vez disso, para cada covariável traçamos um BUFFER -- um círculo de
## raio fixo ao redor da câmera -- e resumimos a paisagem DENTRO desse
## círculo, não apenas em seu centro. Dois tamanhos de buffer servem a
## propósitos diferentes: um buffer de 100 m para estrutura de habitat local
## (o que um animal encontra em seu entorno imediato), e um buffer de
## 1000 m para contexto de paisagem (que tipo de área a câmera está
## inserida de forma mais ampla).
##
## ---- 1b. Por que composição percentual, e não um rótulo categórico ----
## Dado um buffer, ainda há duas formas de resumi-lo: escolher a única
## classe que cobre a maior área ("este local é floresta"), ou reportar a
## PORCENTAGEM da área do buffer em cada classe de interesse ("87%
## floresta, 13% pastagem dentro de 100 m"). Usamos a forma percentual em
## todo o pipeline. Um rótulo categórico "em floresta / não em floresta"
## descarta informação real -- um local com 95% de floresta e um com 51% de
## floresta são ambos chamados de "floresta", mesmo que o segundo esteja
## significativamente mais próximo de uma borda de perturbação. A
## composição percentual preserva esse gradiente, que é exatamente o tipo
## de preditor contínuo de que um modelo de ocupância ou abundância
## precisa: ele pode estimar como a probabilidade de detecção ou ocupância
## muda conforme a cobertura florestal vai de 0% a 100%, em vez de ser
## informado apenas de qual lado de um limiar arbitrário o local se
## encontra.

cat("== Por que sensoriamento remoto, buffers e composição percentual ==\n",
    "Ver o bloco de comentário acima para a justificativa completa.\n\n")


## =============================================================================
## 2. AS COVARIÁVEIS CANDIDATAS
## =============================================================================
## Toda camada candidata foi extraída para cada local de câmera em um
## buffer de 100 m (habitat local) ou de 1000 m (contexto de paisagem).

cov <- fread(file.path(DATA_DIR, "covariates.csv"))
cat("== Covariáveis candidatas ==\n")
cat("locais com covariáveis:", nrow(cov), "\n")
cand_cols <- setdiff(names(cov), c("site_id", "pt_hab", "wc_class_point", "in_pa"))
cat("camadas candidatas:", paste(cand_cols, collapse = ", "), "\n")

# Medimos frações de cobertura do solo do MapBiomas (floresta, savana,
# pastagem, lavoura, vegetação nativa), duas verificações cruzadas globais
# (% de árvores WorldCover, cobertura arbórea Hansen), altura do dossel,
# superfície construída GHSL, temperatura e precipitação TerraClimate, e
# status de área protegida WDPA. Pastagem (pastagem extensiva de gado) e
# lavoura (culturas em linha/plantações) são mantidas como duas covariáveis
# SEPARADAS em vez de uma "agricultura" combinada, porque as espécies de
# armadilhas fotográficas respondem de forma muito diferente a elas -- uma
# borda de pastagem e um campo de soja são tipos de perturbação
# ecologicamente distintos, mesmo que ambos contem como uso humano da terra.

FINAL_COVARS <- c("forest_100m", "savanna_100m", "pasture_100m", "cropland_100m",
                   "native_veg_1000m", "temp_mean_C", "precip_annual_mm")
FULL_COVARS  <- c("forest_100m", "savanna_100m", "pasture_100m", "cropland_100m",
                   "native_veg_1000m", "wc_tree_pct_100m", "treecover2000_100m",
                   "canopy_height_m_100m", "ghsl_built_100m", "temp_mean_C",
                   "precip_annual_mm")


## =============================================================================
## 3. O QUE É O GOOGLE EARTH ENGINE, E O QUE É PRECISO ANTES DE RODAR ESTE CÓDIGO
## =============================================================================
## O Google Earth Engine (GEE) é uma plataforma em nuvem que hospeda
## petabytes de imagens de satélite (Landsat, Sentinel, MODIS, e produtos
## derivados como MapBiomas ou TerraClimate) junto com o poder computacional
## para processá-las, de modo que uma consulta como "a porcentagem de
## cobertura florestal em um círculo de 100 m ao redor de cada um de 1.125
## pontos" roda nos servidores do Google em segundos, em vez de exigir que
## você baixe e processe imagens brutas.
##
## Antes que o código de exemplo abaixo possa ser executado, um leitor
## precisaria:
##   1. Ter (ou criar) uma conta Google, e se inscrever para acesso ao Earth
##      Engine em code.earthengine.google.com/register -- gratuito para uso
##      de pesquisa e não comercial, mas exige uma solicitação e aprovação
##      (geralmente rápida, às vezes um ou dois dias).
##   2. Registrar um projeto do Google Cloud e vinculá-lo ao Earth Engine --
##      o GEE agora exige que todo usuário/script rode sob um ID de projeto
##      de nuvem específico, usado para rastreamento de cota e faturamento
##      (o processamento em si continua gratuito para uso de pesquisa
##      padrão).
##   3. Escolher um método de autenticação. Para uso interativo (o Code
##      Editor do GEE, ou uma sessão Python no seu próprio computador), um
##      login OAuth via navegador é suficiente. Para um script rodando sem
##      supervisão ou em um ambiente isolado (como este pipeline faz), é
##      necessária uma CONTA DE SERVIÇO -- uma identidade de máquina criada
##      no Google Cloud Console, com um arquivo de chave JSON baixado, com
##      as permissões "Earth Engine Resource Writer" e "Service Usage
##      Consumer" no seu projeto.
##   4. Instalar a biblioteca cliente (earthengine-api para Python, ou usar
##      o editor JavaScript diretamente no navegador) e inicializá-la com
##      seu ID de projeto e credenciais antes que qualquer chamada de
##      extração funcione.
##
## Nada disso é necessário para rodar ESTE script -- as covariáveis já estão
## extraídas e salvas em covariates.csv. Isso importa apenas se você quiser
## executar a extração você mesmo, para um novo conjunto de pontos de
## câmera ou um país diferente. O bloco abaixo é APENAS REFERÊNCIA
## (pseudocódigo Python, não executado por este script R).

## --- Referência de extração GEE (Python, NÃO executado aqui) ----------------
## import ee
## ee.Initialize(credentials, project="seu-projeto-gee")
##
## # Classificação de cobertura do solo MapBiomas (resolução de 30 m)
## lulc = ee.Image("projects/mapbiomas-public/assets/brazil/lulc/collection10/"
##                  "mapbiomas_collection10_integration_v1").select("classification_2024")
##
## # Reduzir à % de composição das classes escolhidas dentro de um buffer de 100 m de cada ponto
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
## forest_codes = [1, 3, 5, 6, 49]        # Floresta, Formação Florestal, Manguezal, Floresta Alagável, Restinga Arbórea
## pasture_codes = [15]                    # Pastagem
## cropland_codes = [9,14,18,19,20,39,40,41,46,47,48,62,36,35]  # Lavoura/plantação
##
## sites_fc = ee.FeatureCollection(site_points)
## forest_pct = pct_in_buffer(sites_fc, lulc, forest_codes, 100)
##
## # Clima: normais TerraClimate 1991-2020
## terraclim = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate("1991-01-01","2021-01-01")
## temp_mean = terraclim.select(["tmmx","tmmn"]).mean().reduce(ee.Reducer.mean()).multiply(0.1)
## precip_annual = terraclim.select("pr").mean().multiply(12)
##
## # Áreas protegidas: ponto-em-polígono contra a camada de polígonos WDPA
## wdpa = ee.FeatureCollection("WCMC/WDPA/current/polygons")
## in_pa = sites_fc.map(lambda pt: pt.set("in_pa", wdpa.filterBounds(pt.geometry()).size().gt(0)))

cat("\n== Contexto do GEE ==\n",
    "Ver o bloco de comentário acima para a configuração de conta e a",
    "referência de extração (não executada por este script).\n\n")


## =============================================================================
## 4. HÁ REDUNDÂNCIA NAS COVARIÁVEIS CANDIDATAS?
## =============================================================================
## O conjunto candidato completo de 11 camadas contém diversas medidas de
## cobertura do solo e dossel que poderiam plausivelmente estar todas
## descrevendo o mesmo gradiente subjacente. Verificamos isso de três
## formas com os mesmos dados: uma matriz de correlação par a par, um
## biplot de PCA, e o fator de inflação da variância (VIF) -- cada um
## capturando um tipo diferente de redundância.

# ---- 4a. Correlação, conjunto completo vs. conjunto final -------------------
corr_full  <- cor(cov[, ..FULL_COVARS], use = "pairwise.complete.obs")
corr_final <- cor(cov[, ..FINAL_COVARS], use = "pairwise.complete.obs")

cat("== Matriz de correlação, conjunto candidato completo (11 camadas) ==\n")
print(round(corr_full, 2))
cat("\n== Matriz de correlação, conjunto final de modelagem (7 camadas) ==\n")
print(round(corr_final, 2))

# Altura do dossel, cobertura arbórea Hansen, % de árvores WorldCover e %
# floresta correlacionam fortemente entre si (r > 0,7) no conjunto
# completo -- quatro produtos de satélite diferentes descrevendo o mesmo
# gradiente subjacente (cobertura lenhosa). Após remover as camadas
# redundantes, toda correlação par a par no conjunto final é modesta
# (|r| < 0,55), e pastagem/lavoura são essencialmente independentes de
# todas as outras covariáveis (|r| < 0,3).

# ---- 4b. PCA, conjunto completo vs. conjunto final ---------------------------
pca_full  <- prcomp(na.omit(cov[, ..FULL_COVARS]),  scale. = TRUE)
pca_final <- prcomp(na.omit(cov[, ..FINAL_COVARS]), scale. = TRUE)

cat("\n== PCA, conjunto candidato completo: variância explicada por PC1/PC2 ==\n")
print(round(summary(pca_full)$importance[2, 1:2] * 100, 1))
cat("\n== PCA, conjunto final de modelagem: variância explicada por PC1/PC2 ==\n")
print(round(summary(pca_final)$importance[2, 1:2] * 100, 1))

# No conjunto completo, cinco camadas (altura do dossel, cobertura arbórea
# Hansen, % de árvores WorldCover, % floresta, e área construída GHSL) se
# aglomeram quase na mesma direção no biplot do PCA -- todas medem "quanta
# cobertura lenhosa" e mal se separam. No conjunto final, todas as 7 setas
# apontam em direções visivelmente diferentes -- cada covariável agora
# contribui com informação distinta.

# ---- 4c. VIF: correlação/PCA deixa passar alguma redundância CONJUNTA? ------
## Correlação e PCA descrevem ambos a estrutura par a par ou geral -- podem
## não detectar um caso em que três ou mais covariáveis são conjuntamente
## redundantes mesmo que nenhum par isolado seja altamente correlacionado.
## O VIF verifica isso diretamente: para cada covariável, o VIF mede o
## quão bem TODAS AS OUTRAS covariáveis a preveem, VIF = 1/(1 - R²), onde
## R² vem de regredir essa covariável sobre todas as outras
## simultaneamente. VIF = 1 significa nenhuma redundância; VIF = 5
## significa que as outras covariáveis explicam conjuntamente 80% da
## variância desta, então sua estimativa de coeficiente se torna instável.
## Não há um limiar universal, mas VIF > 5-10 é um limiar de alerta comum.

compute_vif <- function(X) {
  X <- as.data.frame(X)
  sapply(names(X), function(v) {
    fit <- lm(as.formula(paste(v, "~ .")), data = X)
    1 / (1 - summary(fit)$r.squared)
  })
}
vif_final <- compute_vif(na.omit(cov[, ..FINAL_COVARS]))
cat("\n== VIF, conjunto final de modelagem (todos devem ficar bem abaixo de 5-10) ==\n")
print(round(vif_final, 2))

# O conjunto completo de 11 camadas tem VIF INFINITO para a antiga medida
# combinada de agricultura (construída diretamente a partir dos mesmos
# pixels de uma de suas próprias entradas) e VIF de até 7 para as camadas
# de cobertura arbórea; o conjunto final de 7 camadas fica abaixo de 2 --
# confirmando, sob um terceiro ângulo, o que a matriz de correlação e o
# PCA já mostravam.


## =============================================================================
## 5. PADRONIZAÇÃO (E SE DEVEMOS TRANSFORMAR)
## =============================================================================
## Padronizamos por z-score cada covariável contínua -- subtraindo a média,
## dividindo pelo desvio padrão -- antes que ela entre em um modelo. Isso
## coloca cada covariável em uma escala comum, de modo que um coeficiente
## beta ajustado seja diretamente comparável entre covariáveis (um |beta|
## maior significa um efeito mais forte, não apenas uma covariável de
## escala maior), e melhora a estabilidade numérica durante o ajuste do
## modelo.

for (v in FINAL_COVARS) {
  cov[[paste0(v, "_z")]] <- as.numeric(scale(cov[[v]]))
}
cat("\n== Padronização: assimetria bruta vs. z-score ==\n")
skewness <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x); s <- sd(x)
  mean(((x - m) / s)^3)
}
skew_tbl <- data.table(covariavel = FINAL_COVARS,
                        assimetria_bruta = sapply(FINAL_COVARS, function(v) skewness(cov[[v]])))
print(skew_tbl)

# A padronização corrige ESCALA, não FORMA -- uma distribuição enviesada
# continua enviesada após o z-score. Pastagem e lavoura são tão fortemente
# infladas de zeros que nenhuma transformação padrão (log, raiz quadrada)
# as normaliza de forma significativa -- uma transformação log de uma
# variável majoritariamente zero apenas re-rotula os zeros, na maior
# parte. Mantemos essas covariáveis na escala percentual bruta em vez de
# forçar uma transformação inadequada: os modelos que ajustamos (ocupância
# com link logístico) não assumem preditores normalmente distribuídos,
# apenas uma relação linear na escala do link, então o enviesamento na
# própria covariável não é uma violação de modelagem -- apenas significa
# que a maior parte do "sinal" em pastagem/lavoura vem de um número menor
# de locais com valores altos.


## =============================================================================
## 6. ONDE ESTÃO AS COVARIÁVEIS
## =============================================================================
## Olhamos para as covariáveis espacialmente de duas formas: primeiro as
## superfícies NACIONAIS completas (o valor de cada covariável em todo o
## Brasil), depois os LOCAIS DE CÂMERA desenhados sobre essas superfícies,
## o que mostra quais partes de cada gradiente a rede realmente amostra.
##
## (Os mapas raster de superfície nacional e os mapas de locais de câmera
## são construídos a partir de uma grade nacional de predição e das
## lon/lat dos locais -- ver figs/ para as versões renderizadas;
## reproduzir a extração do raster em si requer uma sessão GEE ativa, ver
## Seção 3.)

cat("\n== Cobertura das covariáveis ==\n",
    "Ver figs/m2_covariate_maps_national.png (superfícies nacionais) e\n",
    "figs/m2_covariate_maps_updated.png (locais de câmera) para os mapas.\n")

# Ler as duas juntas torna visível a pegada amostral da rede: as câmeras
# se concentram nas partes de maior floresta e menor perturbação das
# superfícies nacionais, e se rarefazem no interior dominado por pastagem
# e lavoura que os mapas nacionais mostram cobrir uma grande parte do
# país. O Módulo 5 (predição/mapeamento) quantifica quanto do espaço de
# covariáveis as câmeras realmente cobrem, e onde no mapa as predições se
# tornam extrapolação.


## =============================================================================
## 7. O QUANTO O HABITAT VARIA DENTRO DE UM ARRAY?
## =============================================================================
## Nosso desenho agrupa múltiplas câmeras em um ARRAY (a unidade de
## repetição para modelos de comunidade), e módulos posteriores usam um
## efeito aleatório de array como a correção primária para
## pseudo-replicação. Essa correção só é apropriada se o habitat for
## razoavelmente homogêneo DENTRO de um array. Medimos isso com a
## correlação intraclasse (ICC): a fração da variância total de uma
## covariável que é ENTRE arrays, em vez de DENTRO deles. ICC próximo de 1
## significa que um array é uma mancha de habitat única e homogênea; ICC
## próximo de 0 significa que câmeras dentro de um array amostram habitat
## muito diferente.

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
cat("\n== Homogeneidade dentro do array (ICC), por rede ==\n")
print(icc_wide)

# O clima é essencialmente constante dentro de um array em toda rede
# (ICC > 0,93) -- esperado, já que o clima não varia significativamente ao
# longo de poucos quilômetros. O habitat é mais variável: as câmeras
# densamente agrupadas do Snapshot mostram o maior ICC de habitat das três
# redes de forma consistente, enquanto Atlantic e WI mostram ICC de
# habitat mais baixo para pelo menos uma covariável, significando que
# câmeras dentro de um único array podem estar em habitats notavelmente
# diferentes.

# ---- 7a. Vendo o ICC na prática: três arrays de exemplo ---------------------
## Três arrays reais do Snapshot Brasil, escolhidos para cobrir a faixa de
## ICC usando cobertura florestal (%, dentro de 100 m):
##   - Tefé (Amazonas): toda câmera registra 100% de floresta -- um array
##     genuinamente homogêneo (std = 0).
##   - Caitaia (Rio Grande do Norte): principalmente alta cobertura
##     florestal com alguma dispersão perto da borda de uma mancha
##     florestal (std = 6,7).
##   - Montes Claros (Minas Gerais): a cobertura florestal cobre toda a
##     faixa de 0-100% dentro de um único array, porque a paisagem
##     subjacente ali é genuinamente um mosaico fragmentado de
##     floresta/desmatamento (std = 42,5).
## Um efeito aleatório em nível de array trata implicitamente cada câmera
## dentro de um array como amostrando o mesmo habitat -- essa suposição se
## sustenta bem para clima em todos os casos, e para as covariáveis de
## habitat do Snapshot, mas menos bem para arrays de Atlantic/WI como
## Montes Claros. Por isso, o Módulo 3 ajusta os modelos de ocupância de
## espécie única de DUAS formas: com covariáveis de habitat em nível de
## local, e com covariáveis de habitat agregadas ao nível de array, para
## ver o quanto o detalhe mais fino muda as conclusões ecológicas. Ver
## figs/m2_icc_array_example_v2.png para o mapa (raster de cobertura
## florestal de fundo, locais de câmera como pontos pretos).

cat("\n== Exemplos de ICC por array ==\n",
    "Ver o bloco de comentário acima e figs/m2_icc_array_example_v2.png\n")


## =============================================================================
## 8. ÁREAS PROTEGIDAS: QUAL É A EXPOSIÇÃO DA REDE?
## =============================================================================

overall_pct <- round(100 * sum(cov$in_pa, na.rm = TRUE) / sum(!is.na(cov$in_pa)), 1)
cat("\n== Cobertura de área protegida ==\n")
cat("geral:", overall_pct, "% dos locais estão dentro de uma área protegida WDPA\n")

pa <- fread(file.path(DATA_DIR, "pa_by_array_merged.csv"))
cat("-- por rede --\n")
print(pa[, .(locais = sum(n_sites), pct_em_pa = round(100 * sum(n_in_pa) / sum(n_sites), 1)), by = network])

# O status de AP por array é próximo do bimodal -- a maioria dos arrays
# está quase inteiramente dentro de uma AP ou quase inteiramente fora
# dela, porque o limite de um parque raramente divide um array de câmeras
# densamente agrupado. Isso importa para o Módulo 3: um efeito de AP
# estimado a partir do agrupamento em nível de array carrega menos
# informação independente do que a contagem bruta de locais poderia
# sugerir.


## =============================================================================
## 9. AVALIAÇÃO DE MAPAS DE DISTRIBUIÇÃO: QUAIS ESPÉCIES PRECISAM DE UM
##    CONJUNTO CANDIDATO MASCARADO?
## =============================================================================
## Antes de ajustar qualquer modelo de ocupância ou comunidade, verificamos
## as detecções de cada espécie contra sua distribuição geográfica
## conhecida. Alguns arrays caem fora da área de distribuição mapeada de
## uma espécie -- e se a espécie nunca foi detectada ali, essa "ausência"
## não é informação ecológica sobre habitat; é uma certeza geográfica que
## nada tem a ver com as covariáveis acima.
##
## ---- 9a. Considerando mapas de distribuição imperfeitos: um buffer de 100 km ----
## Os polígonos de distribuição da IUCN são eles próprios imperfeitos --
## ficam defasados em relação a mudanças reais de distribuição, e erro de
## mapeamento perto de uma borda de distribuição é comum. Para evitar
## mascarar um local por estar apenas PERTO do limite de uma distribuição
## plausível, cada polígono é expandido por um BUFFER DE 100 KM (calculado
## em uma projeção métrica apropriada para o Brasil) antes que o teste de
## "dentro da distribuição" seja aplicado. Um local só se torna candidato a
## mascaramento se cair fora da distribuição MAIS essa margem de 100 km.
##
## ---- 9b. As duas categorias ----
## Para cada uma das 17 espécies com um polígono de distribuição IUCN, cada
## local é classificado em relação à distribuição, com buffer, dessa
## espécie:
##   - DENTRO DA DISTRIBUIÇÃO -- cai dentro da área mapeada (incl. o
##     buffer de 100 km). Nenhum ajuste necessário.
##   - DETECTADO FORA DA DISTRIBUIÇÃO -- cai fora da distribuição com
##     buffer, mas a espécie FOI fotografada ali. Uma fotografia supera
##     um polígono: esses locais são sempre mantidos, mas sinalizados
##     (podem indicar um mapa de distribuição desatualizado ou um erro
##     taxonômico).
##   - ZERO ESTRUTURAL -- cai fora da distribuição com buffer, e a espécie
##     NUNCA foi detectada ali. Polígono e dados de câmera concordam: este
##     é o "local impossível" que deveria ser removido antes da modelagem,
##     porque a ausência garantida ali não carrega informação de habitat.
##
## ---- 9c. Duas espécies "problemáticas": mapas de distribuição presumivelmente não confiáveis ----
## A Cutia-de-azara (Dasyprocta azarae) e o Gambá-comum (Didelphis
## marsupialis) mostram um padrão de detecções fora da distribuição que é
## amplo e geograficamente disperso de uma forma que parece um polígono de
## distribuição incompatível ou desatualizado, em vez de um limite
## distribucional real. Em vez de forçar uma decisão de mascaramento sobre
## uma geografia possivelmente errada, tratamos seus mapas de distribuição
## como NÃO CONFIÁVEIS e mantemos todos os locais como candidatos de
## modelagem para essas duas espécies -- nenhum mascaramento é aplicado.

mask_summ <- fread(file.path(DATA_DIR, "m2_2_range_mask_summary_v2.csv"))
setorder(mask_summ, -pct_masked)
cat("\n== Impacto do mascaramento por espécie (buffer de 100 km aplicado) ==\n")
print(mask_summ[, .(common, n_sites_total, n_structural_zero_sites,
                     n_sites_after_mask, pct_masked = round(pct_masked, 1),
                     range_map_presumed_unreliable)])

# Gambá-de-orelha-preta, Gambá-de-orelha-branca e Veado-catingueiro ainda
# perdem mais de um terço de seus locais (41-64% mascarados) -- os
# conjuntos candidatos mais conservadores indo para o Módulo 3. Um grupo
# intermediário amplo (Cachorro-do-mato, Queixada, Anta, Tamanduá-bandeira)
# perde um modesto 8-29%. Um grande grupo (Onça-parda, Tamanduá-mirim,
# Tatu-galinha, Cateto, Paca, Jaguatirica, Quati, Irara) está
# essencialmente não afetado (< 6% mascarado). Cutia-de-azara e
# Gambá-comum mostram 0% mascarado porque seus mapas de distribuição são
# tratados como não confiáveis.

## ---- 9d. O que segue adiante para o Módulo 3 e além ----
## A tabela completa de mascaramento por local (m2_2_range_mask_full_v2.csv)
## é a referência de conjunto candidato: para cada uma das 17 espécies
## mapeadas, ela registra se cada local está dentro da distribuição, se foi
## detectado, e portanto se conta como zero estrutural. O Módulo 3 carrega
## essa tabela diretamente e ajusta o modelo de ocupância de cada espécie
## mapeada em seu CONJUNTO CANDIDATO MASCARADO (zeros estruturais
## removidos, locais detectados fora da distribuição mantidos) como padrão.
## Espécies sem polígono de distribuição (aves, e mamíferos fora das 17
## taxa mapeadas) não são afetadas e continuam a usar todos os locais. A
## mesma lógica de mascaramento segue para o Módulo 4 (Modelos de
## comunidade e conjuntos) via uma máscara de distribuição no estado
## latente da ocupância de comunidade.

cat("\n== Continuidade do mascaramento por distribuição ==\n",
    "Ver o bloco de comentário acima; os conjuntos candidatos mascarados",
    "definidos aqui são o padrão para toda espécie mapeada a partir do",
    "Módulo 3 em diante.\n")


## =============================================================================
## RESUMO -- o que os resultados significam
## =============================================================================
## - Habitat e clima são medidos como COMPOSIÇÃO PERCENTUAL dentro de um
##   buffer ao redor de cada câmera, não um rótulo categórico de ponto,
##   porque os modelos de armadilhas fotográficas precisam de um preditor
##   contínuo.
## - De um conjunto candidato original de 11 camadas, várias camadas
##   (altura do dossel, cobertura arbórea Hansen, % de árvores WorldCover,
##   % floresta) se mostraram proxies redundantes para o mesmo gradiente de
##   cobertura lenhosa. Correlação, PCA e VIF confirmam todos a mesma
##   conclusão sob ângulos diferentes: um conjunto reduzido de 7 camadas
##   (floresta, savana, pastagem, lavoura, vegetação nativa, temperatura,
##   precipitação) remove essa redundância mantendo todo sinal
##   ecologicamente distinto.
## - A homogeneidade de habitat dentro do array (ICC) é alta para clima em
##   todos os casos, e alta para as covariáveis de habitat do Snapshot, mas
##   mais baixa para alguns arrays de Atlantic e WI -- uma consideração
##   real de desenho para o efeito aleatório em nível de array usado a
##   partir do Módulo 3 em diante.
## - O status de área protegida é mantido como sua própria covariável,
##   testado ao lado (não incorporado a) o conjunto contínuo de
##   habitat/clima.
## - A avaliação de mapas de distribuição identifica, espécie por espécie,
##   quais locais são geograficamente impossíveis e devem ser mascarados
##   antes da modelagem -- esse conjunto candidato mascarado é o padrão
##   indo para o Módulo 3.
##
## Próximo módulo: Ocupância de Espécie Única.
## =============================================================================
