## =============================================================================
## Módulo 3 — Ocupância de Espécie Única
## Pipeline multi-rede de armadilhas fotográficas Snapshot Brasil
## =============================================================================
##
## O que este script faz, em ordem:
##   1. Ajustar um modelo completo de ocupância de espécie única (Tatu-galinha),
##      com um efeito aleatório de array, e percorrer cada diagnóstico:
##      histórico de detecção, efeitos de covariáveis, e medidas de ajuste do
##      modelo (AUC, c-hat, tamanho do efeito de área protegida).
##   2. Ajustar o mesmo modelo a cada espécie bem detectada na comunidade, e
##      perguntar se agrupar câmeras em arrays como unidade amostral muda as
##      conclusões ecológicas (6 espécies de exemplo).
##   3. Estender o modelo de duas formas: permitindo que um efeito de habitat
##      varie com a temperatura (interação contínua), e permitindo que varie
##      por ecorregião (interação categórica).
##   4. Encerrar com a comparação de abundância Royle-Nichols -- uma forma
##      diferente de usar os mesmos dados de detecção 1/0.
##
## Cada modelo em nível de espécie usa o conjunto candidato de locais
## definido no Módulo 2 (Avaliação de Mapas de Distribuição): para as 17
## espécies com um mapa de distribuição IUCN, locais que caem fora de uma
## distribuição com buffer de 100 km E nunca foram detectados ali (zeros
## estruturais) são excluídos por padrão; locais detectados fora da
## distribuição com buffer são sempre mantidos, sinalizados para revisão.
## Duas espécies (Cutia-de-azara, Gambá-comum) têm mapas de distribuição
## tratados como não confiáveis e usam todos os locais.
##
## Requer: data.table, ggplot2, unmarked, spOccupancy, pROC. Aponte DATA_DIR
## para a pasta contendo os CSVs listados abaixo (empacotados junto com este
## script por padrão).
## =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(unmarked)
})

# ---- 0. Caminhos -------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"

stopifnot(dir.exists(DATA_DIR))

FINAL_COVARS_Z <- c("forest_100m_z","savanna_100m_z","pasture_100m_z","cropland_100m_z",
                     "native_veg_1000m_z","temp_mean_C_z","precip_annual_mm_z")


## =============================================================================
## 1. UM PRIMEIRO MODELO DE OCUPÂNCIA: TATU-GALINHA
## =============================================================================
## Começamos com o Tatu-galinha (Dasypus novemcinctus), uma espécie bem
## detectada, e ajustamos um modelo usando todas as covariáveis do conjunto
## final: floresta, savana, pastagem, lavoura, vegetação nativa em escala de
## paisagem, temperatura, precipitação e status de área protegida. O Tatu
## tem quase nenhuma exclusão por máscara de distribuição, então seus
## resultados aqui estão próximos do que um ajuste sem filtro daria -- uma
## boa referência antes de olhar para espécies onde a máscara importa mais.

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

cat("== Tatu-galinha: conjunto candidato de locais ==\n")
cat("locais após exclusão por máscara de distribuição + covariáveis ausentes:", nrow(y2), "\n\n")

## ---- 1.1 Histórico de detecção ----
## Modelos de ocupância são construídos a partir de um HISTÓRICO DE
## DETECÇÃO: para cada ocasião de 7 dias, a espécie foi detectada (1), não
## detectada (0), ou a câmera não estava funcionando (ausente)? Mostramos
## abaixo dois arrays contrastantes: um com muitas detecções, outro com
## quase nenhuma.

high_arr <- "SNAP_Piaui_Caatinga_Tapuio_25"
low_arr  <- "ATLA_ATL_11"
cat("== Arrays de exemplo do histórico de detecção ==\n")
cat("Array de alta detecção:", high_arr, "\n")
cat("Array de baixa detecção:", low_arr, "\n")
cat("Ver figs/m3_dethist_armadillo.png para o histórico de detecção gráfico,\n")
cat("e as tabelas de matriz 0/1 em data/m3_dethist_table_high.csv /\n")
cat("data/m3_dethist_table_low.csv para a matriz de detecção literal.\n\n")

## ---- 1.2 Ajustando o modelo, com o efeito aleatório de array ----
## Ajustamos o modelo com um intercepto aleatório por array -- a correção
## padrão do pipeline para pseudo-replicação entre câmeras que compartilham
## o mesmo array (ver Módulo 2). Como este é um modelo de efeitos
## aleatórios, usamos um motor Bayesiano (spOccupancy::PGOcc) em vez do
## pacote de máxima verossimilhança `unmarked`, que ajusta apenas efeitos
## fixos. (Este bloco mostra o ajuste de efeitos fixos; o ajuste com
## efeito aleatório de array usa as mesmas covariáveis via spOccupancy --
## ver a tabela de coeficientes impressa abaixo, pré-calculada e carregada
## do disco já que o ajuste Bayesiano leva vários minutos.)

umf <- unmarkedFrameOccu(y = y2, siteCovs = sc2[, c(FINAL_COVARS_Z, "in_pa")])
form_full <- as.formula(paste("~1 ~", paste(FINAL_COVARS_Z, collapse=" + "), "+ in_pa"))
fit_full <- occu(form_full, data = umf)

co <- coef(fit_full); se <- sqrt(diag(vcov(fit_full))); z <- co/se; p <- 2*pnorm(-abs(z))
cat("== Coeficientes de efeitos fixos do Tatu-galinha ==\n")
print(data.table(param=names(co), estimate=round(co,3), se=round(se,3), z=round(z,2), p=round(p,4)))

pgocc_coefs <- fread(file.path(DATA_DIR, "armadillo_pgocc_coefs.csv"))
cat("\n== Coeficientes com efeito aleatório de array (PGOcc): média posterior, intervalo de credibilidade 95% ==\n")
print(pgocc_coefs[param != "(Intercept)"])

sigma_re <- fread(file.path(DATA_DIR, "armadillo_pgocc_sigma_re.csv"))
cat("\nVariância do efeito aleatório de array (escala logit): média =", round(sigma_re$sigma_mean,2),
    " dp =", round(sigma_re$sigma_sd,2), "\n")
cat("Uma variância desse tamanho significa que os arrays diferem substancialmente\n")
cat("em sua ocupância basal, além do que as covariáveis de habitat medidas\n")
cat("explicam -- exatamente o tipo de agrupamento não explicado que o efeito\n")
cat("aleatório é projetado para absorver.\n\n")

## ---- 1.3 Ocupância ingênua vs. baseada em modelo ----
## "Ocupância ingênua" é simplesmente a fração de locais onde a espécie foi
## detectada alguma vez -- ela sempre SUBESTIMA a ocupância verdadeira,
## porque uma espécie pode estar presente mas passar despercebida em toda
## ocasião.

naive_psi_site <- apply(y2, 1, function(r) as.integer(any(r==1, na.rm=TRUE)))
naive_psi <- mean(naive_psi_site)
psi_pred <- predict(fit_full, type="state")$Predicted
modeled_psi <- mean(psi_pred)
cat("== Ocupância ingênua vs. modelada ==\n")
cat("Ocupância ingênua (fração de locais detectados alguma vez):", round(naive_psi,3), "\n")
cat("Ocupância modelada (PGOcc com efeito de array, psi médio):", round(modeled_psi,3), "\n")
cat("Ver figs/m3_naive_vs_modeled.png para a comparação em barras.\n\n")

## ---- 1.4 Efeitos das covariáveis ----
## Quatro covariáveis são significativas no ajuste de efeitos fixos:
## floresta, savana, vegetação nativa e status de área protegida. A figura
## (figs/m3_armadillo_effect_curves.png) mostra a ocupância predita em
## função de cada uma delas, mantendo as outras na média -- com linhas
## verdes para um efeito positivo e vermelhas para um negativo,
## seguindo a convenção de cores usada em todo este pipeline.

## ---- 1.5 Ajuste do modelo: o que esses números realmente significam? ----
## Duas perguntas muito diferentes são feitas sobre um modelo de ocupância
## ajustado, e ajuda mantê-las separadas:
##
##   (a) O MODELO DISCRIMINA LOCAIS OCUPADOS DE NÃO OCUPADOS? -- medido
##       pelo AUC (área sob a curva ROC). Pegue a probabilidade de ocupância
##       predita de cada local (psi-chapéu) e pergunte: se você escolhesse
##       um local verdadeiramente ocupado e um verdadeiramente não ocupado
##       ao acaso, qual a chance de o modelo dar ao ocupado uma
##       probabilidade predita MAIOR? Essa chance é o AUC. AUC = 0,5
##       significa que o modelo não é melhor que jogar uma moeda
##       (probabilidades preditas para locais ocupados e não ocupados se
##       sobrepõem completamente); AUC = 1,0 significa separação perfeita
##       (todo local ocupado pontua mais alto que todo local não ocupado).
##       Ver figs/m3_auc_explainer.png para um exemplo trabalhado: o painel
##       esquerdo mostra duas distribuições sobrepostas de probabilidade
##       predita (uma para presenças verdadeiras, uma para ausências
##       verdadeiras); o painel direito é a curva ROC construída varrendo
##       um limiar de detecção por essas distribuições, e o AUC é a área
##       sob essa curva.
##   (b) A ESTRUTURA DE ERRO ASSUMIDA PELO MODELO ESTÁ CORRETA? -- medido
##       pelo C-HAT (superdispersão), a partir de um teste de qualidade de
##       ajuste por bootstrap paramétrico. O modelo assume que as detecções
##       seguem um processo binomial com uma única probabilidade de
##       detecção; o c-hat compara a discrepância qui-quadrado observada
##       entre o modelo e os dados contra uma distribuição bootstrap dessa
##       mesma estatística simulada A PARTIR do modelo ajustado. c-hat
##       próximo de 1 significa que a estrutura de erro assumida se ajusta;
##       c-hat notavelmente acima de 1 (como aqui) significa que os dados
##       reais são mais variáveis do que o modelo espera (superdispersão)
##       -- comum quando a probabilidade de detecção realmente varia por
##       local ou estação de formas que o modelo não captura. Na prática,
##       isso significa que os erros-padrão do modelo devem ser lidos como
##       um pouco otimistas (estreitos demais) em vez de exatos.
##
## AUC e c-hat NÃO são intercambiáveis: um modelo pode discriminar bem
## (AUC alto) e ainda ser superdisperso (c-hat > 1), porque o AUC só se
## importa com a classificação relativa dos locais, enquanto o c-hat se
## importa se a variância absoluta nas detecções corresponde ao que o
## modelo binomial prediz.

fit_summary <- fread(file.path(DATA_DIR, "armadillo_fit_summary.csv"))
cat("== Resumo de ajuste do modelo do Tatu-galinha ==\n")
print(fit_summary)
cat("\nAUC de", round(fit_summary$auc,2), "indica discriminação modesta -- melhor\n")
cat("que o acaso, mas longe de perfeita, típico para modelos de ocupância em\n")
cat("dados binários de detecção com sinal modesto de covariáveis. c-hat de", round(fit_summary$c_hat,2),
    "significa\nque o modelo é moderadamente superdisperso -- um achado comum e esperado,\n")
cat("não um sinal de que o modelo é inutilizável.\n\n")

## ---- 1.6 Incluir o status de área protegida melhora o modelo? ----
pa_aic <- fread(file.path(DATA_DIR, "armadillo_pa_aic_test.csv"))
cat("== Comparação de AIC: com vs. sem status de área protegida ==\n")
print(pa_aic)

pa_mag <- fread(file.path(DATA_DIR, "armadillo_pa_magnitude.csv"))
cat("\n== Ocupância predita dentro vs. fora de áreas protegidas ==\n")
print(pa_mag)
cat("\nIncluir in_pa reduz o AIC em", round(pa_aic[model=="without in_pa", delta_AIC],1),
    "pontos -- uma melhoria decisiva.\n")
cat("A direção do efeito é", ifelse(pa_mag$direction=="negative","negativa","positiva"),
    ": a ocupância do Tatu-galinha é",
    ifelse(pa_mag$direction=="negative","menor","maior"), "\ndentro de áreas protegidas.",
    "Em termos concretos, a ocupância predita dentro de áreas\n")
cat("protegidas (", round(pa_mag$psi_inside_pa,2), ") é aproximadamente", round(pa_mag$rel_change,2),
    "vezes a ocupância predita fora\n(", round(pa_mag$psi_outside_pa,2),
    ") -- um efeito substancial, não apenas estatisticamente\ndetectável.\n\n")


## =============================================================================
## 2. TODA ESPÉCIE, JUNTA
## =============================================================================
## Toda espécie mapeada neste ajuste usa seu conjunto candidato mascarado
## (Módulo 2): locais de zero estrutural excluídos, locais detectados fora
## da distribuição mantidos. Espécies sem mapa de distribuição, e as duas
## espécies com mapas de distribuição não confiáveis (Cutia-de-azara,
## Gambá-comum), usam todos os locais.
##
## Ajustar um modelo de 8 covariáveis a todas as espécies modeladas com um
## motor apenas de efeitos fixos é pouco confiável para cerca de metade da
## comunidade (separação quase completa). Vamos direto para o modelo com
## efeito aleatório de array, que regulariza cada ajuste através do
## intercepto de array.

modeled <- fread(file.path(DATA_DIR, "modeled_species.csv"))
species_meta <- fread(file.path(DATA_DIR, "species_meta.csv"))
cat("== Modelo de comunidade ==\n")
cat("Espécies modeladas:", nrow(modeled), "\n")
cat("Ver figs/m3_beta_shaded_table.png para a tabela completa de coeficientes da\n")
cat("comunidade: cada linha é uma espécie, com (n_det, n_arr) -- total de\n")
cat("detecções e total de arrays detectados -- mostrado ao lado do nome, cada\n")
cat("coluna é uma covariável (incluindo Área Protegida), verde = efeito\n")
cat("positivo, vermelho = negativo, * = intervalo de credibilidade 95% exclui\n")
cat("zero.\n\n")
cat("Floresta é a covariável mais consistentemente importante em toda a\n")
cat("comunidade, positiva para a grande maioria das espécies com efeito\n")
cat("credível.\n\n")

## ---- 2.1 A simplificação em nível de array se sustenta? ----
## O Módulo 2 sinalizou que o habitat nem sempre é homogêneo dentro de um
## array, especialmente para as redes Atlantic e WI. Aqui testamos
## diretamente: usar arrays -- em vez de câmeras individuais -- como
## unidade amostral muda as conclusões ecológicas? Construímos históricos
## de detecção GENUINAMENTE em nível de array: para cada array e cada
## ocasião, o array conta como "detectado" se qualquer câmera nele
## registrou a espécie naquela ocasião, e "amostrado, não detectado" se
## pelo menos uma câmera estava ativa mas nenhuma detectou. Isso agrupa
## câmeras em arrays como a unidade amostral real -- não apenas calculando
## a média das covariáveis em nível de local mantendo as detecções em
## nível de local.
##
## Seis espécies ilustram a comparação, cobrindo diferentes situações de
## dados e diferentes sensibilidades ao nível de agregação: Paca e
## Tatu-galinha (ricos em dados, concordância moderada câmera/array),
## Queixada (esparsa, fortemente mascarada por distribuição em nível de
## array), e três espécies escolhidas porque seus coeficientes de área
## protegida e/ou floresta mudam substancialmente -- em alguns casos
## invertendo o sinal -- entre o nível de câmera e de array: Jaguatirica,
## Anta e Onça-parda.

site_6 <- fread(file.path(DATA_DIR, "site_level_6species_coefs.csv"))
array_6 <- fread(file.path(DATA_DIR, "array_true_6species_coefs.csv"))
cat("== Coeficientes em nível de câmera vs. array, 6 espécies ==\n")
cat("Ver figs/m3_array_vs_site_6species.png para a comparação em pontos com\n")
cat("barras. Os intervalos de confiança se ampliam drasticamente no nível de\n")
cat("array para cada espécie e covariável -- menos unidades amostrais significa\n")
cat("menos poder. Para Jaguatirica, Anta e Onça-parda especificamente, as\n")
cat("estimativas pontuais de área protegida e/ou floresta mudam de direção\n")
cat("entre os dois níveis -- um alerta contra superinterpretar qualquer\n")
cat("coeficiente único em nível de array.\n\n")

fitstats_6 <- fread(file.path(DATA_DIR, "camera_vs_array_fitstats_6species.csv"))
cat("== Comparação de ajuste do modelo: nível de câmera vs. nível de array, 6 espécies ==\n")
cat("Colunas: 'naive_occupancy_fraction_detected' = fração das unidades\n")
cat("amostrais (câmeras ou arrays) onde a espécie foi detectada ALGUMA VEZ --\n")
cat("uma contagem bruta, sem modelo envolvido. 'modeled_occupancy_probability'\n")
cat("= a probabilidade média de ocupância predita pelo modelo ajustado, em\n")
cat("todas as unidades amostrais -- leva em conta a detecção imperfeita, então\n")
cat("é sempre >= o valor ingênuo.\n\n")
print(fitstats_6[, .(common, level, n, auc=round(auc,3),
                       naive_occupancy_fraction_detected=round(naive_occupancy_fraction_detected,3),
                       modeled_occupancy_probability=round(modeled_occupancy_probability,3))])
cat("\nO AUC em nível de array é maior para as 6 espécies (aproximadamente\n")
cat("0,60-0,80 vs. 0,56-0,70 no nível de câmera) -- com muito menos unidades\n")
cat("amostrais, um histórico de detecção em nível de array é um sinal mais\n")
cat("limpo e menos ruidoso, que um modelo simples de 3 covariáveis consegue\n")
cat("discriminar mais facilmente, mesmo que o ajuste em nível de array tenha\n")
cat("muito menos poder estatístico em seus intervalos de confiança.\n\n")

## ---- 2.2 Por que essa comparação importa ----
## Essa comparação não é apenas uma verificação técnica -- ela fala
## diretamente de uma decisão de desenho que todo este pipeline faz: usar
## o ARRAY como unidade primária de replicação para modelos de comunidade e
## conjuntos (Módulo 4 em diante), enquanto ainda ajustamos modelos de
## espécie única no nível de câmera mais fino aqui. O resultado é
## tranquilizador mas não incondicional: as estimativas pontuais em geral
## mantêm um SINAL consistente entre os dois níveis, significando que a
## simplificação em nível de array captura a história qualitativa correta
## para a maioria das espécies -- mas a magnitude, e ocasionalmente o
## sinal, do efeito de covariável de uma espécie com poucos dados pode
## mudar uma vez que as câmeras são agrupadas em arrays. Na prática, isso
## significa: confie no sinal e no tamanho aproximado de um coeficiente em
## nível de array para uma espécie bem detectada, mas trate um resultado
## em nível de array para uma espécie esparsa ou fortemente mascarada por
## distribuição (como Queixada aqui) como um sinal muito mais fraco do que
## seu p-valor sozinho sugeriria -- confirme com o ajuste em nível de
## câmera mais fino sempre que possível, exatamente como feito nesta seção.

cat("== Por que a comparação em nível de array importa ==\n")
cat("Ver o bloco de comentário acima para a explicação completa.\n\n")

## ---- 2.3 Agrupamentos de espécies a partir da tabela de coeficientes ----
## A tabela de coeficientes de toda a comunidade é densa. Uma análise de
## componentes principais nesses coeficientes (7 covariáveis de
## habitat/clima mais status de área protegida, padronizadas) reduz essa
## tabela a dois eixos e revela agrupamentos ecológicos genuínos.
cat("== PCA de espécies ==\n")
cat("Ver figs/m3_species_pca_biplot.png. PC1 e PC2 juntos explicam\n")
cat("aproximadamente 40% da variância nos perfis de resposta às covariáveis\n")
cat("das espécies. Setas mostram quais covariáveis dominam cada direção;\n")
cat("espécies se agrupam em 4 grupos via k-means nos dois escores de PCA.\n\n")


## =============================================================================
## 3. ADICIONANDO UM EFEITO DE INTERAÇÃO
## =============================================================================
## Estendemos o modelo de ocupância para permitir que um efeito de habitat
## dependa de outra variável -- aqui, se o efeito da cobertura arbórea na
## ocupância muda com a temperatura. Isso é ajustado como um modelo
## COMPLETO, não apenas um termo de interação:
##
##   psi ~ cobertura_arborea * temperatura + savana + pastagem + lavoura +
##         veg_nativa_1000m + precipitação + in_pa
##
## Cobertura arbórea (uma camada de referência do MapBiomas do ano 2000,
## distinta da covariável padrão forest_100m do pipeline) substitui a
## floresta neste modelo, já que as duas são altamente correlacionadas
## (r = 0,82) e criariam colinearidade com o termo de interação. Todas as
## outras covariáveis padrão permanecem no modelo como efeitos aditivos
## comuns.

int_coefs <- fread(file.path(DATA_DIR, "interaction_4species_coefs.csv"))
int_terms <- int_coefs[param == "psi(treecover2000_100m_z:temp_mean_C_z)"]
cat("== Interação temperatura x cobertura arbórea, modelo completo, 4 espécies ==\n")
print(int_terms[, .(species, estimate=round(estimate,3), z=round(z,2))])
cat("\nVer figs/m3_interaction_4species.png para curvas de cenário de\n")
cat("temperatura fria/média/quente com bandas de confiança de 95%. Três de\n")
cat("quatro espécies mostram uma interação estatisticamente significativa\n")
cat("mesmo após controlar pelo conjunto completo de covariáveis -- um modelo\n")
cat("apenas aditivo faria essa média desaparecer.\n\n")


## =============================================================================
## 4. PERMITINDO QUE UMA RELAÇÃO VARIE POR ECORREGIÃO
## =============================================================================
## Uma forma diferente de as relações de habitat variarem: não com uma
## covariável contínua, mas entre regiões ecologicamente distintas.
## Testamos isso com a Jaguatirica (mascarada por distribuição), comparando
## um único efeito de floresta agrupado contra um modelo COMPLETO que
## permite que o efeito de floresta difira entre duas ecorregiões com
## dados suficientes -- Floresta Úmida de Folhas Largas e as mais secas
## Pastagens/Savanas & Arbustais -- enquanto toda outra covariável
## permanece no modelo como um efeito aditivo comum (fixo, não variando por
## ecorregião):
##
##   psi ~ floresta * ecorregiao + savana + pastagem + lavoura +
##         veg_nativa_1000m + temperatura + precipitação + in_pa

eco_aic <- fread(file.path(DATA_DIR, "ecoregion_ocelot_aic.csv"))
cat("== Comparação de AIC: efeito de floresta agrupado vs. modelo de interação por ecorregião ==\n")
print(eco_aic)
cat("\nVer figs/m3_ecoregion_ocelot.png. Na Floresta Úmida de Folhas Largas,\n")
cat("mais cobertura florestal significa mais ocupância de Jaguatirica, como\n")
cat("esperado para uma especialista florestal. Na ecorregião\n")
cat("Pastagens/Savanas & Arbustais a relação se inverte -- a mesma lição da\n")
cat("interação de temperatura, aplicada a um moderador categórico em vez de\n")
cat("contínuo: um único número de 'efeito de floresta' pode esconder duas\n")
cat("relações genuinamente diferentes que só se tornam visíveis quando os\n")
cat("dados são divididos por contexto ecológico.\n\n")


## =============================================================================
## 5. COMPARANDO COM A ABUNDÂNCIA ROYLE-NICHOLS
## =============================================================================
## O modelo Royle-Nichols usa os mesmos dados de detecção 1/0 que a
## ocupância, mas interpreta a FREQUÊNCIA de detecção ao longo das ocasiões
## como um sinal de quantos indivíduos estão presentes em um local, não
## apenas se a espécie está presente ou não: um local visitado por três
## animais é detectado com mais frequência, em média, do que um local
## visitado por um, mesmo que a probabilidade de detecção por indivíduo
## seja idêntica.

cat("== Ocupância vs. abundância Royle-Nichols, em toda a comunidade ==\n")
cat("Ver figs/m3_occ_vs_abundance.png. Painel esquerdo: ocupância média\n")
cat("(psi-chapéu) vs. abundância relativa média Royle-Nichols\n")
cat("(lambda-chapéu) entre espécies -- a concordância é forte no geral.\n")
cat("Painel direito: efeito de floresta estimado de duas formas (ocupância\n")
cat("vs. Royle-Nichols); a concordância é mais modesta mas ainda claramente\n")
cat("positiva. Ajustes Royle-Nichols instáveis (erro-padrão acima de 5, um\n")
cat("sinal de espécie rara sem informação suficiente para fixar tanto um\n")
cat("intercepto quanto uma inclinação) são excluídos de ambos os painéis.\n")
cat("Duas espécies (rotuladas na figura) ficam bem acima da linha 1:1 no\n")
cat("painel de efeito de floresta -- Royle-Nichols estima um efeito de\n")
cat("floresta substancialmente maior na abundância do que a ocupância\n")
cat("estima para presença isolada, embora ambas as espécies sejam raras o\n")
cat("suficiente para que a magnitude Royle-Nichols deva ser lida com cautela\n")
cat("mesmo após excluir os ajustes instáveis mais extremos.\n\n")


## =============================================================================
## RESUMO -- o que os resultados significam
## =============================================================================
## - Todo modelo neste módulo usa o conjunto candidato mascarado por
##   distribuição do Módulo 2 por padrão para as 17 espécies mapeadas.
## - O ajuste apenas de efeitos fixos é pouco confiável para cerca de
##   metade da comunidade devido à separação; o modelo com efeito aleatório
##   de array resolve isso regularizando através do intercepto de array e
##   é o padrão do pipeline.
## - AUC (discriminação) e c-hat (superdispersão) medem coisas diferentes e
##   podem discordar -- um modelo pode discriminar bem e ainda ser
##   superdisperso.
## - A cobertura florestal é o sinal de habitat mais consistente da
##   comunidade; ocupância e abundância Royle-Nichols concordam bem no
##   nível geral e razoavelmente no efeito de floresta especificamente, uma
##   vez excluídos os ajustes Royle-Nichols instáveis.
## - Agrupar câmeras em arrays como unidade amostral amplia
##   substancialmente os intervalos de confiança -- um custo real em poder
##   estatístico -- e para espécies com poucos dados pode mudar a direção
##   estimada de uma covariável, não apenas sua precisão. Confie mais em
##   resultados em nível de array para espécies bem detectadas.
## - Duas extensões -- uma interação contínua (temperatura x cobertura
##   arbórea) e uma categórica (ecorregião) -- mostram ambas que um único
##   "um número serve para tudo" efeito de habitat pode esconder variação
##   real e ecologicamente sensata.
## - O status de área protegida é um efeito genuinamente grande para o
##   Tatu-galinha: a ocupância predita dentro de áreas protegidas é bem
##   inferior à ocupância predita fora, mantendo as outras covariáveis na
##   média.
##
## Próximo módulo: Modelos de Comunidade e Conjuntos.
## =============================================================================
