## =============================================================================
## Módulo 1 — Fundação e Históricos de Detecção
## Pipeline de análise de armadilhas fotográficas do Snapshot Brasil
## =============================================================================
##
## O que este script faz, em ordem:
##   1. Carregar a própria exportação bruta do Snapshot Brasil e descrever sua
##      estrutura (instalações, placenames, sequências, arrays) e o tipo e
##      altura das armadilhas fotográficas.
##   2. Resumir a própria comunidade de espécies do Snapshot Brasil, antes de
##      qualquer combinação.
##   3. Quantificar o custo do limiar de independência de 60 minutos (vs. o
##      padrão convencional de 1 minuto) usando os próprios registros de tempo
##      do Snapshot Brasil.
##   4. Carregar o conjunto de dados COMBINADO final (três redes: Snapshot
##      Brasil, Mata Atlântica, projetos públicos WI) e descrever como foi
##      montado.
##   5. Explicar o pipeline instalação -> sítio -> array (por que trocas de
##      cartão e instalações sobrepostas precisam de tratamento cuidadoso).
##   6. Construir um histórico de detecção para uma espécie/array exemplo e
##      mostrar tanto a forma gráfica quanto a matriz literal 0/1.
##   7. Discutir a escolha da janela de ocasião (por que 7 dias).
##   8. Resumir a comunidade de mamíferos combinada (taxa de detecção vs.
##      ocupação ingênua).
##
## Restrito a MAMÍFEROS durante todo o processo, porque a compilação da Mata
## Atlântica registrou apenas mamíferos — manter aves tornaria toda ave
## falsamente "ausente" em todo array da Mata Atlântica (um artefato de
## dados, não ecologia).
##
## Requer: data.table, ggplot2. Aponte DATA_DIR para a pasta contendo os CSVs
## listados abaixo (empacotados junto a este script por padrão).
## =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})

# ---- 0. Caminhos -------------------------------------------------------------
# Edite DATA_DIR se mover este script para longe de sua pasta data/ empacotada.
DATA_DIR <- "data"
FIGS_DIR <- "figs"   # usado apenas se quiser resalvar figuras a partir deste script

stopifnot(dir.exists(DATA_DIR))


## =============================================================================
## 1. SNAPSHOT BRASIL — estrutura e características (sua própria exportação bruta)
## =============================================================================
## O Wildlife Insights (a plataforma sobre a qual o Snapshot Brasil roda)
## organiza os dados assim:
##   - uma INSTALAÇÃO (deployment): uma câmera funcionando continuamente em um
##     PLACENAME por um período. Reparar uma câmera no meio da temporada
##     (troca de bateria/cartão) cria uma segunda linha de instalação no mesmo
##     placename — ainda uma única câmera física.
##   - uma SEQUÊNCIA: uma rajada de fotos agrupada como um único evento de
##     disparo, com uma identificação de espécie.
##   - um ARRAY: a unidade básica de amostragem do Snapshot Brasil — um grupo
##     de estações de câmera espacialmente agrupadas (uma "grade"/"subprojeto")
##     instaladas em conjunto.

snap_dep <- fread(file.path(DATA_DIR, "snapshot_brasil_deployments.csv"))
snap_seq <- fread(file.path(DATA_DIR, "snapshot_brasil_sequences.csv"))

cat("== Estrutura do Snapshot Brasil ==\n")
cat("instalações:", nrow(snap_dep),
    "| placenames distintos:", uniqueN(snap_dep$placename),
    "| arrays (subprojetos):", uniqueN(snap_dep$subproject_name), "\n")
cat("total de câmeras-dia (instalações brutas):",
    sum(as.integer(as.Date(snap_dep$end_date) - as.Date(snap_dep$start_date)), na.rm = TRUE), "\n")
cat("total de sequências registradas:", nrow(snap_seq), "\n")

# NOTA: data.table::fread lê células CSV em branco como strings vazias (""),
# NÃO como NA. Um filtro como `!is.na(genus)` mantém silenciosamente linhas
# com gênero em branco (sequências não identificadas, ex.: rótulos de nível
# de família como "Rodent") — sempre filtre também por `genus != ""` ao
# excluir registros não identificados.

# ---- 1a. Tipo e altura da armadilha fotográfica ------------------------------
# Cada instalação registra feature_type (para onde a câmera estava apontada)
# e sensor_height (a que altura foi montada). Isso importa ecologicamente: uma
# câmera em uma toca/fonte de água amostra diferente de um ponto de grade
# aleatorizado, e a altura de montagem restringe quais espécies/comportamentos
# são capturados.

height_tbl <- snap_dep[, .N, by = sensor_height][order(-N)]
cat("\n-- distribuição de altura do sensor --\n"); print(height_tbl)

feature_tbl <- snap_dep[, .N, by = feature_type][order(-N)]
cat("\n-- distribuição de tipo de característica --\n"); print(feature_tbl)

# Verificamos especificamente um subconjunto de câmeras em altura de dossel
# (instalações arbóreas precisariam de tratamento separado — comunidade
# diferente, e a suposição de "ativo" para um mamífero terrestre deixa de
# valer ali). NENHUM EXISTE nos dados atuais: toda altura registrada é ao
# nível do solo (altura do joelho ~40-50cm, altura do peito ~1-1,3m, ou
# "Outra" resolvendo em ~40cm ou ~1m acima do solo no campo de texto livre).
# Nenhuma exclusão por dossel é necessária.
cat("\nVerificação de câmera em dossel/arbórea: categorias de sensor_height presentes são",
    paste(unique(snap_dep$sensor_height), collapse = ", "),
    "— todas ao nível do solo. Não existe conjunto de dossel nestes dados.\n")


## =============================================================================
## 2. SNAPSHOT BRASIL — sua própria comunidade de espécies (antes de qualquer combinação)
## =============================================================================

snap_mammals <- snap_seq[class == "Mammalia" & genus != "Homo" &
                            genus != "" & species != ""]
snap_mammals[, sci := paste(genus, species)]

cat("\n== Comunidade própria do Snapshot Brasil ==\n")
cat("táxons de mamíferos distintos (gênero+espécie bruto, antes da harmonização):",
    uniqueN(snap_mammals$sci), "\n")
cat("total de sequências de mamíferos identificadas:", nrow(snap_mammals), "\n")

top_species <- snap_mammals[, .N, by = .(sci, common_name)][order(-N)][1:10]
cat("\n-- dez táxons de mamíferos mais registrados --\n"); print(top_species)


## =============================================================================
## 3. O INTERVALO DE INDEPENDÊNCIA — 1 minuto vs. 60 minutos
## =============================================================================
## Todo conjunto de dados de armadilhas fotográficas precisa de uma regra para
## transformar um fluxo de fotos em "eventos" de detecção discretos e
## independentes (caso contrário, um animal que permanece diante da câmera, ou
## um grupo passando junto, seria contado muitas vezes). A convenção usual é
## 1 MINUTO. Este pipeline usa 60 MINUTOS em vez disso, por consistência com a
## compilação da Mata Atlântica (cujos estudos-fonte não foram registrados com
## precisão suficiente para sustentar uma regra de 1 minuto). Essa mudança não
## é gratuita — quantificamos exatamente o que ela custa usando os próprios
## registros de tempo do Snapshot Brasil (a única rede neste pipeline com
## tempos de sequência de granularidade fina).

snap_seq[, start_time := as.POSIXct(start_time)]
snap_seq[, end_time   := as.POSIXct(end_time)]

mammals_ts <- snap_seq[class == "Mammalia" & genus != "Homo" &
                          genus != "" & species != ""]
mammals_ts[, sci := paste(genus, species)]
setorder(mammals_ts, deployment_id, sci, start_time)

# Lacuna (minutos) até a sequência anterior da MESMA espécie na MESMA
# instalação. A primeira sequência por instalação/espécie não tem sequência
# anterior, então é sempre contada como independente (lacuna = Inf).
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

cat("\n== Intervalo de independência: 1-min vs 60-min (Snapshot Brasil) ==\n")
cat("total de sequências de mamíferos:", n_total, "\n")
cat("independentes @ 1-min: ", n_indep_1min, "\n")
cat("independentes @ 60-min:", n_indep_60min, "\n")
cat("detecções colapsadas (1min -> 60min):", n_collapsed,
    sprintf(" (%.1f%%)\n", pct_collapsed))

# Detalhamento por espécie: quais espécies perdem mais detecções com a janela
# mais larga? (Restrito a espécies com >=20 detecções independentes a 1-min,
# para uma porcentagem estável.)
by_species <- mammals_ts[, .(n_indep_1min = sum(indep_1min),
                              n_indep_60min = sum(indep_60min)),
                          by = .(sci, common_name)]
by_species[, n_collapsed := n_indep_1min - n_indep_60min]
by_species[, pct_collapsed := 100 * n_collapsed / n_indep_1min]
by_species_top <- by_species[n_indep_1min >= 20][order(-pct_collapsed)][1:12]

cat("\n-- espécies mais afetadas pelo limiar de 60 minutos --\n")
print(by_species_top[, .(common_name, sci, n_indep_1min, n_indep_60min,
                          n_collapsed, pct_collapsed = round(pct_collapsed, 1))])

# O padrão é ecológico, não arbitrário: espécies gregárias ou que permanecem
# no local (mocó, gado, quati, queixadas) perdem mais detecções, porque um
# grupo ou um animal que permanece dispara capturas repetidas dentro da mesma
# hora. Espécies solitárias e de amplo deslocamento são comparativamente
# pouco afetadas.


## =============================================================================
## 4. COMBINANDO TRÊS REDES EM UM ÚNICO CONJUNTO DE DADOS
## =============================================================================
## A análise completa combina:
##   - Snapshot Brasil    : o levantamento nacional coordenado de 2025 (Seção 1-3)
##   - Mata Atlântica     : uma compilação de estudos da Mata Atlântica, 2016-2020
##   - Projetos públicos WI: projetos públicos do Wildlife Insights no Brasil, 2010-2025
## Todos em formato Wildlife Insights.

dep <- fread(file.path(DATA_DIR, "final_deployments.csv"))
det <- fread(file.path(DATA_DIR, "detections.csv"))   # apenas mamíferos

cat("\n== Conjunto de dados combinado ==\n")
cat("redes:", uniqueN(dep$network),
    "| arrays:", uniqueN(dep$array_id),
    "| sítios:", uniqueN(dep$site_id),
    "| instalações:", nrow(dep),
    "| câmeras-dia:", format(sum(dep$camera_days, na.rm = TRUE), big.mark = ","), "\n")
cat("detecções de mamíferos:", format(nrow(det), big.mark = ","),
    "| espécies de mamíferos (todas as identificações):", uniqueN(det$sci_mdd), "\n")

# ---- 4a. Como os dados da Mata Atlântica e WI foram selecionados ------------
# O Snapshot Brasil é um levantamento propositalmente construído com desenho
# consistente. As outras duas redes são compilações de estudos independentes,
# então precisaram de harmonização ao mesmo padrão antes da combinação. As
# MESMAS regras foram aplicadas a ambas:
#   - Excluir câmeras com isca (apenas Mata Atlântica).
#   - Intervalo de independência de 60 minutos, aplicado consistentemente
#     entre redes (a Seção 3 quantifica o que isso custa).
#   - Agrupar câmeras em ARRAYS usando a regra de desenho do Snapshot: câmeras
#     a menos de 200m reduzidas a uma; ligadas em um array se a menos de 5km
#     de uma vizinha; array mantido apenas se >= 10 câmeras.
#   - Uma temporada por array — para as redes multi-anuais Mata
#     Atlântica/WI, cada array é fixado em seu único ano mais bem amostrado
#     (Snapshot já é temporada única e está isento).
#   - Precisão de coordenadas — sítios com coordenadas arredondadas a ~1km
#     foram removidos (suas covariáveis ambientais em buffer de 100m estariam
#     erradas).
# O resultado: "sítio" e "array" significam a mesma coisa independentemente
# da rede.

# ---- 4b. Uma escolha de consistência: apenas mamíferos ----------------------
# A compilação da Mata Atlântica registrou APENAS mamíferos. Manter aves na
# análise tornaria toda ave falsamente AUSENTE em todo array da Mata
# Atlântica — um artefato de dados, não ecologia. Restringimos a análise
# inteira a mamíferos para manter toda espécie comparável entre as três
# redes.

# ---- 4c. Taxonomia: harmonizando nomes entre redes --------------------------
# As três redes foram compiladas em épocas diferentes com convenções
# taxonômicas diferentes. Todo nome é reconciliado ao Mammal Diversity
# Database (MDD) v2.4, o padrão global atual.

tax <- fread(file.path(DATA_DIR, "taxonomy_changes.csv"))
cat("\n-- atualizações de nomes taxonômicos aplicadas (MDD v2.4) --\n")
print(tax[, .(data_name, mdd_accepted, mdd_common, ambiguous, sources)])


## =============================================================================
## 5. DE INSTALAÇÕES A SÍTIOS, ARRAYS E ESFORÇO
## =============================================================================
## Uma INSTALAÇÃO é uma câmera funcionando por um período contínuo. A
## distinção entre uma linha bruta de instalação e um SÍTIO analítico importa
## por duas razões concretas, ambas das quais poderiam distorcer um modelo:
##
##   (a) TROCAS DE CARTÃO — reparar uma câmera no meio da temporada cria uma
##       segunda linha de instalação no mesmo local. Tratado ingenuamente,
##       isso dobra a contribuição do sítio para a amostra E pode fabricar
##       falsas ocasiões "0" (amostrado, não detectado) através da lacuna de
##       manutenção se o histórico de detecção for construído preenchendo o
##       intervalo COMPLETO entre a primeira e a última data de instalação, em
##       vez da janela ativa real de cada instalação. (Este foi um bug real
##       encontrado e corrigido durante o desenvolvimento deste pipeline:
##       preencher min(início)-max(fim) como se a câmera tivesse funcionado
##       continuamente fabricou ~1.200 dias "0" falsos em 35 sítios
##       multi-instalação, principalmente na rede Mata Atlântica — pior caso,
##       171 de 181 dias do intervalo foram fabricados em um sítio.)
##       Concatenamos instalações correspondentes em um único SÍTIO
##       (site_id), e respeitamos a janela ativa de cada instalação ao
##       construir históricos de detecção, de modo que lacunas entre visitas
##       de manutenção sejam corretamente tratadas como "sem esforço", não
##       "amostrado, nada detectado."
##   (b) INSTALAÇÕES SOBREPOSTAS — duas instalações no mesmo local que se
##       sobrepõem no tempo (um erro de backend, ex.: entrada duplicada).
##       Excluímos ambas e as sinalizamos, em vez de adivinhar qual está
##       correta.
##
## O ARRAY — não a câmera individual — é a unidade de replicação que todo
## modelo posterior (efeitos aleatórios de ocupação, modelos comunitários)
## usa. Um array contendo silenciosamente sítios duplicados ou inflados por
## lacunas enviesaria toda estimativa construída sobre ele.

eff <- dep[, .(camera_days = sum(camera_days, na.rm = TRUE), n_deploy = .N),
           by = .(network, array_id, site_id)]
cat("\n== Resumo de esforço ==\n")
cat("sítios distintos:", uniqueN(eff$site_id),
    "| câmeras-dia mediana por sítio:", round(median(eff$camera_days, na.rm = TRUE)), "\n")
cat("-- esforço por rede --\n")
print(eff[, .(sitios = uniqueN(site_id), camera_days = sum(camera_days)), by = network])


## =============================================================================
## 6. HISTÓRICOS DE DETECÇÃO
## =============================================================================
## Modelos de ocupação não usam fotos brutas. Eles usam um HISTÓRICO DE
## DETECÇÃO: para cada sítio, o período de amostragem é dividido em OCASIÕES
## iguais (janelas de 7 dias aqui), e cada ocasião é marcada 1 (espécie
## detectada ao menos uma vez) ou 0 (não detectada mas a câmera estava ativa).
## Essa sequência 1/0/1/1/0... é a matéria-prima de um modelo de ocupação.
##
## Exemplo: cutia-preta no array de Tefé (Amazonas) — uma única temporada
## contínua de 2025 com 30 câmeras e sem lacunas de manutenção, escolhida
## porque o quadro é direto de ler (sem lacunas atípicas de dados no meio da
## temporada).

dh_tbl <- fread(file.path(DATA_DIR, "m1_dethist_table.csv"), colClasses = "character")
cat("\n== Exemplo de histórico de detecção: cutia-preta, array de Tefé ==\n")
cat("Matriz literal 0/1 (primeiras 15 de 20 ocasiões, todas as 30 câmeras;",
    "em branco = câmera ainda não/não mais ativa):\n")
print(dh_tbl)

# Veja figs/m1_dethist_new_example.png para a versão gráfica desta mesma
# matriz (cinza = ativo/sem detecção, vermelho = detectado, branco = não
# ativo). Vinte de trinta câmeras detectaram a espécie ao menos uma vez; uma
# câmera mostra uma sequência incomumente longa de detecções consecutivas —
# um sinal comportamental real (um indivíduo residente usando aquele local
# repetidamente), não um artefato.


## =============================================================================
## 7. ESCOLHENDO A JANELA DE OCASIÃO
## =============================================================================
## A duração da ocasião é uma ESCOLHA DE MODELAGEM, não uma propriedade dos
## dados. Ela equilibra duas coisas:
##   - Curta demais (ex.: 1 dia): a maioria das ocasiões é 0 simplesmente
##     porque o animal não estava ali NAQUELE DIA, mesmo que use o sítio. A
##     probabilidade de detecção por ocasião é muito baixa e a ocupação é
##     subestimada.
##   - Longa demais (ex.: 30 dias): quase toda ocasião se torna 1 para
##     espécies comuns, e o modelo perde a informação que distingue sítios, e
##     instalações curtas produzem poucas ocasiões para sequer estimar a
##     detecção.
##
## Espécies diferentes podem justificar janelas diferentes (um carnívoro de
## amplo deslocamento e baixa densidade visitando a cada poucas semanas
## precisa de janela mais longa que um roedor residente fotografado toda
## noite), mas uma única janela comum mantém toda espécie em pé de igualdade
## em um pipeline multiespécies.
##
## ADOTAMOS UMA JANELA DE OCASIÃO DE 7 DIAS como padrão do pipeline: curta o
## suficiente para que espécies comuns ainda variem entre ocasiões, longa o
## suficiente para que as instalações curtas comuns em um conjunto de dados
## multi-rede ainda produzam várias ocasiões.

cat("\n== Janela de ocasião ==\nPadrão do pipeline: ocasiões de 7 dias ",
    "(veja o bloco de comentário acima para o raciocínio).\n")


## =============================================================================
## 8. COMO É A COMUNIDADE COMBINADA
## =============================================================================
## Um primeiro olhar sobre a comunidade de mamíferos nas três redes
## combinadas. Para cada espécie calculamos dois resumos simples:
##   - Taxa de detecção : detecções por 100 câmeras-noite (sinal bruto)
##   - Ocupação ingênua : % de sítios onde a espécie foi detectada ao menos
##                         uma vez (não corrigida para detecção imperfeita)

modeled <- fread(file.path(DATA_DIR, "m1_modeled_species.csv"))
cat("\n== Comunidade combinada ==\n")
cat("espécies modeladas (>=30 detecções, >=8 sítios):", nrow(modeled), "\n")
cat("(veja figs/m1_xmas_mammals.png para a figura completa de",
    "taxa de detecção vs. ocupação ingênua)\n")

# A ocupação ingênua sempre SUBESTIMA a ocupação verdadeira — uma espécie
# pode estar presente mas não detectada em nenhuma visita. Corrigir isso é o
# trabalho do Módulo 3 (modelos de ocupação de espécie única).


## =============================================================================
## RESUMO — o que os resultados significam
## =============================================================================
## - A estrutura própria do Snapshot Brasil (instalações, placenames,
##   sequências, arrays) e o tipo/altura de suas armadilhas fotográficas são
##   todas instalações padrão ao nível do solo, sem conjunto em dossel a
##   excluir.
## - O limiar de independência de 60 minutos, adotado por consistência com a
##   rede da Mata Atlântica, tem um custo real e mensurável concentrado em
##   espécies gregárias e que permanecem no local — uma troca que vale a pena
##   conhecer, não uma suposição oculta.
## - Temos um único conjunto de dados harmonizado e APENAS DE MAMÍFEROS,
##   oriundo de três redes, abrangendo um amplo gradiente ambiental pelo
##   Brasil.
## - As unidades analíticas (sítios, arrays, históricos de detecção)
##   significam a mesma coisa nas três fontes — acertar a etapa de
##   instalação-para-sítio (respeitando janelas reais de câmera ativa, não o
##   intervalo total entre visitas de manutenção) importa porque todo modelo
##   posterior trata o array como sua unidade de replicação.
## - Uma janela de ocasião de 7 dias é o padrão do pipeline, escolhida para
##   equilibrar informação contra as instalações curtas nos dados combinados.
##
## Próximo módulo: Covariáveis e Resumo da Amostragem — extração de camadas
## de hábitat, clima, perturbação e áreas protegidas, e triagem de
## colinearidade.
## =============================================================================
