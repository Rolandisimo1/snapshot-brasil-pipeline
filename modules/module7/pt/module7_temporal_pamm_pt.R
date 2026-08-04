## =============================================================================
## Módulo 7 — Interações Temporais: Sobreposição Diária & Tempo desde o Evento (PAMM)
## Pipeline multi-rede de armadilhas fotográficas Snapshot Brasil
## =============================================================================
##
## AUTOCONTIDO: este script reconstrói todas as tabelas e figuras do relatório
## renderizado a partir de tabelas de resumo/curvas pré-computadas, sem
## reexecutar os kernels de atividade computacionalmente pesados (activity::fitact,
## 1000 bootstraps cada) ou os ajustes de risco PAMM (GAMs exponenciais por partes
## do mgcv) -- exatamente seguindo a convenção do próprio relatório ("mostrado
## como código de referência, não reexecutado"). Coloque este script na mesma
## pasta dos arquivos de dados listados abaixo e execute-o do início ao fim.
##
## Arquivos de entrada necessários (em DATA_DIR, padrão "data" ao lado deste script):
##   temporal_feasibility.csv        -- contagens de detecções com timestamp por
##                                       espécie e aprovação/reprovação da barra de
##                                       sobreposição diária
##   hazard_feasibility_verdict.csv  -- pareamento do melhor respondente bem
##                                       distribuído por driver e veredito de
##                                       viabilidade do modelo de risco
##   diel_overlap.csv                -- coeficientes de sobreposição Dhat4 + testes
##                                       de diferença por bootstrap para os pares de
##                                       espécies focais
##   diel_density_curves.csv         -- curvas de densidade de atividade (0-24h) para
##                                       as figuras de sobreposição diária cão/humano/felídeo
##   allspecies_density.csv          -- curvas de densidade de atividade para todas as
##                                       44 espécies modeladas (figuras de pilha de
##                                       atividade da comunidade)
##   activity_group_key.csv          -- pertencimento de espécies aos 4 grupos de
##                                       habitat da PCA
##   allspecies_meta.csv             -- espécie -> ordem + grupo + n, metadados
##                                       prontos para a pilha de atividade geral
##   biome_meta.csv                  -- espécie -> bioma + ordem + n, para a pilha
##                                       de atividade por bioma
##   jaguar_activity_curve.csv       -- curva de densidade de atividade da própria
##                                       onça-pintada (curva de referência sobreposta
##                                       na figura 2x2 do cateto)
##   peccary_jaguar_test.csv         -- teste de sobreposição diária do Cateto,
##                                       locais com onça-pintada presente vs ausente
##   peccary_temp_test.csv           -- mesmo teste, locais quentes vs frios (controle
##                                       da hipótese de calor)
##   peccary_temp_curves.csv         -- densidade de atividade do Cateto, locais
##                                       quentes vs frios
##   peccary_2x2_curves.csv          -- densidade de atividade do Cateto, desenho
##                                       2x2 onça-pintada x temperatura
##   peccary_2x2_meta.csv            -- n/locais/fração noturna por célula do
##                                       desenho 2x2
##   human_nocturnality.csv          -- mudança na fração noturna por espécie,
##                                       locais de baixa vs alta pressão humana
##   seasonal_shift.csv              -- mudança na fração noturna por espécie,
##                                       estação chuvosa vs seca
##   lunar_effect.csv                -- mudança de atividade por espécie, noites de
##                                       lua escura vs lua cheia
##   pamm_results.csv                -- coeficientes de risco PAMM agrupados + por
##                                       subprojeto (tsh_coef, HR/semana, p) para
##                                       cada par driver-respondente focal
##   influence_check.csv             -- checagem de influência por câmera (remover
##                                       as 2 câmeras principais, reajustar) para os
##                                       dois sinais PAMM significativos em nível de
##                                       subprojeto
##   coati_piracicaba_percamera.csv  -- detalhamento da taxa de detecção por câmera
##                                       para o sinal de quati de Piracicaba
##
## Para de fato REAJUSTAR os kernels de atividade ou os modelos de risco PAMM a
## partir das detecções brutas (horas, não segundos -- mostrado apenas como
## referência), defina RUN_LIVE abaixo como TRUE e veja mod7_code/*.R para os
## scripts exatos de ajuste; isso requer os pacotes R `activity`, `overlap` e
## `mgcv` além das tabelas brutas de dados de eventos pareados (PED) por par,
## que não são disponibilizadas aqui, já que os ajustes não se destinam a ser
## reexecutados rotineiramente.
##
## Requer: data.table, ggplot2, ggrepel, patchwork.

# NOTA DE LOCALE: para que os acentos (á, ã, ç, õ) apareçam corretamente nos
# graficos gerados, execute este script com um locale UTF-8, ex.:
#   Sys.setlocale("LC_CTYPE", "pt_BR.UTF-8")   # ou defina LC_ALL=pt_BR.UTF-8 no shell antes de chamar Rscript
# Em locale "C" (o padrao de muitos containers), os caracteres acentuados nos
# titulos/rotulos das figuras aparecem corrompidos (ex.: "di..ria" em vez de
# "diaria"), embora o texto no arquivo-fonte esteja em UTF-8 correto.
tryCatch(Sys.setlocale("LC_CTYPE", "pt_BR.UTF-8"), warning = function(w) invisible(NULL))

RUN_LIVE <- FALSE   # TRUE reajusta tudo a partir dos dados brutos (horas) via mod7_code/*.R;
                    # FALSE (padrão) carrega as tabelas pré-computadas acima e reconstrói
                    # cada figura/tabela a partir delas, seguindo a própria convenção do relatório.

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---- 0. Caminhos ---------------------------------------------------------------
DATA_DIR <- "data"
FIGS_DIR <- "figs"
stopifnot(dir.exists(DATA_DIR))
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR)

# verde = positivo/alto, vermelho = negativo/baixo; um punhado consistente de tons
# usado em todas as figuras deste script.
COL_POS   <- "#1a9850"
COL_NEG   <- "#d73027"
COL_NEUTRAL <- "#4575b4"
COL_DOG   <- "#6a3d9a"
COL_HUMAN <- "#e08214"

read_d <- function(f) fread(file.path(DATA_DIR, f))


## =============================================================================
## PASSO 1: VIABILIDADE -- essa análise pode honestamente ser executada?
## =============================================================================
## Antes de ajustar qualquer coisa: conte os eventos qualificados por espécie, e se
## um modelo dependesse de um punhado de eventos em uma ou duas câmeras, não o ajuste.
## Reportar *por que* uma análise não foi executada é em si um resultado legítimo.

feas <- read_d("temporal_feasibility.csv")
vd   <- read_d("hazard_feasibility_verdict.csv")
cat(sprintf("detecções com timestamp na tabela de espécies focais: %s\n", format(sum(feas$n_timestamped), big.mark=",")))
cat(sprintf("espécies que ultrapassam o limite de sobreposição diária (>=50 eventos): %d\n\n", sum(feas$n_timestamped >= 50)))
cat("Viabilidade do modelo de risco por driver:\n")
print(vd[, .(driver, total_timestamped, best_responder, best_events, best_locations, hazard_feasible)])

# ---- figura de viabilidade: barra de sobreposição diária (todas as espécies) + barra de risco (3 drivers) ----
feas_top <- feas[order(-n_timestamped)][1:10]
feas_top$species <- factor(feas_top$species, levels=rev(feas_top$species))
pA <- ggplot(feas_top, aes(x=species, y=n_timestamped)) +
  geom_col(fill=COL_NEUTRAL) + coord_flip() +
  geom_hline(yintercept=50, color=COL_NEG, linetype="dashed") +
  geom_text(aes(label=n_timestamped), hjust=-0.15, size=3) +
  labs(x=NULL, y="detecções com timestamp (rede inteira)",
       title="A. Viabilidade de sobreposição diária") +
  theme_minimal(base_size=10) + expand_limits(y=max(feas_top$n_timestamped)*1.12)

vd$fill_col <- ifelse(vd$hazard_feasible, COL_POS, COL_NEG)
pB <- ggplot(vd, aes(x=driver, y=best_events, fill=fill_col)) +
  geom_col() + scale_fill_identity() +
  geom_text(aes(label=sprintf("%d ev / %d loc", best_events, best_locations)), vjust=-0.4, size=3) +
  labs(x=NULL, y="melhor par de eventos driver->nativa bem distribuídos",
       title="B. Viabilidade do modelo de risco (melhor respondente por driver)") +
  theme_minimal(base_size=10) + expand_limits(y=max(vd$best_events)*1.15)

p_feas <- pA + pB
ggsave(file.path(FIGS_DIR, "temporal_feasibility.png"), p_feas, width=13, height=6, dpi=150)
cat("\nSalvo figs/temporal_feasibility.png\n")


## =============================================================================
## STEP 2: SOBREPOSIÇÃO DE ATIVIDADE DIÁRIA
## =============================================================================
## Compara QUANDO (tempo solar, 0-24h) cada espécie está ativa, usando o
## coeficiente de sobreposição Dhat4 de Ridout & Linkie (0 = sem sobreposição,
## 1 = mesmo horário). As curvas de densidade aqui já foram precomputadas
## (activity::fitact, 1000 bootstraps, código de referência em
## mod7_code/fit_diel.R); este passo apenas as replota.

ov  <- read_d("diel_overlap.csv")
dc  <- read_d("diel_density_curves.csv")

dog_ov <- ov[sp1 == "Domestic Dog"][order(Dhat4)]
cat("Sobreposição diária (Dhat4) entre o cão e cada nativo (menor = menos sobreposição):\n")
print(dog_ov[, .(native=sp2, n_native=n2, Dhat4, p_differ)])

# ---- sobreposição cão-nativo: três nativos abrangendo a faixa de sobreposição ----
sel_native <- c("Spotted Paca", "Collared Peccary", "Gray Brocket")
dc_long <- melt(dc[, c("hour","Domestic Dog", sel_native), with=FALSE], id.vars="hour",
                 variable.name="species", value.name="density")
p_dog <- ggplot(dc_long, aes(x=hour, y=density, color=species)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=setNames(c(COL_DOG, "#238b45","#fdae61","#3182bd"),
                                      c("Domestic Dog", sel_native))) +
  labs(x="hora (tempo solar)", y="densidade de atividade",
       title="Sobreposição de atividade diária: Cão doméstico vs três nativos") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR, "activity_overlap_dog.png"), p_dog, width=9, height=5.5, dpi=150)
cat("Salvo figs/activity_overlap_dog.png\n")

# ---- Jaguatirica x Onça-parda ----
dc_felid <- melt(dc[, c("hour","Puma","Ocelot"), with=FALSE], id.vars="hour",
                  variable.name="species", value.name="density")
p_felid <- ggplot(dc_felid, aes(x=hour, y=density, color=species)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Puma"="#e6550d","Ocelot"="#31a354")) +
  labs(x="hora (tempo solar)", y="densidade de atividade",
       title="Sobreposição de atividade diária: Jaguatirica x Onça-parda (Dhat4=0.83, p<0.001)") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR, "activity_overlap_felids.png"), p_felid, width=8, height=5.5, dpi=150)
cat("Salvo figs/activity_overlap_felids.png\n")


## =============================================================================
## ETAPA 3: O DIA DE TODA A COMUNIDADE -- atividade empilhada, ponderada por evento
## =============================================================================
## Altura = somada, atividade ponderada por evento entre espécies; cor = ordem
## taxonômica, tonalidade = espécie dentro da ordem. Construído a partir de densidades
## de atividade por espécie pré-calculadas (mod7_data/allspecies_density.csv), ponderadas
## pela contagem de detecção de cada espécie em toda a rede (n) de forma que as espécies
## comuns dominem a pilha visível da mesma forma que dominam o registro de câmeras.

allsp_dens <- read_d("allspecies_density.csv")
allsp_meta <- read_d("allspecies_meta.csv")
group_key  <- read_d("activity_group_key.csv")

order_pal <- c(Artiodactyla="#3182bd", Carnivora="#e6550d", Rodentia="#31a354",
               Didelphimorphia="#756bb1", Cingulata="#636363", Pilosa="#a1d99b",
               Perissodactyla="#fd8d3c", Primates="#fdae6b", Lagomorpha="#bdbdbd")

build_stack <- function(dens_wide, meta, weight_col="n", title="") {
  # dens_wide: matriz hora x espécie de densidade; meta: espécie -> ordem/grupo/n
  long <- melt(dens_wide, id.vars="hour", variable.name="species", value.name="density")
  long[, species := as.character(species)]
  long <- merge(long, meta[, .(species=sci, order, w=get(weight_col))], by="species")
  long[, weighted := density * w]
  agg <- long[, .(activity=sum(weighted)), by=.(hour, order)]
  agg[, order := factor(order, levels=names(order_pal))]
  ggplot(agg, aes(x=hour, y=activity, fill=order)) +
    geom_area(position="stack", alpha=0.85) +
    scale_fill_manual(values=order_pal, drop=FALSE) +
    labs(x="hora (tempo solar)", y="atividade ponderada por evento (somada)", fill="Ordem", title=title) +
    theme_minimal(base_size=11)
}

p_overall <- build_stack(allsp_dens, allsp_meta, title="Atividade diária total, todas as 44 espécies modeladas")
ggsave(file.path(FIGS_DIR, "activity_stack_overall.png"), p_overall, width=10, height=6, dpi=150)
cat("Salvo figs/activity_stack_overall.png\n")

# ---- pelos 4 grupos de hábitat da PCA do Módulo 3 ----
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
cat("Salvo figs/activity_stack_bygroup.png\n")

# ---- tabela-chave: composição de espécies dos quatro grupos ----
group_summary <- group_key[, .N, by=.(group, group_label)][order(group)]
setnames(group_summary, "N", "species_n")
cat("\nComposição de espécies dos quatro grupos de comunidade de hábitat:\n")
print(group_summary)

# ---- por bioma WWF ----
biome_meta <- read_d("biome_meta.csv")
biomes_to_show <- setdiff(unique(biome_meta$biome), "Mangrove")  # poucos eventos, omitido conforme relatório
p_by_biome_list <- list()
for (bm_name in biomes_to_show) {
  bmeta <- biome_meta[biome == bm_name]
  sp_in_biome <- bmeta$sci
  cols_present <- intersect(c("hour", unique(bmeta$common)), names(allsp_dens))
  # colunas de allsp_density são nomes científicos, não nomes comuns -- combinar por sci
  cols_present <- intersect(c("hour", sp_in_biome), names(allsp_dens))
  if (length(cols_present) < 2) next
  p_by_biome_list[[bm_name]] <- build_stack(allsp_dens[, ..cols_present], bmeta, title=bm_name)
}
p_biome <- wrap_plots(p_by_biome_list, ncol=2) + plot_layout(guides="collect")
ggsave(file.path(FIGS_DIR, "activity_stack_bybiome.png"), p_biome, width=13, height=9, dpi=150)
cat("Salvo figs/activity_stack_bybiome.png\n")


## =============================================================================
## PASSO 4: CATETO -- LOCAIS COM ONÇA-PINTADA PRESENTE VS. AUSENTE
## =============================================================================
## O contraste dentro da mesma espécie é o teste honesto para uma resposta
## comportamental (em oposição a uma diferença intrínseca de nicho): a espécie
## se comporta de forma diferente dependendo da presença do seu principal
## predador?

jt <- read_d("peccary_jaguar_test.csv")
cat(sprintf("Sobreposição de atividade diária do cateto por presença de onça-pintada: Dhat4=%.2f (p=%.4f); noturno com presença %.0f%% vs ausência %.0f%%\n",
            jt$Dhat4, jt$p_differ, 100*jt$noct_present, 100*jt$noct_absent))
tt2 <- read_d("peccary_temp_test.csv")
cat(sprintf("Mesmo teste, agora por temperatura: Dhat4=%.2f (p=%.4f); noturno quente %.0f%% vs frio %.0f%%\n",
            tt2$Dhat4, tt2$p_differ, 100*tt2$noct_hot, 100*tt2$noct_cool))

# contraste onça-pintada presente/ausente + quente/frio, lado a lado
pjt_curves <- read_d("peccary_2x2_curves.csv")  # reaproveita a dimensão onça-pintada do delineamento 2x2
jag_curve  <- read_d("jaguar_activity_curve.csv")
ptemp_curves <- read_d("peccary_temp_curves.csv")

d1 <- melt(pjt_curves[, .(hour, `Jaguar present`=jaguar_hot, `Jaguar absent`=nojaguar_hot)],
           id.vars="hour", variable.name="condition", value.name="density")
p_jag <- ggplot(d1, aes(x=hour, y=density, color=condition)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Jaguar present"="#e6550d", "Jaguar absent"="#31a354")) +
  labs(x="hora (tempo solar)", y="densidade de atividade", color=NULL,
       title="Cateto: onça-pintada presente vs. ausente") +
  theme_minimal(base_size=10)

d2 <- melt(ptemp_curves[, .(hour, Hot=hot, Cool=cool)], id.vars="hour",
           variable.name="condition", value.name="density")
p_temp <- ggplot(d2, aes(x=hour, y=density, color=condition)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("Hot"="#e6550d", "Cool"="#3182bd")) +
  labs(x="hora (tempo solar)", y="densidade de atividade", color=NULL,
       title="Cateto: locais quentes vs. frios (controle de temperatura)") +
  theme_minimal(base_size=10)

p_jag_temp <- p_jag + p_temp
ggsave(file.path(FIGS_DIR, "peccary_jaguar_diel_contrast.png"), p_jag_temp, width=13, height=5.5, dpi=150)
cat("Saved figs/peccary_jaguar_diel_contrast.png\n")

# ---- delineamento 2x2: onça-pintada (sim/não) x temperatura (quente/frio) ----
meta2x2 <- read_d("peccary_2x2_meta.csv")
cat("\nTamanhos das células do delineamento 2x2 (n / locais / fração noturna):\n")
print(meta2x2)

# A célula "onça-pintada presente + frio" (44 detecções de cateto, 2 locais) NÃO
# é tratada como uma condição real: esses 2 locais "frios" estão em 23.2-23.4C,
# bem na fronteira do corte quente/frio, e as onças-pintadas neste conjunto de
# dados estão quase totalmente ausentes dos locais frios (apenas 2 das 244
# detecções de onça-pintada em toda a rede). Em vez de plotar uma curva de
# densidade construída a partir dessa amostra escassa e limítrofe -- o que
# produz um pico espúrio, um artefato da suavização por kernel sobre ~44
# pontos --, este painel é deliberadamente deixado em branco com uma anotação
# explicativa, condizente com o próprio delineamento declarado no relatório.
d2x2 <- melt(pjt_curves[, .(hour, nojaguar_cool, nojaguar_hot, jaguar_hot)],
             id.vars=c("hour"), variable.name="cell", value.name="density")
d2x2[, jaguar_status := ifelse(grepl("^nojaguar", cell), "Jaguar absent", "Jaguar present")]
d2x2[, temp := ifelse(grepl("hot$", cell), "Hot", "Cool")]
d2x2[, jaguar_status := factor(jaguar_status, levels=c("Jaguar absent","Jaguar present"))]
d2x2[, temp := factor(temp, levels=c("Cool","Hot"))]

# a própria coluna "jaguar" de jag_curve (um valor de densidade) não pode ter o
# mesmo nome da variável de facetamento -- facet_grid() resolve a associação de
# painel por camada, e uma coluna com o mesmo nome nos dados de outra camada
# (por mais que não tenha relação) é incorporada como se fosse a chave de
# facetamento, explodindo silenciosamente em um painel para cada valor único de
# densidade.
# A curva de referência da onça-pintada só é desenhada nos 3 painéis reais; o
# painel em branco "Onça-pintada presente x Frio" não recebe nenhuma curva,
# condizente com o delineamento original.
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
  label = "onças-pintadas praticamente ausentes dos locais frios\n\napenas 2 das 244 detecções de onça-pintada em toda a rede vêm de locais frios;\nos 2 locais 'frios' de cateto estão em 23.2-23.4C,\nbem na fronteira do corte quente/frio -- não é uma condição real"
)

p_2x2 <- ggplot(d2x2, aes(x=hour, y=density)) +
  geom_line(color=COL_NEUTRAL, linewidth=1) +
  geom_line(data=jag_curve_full, aes(x=hour, y=density_ref), color="black", linetype="dotted", linewidth=0.8) +
  geom_text(data=blank_note, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            size=3, color="grey40", lineheight=0.9) +
  facet_grid(jaguar_status ~ temp, drop=FALSE) +
  labs(x="hora (tempo solar)", y="densidade de atividade do cateto",
       title="Atividade do cateto: presença de onça-pintada x temperatura",
       subtitle="linha pontilhada = curva de atividade da própria onça-pintada (locais de floresta quente, n=242)") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "peccary_jaguar_temp_2x2.png"), p_2x2, width=9, height=8, dpi=150)
cat("Saved figs/peccary_jaguar_temp_2x2.png\n")


## =============================================================================
## PASSO 5: TRÊS EIXOS ADICIONAIS QUE OS DADOS MESCLADOS PODEM CARREGAR
## =============================================================================
## (1) Nocturnalidade induzida pelo ser humano (estilo Gaynor): uma espécie nativa
##     se torna mais noturna sob alta pressão humana?
## (2) Mudança sazonal de atividade: curvas diárias de estação úmida vs seca por espécie.
## (3) Efeito lunar: espécies fortemente noturnas evitam noites de luar brilhante?

# ---- (1) nocturnalidade induzida pelo ser humano ----
hn <- read_d("human_nocturnality.csv")
hn[, sig := p_differ < 0.05]
hn[, common := factor(common, levels=common[order(dnoct)])]
p_hn <- ggplot(hn, aes(x=common, y=dnoct, color=sig)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_segment(aes(xend=common, y=0, yend=dnoct)) +
  geom_point(size=3) +
  scale_color_manual(values=c("TRUE"=COL_NEG, "FALSE"="grey60"), guide="none") +
  coord_flip() +
  labs(x=NULL, y="mudança na fração noturna (pressão humana alta - baixa)",
       title="Nocturnalidade induzida pelo ser humano: 12 espécies com dados suficientes") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "human_nocturnality.png"), p_hn, width=9, height=6, dpi=150)
n_toward_night <- sum(hn$dnoct > 0); n_away <- sum(hn$dnoct < 0)
cat(sprintf("Nocturnalidade induzida pelo ser humano: %d de %d espécies deslocam-se em direção à noite, %d se afastam\n", n_toward_night, nrow(hn), n_away))
cat("Salvo figs/human_nocturnality.png\n")

# ---- (2) mudança sazonal ----
ss <- read_d("seasonal_shift.csv")
ss[, dnoct := noct_wet - noct_dry]
ss[, sig := p_differ < 0.05]
top6 <- ss[sig == TRUE][order(-abs(dnoct))][1:min(6,.N)]
cat(sprintf("\nMudança sazonal: %d de %d espécies mostram diferença significativa entre estação úmida e seca\n", sum(ss$sig), nrow(ss)))
print(top6[, .(common, n_wet, n_dry, Dhat4, dnoct)])
p_seasonal <- ggplot(top6, aes(x=reorder(common, dnoct), y=dnoct)) +
  geom_col(fill=COL_NEUTRAL) + coord_flip() +
  labs(x=NULL, y="mudança na fração noturna (úmida - seca)",
       title="Mudança sazonal de atividade: 6 maiores mudanças significativas") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "seasonal_shift.png"), p_seasonal, width=9, height=5, dpi=150)
cat("Salvo figs/seasonal_shift.png\n")

# ---- (3) efeito lunar ----
# illum_shift = a iluminação média das noites de detecção de cada espécie menos a
# linha de base noturna da comunidade (mean_illum_null); negativo = lunar-fóbico (menos
# detecções em noites brilhantes), positivo = lunar-fílico. A significância é a
# coluna `p` (teste bootstrap de illum_shift contra 0).
le <- read_d("lunar_effect.csv")
sig_le <- le[p < 0.05]
cat(sprintf("\nEfeito lunar: %d de %d espécies fortemente noturnas mostram resposta significativa à lua\n", nrow(sig_le), nrow(le)))
print(sig_le[, .(common, n_night, n_bright, n_dark, illum_shift, p)])
le[, common := factor(common, levels=common[order(illum_shift)])]
p_lunar <- ggplot(le, aes(x=common, y=illum_shift, color=p<0.05)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_segment(aes(xend=common, y=0, yend=illum_shift)) +
  geom_point(size=3) +
  scale_color_manual(values=c("TRUE"=COL_NEG, "FALSE"="grey60"), guide="none") +
  coord_flip() +
  labs(x=NULL, y="mudança na iluminação lunar vs linha de base noturna da comunidade (negativo = lunar-fóbico)",
       title="Efeito da iluminação lunar: 16 espécies fortemente noturnas") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "lunar_effect.png"), p_lunar, width=9, height=7, dpi=150)
cat("Salvo figs/lunar_effect.png\n")


## =============================================================================
## ETAPA 6: MODELOS DE RISCO DE TEMPO DESDE O EVENTO (PAMM)
## =============================================================================
## A taxa de detecção de curto prazo de uma espécie nativa muda nos dias após um
## humano ou cão passar pela mesma câmera? Modelo aditivo por partes misto
## (GAM binomial): y ~ tsh + s(hour,cc) + subproject + s(camera,re). Os coeficientes
## aqui são carregados dos ajustes pré-computados (mod7_code/fit_pamm.R, apenas referência).

res <- read_d("pamm_results.csv")

# ---- humano: agrupado, em toda a rede ----
pooled_h <- res[scope == "pooled" & grepl("Human", pair)]
cat("Risco agrupado de tempo desde o humano (nenhuma nativa mostra efeito significativo em toda a rede):\n")
print(pooled_h[, .(pair, events, HR_per_week, p)])
cat(sprintf("efeitos humanos agrupados significativos (p<0.05): %d de %d\n\n", sum(pooled_h$p < 0.05), nrow(pooled_h)))

pooled_h[, sig := p < 0.05]
p_pamm_human <- ggplot(pooled_h, aes(x=reorder(pair, HR_per_week), y=HR_per_week)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey50") +
  geom_point(size=3, color=COL_NEUTRAL) + coord_flip() +
  labs(x=NULL, y="razão de risco por semana (1 = sem efeito)",
       title="Risco de tempo desde o humano: quatro nativas focais, em toda a rede") +
  theme_minimal(base_size=10)
ggsave(file.path(FIGS_DIR, "pamm_human.png"), p_pamm_human, width=9, height=5, dpi=150)
cat("Salvo figs/pamm_human.png\n")

# ---- Quati em Piracicaba: um sinal local que vale a pena investigar, e sua checagem de influência ----
inf <- read_d("influence_check.csv")
coati_inf <- inf[grepl("Coati", signal)]
cat("\nChecagem de influência por câmera, sinal do quati em Piracicaba:\n")
print(coati_inf[, .(signal, full_p, drop2_p, robust)])

coati_pc <- read_d("coati_piracicaba_percamera.csv")
coati_pc_long <- melt(coati_pc, id.vars=c("camera","dets","bins"),
                       measure.vars=c("rate_early","rate_late"),
                       variable.name="period", value.name="rate")
coati_pc_long[, period := ifelse(period=="rate_early", "primeiros dias após o humano", "posterior")]
p_coati_a <- ggplot(pooled_h[grepl("Coati", pair)] , aes(x=1, y=HR_per_week)) +
  geom_point(size=4, color=COL_NEUTRAL) +
  annotate("text", x=1, y=pooled_h[grepl("Coati",pair)]$HR_per_week, label="agrupado", vjust=-1.5) +
  labs(title="Taxa de detecção do quati\nvs dias desde o humano", x=NULL, y="HR por semana") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_blank())
p_coati_b <- ggplot(coati_pc_long, aes(x=reorder(camera, -dets), y=rate, fill=period)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c("primeiros dias após o humano"=COL_HUMAN, "posterior"="grey60")) +
  labs(x=NULL, y="detecções por 1.000 bins", fill=NULL,
       title="Detalhamento por câmera: concentração do sinal") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=30, hjust=1))
p_coati <- p_coati_a + p_coati_b + plot_layout(widths=c(1,2))
ggsave(file.path(FIGS_DIR, "pamm_coati.png"), p_coati, width=11, height=5.5, dpi=150)
cat("Salvo figs/pamm_coati.png\n")

# ---- Cães e Paca ----
dog_res <- res[grepl("Dog", pair)]
cat("\nModelo de tempo desde o cão para Paca (agrupado e melhor subprojeto amostrado):\n")
print(dog_res[, .(scope, events, HR_per_week, p)])

dog_pooled <- dog_res[scope=="pooled"]
dog_sub    <- dog_res[scope!="pooled"]
p_dog_pamm <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", color="grey50") +
  geom_point(data=dog_pooled, aes(x="pooled", y=HR_per_week), size=4, color=COL_DOG) +
  geom_point(data=dog_sub, aes(x="best subproject\n(dashed line)", y=HR_per_week), size=4, color=COL_DOG, shape=17) +
  labs(x=NULL, y="razão de risco por semana (1 = sem efeito)",
       title="Taxa de detecção da Paca vs dias desde um cão") +
  theme_minimal(base_size=11)
ggsave(file.path(FIGS_DIR, "pamm_dog_paca.png"), p_dog_pamm, width=8, height=5.5, dpi=150)
cat("Salvo figs/pamm_dog_paca.png\n")

cat("\n=== Script do Módulo 7 completo: triagem de viabilidade, sobreposição diária, pilhas de atividade,\n")
cat("    contraste cateto/onça-pintada, 3 eixos estendidos, e resultados PAMM todos reconstruídos. ===\n")
cat("Veja module7_temporal_pamm.qmd para a análise completa narrada e citações.\n")
