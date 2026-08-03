## =============================================================================
## Módulo 6 — Interações entre Espécies
## Pipeline de câmeras-trap multi-rede do Snapshot Brasil
## =============================================================================
##
## AUTOCONTIDO: este script constrói todas as tabelas e figuras do relatório
## renderizado a partir dos arquivos brutos do conjunto de dados combinado, mais um pequeno conjunto de tabelas
## de referência que são, elas próprias, o resultado de um trabalho anterior pontual
## (máscara de distribuição das espécies, correlações residuais do GJAM do Módulo 4, e os
## checkpoints de histórico de detecção/covariável em nível de array do Módulo 5/6's
## agregação de array). Duas classes de modelo estão envolvidas:
##
## 1. Os 30 ajustes de interação mediada pela abundância (modelos N-mixture pareados,
## NIMBLE/MCMC) são cada um custoso (até 300.000 iterações ao longo de
## três rodadas de reintentativa escalonadas; o lote completo leva muitas horas
## sequencialmente, e o sandbox no qual este pipeline foi construído não consegue executar
## mais de um desses ajustes por vez). Por padrão, este script CARREGA a
## tabela de coeficientes pré-computada distribuída junto com ele em vez de
## reajustar; defina RUN_LIVE_PAIRS <- TRUE abaixo para reajustar cada par do
## zero usando fit_pair_array.R (distribuído em mod6_code/), que este
## script irá carregar (source) e chamar uma vez por par.
## 2. A matriz de correlação residual do GJAM (todas as 44 espécies) é o resultado
## do Módulo 4, não reajustada aqui; é carregada como uma tabela de referência.
##
## ENTRADAS (colocar em DATA_DIR, padrão "data/", junto com este script):
## all_pairs_array_level_coefs.csv   -- resumo posterior completo (média/dp/lcl/ucl/Rhat)
## para cada parâmetro de todos os 30 pares em nível de array
## convergence_diagnostics_summary.csv -- Rhat máximo por par e indicador de convergência
## gjam_residual_correlation.csv     -- matriz de correlação residual de 44x44 espécies
## (saída do GJAM do Módulo 4, usada como está)
## species_common_names.csv          -- tabela de referência científico -> nome comum
## pair_topology.csv                 -- as frações de detecção dos 20 pares originais
## e tipo de topologia, para referência
## Nota: este módulo ajusta sobre a tabela de covariáveis em nível de array do Módulo 5, que é, ela própria, as mesmas
## 11 covariáveis do arquivo em nível de câmera, reduzidas ao conjunto de 6 covariáveis
## ARRAY_SET no momento do ajuste -- não um arquivo separadamente reduzido.
## Defina RUN_LIVE_PAIRS <- TRUE (e forneça mod6_code/fit_pair_array.R, além dos
## arquivos pair_inputs/*.rds que ele espera) para reajustar do zero em vez disso.
## =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})

DATA_DIR <- "data"
FIGS_DIR <- "figs"
dir.create(FIGS_DIR, showWarnings = FALSE)

RUN_LIVE_PAIRS <- FALSE   # # TRUE reajusta todos os 30 pares via NIMBLE/MCMC (muitas horas); FALSE carrega os precomputados

## =============================================================================
## PASSO 1: CARREGAR RESULTADOS (pré-computados por padrão; reajustados se RUN_LIVE_PAIRS)
## =============================================================================
if (RUN_LIVE_PAIRS) {
  source(file.path("mod6_code", "fit_pair_array.R"))  # # espera mod6_array/pair_inputs/*.rds
  stop("Live refit path invoked -- see fit_pair_array.R and the module's methods note ",
       "for the full escalating-retry MCMC procedure; this path is not run automatically ",
       "because a single pair can take hours and pairs must be fit strictly sequentially.")
} else {
  all_pairs <- fread(file.path(DATA_DIR, "all_pairs_array_level_coefs.csv"))
  conv <- fread(file.path(DATA_DIR, "convergence_diagnostics_summary.csv"))
  cat("Loaded precomputed array-level fits:", uniqueN(all_pairs$pair_idx, na.rm=TRUE) * 3, "topology groups,",
      "converged:", sum(conv$converged), "of", nrow(conv), "pairs\n")
}

common_names <- fread(file.path(DATA_DIR, "species_common_names.csv"))
name_map <- setNames(common_names$common, common_names$species)
to_common <- function(sci) ifelse(sci %in% names(name_map), name_map[sci], sci)

## =============================================================================
## PASSO 2: O COEFICIENTE DE INTERAÇÃO (gamma0) PARA TODOS OS 30 PARES
## =============================================================================
## gamma0 é o efeito da log-abundância da espécie motriz sobre a log-abundância da espécie respondente
## no modelo N-mixture pareado -- a estimativa real da interação (gamma1,
## gamma2, gamma3 são termos de interação com temperatura testados separadamente; veja o
## texto do módulo para saber por que eles não são a quantidade principal).
gamma0 <- all_pairs[param == "gamma0"]
gamma0 <- merge(gamma0, conv[, .(driver, responder, source, converged)],
                 by = c("driver", "responder", "source"), all.x = TRUE)
gamma0[, sig := (lcl > 0) | (ucl < 0)]
gamma0[, driver_c := to_common(driver)]
gamma0[, responder_c := to_common(responder)]
cat("Significant pairs (95% CI excludes 0):", sum(gamma0$sig), "of", nrow(gamma0), "\n")

sig_pairs <- gamma0[sig == TRUE, .(driver_c, responder_c, source, mean, lcl, ucl)]
fwrite(sig_pairs, "gamma0_significant_pairs.csv")

## =============================================================================
## PASSO 3: FOREST PLOT -- TODOS OS 30 PARES, ORDENADOS
## =============================================================================
gamma0[, label := paste0(driver_c, " \u2192 ", responder_c)]
gamma0_plot <- gamma0[order(mean)]
gamma0_plot[, y := .I]
gamma0_plot[, col := fifelse(!sig, "grey80", fifelse(mean > 0, "#2166ac", "#b2182b"))]

p_forest <- ggplot(gamma0_plot, aes(y = reorder(label, mean))) +
  geom_segment(aes(x = lcl, xend = ucl, yend = reorder(label, mean), color = col), linewidth = 1) +
  geom_point(aes(x = mean, color = col), size = 2) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_color_identity() +
  labs(x = "Interaction coefficient (\u03b30): driver's log-abundance effect on responder's log-abundance",
       y = NULL, title = "Abundance-mediated interaction strength, all 30 array-level pairs") +
  theme_minimal(base_size = 9)
ggsave(file.path(FIGS_DIR, "gamma0_forest_all_pairs.png"), p_forest, width = 9, height = 11, dpi = 150)
cat("Saved figs/gamma0_forest_all_pairs.png\n")

## =============================================================================
## ETAPA 4: A REDE DE INTERAÇÃO MEDIADA PELA ABUNDÂNCIA
## =============================================================================
suppressMessages(library(igraph))
edges_abu <- gamma0[, .(from = driver_c, to = responder_c, sig, sign = sign(mean))]
g_abu <- graph_from_data_frame(edges_abu, directed = TRUE)
E(g_abu)$color <- ifelse(!edges_abu$sig, "grey85", ifelse(edges_abu$sign > 0, "#2166ac", "#b2182b"))
E(g_abu)$width <- ifelse(edges_abu$sig, 2.5, 0.6)
E(g_abu)$lty   <- ifelse(edges_abu$sig, 1, 2)
set.seed(42)
layout_abu <- layout_with_fr(g_abu)

png(file.path(FIGS_DIR, "interaction_network_abundance.png"), width = 1600, height = 1300, res = 140)
par(mar = c(4, 1, 3, 1))
plot(g_abu, layout = layout_abu, vertex.size = 14, vertex.color = "#4575b4",
     vertex.label.cex = 0.75, vertex.label.color = "black", edge.arrow.size = 0.5,
     main = "Abundance-mediated species interaction network \u2014 array level, 30 driver\u2192responder pairs")
legend("bottomleft", legend = c("Significant, positive", "Significant, negative", "Not significant"),
       col = c("#2166ac", "#b2182b", "grey85"), lty = c(1,1,2), lwd = c(2.5,2.5,0.6), bty = "n", cex = 0.8)
dev.off()
cat("Saved figs/interaction_network_abundance.png\n")

## =============================================================================
## PASSO 5: A REDE DE CO-OCORRÊNCIA RESIDUAL DO GJAM
## =============================================================================
gjam_corr <- fread(file.path(DATA_DIR, "gjam_residual_correlation.csv"))
sp_gjam <- gjam_corr$species
mat_gjam <- as.matrix(gjam_corr[, -"species"])
rownames(mat_gjam) <- sp_gjam

THRESH <- 0.4
edges_gjam <- data.table()
for (i in seq_along(sp_gjam)) for (j in seq_len(i-1)) {
  v <- mat_gjam[i, j]
  if (!is.na(v) && abs(v) > THRESH) {
    edges_gjam <- rbind(edges_gjam, data.table(from = to_common(sp_gjam[i]), to = to_common(sp_gjam[j]), r = v))
  }
}
cat("GJAM edges above |r|>", THRESH, ":", nrow(edges_gjam), "of", choose(length(sp_gjam), 2), "possible pairs\n")

g_gjam <- graph_from_data_frame(edges_gjam, directed = FALSE)
E(g_gjam)$width <- abs(edges_gjam$r) * 4
E(g_gjam)$color <- "#2166ac"
set.seed(11)
layout_gjam <- layout_with_fr(g_gjam)

png(file.path(FIGS_DIR, "gjam_network_full.png"), width = 1400, height = 1200, res = 140)
par(mar = c(1, 1, 3, 1))
plot(g_gjam, layout = layout_gjam, vertex.size = 16, vertex.color = "#66c2a5",
     vertex.label.cex = 0.8, vertex.label.color = "black",
     main = paste0("GJAM residual co-occurrence network (|r| > ", THRESH, ", habitat effects removed)"))
dev.off()
cat("Saved figs/gjam_network_full.png\n")

cat("\n=== Module 6 script complete: all figures + tables built from precomputed fits. ===\n")
cat("See module6_species_interactions.qmd for the narrated analysis and citations.\n")
