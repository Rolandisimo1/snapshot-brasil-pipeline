params <-
list(run_live = FALSE)

library(data.table)
DATA <- "mod6_data"
FIGS <- "mod6_figs"
read_d <- function(f) fread(file.path(DATA, f))

knitr::include_graphics(file.path(FIGS, "occ_vs_abundance.png"))

siv_cam <- read_d("species_interaction_values.csv")
ppc_ok_cam <- siv_cam[bayes_p > .05 & bayes_p < .95 & chat < 1.15, .N]
knitr::kable(
  data.frame(
    Quantity = c("Pares ajustados (a-priori)", "Pares ajustados (sinalizados pelo GJAM)",
                 "Total de ajustes (incl. direção reversa)", "Reajustados com ODRE",
                 "Aprovados na verificação preditiva a posteriori"),
    Value = c(12, 8, nrow(siv_cam), siv_cam[used_odre == TRUE, .N], ppc_ok_cam)),
  caption = "Resumo do ajuste em nível de câmera. A maioria dos pares exigiu o efeito aleatório de sobredispersão."
)

knitr::include_graphics(file.path(FIGS, "siv_forest_camera.png"))

knitr::include_graphics(file.path(FIGS, "interaction_network_camera.png"))

knitr::include_graphics(file.path(FIGS, "abundance_response_curves_camera.png"))

knitr::include_graphics(file.path(FIGS, "topdown_bottomup_comparison.png"))

tally <- jsonlite::fromJSON(file.path(DATA, "topdown_bottomup_tally.json"))
al <- jsonlite::fromJSON(file.path(DATA, "amir_luskin_comparison.json"))
knitr::kable(
  data.frame(
    Quantity = c("Testes cross-tróficos", "Top-down suportado", "Bottom-up suportado", "Não suportado"),
    `Neotropical (nível de câmera)` = c(tally$our_cross_n,
                    paste0(tally$our_topdown, " (", tally$our_topdown_pct, "%)"),
                    paste0(tally$our_bottomup, " (", tally$our_bottomup_pct, "%)"),
                    paste0(tally$our_unsupported, " (", tally$our_unsup_pct, "%)")),
    `Sudeste Asiático (Amir & Luskin)` = c(al$preferred_pairs,
                    paste0(al$topdown_supported_n, " (", al$topdown_supported_pct, "%)"),
                    paste0(al$bottomup_supported_n, " (", al$bottomup_supported_pct, "%)"),
                    paste0(al$unsupported_n, " (", al$unsupported_pct, "%)")),
    check.names = FALSE),
  caption = "Balanço top-down/bottom-up em nível de câmera vs. o estudo comparativo do sudeste asiático."
)

conv <- read_d("convergence_diagnostics_summary.csv")
cat("Array-level pairs converged:", sum(conv$converged), "of", nrow(conv), "\n")

knitr::include_graphics(file.path(FIGS, "gamma0_forest_all_pairs.png"))

sig <- read_d("gamma0_significant_pairs.csv")
knitr::kable(sig[, .(Driver=driver_c, Responder=responder_c, Source=source,
                       Mean=round(mean,2), Lower95=round(lcl,2), Upper95=round(ucl,2))],
             caption = "Os 7 pares em nível de array com intervalo de credibilidade de 95% que exclui zero.")

knitr::include_graphics(file.path(FIGS, "abundance_response_curves_array.png"))

knitr::include_graphics(file.path(FIGS, "interaction_network_abundance.png"))

knitr::include_graphics(file.path(FIGS, "camera_vs_array_comparison.png"))

knitr::include_graphics(file.path(FIGS, "camera_vs_array_slopes.png"))

cva <- read_d("camera_vs_array_comparison.csv")
n_pairs <- nrow(cva)
cam_sig <- sum(cva$g0_sig)
arr_sig <- sum(cva$g0_sig_arr)
sign_agree <- sum(cva$sign_agree)
both_sig <- sum(cva$both_sig)
cam_only <- sum(cva$g0_sig & !cva$g0_sig_arr)
arr_only <- sum(!cva$g0_sig & cva$g0_sig_arr)
knitr::kable(
  data.frame(
    Quantity = c("Pares comparados", "Significativo em nível de câmera", "Significativo em nível de array",
                 "Significativo em ambos", "Significativo apenas em câmera", "Significativo apenas em array",
                 "Concordância de sinal (ambas resoluções)"),
    Value = c(n_pairs,
              paste0(cam_sig, "/", n_pairs, " (", round(100*cam_sig/n_pairs), "%)"),
              paste0(arr_sig, "/", n_pairs, " (", round(100*arr_sig/n_pairs), "%)"),
              both_sig, cam_only, arr_only,
              paste0(sign_agree, "/", n_pairs, " (", round(100*sign_agree/n_pairs), "%)"))),
  caption = "Concordância entre nível de câmera e nível de array nos 20 pares compartilhados."
)

knitr::include_graphics(file.path(FIGS, "gjam_network_full.png"))

gcorr <- read_d("gjam_residual_correlation.csv")
cat("Species pairs with |r| > 0.4:", 12, "of", choose(44,2), "possible pairs\n")
