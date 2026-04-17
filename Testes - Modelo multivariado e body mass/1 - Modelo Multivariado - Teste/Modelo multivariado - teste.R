# -------------------------------------------------------------------------
# MODELOS GAM COMPLETOS - Com partição de variância e partial residuals
# Seguindo a metodologia de Coelho et al.
# -------------------------------------------------------------------------

library(mgcv)
library(ggplot2)
library(MetBrewer)
library(dplyr)
library(patchwork)

# =========================================================================
# 1. PREPARAÇÃO DOS DADOS (assumindo que out_wide_SA já existe)
# =========================================================================

# Identificar a coluna de abundância
coluna_abundancia <- names(out_wide_SA)[grep("Abundance|abundance|richness", 
                                             names(out_wide_SA))][1]

# Identificar preditoras (ajuste conforme seus nomes reais)
preditoras <- c()
possiveis_isolation <- c("Isolation_Mean", "Isolation Mean", "isolation_mean")
possiveis_fwd <- c("Frequency_Weighted_Distance", "Frequency Weighted Distance", "FWD")
possiveis_temp <- c("Temperature", "temp", "temperature")
possiveis_precip <- c("Precipitation_Log", "Precipitation Log", "precipitation_log")

preditoras[1] <- possiveis_isolation[possiveis_isolation %in% names(out_wide_SA)][1]
preditoras[2] <- possiveis_fwd[possiveis_fwd %in% names(out_wide_SA)][1]
preditoras[3] <- possiveis_temp[possiveis_temp %in% names(out_wide_SA)][1]
preditoras[4] <- possiveis_precip[possiveis_precip %in% names(out_wide_SA)][1]

# Remover NAs
out_clean <- out_wide_SA %>%
  select(es_SA, all_of(coluna_abundancia), all_of(preditoras)) %>%
  drop_na()

cat("=== CONFIGURAÇÃO ===\n")
cat("Variável resposta:", coluna_abundancia, "\n")
cat("Preditoras:", paste(preditoras, collapse = ", "), "\n")
cat("Total de observações:", nrow(out_clean), "\n")

# =========================================================================
# 2. FUNÇÃO PARA ANÁLISE COMPLETA DE UM SERVIÇO (com partição de variância)
# =========================================================================

analisar_servico_completo <- function(dados, nome_servico) {
  
  cat("\n========================================\n")
  cat("Analisando:", nome_servico, "\n")
  cat("n =", nrow(dados), "\n")
  
  if(nrow(dados) < 10) {
    cat("⚠ Dados insuficientes. Pulando...\n")
    return(NULL)
  }
  
  # -----------------------------------------------------------------------
  # Modelo completo (global)
  # -----------------------------------------------------------------------
  form_completo <- as.formula(paste(coluna_abundancia, "~", 
                                    paste(paste0("s(", preditoras, ", k = 4)"), 
                                          collapse = " + ")))
  
  b <- tryCatch({
    gam(form_completo, data = dados, family = poisson())
  }, error = function(e) {
    cat("Erro no modelo completo:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(b)) return(NULL)
  if(!b$converged) {
    cat("Modelo não convergiu\n")
    return(NULL)
  }
  
  cat("✓ Modelo completo convergiu\n")
  cat("  Deviance explicada:", round(1 - b$deviance/b$null.deviance, 3), "\n")
  
  # -----------------------------------------------------------------------
  # Modelo nulo
  # -----------------------------------------------------------------------
  b0 <- gam(as.formula(paste(coluna_abundancia, "~ 1")), 
            data = dados, family = poisson())
  
  # -----------------------------------------------------------------------
  # Modelos reduzidos (usando sp do modelo completo - igual ao exemplo)
  # -----------------------------------------------------------------------
  
  # Sem Isolation
  form_sem1 <- as.formula(paste(coluna_abundancia, "~", 
                                paste(paste0("s(", preditoras[-1], ", k = 4)"), 
                                      collapse = " + ")))
  b1 <- gam(form_sem1, data = dados, family = poisson(), sp = b$sp[2:4])
  
  # Sem FWD
  form_sem2 <- as.formula(paste(coluna_abundancia, "~", 
                                paste(paste0("s(", preditoras[c(1,3,4)], ", k = 4)"), 
                                      collapse = " + ")))
  b2 <- gam(form_sem2, data = dados, family = poisson(), sp = b$sp[c(1,3,4)])
  
  # Sem Temperature
  form_sem3 <- as.formula(paste(coluna_abundancia, "~", 
                                paste(paste0("s(", preditoras[c(1,2,4)], ", k = 4)"), 
                                      collapse = " + ")))
  b3 <- gam(form_sem3, data = dados, family = poisson(), sp = b$sp[c(1,2,4)])
  
  # Sem Precipitation
  form_sem4 <- as.formula(paste(coluna_abundancia, "~", 
                                paste(paste0("s(", preditoras[1:3], ", k = 4)"), 
                                      collapse = " + ")))
  b4 <- gam(form_sem4, data = dados, family = poisson(), sp = b$sp[1:3])
  
  # Modelos com apenas pares de variáveis (para efeitos conjuntos)
  form_pc <- as.formula(paste(coluna_abundancia, "~", 
                              paste(paste0("s(", preditoras[3:4], ", k = 4)"), 
                                    collapse = " + ")))
  b_pc <- gam(form_pc, data = dados, family = poisson(), sp = b$sp[3:4])
  
  form_geo <- as.formula(paste(coluna_abundancia, "~", 
                               paste(paste0("s(", preditoras[1:2], ", k = 4)"), 
                                     collapse = " + ")))
  b_geo <- gam(form_geo, data = dados, family = poisson(), sp = b$sp[1:2])
  
  # -----------------------------------------------------------------------
  # Cálculo da partição de variância (seguindo Coelho et al.)
  # -----------------------------------------------------------------------
  dev_nulo <- deviance(b0)
  dev_completo <- deviance(b)
  
  # Deviance total explicada
  total <- (dev_nulo - dev_completo) / dev_nulo
  
  # Contribuição única de cada variável
  isolation_alone <- (deviance(b1) - dev_completo) / dev_nulo
  fwd_alone <- (deviance(b2) - dev_completo) / dev_nulo
  temp_alone <- (deviance(b3) - dev_completo) / dev_nulo
  precip_alone <- (deviance(b4) - dev_completo) / dev_nulo
  
  # Efeitos conjuntos (como no exemplo)
  geo_clim_joint <- (deviance(b_pc) - dev_completo) / dev_nulo
  clim_self_joint <- (deviance(b_geo) - dev_completo) / dev_nulo
  
  # Efeito conjunto total
  joint_effect <- total - (isolation_alone + fwd_alone + temp_alone + precip_alone)
  
  # Resultados
  resultado <- data.frame(
    Servico = nome_servico,
    Total_Deviance = total,
    Isolation_Alone = isolation_alone,
    FWD_Alone = fwd_alone,
    Temperature_Alone = temp_alone,
    Precipitation_Alone = precip_alone,
    Joint_Geography_Climate = geo_clim_joint,
    Joint_Climate_Itself = clim_self_joint,
    Joint_Effect = joint_effect,
    N = nrow(dados),
    AIC = AIC(b)
  )
  
  # Mostrar resultados
  cat("\n--- PARTIÇÃO DE VARIÂNCIA ---\n")
  cat("Total deviance explicada:", round(total * 100, 1), "%\n")
  cat("Isolation sozinha:", round(isolation_alone * 100, 1), "%\n")
  cat("FWD sozinha:", round(fwd_alone * 100, 1), "%\n")
  cat("Temperature sozinha:", round(temp_alone * 100, 1), "%\n")
  cat("Precipitation sozinha:", round(precip_alone * 100, 1), "%\n")
  cat("Efeito conjunto:", round(joint_effect * 100, 1), "%\n")
  
  # Retornar modelo e resultados
  return(list(
    resultados = resultado,
    modelo = b,
    dados = dados,
    b1 = b1, b2 = b2, b3 = b3, b4 = b4  # guardar modelos reduzidos para plots
  ))
}

# =========================================================================
# 3. APLICAR PARA TODOS OS SERVIÇOS
# =========================================================================

servicos <- unique(out_clean$es_SA)
resultados_lista <- list()
modelos_lista <- list()

for(i in seq_along(servicos)) {
  dados_servico <- out_clean %>% filter(es_SA == servicos[i])
  analise <- analisar_servico_completo(dados_servico, servicos[i])
  
  if(!is.null(analise)) {
    resultados_lista[[i]] <- analise$resultados
    modelos_lista[[servicos[i]]] <- list(
      modelo = analise$modelo,
      b1 = analise$b1, b2 = analise$b2, b3 = analise$b3, b4 = analise$b4,
      dados = analise$dados
    )
  }
}

# Combinar resultados
resultados_finais <- bind_rows(resultados_lista)

# =========================================================================
# 4. TABELA DE RESULTADOS (como no exemplo)
# =========================================================================

cat("\n\n=========================================================================\n")
cat("RESULTADOS DA PARTIÇÃO DE VARIÂNCIA\n")
cat("=========================================================================\n")

# Formatar tabela
tabela_resultados <- resultados_finais %>%
  mutate(across(where(is.numeric), ~ round(. * 100, 1))) %>%
  arrange(desc(Total_Deviance))

print(tabela_resultados)

# Salvar
if(!dir.exists("Results")) dir.create("Results")
write.csv(tabela_resultados, "Results/Particao_Variancia_Completa.csv", row.names = FALSE)

# =========================================================================
# 5. GRÁFICO DE CONTRIBUIÇÃO DAS VARIÁVEIS (barplot empilhado)
# =========================================================================

# Preparar dados para o gráfico
dados_plot <- resultados_finais %>%
  select(Servico, Isolation_Alone, FWD_Alone, Temperature_Alone, 
         Precipitation_Alone, Joint_Effect) %>%
  pivot_longer(cols = -Servico, names_to = "Componente", values_to = "Proporcao")

# Renomear componentes para o gráfico
dados_plot$Componente <- factor(dados_plot$Componente,
                                levels = c("Isolation_Alone", "FWD_Alone", 
                                           "Temperature_Alone", "Precipitation_Alone",
                                           "Joint_Effect"),
                                labels = c("Isolation", "FWD", "Temperature", 
                                           "Precipitation", "Joint Effect"))

# Gráfico de barras empilhadas
fig_contribuicoes <- ggplot(dados_plot, aes(x = reorder(Servico, Proporcao, sum), 
                                            y = Proporcao * 100, fill = Componente)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = met.brewer("Egypt", 5)) +
  labs(x = "Ecosystem Service", 
       y = "Proportion of Deviance Explained (%)",
       title = "Variance Partitioning",
       subtitle = "Individual and joint contributions of environmental variables",
       fill = "Component") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  coord_flip()

print(fig_contribuicoes)

# Salvar
if(!dir.exists("Figures")) dir.create("Figures")
ggsave("Figures/Variance_Partitioning.jpg", fig_contribuicoes, 
       width = 12, height = 8, dpi = 300)

# =========================================================================
# 6. PARTIAL RESIDUALS PLOTS (para o melhor serviço)
# =========================================================================

# Escolher o melhor serviço (maior deviance explicada)
melhor_servico <- resultados_finais$Servico[which.max(resultados_finais$Total_Deviance)]
cat("\n\n=== GERANDO PARTIAL RESIDUALS PARA O MELHOR SERVIÇO ===\n")
cat("Serviço:", melhor_servico, "\n")
cat("Deviance explicada:", round(max(resultados_finais$Total_Deviance) * 100, 1), "%\n")

# Obter dados e modelos do melhor serviço
melhor_dados <- modelos_lista[[melhor_servico]]$dados
melhor_modelo <- modelos_lista[[melhor_servico]]$modelo
b1 <- modelos_lista[[melhor_servico]]$b1  # sem Isolation
b2 <- modelos_lista[[melhor_servico]]$b2  # sem FWD
b3 <- modelos_lista[[melhor_servico]]$b3  # sem Temperature
b4 <- modelos_lista[[melhor_servico]]$b4  # sem Precipitation

# Função para criar partial residuals plot (estilo Coelho et al.)
plot_partial_residuals <- function(dados, var_x, var_nome, residuals_model, cor) {
  
  plot_data <- data.frame(
    x = dados[[var_x]],
    residuals = residuals(residuals_model, type = "response")
  )
  
  ggplot(plot_data, aes(x = x, y = residuals)) +
    geom_point(size = 3, col = cor, alpha = 0.6) +
    labs(x = var_nome, y = "Partial Residuals") +
    stat_smooth(method = gam, formula = y ~ s(x, k = 4), 
                col = met.brewer("Hiroshige")[7], fill = "gray70", alpha = 0.3) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(color = "black", size = 10),
          plot.title = element_text(hjust = 0.5))
}

# Cores
cores <- met.brewer("Egypt", 4)

# Criar os 4 partial residuals plots
p1 <- plot_partial_residuals(melhor_dados, preditoras[1], "Isolation Mean (km)", b1, cores[1])
p2 <- plot_partial_residuals(melhor_dados, preditoras[2], "Freq. Weighted Distance", b2, cores[2])
p3 <- plot_partial_residuals(melhor_dados, preditoras[3], "Temperature (°C)", b3, cores[3])
p4 <- plot_partial_residuals(melhor_dados, preditoras[4], "Precipitation (sqrt)", b4, cores[4])

# Combinar em uma figura 2x2
fig_partial <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = paste("Partial Residuals -", melhor_servico),
    subtitle = paste("Deviance Explained:", round(max(resultados_finais$Total_Deviance) * 100, 1), "%"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                  plot.subtitle = element_text(hjust = 0.5, size = 12))
  )

print(fig_partial)

# Salvar
ggsave(paste0("Figures/Partial_Residuals_", gsub(" ", "_", melhor_servico), ".jpg"), 
       fig_partial, width = 12, height = 10, dpi = 300)

# =========================================================================
# 7. GRÁFICO DE BARRAS DAS CONTRIBUIÇÕES INDIVIDUAIS
# =========================================================================

# Preparar dados para contribuições individuais
individuais <- resultados_finais %>%
  select(Servico, Isolation_Alone, FWD_Alone, Temperature_Alone, Precipitation_Alone) %>%
  pivot_longer(cols = -Servico, names_to = "Variavel", values_to = "Contribuicao")

individuais$Variavel <- factor(individuais$Variavel,
                               levels = c("Isolation_Alone", "FWD_Alone", 
                                          "Temperature_Alone", "Precipitation_Alone"),
                               labels = c("Isolation", "FWD", "Temperature", "Precipitation"))

fig_individuais <- ggplot(individuais, aes(x = reorder(Servico, Contribuicao, sum), 
                                           y = Contribuicao * 100, fill = Variavel)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = met.brewer("Egypt", 4)) +
  labs(x = "Ecosystem Service", 
       y = "Unique Contribution (%)",
       title = "Individual Contributions of Each Variable",
       fill = "Variable") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_flip()

print(fig_individuais)
ggsave("Figures/Individual_Contributions.jpg", fig_individuais, 
       width = 12, height = 8, dpi = 300)

# =========================================================================
# 8. RESUMO FINAL
# =========================================================================

cat("\n\n=========================================================================\n")
cat("RESUMO FINAL DA ANÁLISE\n")
cat("=========================================================================\n")

cat("\n📊 ESTATÍSTICAS GERAIS:\n")
cat("  - Total de serviços analisados:", nrow(resultados_finais), "\n")
cat("  - Média da deviance explicada:", round(mean(resultados_finais$Total_Deviance) * 100, 1), "%\n")
cat("  - Mediana:", round(median(resultados_finais$Total_Deviance) * 100, 1), "%\n")
cat("  - Mínimo:", round(min(resultados_finais$Total_Deviance) * 100, 1), "%\n")
cat("  - Máximo:", round(max(resultados_finais$Total_Deviance) * 100, 1), "%\n")

cat("\n🏆 TOP 3 SERVIÇOS:\n")
top3 <- resultados_finais %>% 
  arrange(desc(Total_Deviance)) %>% 
  head(3)
for(i in 1:nrow(top3)) {
  cat("  ", i, ".", top3$Servico[i], "-", round(top3$Total_Deviance[i] * 100, 1), "%\n")
}

cat("\n📁 ARQUIVOS GERADOS:\n")
cat("  - Results/Particao_Variancia_Completa.csv\n")
cat("  - Figures/Variance_Partitioning.jpg\n")
cat("  - Figures/Individual_Contributions.jpg\n")
cat("  - Figures/Partial_Residuals_", gsub(" ", "_", melhor_servico), ".jpg\n", sep = "")

cat("\n✅ Análise concluída com sucesso!\n")