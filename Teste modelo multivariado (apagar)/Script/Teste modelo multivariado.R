# -------------------------------------------------------------------------
# MODELOS MULTIVARIADOS GAM (múltiplas variáveis preditoras)
# -------------------------------------------------------------------------

# Primeiro, vamos criar o dataframe wide (formato original, não pivotado)
# Isso é necessário porque seus modelos multivariados usam várias colunas

data_wide_list_SA <- list()

for (i in seq_along(es_names_SA)) {
  res_map_SA <- lets.maplizer.env(res_SA, y = ES_SA[[es_names_SA[i]]], z = ES_SA$species, func = sum)
  out_SA <- lets.envcells(res_SA)
  
  # Retrieving Temperature (column 2), Precipitation (column 3), and Abundance (column 4)
  temp_climate_SA <- tibble(
    es_SA = es_names_SA[i],
    Temperature = res_map_SA$Matrix_env[, 2],
    `Precipitation (Log)` = res_map_SA$Matrix_env[, 3],
    Abundance_SA = res_map_SA$Matrix_env[, 4]
  )
  
  # Retrieving the geographic variables from the out object
  temp_geo_SA <- as_tibble(out_SA[keep_SA, -1])
  
  # Putting everything together and removing NAs
  temp_wide_SA <- bind_cols(temp_climate_SA, temp_geo_SA) %>% drop_na()
  
  # Saving in the list
  data_wide_list_SA[[i]] <- temp_wide_SA
}

# Combining all services into a single data frame
out_wide_SA <- bind_rows(data_wide_list_SA)

# Renomeando colunas para nomes sem espaços (evita problemas)
out_wide_SA <- out_wide_SA %>%
  rename(
    Isolation_Mean = `Isolation (Mean)`,
    Frequency_Weighted_Distance = `Frequency Weighted Distance`,
    Precipitation_Log = `Precipitation (Log)`
  )

# Verificar os nomes das colunas
names(out_wide_SA)

# -------------------------------------------------------------------------
# RODANDO OS MODELOS MULTIVARIADOS PARA CADA SERVIÇO
# -------------------------------------------------------------------------

results_gam_SA <- out_wide_SA %>%
  group_by(es_SA) %>%
  nest() %>%
  mutate(
    # Modelo 1: Nulo (só intercepto)
    model_null = map(data, ~ mgcv::gam(Abundance_SA ~ 1, 
                                       data = .x, 
                                       family = poisson(link = "log"))),
    
    # Modelo 2: Isolation + Temperature + Precipitation
    model_2 = map(data, ~ mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                      s(Temperature) + 
                                      s(Precipitation_Log), 
                                    data = .x, 
                                    family = poisson(link = "log"))),
    
    # Modelo 3: Frequency Weighted Distance + Temperature + Precipitation
    model_3 = map(data, ~ mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                                      s(Temperature) + 
                                      s(Precipitation_Log), 
                                    data = .x, 
                                    family = poisson(link = "log"))),
    
    # Modelo 4: Isolation + Frequency Weighted Distance + Temperature
    model_4 = map(data, ~ mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                      s(Frequency_Weighted_Distance) + 
                                      s(Temperature), 
                                    data = .x, 
                                    family = poisson(link = "log"))),
    
    # Modelo 5: Isolation + Frequency Weighted Distance + Precipitation
    model_5 = map(data, ~ mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                      s(Frequency_Weighted_Distance) + 
                                      s(Precipitation_Log), 
                                    data = .x, 
                                    family = poisson(link = "log"))),
    
    # Modelo 6: Temperature + Precipitation (apenas)
    model_6 = map(data, ~ mgcv::gam(Abundance_SA ~ s(Temperature) + 
                                      s(Precipitation_Log), 
                                    data = .x, 
                                    family = poisson(link = "log"))),
    
    # Modelo 7: Isolation + Frequency Weighted Distance (apenas)
    model_7 = map(data, ~ mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                      s(Frequency_Weighted_Distance), 
                                    data = .x, 
                                    family = poisson(link = "log"))),
    
    # Modelo Original (Frequency apenas)
    model_original = map(data, ~ mgcv::gam(Abundance_SA ~ s(Frequency), 
                                           data = .x, 
                                           family = poisson(link = "log")))
  )

# -------------------------------------------------------------------------
# EXTRAINDO MÉTRICAS DE COMPARAÇÃO (Deviance Explicada e AIC)
# -------------------------------------------------------------------------

results_metrics <- results_gam_SA %>%
  mutate(
    # Deviance explicada (pseudo-R2) para cada modelo
    Deviance_Null = map_dbl(model_null, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_2 = map_dbl(model_2, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_3 = map_dbl(model_3, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_4 = map_dbl(model_4, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_5 = map_dbl(model_5, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_6 = map_dbl(model_6, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_7 = map_dbl(model_7, ~ 1 - (.x$deviance / .x$null.deviance)),
    Deviance_Original = map_dbl(model_original, ~ 1 - (.x$deviance / .x$null.deviance)),
    
    # AIC para cada modelo
    AIC_Null = map_dbl(model_null, AIC),
    AIC_2 = map_dbl(model_2, AIC),
    AIC_3 = map_dbl(model_3, AIC),
    AIC_4 = map_dbl(model_4, AIC),
    AIC_5 = map_dbl(model_5, AIC),
    AIC_6 = map_dbl(model_6, AIC),
    AIC_7 = map_dbl(model_7, AIC),
    AIC_Original = map_dbl(model_original, AIC)
  ) %>%
  select(es_SA, starts_with("Deviance_"), starts_with("AIC_"))

# Visualizar os resultados
print(results_metrics)

# -------------------------------------------------------------------------
# QUAL MODELO É MELHOR PARA CADA SERVIÇO (menor AIC)
# -------------------------------------------------------------------------

# Transformar para formato longo para facilitar a comparação
results_long <- results_metrics %>%
  pivot_longer(
    cols = -es_SA,
    names_to = c(".value", "Modelo"),
    names_pattern = "(Deviance|AIC)_(.*)"
  ) %>%
  mutate(Modelo = factor(Modelo, 
                         levels = c("Null", "2", "3", "4", "5", "6", "7", "Original"),
                         labels = c("Nulo", "Modelo_2", "Modelo_3", "Modelo_4", 
                                    "Modelo_5", "Modelo_6", "Modelo_7", "Original")))

# Para cada serviço, qual modelo tem menor AIC?
melhor_modelo <- results_long %>%
  group_by(es_SA) %>%
  filter(AIC == min(AIC, na.rm = TRUE)) %>%
  select(es_SA, Melhor_Modelo = Modelo, AIC, Deviance)

print(melhor_modelo)

# -------------------------------------------------------------------------
# EXTRAINDO P-VALUES E ESTATÍSTICAS DO MODELO 2 (exemplo)
# -------------------------------------------------------------------------

# Para extrair estatísticas detalhadas de um modelo específico
resultados_detalhados <- results_gam_SA %>%
  mutate(
    summary_model_2 = map(model_2, ~ {
      s <- summary(.x)
      tibble(
        R_sq = s$r.sq,
        Deviance_Explained = s$dev.expl,
        n = length(.x$residuals)
      )
    }),
    
    # Extrair coeficientes/p-valores de cada variável no modelo 2
    coeficientes_2 = map(model_2, ~ {
      coef_summary <- summary(.x)$s.table
      tibble(
        Variavel = rownames(coef_summary),
        edf = coef_summary[, 1],
        Ref_df = coef_summary[, 2],
        F_statistic = coef_summary[, 3],
        p_value = coef_summary[, 4]
      )
    })
  ) %>%
  unnest(summary_model_2, coeficientes_2)

# Ver resultados detalhados do modelo 2
print(resultados_detalhados %>% select(es_SA, Variavel, edf, F_statistic, p_value))

# -------------------------------------------------------------------------
# SALVANDO OS RESULTADOS
# -------------------------------------------------------------------------

# Salvar tabela de comparação de modelos
write_csv(results_metrics, "Results/Comparacao_Modelos_Multivariados.csv")

# Salvar melhor modelo por serviço
write_csv(melhor_modelo, "Results/Melhor_Modelo_por_Servico.csv")

# Salvar estatísticas detalhadas do modelo 2
write_csv(resultados_detalhados, "Results/Estatisticas_Modelo_2.csv")

print("Análise multivariada concluída com sucesso!")

# -------------------------------------------------------------------------
# PARTIÇÃO DE VARIÂNCIA (Decomposition of Deviance)
# Com 4 variáveis: Isolation_Mean, Frequency_Weighted_Distance, 
#                  Temperature, Precipitation_Log
# -------------------------------------------------------------------------

library(MetBrewer)
library(ggplot2)
library(patchwork)

# -------------------------------------------------------------------------
# VERSÃO PARA TODOS OS SERVIÇOS (com 4 variáveis)
# -------------------------------------------------------------------------

# Função para calcular partição de variância para um dataset
compute_variance_partitioning <- function(data) {
  tryCatch({
    # Modelo nulo
    b0 <- mgcv::gam(Abundance_SA ~ 1, data = data, family = poisson(link = "log"))
    
    # Modelo completo (todas as 4 variáveis)
    b_completo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                              s(Frequency_Weighted_Distance) + 
                              s(Temperature) + 
                              s(Precipitation_Log), 
                            data = data, 
                            family = poisson(link = "log"))
    
    # Modelos sem cada variável (para calcular contribuição única)
    b_sem_isolation <- mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                                   s(Temperature) + 
                                   s(Precipitation_Log), 
                                 data = data, 
                                 family = poisson(link = "log"))
    
    b_sem_fwd <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                             s(Temperature) + 
                             s(Precipitation_Log), 
                           data = data, 
                           family = poisson(link = "log"))
    
    b_sem_temp <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                              s(Frequency_Weighted_Distance) + 
                              s(Precipitation_Log), 
                            data = data, 
                            family = poisson(link = "log"))
    
    b_sem_precip <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                s(Frequency_Weighted_Distance) + 
                                s(Temperature), 
                              data = data, 
                              family = poisson(link = "log"))
    
    dev_nulo <- deviance(b0)
    dev_completo <- deviance(b_completo)
    dev_sem_isolation <- deviance(b_sem_isolation)
    dev_sem_fwd <- deviance(b_sem_fwd)
    dev_sem_temp <- deviance(b_sem_temp)
    dev_sem_precip <- deviance(b_sem_precip)
    
    # Deviance total explicada pelo modelo completo
    total_dev <- (dev_nulo - dev_completo) / dev_nulo
    
    # Deviance explicada quando cada variável é removida
    isolation_without_others <- (dev_nulo - dev_sem_isolation) / dev_nulo
    fwd_without_others <- (dev_nulo - dev_sem_fwd) / dev_nulo
    temp_without_others <- (dev_nulo - dev_sem_temp) / dev_nulo
    precip_without_others <- (dev_nulo - dev_sem_precip) / dev_nulo
    
    # Contribuição única de cada variável
    isolation_alone <- total_dev - isolation_without_others
    fwd_alone <- total_dev - fwd_without_others
    temp_alone <- total_dev - temp_without_others
    precip_alone <- total_dev - precip_without_others
    
    # Efeito conjunto (interação entre todas as variáveis)
    joint_effect <- total_dev - (isolation_alone + fwd_alone + temp_alone + precip_alone)
    
    return(tibble(
      Total_Deviance = total_dev,
      Isolation_Alone = isolation_alone,
      FWD_Alone = fwd_alone,
      Temperature_Alone = temp_alone,
      Precipitation_Alone = precip_alone,
      Joint_Effect = joint_effect,
      Unexplained = 1 - total_dev
    ))
  }, error = function(e) {
    return(tibble(
      Total_Deviance = NA,
      Isolation_Alone = NA,
      FWD_Alone = NA,
      Temperature_Alone = NA,
      Precipitation_Alone = NA,
      Joint_Effect = NA,
      Unexplained = NA
    ))
  })
}

# Aplicar para todos os serviços
variance_partitioning_results <- out_wide_SA %>%
  group_by(es_SA) %>%
  nest() %>%
  mutate(partitioning = map(data, compute_variance_partitioning)) %>%
  unnest(partitioning) %>%
  select(-data)

# Visualizar resultados
print(variance_partitioning_results)

# Salvar resultados
write_csv(variance_partitioning_results, "Results/Variance_Partitioning_4vars.csv")

# -------------------------------------------------------------------------
# GRÁFICO DE PARTIÇÃO DE VARIÂNCIA (stacked bar plot)
# -------------------------------------------------------------------------

# Transformar para formato longo
variance_long <- variance_partitioning_results %>%
  pivot_longer(cols = c(Isolation_Alone, FWD_Alone, Temperature_Alone, 
                        Precipitation_Alone, Joint_Effect, Unexplained),
               names_to = "Component",
               values_to = "Proportion") %>%
  mutate(Component = factor(Component,
                            levels = c("Isolation_Alone", "FWD_Alone", 
                                       "Temperature_Alone", "Precipitation_Alone",
                                       "Joint_Effect", "Unexplained"),
                            labels = c("Isolation Alone", 
                                       "Freq. Weighted Distance Alone", 
                                       "Temperature Alone",
                                       "Precipitation Alone",
                                       "Joint Effect", 
                                       "Unexplained")))

# Gráfico de barras empilhadas
ggplot(variance_long, aes(x = es_SA, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = met.brewer("Egypt", 6)) +
  labs(x = "Ecosystem Service", 
       y = "Proportion of Deviance Explained",
       title = "Variance Partitioning: 4 Environmental Variables",
       fill = "Component") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  coord_flip()

ggsave("Figures/Variance_Partitioning_4vars.jpg", width = 12, height = 8)

# -------------------------------------------------------------------------
# PARTIAL RESIDUALS PLOTS (para todas as 4 variáveis)
# -------------------------------------------------------------------------

# Função para criar partial residuals plot
plot_partial_residuals <- function(data, var_name, var_label, color_index) {
  
  if(var_name == "Isolation_Mean") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                          s(Temperature) + 
                          s(Precipitation_Log), 
                        data = data, 
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data$Isolation_Mean
    
  } else if(var_name == "Frequency_Weighted_Distance") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Temperature) + 
                          s(Precipitation_Log), 
                        data = data, 
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data$Frequency_Weighted_Distance
    
  } else if(var_name == "Temperature") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Precipitation_Log), 
                        data = data, 
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data$Temperature
    
  } else if(var_name == "Precipitation_Log") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Temperature), 
                        data = data, 
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data$Precipitation_Log
    
  } else {
    stop("Variável não reconhecida")
  }
  
  # Criar dataframe para plot
  plot_data <- tibble(x = x_var, y = residuals)
  
  # Gráfico
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(size = 2.5, col = met.brewer("Egypt")[color_index], alpha = 0.5) +
    labs(x = var_label, y = "Partial Residuals") +
    stat_smooth(method = gam, 
                formula = y ~ s(x, bs = "cs", k = 4), 
                col = met.brewer("Hiroshige")[7], 
                fill = "gray70",
                alpha = 0.3) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          axis.text = element_text(color = "black", size = 9),
          axis.line = element_line(colour = "black"))
}

# Gerar partial residuals plots para todos os serviços
servicos_unicos <- unique(out_wide_SA$es_SA)

# Criar PDF com todos os partial residuals plots
pdf("Figures/Partial_Residuals_All_Services_4vars.pdf", width = 12, height = 16)

for(servico in servicos_unicos) {
  
  cat("Plotando para:", servico, "\n")
  
  # Filtrar dados do serviço
  data_servico <- out_wide_SA %>% filter(es_SA == servico)
  
  # Criar plots para as 4 variáveis
  p1 <- plot_partial_residuals(data_servico, "Isolation_Mean", "Isolation Mean (km)", 1)
  p2 <- plot_partial_residuals(data_servico, "Frequency_Weighted_Distance", "Frequency Weighted Distance", 2)
  p3 <- plot_partial_residuals(data_servico, "Temperature", "Temperature (°C)", 3)
  p4 <- plot_partial_residuals(data_servico, "Precipitation_Log", "Precipitation (sqrt)", 4)
  
  # Combinar os 4 plots (2x2)
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = servico,
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  print(combined)
}

dev.off()

print("Partial residuals plots salvos em Figures/Partial_Residuals_All_Services_4vars.pdf")

# -------------------------------------------------------------------------
# GRÁFICO DE CONTRIBUIÇÃO INDIVIDUAL (barras lado a lado)
# -------------------------------------------------------------------------

# Para visualizar melhor a contribuição única de cada variável
individual_contributions <- variance_partitioning_results %>%
  select(es_SA, Isolation_Alone, FWD_Alone, Temperature_Alone, Precipitation_Alone) %>%
  pivot_longer(cols = c(Isolation_Alone, FWD_Alone, Temperature_Alone, Precipitation_Alone),
               names_to = "Variable",
               values_to = "Proportion") %>%
  mutate(Variable = factor(Variable,
                           levels = c("Isolation_Alone", "FWD_Alone", 
                                      "Temperature_Alone", "Precipitation_Alone"),
                           labels = c("Isolation", "FWD", "Temperature", "Precipitation")))

ggplot(individual_contributions, aes(x = es_SA, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = met.brewer("Egypt", 4)) +
  labs(x = "Ecosystem Service", 
       y = "Unique Contribution (Proportion of Deviance)",
       title = "Individual Contributions of Each Variable",
       fill = "Variable") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

ggsave("Figures/Individual_Contributions_4vars.jpg", width = 14, height = 7)

# -------------------------------------------------------------------------
# TABELA RESUMO DOS RESULTADOS (ordenada por total deviance)
# -------------------------------------------------------------------------

# Ordenar por total deviance explicada
resultados_ordenados <- variance_partitioning_results %>%
  arrange(desc(Total_Deviance)) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

print(resultados_ordenados)

# Salvar versão formatada
write_csv(resultados_ordenados, "Results/Variance_Partitioning_Ranked_4vars.csv")

# -------------------------------------------------------------------------
# SUMÁRIO ESTATÍSTICO DOS MODELOS
# -------------------------------------------------------------------------

# Função para extrair estatísticas dos modelos completos
extract_model_stats <- function(data) {
  tryCatch({
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Temperature) + 
                          s(Precipitation_Log), 
                        data = data, 
                        family = poisson(link = "log"))
    
    summary_modelo <- summary(modelo)
    
    # Extrair informações das variáveis suavizadas
    smooth_table <- as.data.frame(summary_modelo$s.table)
    
    tibble(
      Isolation_edf = smooth_table[1, 1],
      Isolation_F = smooth_table[1, 3],
      Isolation_p = smooth_table[1, 4],
      FWD_edf = smooth_table[2, 1],
      FWD_F = smooth_table[2, 3],
      FWD_p = smooth_table[2, 4],
      Temperature_edf = smooth_table[3, 1],
      Temperature_F = smooth_table[3, 3],
      Temperature_p = smooth_table[3, 4],
      Precipitation_edf = smooth_table[4, 1],
      Precipitation_F = smooth_table[4, 3],
      Precipitation_p = smooth_table[4, 4],
      R_sq = summary_modelo$r.sq,
      Deviance_Explained = summary_modelo$dev.expl,
      GCV = summary_modelo$sp.criterion
    )
  }, error = function(e) {
    return(tibble(
      Isolation_edf = NA, Isolation_F = NA, Isolation_p = NA,
      FWD_edf = NA, FWD_F = NA, FWD_p = NA,
      Temperature_edf = NA, Temperature_F = NA, Temperature_p = NA,
      Precipitation_edf = NA, Precipitation_F = NA, Precipitation_p = NA,
      R_sq = NA, Deviance_Explained = NA, GCV = NA
    ))
  })
}

# Aplicar para todos os serviços
model_stats <- out_wide_SA %>%
  group_by(es_SA) %>%
  nest() %>%
  mutate(stats = map(data, extract_model_stats)) %>%
  unnest(stats) %>%
  select(-data)

# Visualizar estatísticas
print(model_stats)

# Salvar
write_csv(model_stats, "Results/Model_Statistics_4vars.csv")

print("Análise completa concluída com sucesso!")

# -------------------------------------------------------------------------
# PARTIAL RESIDUALS PLOTS (Adaptado de Coelho et al.)
# Para Isolation_Mean, Frequency_Weighted_Distance, Temperature, Precipitation_Log
# -------------------------------------------------------------------------

library(MetBrewer)
library(ggplot2)
library(patchwork)

# -------------------------------------------------------------------------
# VERSÃO 1: Para um serviço específico (exemplo)
# -------------------------------------------------------------------------

# Escolha um serviço para visualizar (ex: o primeiro da lista)
servico_exemplo <- unique(out_wide_SA$es_SA)[1]
cat("Analisando serviço:", servico_exemplo, "\n")

# Filtrar dados para o serviço escolhido
data_servico <- out_wide_SA %>% filter(es_SA == servico_exemplo)

# Ajustar o modelo completo
modelo_completo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                               s(Frequency_Weighted_Distance) + 
                               s(Temperature) + 
                               s(Precipitation_Log), 
                             data = data_servico, 
                             family = poisson(link = "log"))

# -------------------------------------------------------------------------
# Método 1: Partial residuals calculados manualmente (como no código de Coelho et al.)
# -------------------------------------------------------------------------

# Para Isolation_Mean: residuais do modelo SEM Isolation, plotado contra Isolation
modelo_sem_isolation <- mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                                    s(Temperature) + 
                                    s(Precipitation_Log), 
                                  data = data_servico, 
                                  family = poisson(link = "log"))

# Partial residuals = residuais do modelo sem a variável
partial_residuals_isolation <- residuals(modelo_sem_isolation, type = "response")

# DataFrame para plot
partial_data_isolation <- data.frame(
  x = data_servico$Isolation_Mean,
  y = partial_residuals_isolation
)

# Gráfico partial residuals para Isolation_Mean
PartialResidual_Isolation <- ggplot(partial_data_isolation, aes(x = x, y = y)) +
  geom_point(size = 4, col = met.brewer("Egypt")[1], alpha = 0.7) +
  labs(x = "Isolation Mean (km)", y = "Partial Residuals") +
  stat_smooth(method = gam, 
              formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], 
              fill = "gray70",
              alpha = 0.3) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(colour = "black"))

# Para Frequency_Weighted_Distance
modelo_sem_fwd <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                              s(Temperature) + 
                              s(Precipitation_Log), 
                            data = data_servico, 
                            family = poisson(link = "log"))

partial_residuals_fwd <- residuals(modelo_sem_fwd, type = "response")

partial_data_fwd <- data.frame(
  x = data_servico$Frequency_Weighted_Distance,
  y = partial_residuals_fwd
)

PartialResidual_FWD <- ggplot(partial_data_fwd, aes(x = x, y = y)) +
  geom_point(size = 4, col = met.brewer("Egypt")[2], alpha = 0.7) +
  labs(x = "Frequency Weighted Distance", y = "Partial Residuals") +
  stat_smooth(method = gam, 
              formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], 
              fill = "gray70",
              alpha = 0.3) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(colour = "black"))

# Para Temperature
modelo_sem_temp <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                               s(Frequency_Weighted_Distance) + 
                               s(Precipitation_Log), 
                             data = data_servico, 
                             family = poisson(link = "log"))

partial_residuals_temp <- residuals(modelo_sem_temp, type = "response")

partial_data_temp <- data.frame(
  x = data_servico$Temperature,
  y = partial_residuals_temp
)

PartialResidual_Temp <- ggplot(partial_data_temp, aes(x = x, y = y)) +
  geom_point(size = 4, col = met.brewer("Egypt")[3], alpha = 0.7) +
  labs(x = "Temperature (°C)", y = "Partial Residuals") +
  stat_smooth(method = gam, 
              formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], 
              fill = "gray70",
              alpha = 0.3) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(colour = "black"))

# Para Precipitation_Log
modelo_sem_precip <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                 s(Frequency_Weighted_Distance) + 
                                 s(Temperature), 
                               data = data_servico, 
                               family = poisson(link = "log"))

partial_residuals_precip <- residuals(modelo_sem_precip, type = "response")

partial_data_precip <- data.frame(
  x = data_servico$Precipitation_Log,
  y = partial_residuals_precip
)

PartialResidual_Precip <- ggplot(partial_data_precip, aes(x = x, y = y)) +
  geom_point(size = 4, col = met.brewer("Egypt")[4], alpha = 0.7) +
  labs(x = "Precipitation (sqrt)", y = "Partial Residuals") +
  stat_smooth(method = gam, 
              formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], 
              fill = "gray70",
              alpha = 0.3) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(colour = "black"))

# Combinar os 4 gráficos
combined_partial_residuals <- (PartialResidual_Isolation + PartialResidual_FWD) /
  (PartialResidual_Temp + PartialResidual_Precip) +
  plot_annotation(title = paste("Partial Residuals -", servico_exemplo),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))

# Exibir
print(combined_partial_residuals)

# Salvar
ggsave(paste0("Figures/Partial_Residuals_", servico_exemplo, ".jpg"), 
       combined_partial_residuals, width = 14, height = 10)

# -------------------------------------------------------------------------
# VERSÃO 2: Para TODOS os serviços (loop)
# -------------------------------------------------------------------------

# Criar PDF com todos os partial residuals plots
pdf("Figures/Partial_Residuals_All_Services_Complete.pdf", width = 14, height = 10)

for(servico in unique(out_wide_SA$es_SA)) {
  
  cat("Plotando partial residuals para:", servico, "\n")
  
  # Filtrar dados
  data_servico <- out_wide_SA %>% filter(es_SA == servico)
  
  # Ajustar modelos
  modelo_sem_isolation <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                s(Temperature) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  modelo_sem_fwd <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Temperature) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  modelo_sem_temp <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  modelo_sem_precip <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Temperature), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Gráfico para Isolation
  if(!is.null(modelo_sem_isolation)) {
    p1 <- ggplot(data.frame(x = data_servico$Isolation_Mean, 
                            y = residuals(modelo_sem_isolation, type = "response")), 
                 aes(x = x, y = y)) +
      geom_point(size = 3, col = met.brewer("Egypt")[1], alpha = 0.6) +
      labs(x = "Isolation Mean (km)", y = "Partial Residuals") +
      stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
                  col = met.brewer("Hiroshige")[7], fill = "gray70", alpha = 0.3) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = NA),
            axis.line = element_line(colour = "black"))
  } else { p1 <- plot_spacer() }
  
  # Gráfico para FWD
  if(!is.null(modelo_sem_fwd)) {
    p2 <- ggplot(data.frame(x = data_servico$Frequency_Weighted_Distance, 
                            y = residuals(modelo_sem_fwd, type = "response")), 
                 aes(x = x, y = y)) +
      geom_point(size = 3, col = met.brewer("Egypt")[2], alpha = 0.6) +
      labs(x = "Frequency Weighted Distance", y = "Partial Residuals") +
      stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
                  col = met.brewer("Hiroshige")[7], fill = "gray70", alpha = 0.3) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = NA),
            axis.line = element_line(colour = "black"))
  } else { p2 <- plot_spacer() }
  
  # Gráfico para Temperature
  if(!is.null(modelo_sem_temp)) {
    p3 <- ggplot(data.frame(x = data_servico$Temperature, 
                            y = residuals(modelo_sem_temp, type = "response")), 
                 aes(x = x, y = y)) +
      geom_point(size = 3, col = met.brewer("Egypt")[3], alpha = 0.6) +
      labs(x = "Temperature (°C)", y = "Partial Residuals") +
      stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
                  col = met.brewer("Hiroshige")[7], fill = "gray70", alpha = 0.3) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = NA),
            axis.line = element_line(colour = "black"))
  } else { p3 <- plot_spacer() }
  
  # Gráfico para Precipitation
  if(!is.null(modelo_sem_precip)) {
    p4 <- ggplot(data.frame(x = data_servico$Precipitation_Log, 
                            y = residuals(modelo_sem_precip, type = "response")), 
                 aes(x = x, y = y)) +
      geom_point(size = 3, col = met.brewer("Egypt")[4], alpha = 0.6) +
      labs(x = "Precipitation (sqrt)", y = "Partial Residuals") +
      stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
                  col = met.brewer("Hiroshige")[7], fill = "gray70", alpha = 0.3) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = NA),
            axis.line = element_line(colour = "black"))
  } else { p4 <- plot_spacer() }
  
  # Combinar
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = paste("Partial Residuals -", servico),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  
  print(combined)
}

dev.off()

cat("PDF salvo em: Figures/Partial_Residuals_All_Services_Complete.pdf\n")

# -------------------------------------------------------------------------
# VERSÃO 3: Usando o pacote mgcv para extrair partial residuals diretamente
# -------------------------------------------------------------------------

# Método alternativo usando predict com type = "terms"
# Isso dá os efeitos suavizados de cada variável

# Para um serviço específico
servico_plot <- unique(out_wide_SA$es_SA)[1]
data_plot <- out_wide_SA %>% filter(es_SA == servico_plot)

# Ajustar modelo
modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                      s(Frequency_Weighted_Distance) + 
                      s(Temperature) + 
                      s(Precipitation_Log), 
                    data = data_plot, 
                    family = poisson(link = "log"))

# Extrair efeitos suavizados (smooth effects)
pred_data <- data_plot
pred_matrix <- predict(modelo, type = "terms", se.fit = FALSE)

# Cada coluna representa o efeito suavizado de cada variável
colnames(pred_matrix) <- c("Isolation_effect", "FWD_effect", "Temperature_effect", "Precipitation_effect")

# Adicionar ao dataframe
data_with_effects <- cbind(data_plot, pred_matrix)

# Plot dos efeitos suavizados
p1_effect <- ggplot(data_with_effects, aes(x = Isolation_Mean, y = Isolation_effect)) +
  geom_point(size = 3, col = met.brewer("Egypt")[1], alpha = 0.6) +
  labs(x = "Isolation Mean (km)", y = "Smooth Effect") +
  stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw()

p2_effect <- ggplot(data_with_effects, aes(x = Frequency_Weighted_Distance, y = FWD_effect)) +
  geom_point(size = 3, col = met.brewer("Egypt")[2], alpha = 0.6) +
  labs(x = "Frequency Weighted Distance", y = "Smooth Effect") +
  stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw()

p3_effect <- ggplot(data_with_effects, aes(x = Temperature, y = Temperature_effect)) +
  geom_point(size = 3, col = met.brewer("Egypt")[3], alpha = 0.6) +
  labs(x = "Temperature (°C)", y = "Smooth Effect") +
  stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw()

p4_effect <- ggplot(data_with_effects, aes(x = Precipitation_Log, y = Precipitation_effect)) +
  geom_point(size = 3, col = met.brewer("Egypt")[4], alpha = 0.6) +
  labs(x = "Precipitation (sqrt)", y = "Smooth Effect") +
  stat_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4), 
              col = met.brewer("Hiroshige")[7], se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw()

combined_effects <- (p1_effect + p2_effect) / (p3_effect + p4_effect) +
  plot_annotation(title = paste("Smooth Effects -", servico_plot))

ggsave(paste0("Figures/Smooth_Effects_", servico_plot, ".jpg"), 
       combined_effects, width = 14, height = 10)

print("Análise de partial residuals concluída!")

# -------------------------------------------------------------------------
# MODELOS GAM COM NEGATIVE BINOMIAL (para serviços com superdispersão)
# -------------------------------------------------------------------------

library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)

# Identificar serviços com superdispersão (Dispersion Ratio > 1.5)
servicos_overdisp <- c("cultural_service_charismatic_species", 
                       "total", 
                       "disease_sentinelling", 
                       "seed_dispersal", 
                       "pest_and_disease_control", 
                       "pollination")

# Função para rodar modelos com diferentes famílias
compare_families <- function(data_servico, servico_nome) {
  
  cat("Processando:", servico_nome, "\n")
  
  # Modelo Poisson
  model_poisson <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Temperature) + 
                s(Precipitation_Log),
              data = data_servico,
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Modelo Negative Binomial
  model_nb <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Temperature) + 
                s(Precipitation_Log),
              data = data_servico,
              family = nb())
  }, error = function(e) NULL)
  
  # Comparação
  resultado <- tibble(
    es_SA = servico_nome,
    Poisson_AIC = ifelse(is.null(model_poisson), NA, AIC(model_poisson)),
    NB_AIC = ifelse(is.null(model_nb), NA, AIC(model_nb)),
    Poisson_Deviance = ifelse(is.null(model_poisson), NA, 
                              1 - (model_poisson$deviance / model_poisson$null.deviance)),
    NB_Deviance = ifelse(is.null(model_nb), NA,
                         1 - (model_nb$deviance / model_nb$null.deviance)),
    Melhor_Modelo = case_when(
      NB_AIC < Poisson_AIC ~ "Negative Binomial",
      TRUE ~ "Poisson"
    )
  )
  
  # Salvar modelos se necessário
  if(!is.null(model_nb)) {
    assign(paste0("model_", gsub(" ", "_", servico_nome), "_NB"), 
           model_nb, envir = .GlobalEnv)
  }
  
  return(resultado)
}

# Aplicar para serviços com superdispersão
comparacao_familias <- map2_dfr(
  out_wide_SA %>% filter(es_SA %in% servicos_overdisp) %>% group_by(es_SA) %>% nest() %>% pull(data),
  servicos_overdisp,
  ~ compare_families(.x, .y)
)

print("Comparação Poisson vs Negative Binomial:")
print(comparacao_familias)

# Salvar comparação
write_csv(comparacao_familias, "Results/Poisson_vs_NB_Comparison.csv")

# -------------------------------------------------------------------------
# RODAR MODELOS FINAIS COM A MELHOR FAMÍLIA PARA CADA SERVIÇO
# -------------------------------------------------------------------------

# Função para rodar modelo final com a família apropriada
run_final_model <- function(data_servico, servico_nome, use_nb = FALSE) {
  
  if(use_nb) {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Temperature) + 
                          s(Precipitation_Log),
                        data = data_servico,
                        family = nb())
    familia_usada <- "Negative Binomial"
  } else {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Temperature) + 
                          s(Precipitation_Log),
                        data = data_servico,
                        family = poisson(link = "log"))
    familia_usada <- "Poisson"
  }
  
  return(list(modelo = modelo, familia = familia_usada))
}

# Definir quais serviços usam NB
servicos_nb <- comparacao_familias %>% 
  filter(Melhor_Modelo == "Negative Binomial") %>% 
  pull(es_SA)

# Rodar modelos finais
modelos_finais <- list()

for(servico in unique(out_wide_SA$es_SA)) {
  data_servico <- out_wide_SA %>% filter(es_SA == servico)
  use_nb <- servico %in% servicos_nb
  
  resultado <- run_final_model(data_servico, servico, use_nb)
  modelos_finais[[servico]] <- resultado$modelo
  
  cat(servico, "->", resultado$familia, "\n")
}

# -------------------------------------------------------------------------
# FIGURAS DE ALTA QUALIDADE PARA PUBLICAÇÃO
# -------------------------------------------------------------------------

library(ggplot2)
library(MetBrewer)
library(patchwork)
library(gridExtra)

# Definir tema para publicação
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(colour = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = base_size + 4, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      legend.position = "bottom",
      legend.key = element_rect(fill = "white"),
      strip.background = element_rect(fill = "gray95", colour = "black"),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

# -------------------------------------------------------------------------
# FIGURA 1: Deviance Explicada por Serviço (Barras)
# -------------------------------------------------------------------------

fig1_deviance <- tabela_geral_formatada %>%
  mutate(es_SA = factor(es_SA, levels = es_SA[order(Deviance_Explained)])) %>%
  ggplot(aes(x = es_SA, y = Deviance_Explained, fill = Quality)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(round(Deviance_Explained, 1), "%")), 
            hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("Excellent" = "#2E8B57")) +
  labs(x = "", 
       y = "Deviance Explained (%)",
       title = "Model Performance by Ecosystem Service",
       fill = "Quality") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip(ylim = c(0, 100)) +
  geom_vline(xintercept = c(30, 50, 70), linetype = "dashed", alpha = 0.3)

ggsave("Figures/Figure1_Deviance_Explained.jpg", fig1_deviance, 
       width = 10, height = 8, dpi = 300)

# -------------------------------------------------------------------------
# FIGURA 2: Importância das Variáveis (Heatmap de p-valores)
# -------------------------------------------------------------------------

# Preparar dados para heatmap
heatmap_data <- tabela_smooth_formatada %>%
  mutate(
    log_p = -log10(p_value + 1e-10),
    Significance_plot = case_when(
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      TRUE ~ "p > 0.05"
    ),
    Variable_clean = case_when(
      Variable == "Isolation_Mean" ~ "Isolation",
      Variable == "Frequency_Weighted_Distance" ~ "FWD",
      Variable == "Temperature" ~ "Temperature",
      Variable == "Precipitation_Log" ~ "Precipitation"
    ),
    es_SA_clean = case_when(
      es_SA == "cultural_service_charismatic_species" ~ "Cultural",
      es_SA == "nutrient_transporting_horizontal_and_vertical" ~ "Nutrient",
      es_SA == "pest_and_disease_control" ~ "Pest control",
      es_SA == "top_down_regulation" ~ "Top-down",
      es_SA == "ecosystem_engineering" ~ "Ecosystem eng.",
      es_SA == "disease_sentinelling" ~ "Disease sent.",
      es_SA == "carrion_control" ~ "Carrion",
      es_SA == "rodent_control" ~ "Rodent",
      es_SA == "seed_dispersal" ~ "Seed dispersal",
      es_SA == "pollination" ~ "Pollination",
      es_SA == "ecotourism" ~ "Ecotourism",
      TRUE ~ es_SA
    )
  )

fig2_heatmap <- ggplot(heatmap_data, aes(x = Variable_clean, y = es_SA_clean, fill = log_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Significance_plot), size = 2.5) +
  scale_fill_gradientn(colors = met.brewer("Hiroshige", 5), 
                       name = "-log10(p-value)",
                       limits = c(0, max(heatmap_data$log_p, na.rm = TRUE))) +
  labs(x = "Environmental Variable", 
       y = "Ecosystem Service",
       title = "Variable Significance Across Ecosystem Services") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        legend.position = "right")

ggsave("Figures/Figure2_Variable_Significance_Heatmap.jpg", fig2_heatmap,
       width = 12, height = 8, dpi = 300)

# -------------------------------------------------------------------------
# FIGURA 3: Partial Residuals para os Melhores Serviços
# -------------------------------------------------------------------------

# Função para criar partial residuals plot de alta qualidade
plot_partial_residuals_pub <- function(data_servico, servico_nome, var_name, var_label, color_palette) {
  
  # Ajustar modelo sem a variável de interesse
  if(var_name == "Isolation_Mean") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                          s(Temperature) + 
                          s(Precipitation_Log),
                        data = data_servico,
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data_servico$Isolation_Mean
  } else if(var_name == "Frequency_Weighted_Distance") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Temperature) + 
                          s(Precipitation_Log),
                        data = data_servico,
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data_servico$Frequency_Weighted_Distance
  } else if(var_name == "Temperature") {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Precipitation_Log),
                        data = data_servico,
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data_servico$Temperature
  } else {
    modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                          s(Frequency_Weighted_Distance) + 
                          s(Temperature),
                        data = data_servico,
                        family = poisson(link = "log"))
    residuals <- residuals(modelo, type = "response")
    x_var <- data_servico$Precipitation_Log
  }
  
  plot_data <- data.frame(x = x_var, y = residuals)
  
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, size = 2, color = color_palette[1]) +
    geom_smooth(method = gam, formula = y ~ s(x, bs = "cs", k = 4),
                color = color_palette[2], fill = color_palette[3], alpha = 0.3) +
    labs(x = var_label, y = "Partial Residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_publication(base_size = 10) +
    theme(plot.margin = margin(5, 5, 5, 5))
}

# Selecionar top 4 serviços para figuras
top_services <- c("cultural_service_charismatic_species", 
                  "nutrient_transporting_horizontal_and_vertical",
                  "seed_dispersal",
                  "ecotourism")

cores_met <- met.brewer("Egypt", 4)

# Criar figura combinada para o melhor serviço
melhor_servico <- top_services[1]
data_melhor <- out_wide_SA %>% filter(es_SA == melhor_servico)

p1 <- plot_partial_residuals_pub(data_melhor, melhor_servico, "Isolation_Mean", 
                                 "Isolation Mean (km)", cores_met[c(1,3,2)])
p2 <- plot_partial_residuals_pub(data_melhor, melhor_servico, "Frequency_Weighted_Distance",
                                 "Freq. Weighted Distance", cores_met[c(2,4,1)])
p3 <- plot_partial_residuals_pub(data_melhor, melhor_servico, "Temperature",
                                 "Temperature (°C)", cores_met[c(3,1,4)])
p4 <- plot_partial_residuals_pub(data_melhor, melhor_servico, "Precipitation_Log",
                                 "Precipitation (sqrt)", cores_met[c(4,2,3)])

fig3_partial <- (p1 + p2) / (p3 + p4) +
  plot_annotation(title = paste("Partial Residuals -", melhor_servico),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

ggsave("Figures/Figure3_Partial_Residuals_Best_Service.jpg", fig3_partial,
       width = 10, height = 8, dpi = 300)

# -------------------------------------------------------------------------
# FIGURA 4: Smooth Effects (curvas GAM)
# -------------------------------------------------------------------------

# Função para plotar smooth effects
plot_smooth_effect <- function(model, var_name, var_label, color) {
  
  summary_model <- summary(model)
  
  # Extrair smooth effect
  if(var_name == "Isolation_Mean") {
    smooth_data <- predict(model, type = "terms", se.fit = TRUE)
    smooth_values <- smooth_data$fit[,1]
    smooth_se <- smooth_data$se.fit[,1]
    x_var <- model$model$Isolation_Mean
  } else if(var_name == "Frequency_Weighted_Distance") {
    smooth_data <- predict(model, type = "terms", se.fit = TRUE)
    smooth_values <- smooth_data$fit[,2]
    smooth_se <- smooth_data$se.fit[,2]
    x_var <- model$model$Frequency_Weighted_Distance
  } else if(var_name == "Temperature") {
    smooth_data <- predict(model, type = "terms", se.fit = TRUE)
    smooth_values <- smooth_data$fit[,3]
    smooth_se <- smooth_data$se.fit[,3]
    x_var <- model$model$Temperature
  } else {
    smooth_data <- predict(model, type = "terms", se.fit = TRUE)
    smooth_values <- smooth_data$fit[,4]
    smooth_se <- smooth_data$se.fit[,4]
    x_var <- model$model$Precipitation_Log
  }
  
  plot_data <- data.frame(
    x = x_var,
    y = smooth_values,
    ymin = smooth_values - 1.96 * smooth_se,
    ymax = smooth_values + 1.96 * smooth_se
  )
  
  # Ordenar para suavizar a linha
  plot_data <- plot_data[order(plot_data$x),]
  
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = color, alpha = 0.3) +
    geom_line(color = color, size = 1.2) +
    labs(x = var_label, y = "Smooth Effect (s)") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_publication(base_size = 10)
}

# Criar figura de smooth effects para o melhor serviço
modelo_melhor <- modelos_finais[[melhor_servico]]

s1 <- plot_smooth_effect(modelo_melhor, "Isolation_Mean", "Isolation Mean (km)", "#1B9E77")
s2 <- plot_smooth_effect(modelo_melhor, "Frequency_Weighted_Distance", "FWD", "#D95F02")
s3 <- plot_smooth_effect(modelo_melhor, "Temperature", "Temperature (°C)", "#7570B3")
s4 <- plot_smooth_effect(modelo_melhor, "Precipitation_Log", "Precipitation (sqrt)", "#E7298A")

fig4_smooth <- (s1 + s2) / (s3 + s4) +
  plot_annotation(title = paste("GAM Smooth Effects -", melhor_servico),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

ggsave("Figures/Figure4_Smooth_Effects_Best_Service.jpg", fig4_smooth,
       width = 10, height = 8, dpi = 300)

# -------------------------------------------------------------------------
# TABELA COMPLETA DE RESULTADOS DOS MODELOS GAM
# Avaliação da qualidade dos modelos
# -------------------------------------------------------------------------

library(mgcv)
library(dplyr)
library(purrr)
library(tidyr)
library(broom)

# -------------------------------------------------------------------------
# FUNÇÃO PARA EXTRAIR TODAS AS ESTATÍSTICAS DE UM MODELO
# -------------------------------------------------------------------------

extract_all_gam_stats <- function(model, data, var_names) {
  tryCatch({
    summary_model <- summary(model)
    
    # Estatísticas gerais do modelo
    general_stats <- tibble(
      Deviance_Explained = summary_model$dev.expl,
      R_sq = summary_model$r.sq,
      AIC = AIC(model),
      BIC = BIC(model),
      GCV = summary_model$sp.criterion,
      n = length(model$residuals),
      Null_Deviance = model$null.deviance,
      Residual_Deviance = model$deviance
    )
    
    # Estatísticas para cada variável smooth
    if(length(summary_model$s.table) > 0) {
      smooth_table <- as.data.frame(summary_model$s.table)
      smooth_stats <- data.frame(
        Variable = var_names,
        edf = smooth_table[, 1],
        Ref_df = smooth_table[, 2],
        F_statistic = smooth_table[, 3],
        p_value = smooth_table[, 4],
        stringsAsFactors = FALSE
      )
    } else {
      smooth_stats <- data.frame(
        Variable = var_names,
        edf = NA,
        Ref_df = NA,
        F_statistic = NA,
        p_value = NA
      )
    }
    
    # Teste de superdispersão (para Poisson)
    pearson_residuals <- residuals(model, type = "pearson")
    dispersion_ratio <- sum(pearson_residuals^2) / (length(pearson_residuals) - length(coef(model)))
    
    general_stats$Dispersion_Ratio <- dispersion_ratio
    
    # Avaliação da superdispersão
    if(dispersion_ratio > 1.5) {
      general_stats$Overdispersion_Warning <- "Yes - consider Negative Binomial"
    } else {
      general_stats$Overdispersion_Warning <- "No"
    }
    
    return(list(general = general_stats, smooth = smooth_stats))
    
  }, error = function(e) {
    return(list(general = tibble(Error = e$message), smooth = NULL))
  })
}

# -------------------------------------------------------------------------
# RODAR MODELOS E EXTRAIR ESTATÍSTICAS PARA CADA SERVIÇO
# -------------------------------------------------------------------------

# Lista de serviços únicos
servicos <- unique(out_wide_SA$es_SA)

# Lista para armazenar resultados
resultados_completos <- list()

for(servico in servicos) {
  
  cat("Processando:", servico, "\n")
  
  # Filtrar dados
  data_servico <- out_wide_SA %>% filter(es_SA == servico)
  
  # Modelo completo (4 variáveis)
  modelo <- mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                        s(Frequency_Weighted_Distance) + 
                        s(Temperature) + 
                        s(Precipitation_Log), 
                      data = data_servico, 
                      family = poisson(link = "log"))
  
  # Nomes das variáveis
  var_nomes <- c("Isolation_Mean", "Frequency_Weighted_Distance", 
                 "Temperature", "Precipitation_Log")
  
  # Extrair estatísticas
  stats <- extract_all_gam_stats(modelo, data_servico, var_nomes)
  
  # Combinar resultados
  resultados_completos[[servico]] <- list(
    general = stats$general %>% mutate(es_SA = servico, .before = 1),
    smooth = stats$smooth %>% mutate(es_SA = servico, .before = 1)
  )
}

# -------------------------------------------------------------------------
# TABELA 1: ESTATÍSTICAS GERAIS DOS MODELOS
# -------------------------------------------------------------------------

tabela_geral <- bind_rows(lapply(resultados_completos, function(x) x$general))

# Formatar números
tabela_geral_formatada <- tabela_geral %>%
  mutate(
    Deviance_Explained = round(Deviance_Explained * 100, 2),
    R_sq = round(R_sq, 3),
    AIC = round(AIC, 1),
    BIC = round(BIC, 1),
    GCV = round(GCV, 4),
    Dispersion_Ratio = round(Dispersion_Ratio, 2)
  ) %>%
  arrange(desc(Deviance_Explained))

# Adicionar classificação de qualidade
tabela_geral_formatada <- tabela_geral_formatada %>%
  mutate(
    Quality = case_when(
      Deviance_Explained >= 70 ~ "Excellent",
      Deviance_Explained >= 50 ~ "Good",
      Deviance_Explained >= 30 ~ "Moderate",
      Deviance_Explained >= 10 ~ "Poor",
      TRUE ~ "Very Poor"
    )
  )

# Visualizar
print(tabela_geral_formatada)

# Salvar
write_csv(tabela_geral_formatada, "Results/Model_Quality_Statistics.csv")

# -------------------------------------------------------------------------
# TABELA 2: ESTATÍSTICAS DAS VARIÁVEIS (edf, F, p-valor)
# -------------------------------------------------------------------------

tabela_smooth <- bind_rows(lapply(resultados_completos, function(x) x$smooth))

# Formatar
tabela_smooth_formatada <- tabela_smooth %>%
  mutate(
    edf = round(edf, 2),
    F_statistic = round(F_statistic, 3),
    p_value = round(p_value, 4),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ "ns"
    ),
    Effect_Type = case_when(
      edf < 1.5 ~ "Linear",
      edf < 3 ~ "Weakly Non-linear",
      TRUE ~ "Strongly Non-linear"
    )
  ) %>%
  arrange(es_SA, desc(p_value))

# Visualizar
print(tabela_smooth_formatada)

# Salvar
write_csv(tabela_smooth_formatada, "Results/Variable_Statistics.csv")

# -------------------------------------------------------------------------
# TABELA 3: COMPARAÇÃO ENTRE MODELOS (com AIC e Deviance)
# -------------------------------------------------------------------------

# Função para comparar diferentes modelos para o mesmo serviço
compare_models_for_service <- function(data_servico) {
  
  # Modelo completo (4 variáveis)
  model_full <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Temperature) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Modelo sem Isolation
  model_no_isolation <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Frequency_Weighted_Distance) + 
                s(Temperature) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Modelo sem FWD
  model_no_fwd <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Temperature) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Modelo sem Temperature
  model_no_temp <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Precipitation_Log), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Modelo sem Precipitation
  model_no_precip <- tryCatch({
    mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                s(Frequency_Weighted_Distance) + 
                s(Temperature), 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Modelo nulo
  model_null <- tryCatch({
    mgcv::gam(Abundance_SA ~ 1, 
              data = data_servico, 
              family = poisson(link = "log"))
  }, error = function(e) NULL)
  
  # Extrair métricas
  result <- tibble(
    Model = c("Full", "No Isolation", "No FWD", "No Temp", "No Precip", "Null"),
    AIC = c(
      ifelse(is.null(model_full), NA, AIC(model_full)),
      ifelse(is.null(model_no_isolation), NA, AIC(model_no_isolation)),
      ifelse(is.null(model_no_fwd), NA, AIC(model_no_fwd)),
      ifelse(is.null(model_no_temp), NA, AIC(model_no_temp)),
      ifelse(is.null(model_no_precip), NA, AIC(model_no_precip)),
      ifelse(is.null(model_null), NA, AIC(model_null))
    ),
    Deviance_Explained = c(
      ifelse(is.null(model_full), NA, 1 - (model_full$deviance / model_full$null.deviance)),
      ifelse(is.null(model_no_isolation), NA, 1 - (model_no_isolation$deviance / model_no_isolation$null.deviance)),
      ifelse(is.null(model_no_fwd), NA, 1 - (model_no_fwd$deviance / model_no_fwd$null.deviance)),
      ifelse(is.null(model_no_temp), NA, 1 - (model_no_temp$deviance / model_no_temp$null.deviance)),
      ifelse(is.null(model_no_precip), NA, 1 - (model_no_precip$deviance / model_no_precip$null.deviance)),
      ifelse(is.null(model_null), NA, 0)
    )
  )
  
  return(result)
}

# Aplicar para todos os serviços
comparacao_modelos <- out_wide_SA %>%
  group_by(es_SA) %>%
  nest() %>%
  mutate(comparison = map(data, compare_models_for_service)) %>%
  unnest(comparison) %>%
  select(-data)

# Formatar
comparacao_modelos_formatada <- comparacao_modelos %>%
  mutate(
    AIC = round(AIC, 1),
    Deviance_Explained = round(Deviance_Explained * 100, 2)
  ) %>%
  arrange(es_SA, AIC)

# Visualizar
print(comparacao_modelos_formatada)

# Salvar
write_csv(comparacao_modelos_formatada, "Results/Model_Comparisons.csv")

# -------------------------------------------------------------------------
# TABELA 4: DIAGNÓSTICO DOS RESÍDUOS
# -------------------------------------------------------------------------

# Função para diagnóstico de resíduos
residual_diagnostics <- function(model, data_servico) {
  tryCatch({
    # Resíduos
    residuals_deviance <- residuals(model, type = "deviance")
    residuals_pearson <- residuals(model, type = "pearson")
    residuals_response <- residuals(model, type = "response")
    
    # Teste de normalidade dos resíduos deviance
    shapiro_test <- tryCatch(shapiro.test(sample(residuals_deviance, min(5000, length(residuals_deviance)))),
                             error = function(e) list(p.value = NA))
    
    # Teste de autocorrelação (se dados têm ordem temporal)
    # library(lmtest)
    # dw_test <- dwtest(model)
    
    # Sumarizar
    tibble(
      Mean_Residual = mean(residuals_deviance, na.rm = TRUE),
      SD_Residual = sd(residuals_deviance, na.rm = TRUE),
      Skewness = moments::skewness(residuals_deviance, na.rm = TRUE),
      Kurtosis = moments::kurtosis(residuals_deviance, na.rm = TRUE),
      Shapiro_p = shapiro_test$p.value,
      Prop_Residuals_Within_2SD = mean(abs(residuals_deviance) < 2, na.rm = TRUE),
      Max_Residual = max(abs(residuals_deviance), na.rm = TRUE)
    )
  }, error = function(e) {
    tibble(
      Mean_Residual = NA, SD_Residual = NA, Skewness = NA,
      Kurtosis = NA, Shapiro_p = NA, Prop_Residuals_Within_2SD = NA,
      Max_Residual = NA
    )
  })
}

# Aplicar para todos os serviços
diagnosticos_residuos <- out_wide_SA %>%
  group_by(es_SA) %>%
  nest() %>%
  mutate(
    model = map(data, ~ mgcv::gam(Abundance_SA ~ s(Isolation_Mean) + 
                                    s(Frequency_Weighted_Distance) + 
                                    s(Temperature) + 
                                    s(Precipitation_Log), 
                                  data = .x, 
                                  family = poisson(link = "log"))),
    diagnostics = map2(model, data, residual_diagnostics)
  ) %>%
  unnest(diagnostics) %>%
  select(-data, -model)

# Formatar
diagnosticos_formatados <- diagnosticos_residuos %>%
  mutate(
    Mean_Residual = round(Mean_Residual, 3),
    SD_Residual = round(SD_Residual, 3),
    Skewness = round(Skewness, 3),
    Kurtosis = round(Kurtosis, 3),
    Shapiro_p = round(Shapiro_p, 4),
    Prop_Residuals_Within_2SD = round(Prop_Residuals_Within_2SD * 100, 1),
    Residual_Quality = case_when(
      abs(Mean_Residual) < 0.1 & Prop_Residuals_Within_2SD > 95 ~ "Good",
      abs(Mean_Residual) < 0.5 & Prop_Residuals_Within_2SD > 90 ~ "Acceptable",
      TRUE ~ "Poor"
    )
  )

# Visualizar
print(diagnosticos_formatados)

# Salvar
write_csv(diagnosticos_formatados, "Results/Residual_Diagnostics.csv")

# -------------------------------------------------------------------------
# TABELA RESUMO FINAL (todas as métricas em uma tabela)
# -------------------------------------------------------------------------

# Juntar todas as informações
tabela_final <- tabela_geral_formatada %>%
  left_join(
    tabela_smooth_formatada %>%
      group_by(es_SA) %>%
      summarise(
        Significant_Variables = paste(Variable[Significance %in% c("***", "**", "*")], collapse = ", "),
        NonLinear_Variables = paste(Variable[Effect_Type == "Strongly Non-linear"], collapse = ", ")
      ),
    by = "es_SA"
  ) %>%
  left_join(
    diagnosticos_formatados %>% select(es_SA, Residual_Quality, Prop_Residuals_Within_2SD),
    by = "es_SA"
  )

# Reordenar colunas
tabela_final <- tabela_final %>%
  select(
    es_SA,
    Quality,
    Deviance_Explained,
    R_sq,
    AIC,
    GCV,
    Dispersion_Ratio,
    Significant_Variables,
    NonLinear_Variables,
    Residual_Quality,
    Prop_Residuals_Within_2SD,
    n
  )

# Visualizar tabela final
print(tabela_final)

# Salvar
write_csv(tabela_final, "Results/Complete_Model_Evaluation.csv")

# -------------------------------------------------------------------------
# RESULTADO GLOBAL (média entre todos os serviços)
# -------------------------------------------------------------------------

resultado_global <- tabela_geral_formatada %>%
  summarise(
    Metric = c("Mean Deviance Explained (%)", "Mean R²", "Mean AIC", "Mean GCV", "Mean Dispersion Ratio"),
    Value = c(
      mean(Deviance_Explained, na.rm = TRUE),
      mean(R_sq, na.rm = TRUE),
      mean(AIC, na.rm = TRUE),
      mean(GCV, na.rm = TRUE),
      mean(Dispersion_Ratio, na.rm = TRUE)
    )
  ) %>%
  mutate(Value = round(Value, 2))

print("=== RESULTADO GLOBAL DOS MODELOS ===")
print(resultado_global)

# -------------------------------------------------------------------------
# GRÁFICO: QUALIDADE DOS MODELOS POR SERVIÇO
# -------------------------------------------------------------------------

library(ggplot2)

# Gráfico de barras da deviance explicada
ggplot(tabela_geral_formatada, aes(x = reorder(es_SA, Deviance_Explained), 
                                   y = Deviance_Explained, 
                                   fill = Quality)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Excellent" = "#2E8B57",
                               "Good" = "#3CB371",
                               "Moderate" = "#FFD700",
                               "Poor" = "#CD853F",
                               "Very Poor" = "#CD5C5C")) +
  labs(x = "Ecosystem Service", 
       y = "Deviance Explained (%)",
       title = "Model Quality by Ecosystem Service",
       fill = "Quality") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  geom_hline(yintercept = c(30, 50, 70), linetype = "dashed", alpha = 0.5) +
  coord_flip()

ggsave("Figures/Model_Quality_by_Service.jpg", width = 12, height = 8)

# -------------------------------------------------------------------------
# GRÁFICO: IMPORTÂNCIA DAS VARIÁVEIS (p-valores)
# -------------------------------------------------------------------------

ggplot(tabela_smooth_formatada, aes(x = reorder(Variable, -log10(p_value + 1e-10)), 
                                    y = -log10(p_value + 1e-10), 
                                    fill = Significance)) +
  geom_bar(stat = "identity") +
  facet_wrap(~es_SA, scales = "free_x") +
  scale_fill_manual(values = c("***" = "red", 
                               "**" = "orange", 
                               "*" = "yellow", 
                               "." = "lightblue", 
                               "ns" = "gray")) +
  labs(x = "Variable", 
       y = "-log10(p-value)",
       title = "Variable Significance by Ecosystem Service",
       fill = "Significance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8))

ggsave("Figures/Variable_Significance.jpg", width = 15, height = 12)

# -------------------------------------------------------------------------
# MENSAGEM FINAL
# -------------------------------------------------------------------------

cat("\n========================================\n")
cat("TABELAS GERADAS:\n")
cat("1. Results/Model_Quality_Statistics.csv - Qualidade geral dos modelos\n")
cat("2. Results/Variable_Statistics.csv - Estatísticas de cada variável\n")
cat("3. Results/Model_Comparisons.csv - Comparação entre modelos\n")
cat("4. Results/Residual_Diagnostics.csv - Diagnóstico dos resíduos\n")
cat("5. Results/Complete_Model_Evaluation.csv - Tabela final completa\n")
cat("========================================\n")

# Imprimir resumo
cat("\n=== RESUMO DA QUALIDADE DOS MODELOS ===\n")
print(tabela_geral_formatada %>% 
        group_by(Quality) %>% 
        summarise(Count = n()) %>%
        arrange(desc(Count)))