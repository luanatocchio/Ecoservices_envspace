# ANÁLISE DE MASSA CORPORAL (LOG-TRANSFORMADA) NO ESPAÇO AMBIENTAL
library(mgcv)
library(tidyverse)
library(letsR)
library(terra)
library(geodata)
library(sf)
library(broom)

# 1. CARREGAR E PREPARAR OS DADOS 

bio_data <- worldclim_global(var = "bio", res = 10, 
                             path = "Data/Primary Data/worldclim")
clim <- terra::subset(bio_data, c(1, 12))

SA <- st_read("Data/Primary Data/SA/Lim_america_do_sul_2021.shp")
SA <- st_transform(SA, crs(clim))

clim_SA0 <- terra::crop(clim, SA)
clim_SA <- terra::mask(clim_SA0, SA)

v <- values(clim_SA)[, 2]
logv <- sqrt(v)
values(clim_SA)[, 2] <- ifelse(is.infinite(logv), NA, logv)

body_mass_data <- read.csv("Data/Primary Data/body_mass_only.csv")
print(head(body_mass_data))
vale_species <- read.csv("Data/Primary Data/vale_species_last.csv")[[1]]
pam_mammals_SA <- lets.load(file = "Data/Secondary Data/pam_SA.RData") 
pam_mammals_sub_SA <- lets.subsetPAM(x = pam_mammals_SA, names = vale_species, remove.cells = FALSE)
ES_SA <- read.csv("Data/Primary Data/ES_11_TraitsNOVO.csv")
res_SA <- lets.envpam(pam_mammals_sub_SA, envs_SA, n_bins = 30)
envs_SA <- lets.addvar(pam_mammals_sub_SA, clim_SA, onlyvar = TRUE)
vale_species <- read.csv("Data/Primary Data/vale_species_last.csv")[[1]]
pam_mammals_sub_SA <- lets.subsetPAM(x = pam_mammals_SA, names = vale_species, remove.cells = FALSE)

# Detectar o separador automaticamente
primeira_linha <- as.character(body_mass_data[1, 1])
if(grepl(";", primeira_linha)) {
  separador <- ";"
} else if(grepl(",", primeira_linha)) {
  separador <- ","
} else if(grepl(" ", primeira_linha)) {
  separador <- " "
} else if(grepl("\\.", primeira_linha)) {
  separador <- "\\."
} else {
  separador <- NULL
}

if(!is.null(separador)) {
  body_mass_separated <- body_mass_data %>%
    separate(1, into = c("species", "body_mass"), sep = separador) %>%
    mutate(body_mass = as.numeric(body_mass))
} else {
  body_mass_separated <- body_mass_data %>%
    mutate(species = as.character(.[,1]),
           body_mass = as.numeric(str_extract(.[,1], "[0-9]+(\\.[0-9]+)?")))
}

body_mass_separated <- body_mass_separated %>%
  filter(!is.na(species), !is.na(body_mass), body_mass > 0) %>%
  mutate(body_mass_log10 = log10(body_mass))

# 2. FILTRAR ESPÉCIES PRESENTES NO PAM
species_in_pam <- ES_SA$species
BM_filtered <- body_mass_separated %>% filter(species %in% species_in_pam)
if(nrow(BM_filtered) == 0) stop("ERRO: Nenhuma espécie com dados de massa corporal encontrada no PAM!")

# 3. CRIAR DATAFRAMES PARA MÉDIA E SOMA
BM_mean <- BM_filtered %>% select(species, body_mass = body_mass_log10)
BM_sum_raw <- BM_filtered %>% select(species, body_mass_raw = body_mass)

# 4. MAPEAR MASSA CORPORAL NO ESPAÇO AMBIENTAL
media_com_na <- function(x) mean(x, na.rm = TRUE)
soma_com_na <- function(x) sum(x, na.rm = TRUE)

res_map_BM_mean <- lets.maplizer.env(x = res_SA, y = BM_mean$body_mass, z = BM_mean$species, func = media_com_na)
res_map_BM_sum <- lets.maplizer.env(x = res_SA, y = BM_sum_raw$body_mass_raw, z = BM_sum_raw$species, func = soma_com_na)

# 5. EXTRAIR DADOS
keep_BM <- res_SA$Presence_and_Absence_Matrix_env[, 1]
out_BM_cells <- lets.envcells(res_SA)

abund_env_BM_mean <- res_map_BM_mean$Matrix_env[, -(1:3), drop = FALSE]
abund_env_BM_sum_raw <- res_map_BM_sum$Matrix_env[, -(1:3), drop = FALSE]
abund_env_BM_sum_log10 <- log10(abund_env_BM_sum_raw)

# 6. PREPARAR TABELAS DE DADOS
preparar_dados <- function(abund_vec, nome_metrica) {
  tibble(metric = nome_metrica, Body_Mass = abund_vec[, 1], out_BM_cells[keep_BM, -1], .name_repair = "unique") %>%
    drop_na() %>%
    select(metric, Body_Mass, Frequency, `Isolation (Mean)`, `Frequency Weighted Distance`) %>%
    rename(Isolation_Mean = `Isolation (Mean)`, FWD = `Frequency Weighted Distance`) %>%
    pivot_longer(cols = c(Frequency, Isolation_Mean, FWD), names_to = "Variable", values_to = "Value")
}

data_plot_BM_mean <- preparar_dados(abund_env_BM_mean, "body_mass_log10_mean")
data_plot_BM_sum <- preparar_dados(abund_env_BM_sum_log10, "body_mass_log10_sum")
data_plot_BM_combined <- bind_rows(data_plot_BM_mean, data_plot_BM_sum)

# 7. RODAR MODELOS GAM
gam_results_mean <- data_plot_BM_mean %>%
  group_by(Variable) %>%
  summarise(model = list(mgcv::gam(Body_Mass ~ s(Value), data = pick(everything()))), .groups = "drop") %>%
  mutate(Metric = "Mean")

gam_results_sum <- data_plot_BM_sum %>%
  group_by(Variable) %>%
  summarise(model = list(mgcv::gam(Body_Mass ~ s(Value), data = pick(everything()))), .groups = "drop") %>%
  mutate(Metric = "Sum")

# 8. EXTRAIR RESULTADOS
extrair_resultados <- function(gam_results_df) {
  gam_results_df %>%
    mutate(
      summary_stats = map(model, function(m) {
        s <- summary(m)
        tibble(Deviance_Explained = ifelse(!is.null(s$dev.expl), s$dev.expl, NA),
               AIC = AIC(m), R_sq = ifelse(!is.null(s$r.sq), s$r.sq, NA))
      }),
      coefficients = map(model, function(m) {
        coef_summary <- summary(m)$s.table
        if(length(coef_summary) > 0) {
          tibble(edf = coef_summary[, 1], Ref_df = coef_summary[, 2], 
                 F_statistic = coef_summary[, 3], p_value = coef_summary[, 4])
        } else {
          tibble(edf = NA, Ref_df = NA, F_statistic = NA, p_value = NA)
        }
      })
    ) %>%
    unnest(c(summary_stats, coefficients), keep_empty = TRUE) %>%
    select(Metric, Variable, edf, Ref_df, F_statistic, p_value, Deviance_Explained, AIC, R_sq) %>%
    mutate(p_value = round(p_value, 4), Deviance_Explained = round(Deviance_Explained * 100, 1),
           AIC = round(AIC, 1), R_sq = round(R_sq, 3), edf = round(edf, 2),
           F_statistic = round(F_statistic, 3),
           Significance = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                                    p_value < 0.05 ~ "*", is.na(p_value) ~ "NA", TRUE ~ "ns"))
}

table_mean <- extrair_resultados(gam_results_mean)
table_sum <- extrair_resultados(gam_results_sum)
table_BM <- bind_rows(table_mean, table_sum)

# 9. RESULTADOS
print(table_BM)

# 10. GRÁFICO COMPARATIVO (Mean vs Sum)
p_comparison <- data_plot_BM_combined %>%
  ggplot(aes(y = Body_Mass, x = Value, color = metric)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, alpha = 0.3) +
  geom_point(alpha = 0.2, size = 0.8) +
  facet_grid(metric ~ Variable, scales = "free") +
  scale_color_manual(values = c("body_mass_log10_mean" = "darkblue", "body_mass_log10_sum" = "darkred")) +
  theme_bw() +
  labs(title = "Body Mass (log10-transformed): Mean vs Sum", y = "log10(Body Mass)", x = "Variable Value", color = "Metric") +
  theme(strip.background = element_rect(fill = "gray90"), strip.text = element_text(face = "bold", size = 10), legend.position = "bottom")

print(p_comparison)
ggsave("Figures/Body_Mass_Mean_vs_Sum_Comparison.jpg", p_comparison, width = 12, height = 10, dpi = 300)

# 11. DIAGNÓSTICO DE RESÍDUOS
pdf("Figures/Residuals_Diagnosis_Body_Mass_Log10.pdf", width = 10, height = 14)
par(mfrow = c(2, 2))
for(i in 1:nrow(gam_results_mean)) {
  if(!is.null(gam_results_mean$model[[i]])) {
    mgcv::gam.check(gam_results_mean$model[[i]])
    mtext(paste("MEAN - Variable:", gam_results_mean$Variable[i]), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 0.8)
  }
}
for(i in 1:nrow(gam_results_sum)) {
  if(!is.null(gam_results_sum$model[[i]])) {
    mgcv::gam.check(gam_results_sum$model[[i]])
    mtext(paste("SUM - Variable:", gam_results_sum$Variable[i]), side = 3, line = -1.5, outer = TRUE, font = 2, cex = 0.8)
  }
}
dev.off()

# 12. MAPAS AMBIENTAIS
dir.create("Figures/Maps_env_BodyMass", recursive = TRUE, showWarnings = FALSE)

tryCatch({
  png("Figures/Maps_env_BodyMass/Body_Mass_Mean_Environmental_Space.png", width = 3000, height = 1500, res = 300)
  lets.plot.envpam(res_map_BM_mean, world = TRUE)
  title("Body Mass (Mean) in Environmental Space", line = 2, cex.main = 1.5)
  dev.off()
}, error = function(e) cat("Erro no mapa da média:", e$message, "\n"))

tryCatch({
  png("Figures/Maps_env_BodyMass/Body_Mass_Sum_Environmental_Space.png", width = 3000, height = 1500, res = 300)
  lets.plot.envpam(res_map_BM_sum, world = TRUE)
  title("Body Mass (Sum) in Environmental Space", line = 2, cex.main = 1.5)
  dev.off()
}, error = function(e) cat("Erro no mapa da soma:", e$message, "\n"))

View(table_BM)


