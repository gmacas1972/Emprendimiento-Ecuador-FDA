# =============================================================================
# SCRIPT MAESTRO: ANÁLISIS DE DATOS FUNCIONALES Y VALIDACIÓN PARA WOS/SCOPUS
# Emprendedores (P42=6) por Provincia y Área (Urbano/Rural) | 2021-2025
# =============================================================================

# ---------- 0. CONFIGURACIÓN INICIAL ----------
paquetes <- c("fda", "fda.usc", "tidyverse", "haven", "srvyr", 
              "NbClust", "fpc", "pvclust", "mclust", "cluster", 
              "pheatmap", "gridExtra")
nuevos <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(nuevos)) install.packages(nuevos, dependencies = TRUE)

library(fda); library(fda.usc); library(tidyverse); library(haven)
library(srvyr); library(NbClust); library(fpc); library(pvclust)
library(mclust); library(cluster); library(pheatmap); library(gridExtra)

dir.create("Resultados_Finales", showWarnings = FALSE)
cat(">>> Iniciando análisis maestro...\n")

# ---------- 1. CARGA DE DATOS Y CREACIÓN DE VARIABLES ----------
cache_file <- "df_enemdu_master.rds"

if (file.exists(cache_file)) {
  df <- readRDS(cache_file)
} else {
  archivos <- paste0("D:/temp_enemdu_", 2021:2025, ".dta")
  vars_needed <- c("prov", "periodo", "area", "upm", "estrato", "fexp", "condact",
                   "p02", "p03", "p42", "p10a", "p10b", "p49", "p05a", "p05b", 
                   "p51a", "p51b", "p51c", "ingrl", "p63", "p64b", "p65", "p66", 
                   "p67", "p68b", "p69", "p70b", "p71b", "p72b", "p73b", "p74b", "p76", "p78")
  
  leer_anio <- function(path) {
    read_dta(path, col_select = any_of(vars_needed)) %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  }
  
  df <- map_dfr(archivos, leer_anio)
  
  prov_lookup <- c("1"="Azuay","2"="Bolívar","3"="Cañar","4"="Carchi","5"="Cotopaxi",
                   "6"="Chimborazo","7"="El Oro","8"="Esmeraldas","9"="Guayas","10"="Imbabura",
                   "11"="Loja","12"="Los Ríos","13"="Manabí","14"="Morona Santiago","15"="Napo",
                   "16"="Pastaza","17"="Pichincha","18"="Tungurahua","19"="Zamora Chinchipe",
                   "20"="Galápagos","21"="Sucumbíos","22"="Orellana","23"="Santo Domingo","24"="Santa Elena")
  
  df <- df %>%
    mutate(
      anio = periodo %/% 100,
      trimestre = ceiling((periodo %% 100) / 3),
      prov_nombre = dplyr::recode(as.character(prov), !!!prov_lookup, .default = "Desconocida"),
      area_bin = if_else(area == 2, 1, 0) # 0 = Urbano, 1 = Rural
    )
  saveRDS(df, cache_file)
}

# Construcción de variables robusta
df <- df %>%
  mutate(
    across(c(p63, p64b, p65, p66, p67, p68b, p69, p70b, p71b, p72b, p73b, p74b, p76, p78),
           ~ if_else(is.na(.x) | .x == 999999 | .x < 0 | (!is.na(ingrl) & ingrl == -1), 0, .x)),
    ingtotal = p63 + p64b - p65 + p66 - p67 + p68b + p69 + p70b + p71b + p72b + p73b + p74b + p76 + p78,
    ingtotal = if_else(is.na(ingtotal), 0, ingtotal),
    yeduc = case_when(
      p10a == 1 & p10b == 0 ~ 0, p10a == 2 & p10b == 0 ~ 0, p10a == 2 & p10b == 1 ~ 2,
      p10a == 2 & p10b == 2 ~ 4, p10a == 2 & p10b == 3 ~ 6, p10a == 2 & p10b == 4 ~ 7,
      p10a == 3 ~ 1, p10a == 4 ~ 1 + p10b, p10a == 5 ~ p10b, p10a == 6 ~ 7 + p10b,
      p10a == 7 ~ 10 + p10b, p10a == 8 ~ 13 + p10b, p10a == 9 ~ 13 + p10b, p10a == 10 ~ 18 + p10b,
      TRUE ~ 0
    ),
    informal = as.integer(p49 %in% c(2, 3))
  )

# ---------- 2. AGREGACIÓN MUESTRAL (SIN GALÁPAGOS) ----------
res_area <- df %>%
  filter(p42 == 6, provincia != 20) %>% 
  as_survey_design(ids = upm, strata = estrato, weights = fexp, nest = TRUE) %>%
  group_by(provincia, prov_nombre, area_bin, anio, trimestre) %>%
  summarise(
    ingreso_med   = survey_median(ingtotal, na.rm = TRUE, vartype = NULL),
    yeduc_prom    = survey_mean(yeduc, na.rm = TRUE, vartype = NULL),
    prop_informal = survey_mean(informal, na.rm = TRUE, vartype = NULL),
    prop_mujeres  = survey_mean(p02 == 2, na.rm = TRUE, vartype = NULL),
    horas_prom    = survey_mean(p51a + p51b + p51c, na.rm = TRUE, vartype = NULL),
    .groups = "drop"
  )

# ---------- 3. FUNCIONES FDA Y CLUSTERING ----------
construir_fdata_area <- function(datos, variable) {
  mat_ancho <- datos %>%
    filter(!is.na(trimestre), !is.na(.data[[variable]])) %>%
    mutate(per = anio + (trimestre - 1) / 4) %>%
    dplyr::select(provincia, prov_nombre, area_bin, per, value = all_of(variable)) %>%
    pivot_wider(id_cols = c(provincia, prov_nombre, area_bin), names_from = per, values_from = value) %>%
    arrange(provincia, area_bin)
  
  nombres_curva <- paste0(mat_ancho$prov_nombre, ifelse(mat_ancho$area_bin == 0, " (Urb)", " (Rur)"))
  mat <- as.matrix(mat_ancho[, -(1:3)])
  rownames(mat) <- nombres_curva
  list(fdata = fdata(mat, argvals = as.numeric(colnames(mat_ancho)[-(1:3)])), 
       mat = mat, prov_nom = nombres_curva)
}

analisis_fda_completo <- function(obj, titulo, k, exportar_grafico = TRUE) {
  fd <- obj$fdata
  out <- outliers.depth.pond(fd, nb = 200, suav = 0.05, trim = 0.1)
  idx_atip <- which(obj$prov_nom %in% as.character(out$outliers))
  
  if(length(idx_atip) > 0) {
    fd_cl <- fd[-idx_atip, ]; prov_cl <- obj$prov_nom[-idx_atip]
  } else {
    fd_cl <- fd; prov_cl <- obj$prov_nom
  }
  
  D <- metric.lp(fd_cl)
  hc <- hclust(as.dist(D), method = "ward.D2")
  cl <- cutree(hc, k = k)
  
  if(exportar_grafico) {
    png(paste0("Resultados_Finales/Trayectorias_", gsub(" ", "_", titulo), ".png"), width = 1000, height = 700, res = 120)
    plot(fd_cl, col = cl, lwd = 2, main = titulo, ylab = "Valor", xlab = "Año")
    legend("topleft", legend = paste("Cluster", sort(unique(cl))), col = sort(unique(cl)), lwd = 2, bty = "n")
    dev.off()
  }
  invisible(list(cluster = cl, prov_cluster = prov_cl, atipicas = obj$prov_nom[idx_atip], fdata = fd_cl, D = D))
}

# ---------- 4. EJECUCIÓN DE ANÁLISIS Y VALIDACIÓN ----------
variables <- list(
  list(var = "ingreso_med",  k = 3, nombre = "Ingreso Total"),
  list(var = "yeduc_prom",   k = 4, nombre = "Años de Educacion"),
  list(var = "prop_informal",k = 2, nombre = "Informalidad"),
  list(var = "prop_mujeres", k = 3, nombre = "Proporcion Mujeres"),
  list(var = "horas_prom",   k = 3, nombre = "Horas Promedio")
)

resultados_clusters <- list()
tabla_robustez <- data.frame()

for (v in variables) {
  cat("\n======================================================\n")
  cat("PROCESANDO Y VALIDANDO:", v$nombre, "\n")
  
  obj <- construir_fdata_area(res_area, v$var)
  res <- analisis_fda_completo(obj, v$nombre, k = v$k)
  resultados_clusters[[v$var]] <- res
  
  # PRUEBAS DE VALIDACIÓN (Solo para Ingreso como prueba principal, y Jaccard para todas)
  set.seed(123)
  boot <- clusterboot(as.dist(res$D), B = 500, clustermethod = disthclustCBI, 
                      method = "ward.D2", k = v$k, count = FALSE)
  jaccard_prom <- mean(boot$bootmean)
  
  anova <- fanova.onefactor(res$fdata, as.factor(res$cluster), nboot = 500)
  
  cat("Jaccard promedio:", round(jaccard_prom, 3), "| ANOVA p-valor:", anova$pvalue, "\n")
  
  tabla_robustez <- rbind(tabla_robustez, data.frame(
    Dimension = v$nombre,
    Clusters_k = v$k,
    Jaccard_Medio = round(jaccard_prom, 3),
    ANOVA_p_valor = anova$pvalue
  ))
}

# ---------- 5. PRUEBAS AVANZADAS DE SENSIBILIDAD (SOLO INGRESO) ----------
cat("\n>>> Ejecutando pruebas avanzadas para Ingreso (NbClust, pvclust, ARI)...\n")
obj_ingreso <- construir_fdata_area(res_area, "ingreso_med")

# 5.1 NbClust
nb <- NbClust(data = obj_ingreso$mat, distance = "euclidean", min.nc = 2, max.nc = 6, method = "ward.D2", index = "all")
write.csv(as.data.frame(nb$Best.nc), "Resultados_Finales/Validacion_NbClust_Ingreso.csv")

# 5.2 Pvclust (Guardar Gráfico)
png("Resultados_Finales/Validacion_Pvclust_Ingreso.png", width = 1200, height = 800, res = 120)
pv <- pvclust(t(obj_ingreso$mat), method.dist = "euclidean", method.hclust = "ward.D2", nboot = 500)
plot(pv, hang = -1, cex = 0.6, main = "Dendrograma Pvclust (Valores AU/BP) - Ingreso")
dev.off()

# 5.3 K-means ARI
fd_cl_ing <- resultados_clusters[["ingreso_med"]]$fdata
km <- kmeans.fd(fd_cl_ing, ncl = 3, draw = FALSE)
ari <- adjustedRandIndex(resultados_clusters[["ingreso_med"]]$cluster, km$cluster)
cat("Índice de Rand Ajustado (ARI Jerárquico vs K-means):", round(ari, 3), "\n")

# 5.4 Gráfico de Silueta
png("Resultados_Finales/Validacion_Silueta_Ingreso.png", width = 800, height = 600, res = 120)
sil <- silhouette(resultados_clusters[["ingreso_med"]]$cluster, as.dist(resultados_clusters[["ingreso_med"]]$D))
plot(sil, col = 1:3, main = "Gráfico de Silueta (k=3) - Ingreso")
dev.off()

# ---------- 6. MATRIZ FINAL Y TABLAS CRUZADAS ----------
cat("\n>>> Construyendo matriz de tipología y tablas cruzadas...\n")
todas_curvas <- resultados_clusters[["ingreso_med"]]$prov_cluster
matriz_final <- data.frame(Curva = todas_curvas, stringsAsFactors = FALSE)

for (v in variables) {
  res <- resultados_clusters[[v$var]]
  col_name <- v$nombre
  matriz_final[[col_name]] <- res$cluster[match(matriz_final$Curva, res$prov_cluster)]
}

write.csv(matriz_final, "Resultados_Finales/Matriz_Tipologia_Final.csv", row.names = FALSE)
write.csv(tabla_robustez, "Resultados_Finales/Tabla_Robustez_Metodologica.csv", row.names = FALSE)

# Heatmap de la tipología final
mat_heat <- as.matrix(matriz_final[, -1])
rownames(mat_heat) <- matriz_final$Curva
png("Resultados_Finales/Heatmap_Tipologia.png", width = 900, height = 1000, res = 120)
pheatmap(mat_heat, cluster_rows = TRUE, cluster_cols = FALSE, 
         na_col = "white", main = "Tipología Provincial (Urbano/Rural)",
         color = colorRampPalette(c("white", "skyblue", "orange", "darkred"))(4))
dev.off()

# Imprimir tablas cruzadas clave
cat("\n=== TABLA CRUZADA: INGRESO vs INFORMALIDAD ===\n")
print(table(matriz_final$`Ingreso Total`, matriz_final$Informalidad, useNA = "ifany"))

cat("\n=== TABLA CRUZADA: INGRESO vs EDUCACIÓN ===\n")
print(table(matriz_final$`Ingreso Total`, matriz_final$`Años de Educacion`, useNA = "ifany"))

cat("\n=========================================================\n")
cat("¡ANÁLISIS MAESTRO COMPLETADO CON ÉXITO!\n")
cat("Revisa la carpeta 'Resultados_Finales' para encontrar:\n")
cat(" - 5 gráficos de trayectorias (.png)\n")
cat(" - 1 heatmap de tipología (.png)\n")
cat(" - 2 gráficos de validación (Pvclust y Silueta)\n")
cat(" - 3 tablas de datos (.csv): Matriz, Robustez y NbClust.\n")
cat("=========================================================\n")