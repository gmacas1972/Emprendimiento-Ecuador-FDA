# =============================================================================
# SCRIPT AUTOCONTENIDO: GENERAR GRÁFICAS CORREGIDAS DESDE CERO
# =============================================================================
setwd("C:/Users/GUIDO MACAS/OneDrive/Documentos/Emprendimiento_Ecuador_FDA")
# 1. Cargar librerías necesarias
library(tidyverse)
library(fda)
library(fda.usc)
library(pheatmap)
library(pvclust)
library(cluster)

# 2. Leer los datos que ya procesaste y guardaste en tu carpeta
cat(">>> Cargando datos...\n")
res_area <- read.csv("outputs/tables/aggregated_data_46curves.csv")
matriz_final <- read.csv("outputs/tables/typology_matrix_46curves.csv")

# 3. Función para recrear las curvas funcionales
construir_fdata_area <- function(datos, variable) {
  mat_ancho <- datos %>%
    filter(!is.na(quarter), !is.na(.data[[variable]])) %>%
    mutate(per = year + (quarter - 1) / 4) %>%
    dplyr::select(province, prov_name, area_bin, per, value = all_of(variable)) %>%
    pivot_wider(id_cols = c(province, prov_name, area_bin), names_from = per, values_from = value) %>%
    arrange(province, area_bin)
  
  curve_names <- paste0(mat_ancho$prov_name, ifelse(mat_ancho$area_bin == 0, " (Urb)", " (Rur)"))
  mat <- as.matrix(mat_ancho[, -(1:3)])
  rownames(mat) <- curve_names
  list(fdata = fdata(mat, argvals = as.numeric(colnames(mat_ancho)[-(1:3)])), 
       mat = mat, curve_names = curve_names)
}

# =============================================================================
# GRÁFICA 1: TRAYECTORIAS DE INFORMALIDAD
# =============================================================================
cat(">>> Generando gráfico de trayectorias de informalidad...\n")
obj_informal <- construir_fdata_area(res_area, "informality_rate")
D_inf <- metric.lp(obj_informal$fdata)
cl_inf <- cutree(hclust(as.dist(D_inf), method = "ward.D2"), k = 2)

mat_inf <- obj_informal$mat
mean_c1 <- colMeans(mat_inf[cl_inf == 1, , drop = FALSE])
mean_c2 <- colMeans(mat_inf[cl_inf == 2, , drop = FALSE])
argvals <- as.numeric(colnames(mat_inf))

png("outputs/figures/fig_trajectories_informality.png", width = 1000, height = 700, res = 150)
plot(argvals, mean_c1, type = "l", lwd = 3, col = "black", ylim = c(0, 1),
     xlab = "Year", ylab = "Informality Rate", 
     main = "Informality Trajectories by Cluster")
lines(argvals, mean_c2, lwd = 3, col = "red")
legend("bottomleft", legend = c("Cluster 1 (Low Informality)", "Cluster 2 (High Informality)"), 
       col = c("black", "red"), lwd = 3, bty = "n")
grid()
dev.off()

# =============================================================================
# GRÁFICA 2: HEATMAP CORREGIDO
# =============================================================================
cat(">>> Generando Heatmap corregido...\n")
mat_heat <- as.matrix(matriz_final[, -1])
rownames(mat_heat) <- matriz_final$Curva

png("outputs/figures/fig_heatmap_fixed.png", width = 1200, height = 1400, res = 150)
pheatmap(mat_heat, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         angle_col = 45, 
         main = "Functional Typology (Province-Area)",
         color = colorRampPalette(c("white", "skyblue", "orange", "darkred"))(4),
         fontsize = 10,
         border_color = "grey60")
dev.off()

# =============================================================================
# GRÁFICA 3: DENDROGRAMA PVCLUST
# =============================================================================
cat(">>> Generando Dendrograma Pvclust (esto tarda un poco)...\n")
obj_ingreso <- construir_fdata_area(res_area, "income_median")

png("outputs/figures/fig_pvclust_fixed.png", width = 1400, height = 900, res = 150)
par(mar = c(10, 4, 4, 2) + 0.1) 
set.seed(123)
pv <- pvclust(t(obj_ingreso$mat), method.dist = "euclidean", method.hclust = "ward.D2", nboot = 500)
plot(pv, hang = -1, cex = 0.55, main = "Pvclust Dendrogram - Income (AU Values)")
dev.off()

# =============================================================================
# GRÁFICA 4: GRÁFICO DE SILUETA
# =============================================================================
cat(">>> Generando Gráfico de Silueta...\n")
D_ingreso <- metric.lp(obj_ingreso$fdata)
cl_ingreso <- cutree(hclust(as.dist(D_ingreso), method = "ward.D2"), k = 3)

png("outputs/figures/fig_silhouette_fixed.png", width = 900, height = 700, res = 150)
sil <- silhouette(cl_ingreso, as.dist(D_ingreso))
plot(sil, col = 1:3, main = "Silhouette Plot (k=3) - Income", cex.names = 0.8)
dev.off()

cat("\n=========================================================\n")
cat("¡GRÁFICAS CORREGIDAS GUARDADAS CON ÉXITO!\n")
cat("Ve a la carpeta 'outputs/figures/' y verás los 4 archivos con el sufijo '_fixed' o '_informality'.\n")
cat("=========================================================\n")
