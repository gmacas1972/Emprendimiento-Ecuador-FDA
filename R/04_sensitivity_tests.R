
# =============================================================================
# SCRIPT 04: SENSITIVITY TESTS (Smoothing & Outliers)
# =============================================================================
library(tidyverse)
library(fda)
library(fda.usc)
library(mclust)

# Asumimos que 'res_area' y 'construir_fdata_area' están cargados en memoria.
# Si no, debes cargar el caché y definir la función nuevamente aquí.

cat(">>> Running Sensitivity Test 1: Smoothing (Basis Dimension)...\n")
mat_income <- construir_fdata_area(res_area, "income_median")$mat
argvals <- as.numeric(colnames(mat_income))

# Bases: Baja (5) y Alta (15)
basis_low <- create.bspline.basis(rangeval = range(argvals), nbasis = 5)
basis_high <- create.bspline.basis(rangeval = range(argvals), nbasis = 15)

fd_low <- smooth.basis(argvals, t(mat_income), basis_low)$y
fd_high <- smooth.basis(argvals, t(mat_income), basis_high)$y

D_low <- metric.lp(fdata(t(fd_low), argvals = argvals))
D_high <- metric.lp(fdata(t(fd_high), argvals = argvals))

cl_low <- cutree(hclust(as.dist(D_low), method = "ward.D2"), k = 3)
cl_high <- cutree(hclust(as.dist(D_high), method = "ward.D2"), k = 3)

# Original (46 curvas)
D_orig <- metric.lp(fdata(mat_income, argvals = argvals))
cl_original_46 <- cutree(hclust(as.dist(D_orig), method = "ward.D2"), k = 3)

# Calcular ARI
ari_low <- adjustedRandIndex(cl_original_46, cl_low)
ari_high <- adjustedRandIndex(cl_original_46, cl_high)

# Guardar Tabla de Suavizado
resultado_sensibilidad <- data.frame(
  Test = c("Smoothing Low (5 bases)", "Smoothing High (15 bases)"),
  ARI = c(ari_low, ari_high)
)
write.csv(resultado_sensibilidad, "outputs/tables/sensitivity_smoothing_test.csv", row.names = FALSE)


cat(">>> Running Sensitivity Test 2: Outliers...\n")
# Si no tienes el objeto 'res_income' de 43 curvas, recalcúlalo rápido:
if(!exists("res_income")) {
  fd_orig <- construir_fdata_area(res_area, "income_median")$fdata
  out <- outliers.depth.pond(fd_orig, nb = 200, suav = 0.05, trim = 0.1)
  idx_atip <- which(construir_fdata_area(res_area, "income_median")$curve_names %in% as.character(out$outliers))
  fd_cl <- fd_orig[-idx_atip, ]
  D_cl <- metric.lp(fd_cl)
  cl_43 <- cutree(hclust(as.dist(D_cl), method = "ward.D2"), k = 3)
  prov_cl_43 <- construir_fdata_area(res_area, "income_median")$curve_names[-idx_atip]
} else {
  cl_43 <- res_income$cluster
  prov_cl_43 <- res_income$prov_cluster
}

# Comparar 46 vs 43
df_comparacion <- data.frame(
  Curva = names(cl_original_46),
  Cluster_All46 = as.numeric(cl_original_46)
)
df_comparacion$Cluster_Orig43 <- cl_43[match(df_comparacion$Curva, prov_cl_43)]

# Guardar Tabla de Contingencia de Atípicos
tabla_outliers <- as.data.frame(table(df_comparacion$Cluster_All46, df_comparacion$Cluster_Orig43, useNA = "ifany"))
colnames(tabla_outliers) <- c("Cluster_46", "Cluster_43", "Freq")
write.csv(tabla_outliers, "outputs/tables/sensitivity_outliers_contingency.csv", row.names = FALSE)

cat(">>> All sensitivity tables saved in outputs/tables/\n")





