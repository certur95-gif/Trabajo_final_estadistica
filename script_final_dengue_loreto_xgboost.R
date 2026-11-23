###############################################################
# SCRIPT FINAL DE ANALISIS ESTADISTICO - TESIS DENGUE LORETO
# Titulo: Prediccion de casos semanales de dengue en Loreto
# mediante ensamble boosting con variables climaticas y rezagos
# temporales, 2017-2025
#
# Autor del analisis: ajustar nombres segun tesis
# Software: R / RStudio
# Fecha de ejecucion: actualizar al ejecutar
###############################################################

###############################################################
# 0. LIMPIEZA DEL ENTORNO Y PAQUETES
###############################################################

rm(list = ls())
gc()

set.seed(123)
options(scipen = 999)

paquetes <- c(
  "readxl", "dplyr", "tidyr", "janitor", "stringr", "purrr",
  "tibble", "ggplot2", "slider", "xgboost", "Matrix", "vip", "scales",
  "openxlsx"
)
instalar_si_falta <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(paquetes, instalar_si_falta))

# Paquete opcional para optimizacion bayesiana.
# Si no instala o falla, el script usa busqueda aleatoria reproducible.
USAR_OPT_BAYESIANA <- TRUE
if (USAR_OPT_BAYESIANA && !requireNamespace("ParBayesianOptimization", quietly = TRUE)) {
  try(install.packages("ParBayesianOptimization", dependencies = TRUE), silent = TRUE)
}

###############################################################
# 1. CONFIGURACION DE RUTAS
###############################################################

# Recomendacion: coloque este script en la misma carpeta del Excel.
# Si el archivo no se encuentra, RStudio abrira una ventana para seleccionarlo.

archivo_excel <- "clima_departamento_2017_2025_tesis.xlsx"

if (!file.exists(archivo_excel)) {
  message("No se encontro el Excel en el directorio actual. Seleccione el archivo manualmente.")
  archivo_excel <- file.choose()
}

carpeta_salida <- "salidas_tesis_dengue_loreto"
if (!dir.exists(carpeta_salida)) dir.create(carpeta_salida)

###############################################################
# 2. FUNCIONES AUXILIARES
###############################################################

convertir_numero <- function(x) {
  x_chr <- as.character(x)
  x_chr <- stringr::str_replace_all(x_chr, " ", "")
  x_chr <- stringr::str_replace_all(x_chr, ",", ".")
  suppressWarnings(as.numeric(x_chr))
}

calcular_metricas <- function(obs, pred, y_entrenamiento, m = 52) {
  obs <- as.numeric(obs)
  pred <- as.numeric(pred)

  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae  <- mean(abs(obs - pred), na.rm = TRUE)

  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- ifelse(ss_tot == 0, NA_real_, 1 - ss_res / ss_tot)

  mape <- mean(abs((obs - pred) / ifelse(obs == 0, NA, obs)), na.rm = TRUE) * 100

  if (length(y_entrenamiento) > m) {
    denominador_mase <- mean(abs(y_entrenamiento[(m + 1):length(y_entrenamiento)] -
                                    y_entrenamiento[1:(length(y_entrenamiento) - m)]), na.rm = TRUE)
  } else {
    denominador_mase <- NA_real_
  }

  mase <- ifelse(is.na(denominador_mase) || denominador_mase == 0, NA_real_, mae / denominador_mase)

  tibble::tibble(
    RMSE = rmse,
    MAE = mae,
    R2 = r2,
    MAPE = mape,
    MASE = mase
  )
}

remover_predictores_redundantes <- function(datos, predictores, corte = 0.95) {
  predictores <- predictores[predictores %in% names(datos)]
  predictores <- predictores[purrr::map_lgl(datos[predictores], is.numeric)]

  if (length(predictores) <= 2) return(predictores)

  matriz_cor <- suppressWarnings(cor(datos[predictores], use = "pairwise.complete.obs"))
  matriz_cor[is.na(matriz_cor)] <- 0

  eliminar <- character(0)

  repeat {
    matriz_abs <- abs(matriz_cor)
    matriz_abs[lower.tri(matriz_abs, diag = TRUE)] <- 0
    max_cor <- max(matriz_abs, na.rm = TRUE)

    if (max_cor < corte || is.infinite(max_cor)) break

    pos <- which(matriz_abs == max_cor, arr.ind = TRUE)[1, ]
    var1 <- rownames(matriz_abs)[pos[1]]
    var2 <- colnames(matriz_abs)[pos[2]]

    media1 <- mean(abs(matriz_cor[var1, setdiff(colnames(matriz_cor), var1)]), na.rm = TRUE)
    media2 <- mean(abs(matriz_cor[var2, setdiff(colnames(matriz_cor), var2)]), na.rm = TRUE)

    var_eliminar <- ifelse(media1 >= media2, var1, var2)
    eliminar <- c(eliminar, var_eliminar)

    mantener <- setdiff(colnames(matriz_cor), var_eliminar)
    if (length(mantener) <= 2) break
    matriz_cor <- matriz_cor[mantener, mantener, drop = FALSE]
  }

  setdiff(predictores, unique(eliminar))
}

crear_ventanas_rolling <- function(datos, initial = 208, assess = 26, skip = 13) {
  n <- nrow(datos)
  ventanas <- list()
  i <- 1
  fin_train <- initial

  while ((fin_train + assess) <= n) {
    ventanas[[i]] <- list(
      train = 1:fin_train,
      valid = (fin_train + 1):(fin_train + assess)
    )
    fin_train <- fin_train + skip
    i <- i + 1
  }

  ventanas
}

imputar_con_mediana <- function(train_df, new_df, predictores) {
  medianas <- purrr::map_dbl(predictores, function(v) {
    med <- median(train_df[[v]], na.rm = TRUE)
    ifelse(is.na(med), 0, med)
  })
  names(medianas) <- predictores

  for (v in predictores) {
    train_df[[v]][is.na(train_df[[v]])] <- medianas[[v]]
    new_df[[v]][is.na(new_df[[v]])] <- medianas[[v]]
  }

  list(train = train_df, new = new_df, medianas = medianas)
}

entrenar_predecir_xgb <- function(train_df, new_df, predictores, parametros) {
  imp <- imputar_con_mediana(train_df, new_df, predictores)
  train_df <- imp$train
  new_df <- imp$new

  x_train <- as.matrix(train_df[, predictores, drop = FALSE])
  x_new   <- as.matrix(new_df[, predictores, drop = FALSE])

  dtrain <- xgboost::xgb.DMatrix(data = x_train, label = log1p(train_df$casos))
  dnew   <- xgboost::xgb.DMatrix(data = x_new)

  params_xgb <- list(
    objective = "reg:squarederror",
    eval_metric = "mae",
    max_depth = as.integer(parametros$max_depth),
    eta = parametros$eta,
    subsample = parametros$subsample,
    colsample_bytree = parametros$colsample_bytree,
    min_child_weight = parametros$min_child_weight,
    lambda = parametros$lambda,
    alpha = parametros$alpha
  )

  modelo <- xgboost::xgb.train(
    params = params_xgb,
    data = dtrain,
    nrounds = as.integer(parametros$nrounds),
    verbose = 0
  )

  pred <- predict(modelo, dnew)
  pred <- pmax(0, expm1(pred))

  list(modelo = modelo, pred = pred)
}

evaluar_cv_xgb <- function(datos_cv, ventanas, predictores, parametros, metrica = "MAE") {
  resultados <- purrr::map_dfr(ventanas, function(v) {
    train_fold <- datos_cv[v$train, , drop = FALSE]
    valid_fold <- datos_cv[v$valid, , drop = FALSE]

    ajuste <- entrenar_predecir_xgb(train_fold, valid_fold, predictores, parametros)

    calcular_metricas(
      obs = valid_fold$casos,
      pred = ajuste$pred,
      y_entrenamiento = train_fold$casos,
      m = 52
    )
  })

  mean(resultados[[metrica]], na.rm = TRUE)
}

tunear_xgb <- function(datos_cv, ventanas, predictores, nombre_modelo,
                       n_init = 6, n_iter = 12, n_grid_fallback = 30) {

  message("\nTuning del modelo: ", nombre_modelo)

  usar_bayes <- USAR_OPT_BAYESIANA && requireNamespace("ParBayesianOptimization", quietly = TRUE)

  if (usar_bayes) {
    resultado_bayes <- tryCatch({
      fun_objetivo <- function(max_depth, eta, subsample, colsample_bytree,
                               min_child_weight, lambda, alpha, nrounds) {
        parametros <- list(
          max_depth = as.integer(round(max_depth)),
          eta = eta,
          subsample = subsample,
          colsample_bytree = colsample_bytree,
          min_child_weight = min_child_weight,
          lambda = lambda,
          alpha = alpha,
          nrounds = as.integer(round(nrounds))
        )

        mae_cv <- evaluar_cv_xgb(datos_cv, ventanas, predictores, parametros, metrica = "MAE")
        list(Score = -mae_cv, Pred = 0)
      }

      opt <- ParBayesianOptimization::bayesOpt(
        FUN = fun_objetivo,
        bounds = list(
          max_depth = c(2L, 8L),
          eta = c(0.01, 0.25),
          subsample = c(0.60, 1.00),
          colsample_bytree = c(0.60, 1.00),
          min_child_weight = c(1, 20),
          lambda = c(0, 10),
          alpha = c(0, 5),
          nrounds = c(100L, 700L)
        ),
        initPoints = n_init,
        iters.n = n_iter,
        acq = "ucb",
        verbose = 0
      )

      mejor <- ParBayesianOptimization::getBestPars(opt)
      mejor <- as.list(mejor)
      mejor$max_depth <- as.integer(round(mejor$max_depth))
      mejor$nrounds <- as.integer(round(mejor$nrounds))
      mejor$metodo_tuning <- "bayesiano"
      mejor$mae_cv <- -max(opt$scoreSummary$Score, na.rm = TRUE)
      mejor
    }, error = function(e) {
      message("La optimizacion bayesiana fallo. Se usara busqueda aleatoria. Motivo: ", e$message)
      NULL
    })

    if (!is.null(resultado_bayes)) return(resultado_bayes)
  }

  grid <- tibble::tibble(
    max_depth = sample(2:8, n_grid_fallback, replace = TRUE),
    eta = runif(n_grid_fallback, 0.01, 0.25),
    subsample = runif(n_grid_fallback, 0.60, 1.00),
    colsample_bytree = runif(n_grid_fallback, 0.60, 1.00),
    min_child_weight = runif(n_grid_fallback, 1, 20),
    lambda = runif(n_grid_fallback, 0, 10),
    alpha = runif(n_grid_fallback, 0, 5),
    nrounds = sample(seq(100, 700, by = 50), n_grid_fallback, replace = TRUE)
  )

  grid$mae_cv <- purrr::pmap_dbl(grid, function(max_depth, eta, subsample, colsample_bytree,
                                                min_child_weight, lambda, alpha, nrounds) {
    parametros <- list(
      max_depth = as.integer(max_depth),
      eta = eta,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      lambda = lambda,
      alpha = alpha,
      nrounds = as.integer(nrounds)
    )
    evaluar_cv_xgb(datos_cv, ventanas, predictores, parametros, metrica = "MAE")
  })

  mejor <- grid %>% arrange(mae_cv) %>% slice(1) %>% as.list()
  mejor$metodo_tuning <- "busqueda_aleatoria"
  mejor
}

predecir_modelo_final <- function(train_df, eval_df, predictores, parametros, nombre_modelo) {
  ajuste <- entrenar_predecir_xgb(train_df, eval_df, predictores, parametros)

  tibble::tibble(
    modelo = nombre_modelo,
    ano = eval_df$ano,
    semana = eval_df$semana,
    t = eval_df$t,
    observado = eval_df$casos,
    predicho = ajuste$pred,
    residuo = observado - predicho,
    abs_error = abs(residuo)
  )
}

###############################################################
# 3. IMPORTACION Y DEPURACION DE LA BASE ORIGINAL
###############################################################

base_original <- readxl::read_excel(archivo_excel, sheet = "Provincia") %>%
  janitor::clean_names()

# Validacion de columnas obligatorias
columnas_requeridas <- c(
  "departamento", "provincia", "tmean", "tmax", "tmin", "humr", "ptot",
  "ano", "semana", "tmean_clima", "tmax_clima", "tmin_clima",
  "humr_clima", "ptot_clima", "tipo_caso", "casos"
)

faltantes <- setdiff(columnas_requeridas, names(base_original))
if (length(faltantes) > 0) {
  stop("Faltan columnas obligatorias en el Excel: ", paste(faltantes, collapse = ", "))
}

base_limpia <- base_original %>%
  mutate(
    departamento = stringr::str_to_upper(stringr::str_squish(as.character(departamento))),
    provincia = stringr::str_to_upper(stringr::str_squish(as.character(provincia))),
    tipo_caso = stringr::str_to_upper(stringr::str_squish(as.character(tipo_caso))),
    ano = as.integer(ano),
    semana = as.integer(semana),
    across(c(tmean, tmax, tmin, humr, ptot,
             tmean_clima, tmax_clima, tmin_clima, humr_clima, ptot_clima, casos), convertir_numero)
  ) %>%
  filter(departamento == "LORETO", ano >= 2017, ano <= 2025, semana >= 1, semana <= 53)

# Control inicial de estructura
control_inicial <- tibble::tibble(
  indicador = c(
    "Filas originales",
    "Filas filtradas Loreto 2017-2025",
    "Anios incluidos",
    "Provincias incluidas",
    "Semanas minimas",
    "Semanas maximas",
    "Filas con Tipo caso TOTAL",
    "Casos perdidos o no numericos"
  ),
  valor = c(
    nrow(base_original),
    nrow(base_limpia),
    paste(sort(unique(base_limpia$ano)), collapse = ", "),
    length(unique(base_limpia$provincia)),
    min(base_limpia$semana, na.rm = TRUE),
    max(base_limpia$semana, na.rm = TRUE),
    sum(base_limpia$tipo_caso == "TOTAL", na.rm = TRUE),
    sum(is.na(base_limpia$casos))
  )
)

write.csv(control_inicial, file.path(carpeta_salida, "00_control_inicial.csv"), row.names = FALSE)

###############################################################
# 4. CONSOLIDACION DE PROVINCIA-SEMANA A DEPARTAMENTO-SEMANA
###############################################################

variables_clima <- c("tmean", "tmax", "tmin", "humr", "ptot",
                     "tmean_clima", "tmax_clima", "tmin_clima", "humr_clima", "ptot_clima")

# Se usa distinct() para evitar duplicar clima si existieran varias filas por tipo de caso.
clima_departamental <- base_limpia %>%
  select(departamento, provincia, ano, semana, all_of(variables_clima)) %>%
  distinct() %>%
  group_by(ano, semana) %>%
  summarise(
    across(all_of(variables_clima), ~ mean(.x, na.rm = TRUE)),
    n_provincias_clima = n_distinct(provincia),
    .groups = "drop"
  )

casos_departamentales <- base_limpia %>%
  filter(tipo_caso == "TOTAL") %>%
  group_by(ano, semana) %>%
  summarise(
    casos = sum(casos, na.rm = TRUE),
    n_provincias_con_casos = n_distinct(provincia),
    .groups = "drop"
  )

serie_departamental <- clima_departamental %>%
  left_join(casos_departamentales, by = c("ano", "semana")) %>%
  mutate(
    casos = ifelse(is.na(casos), 0, casos),
    n_provincias_con_casos = ifelse(is.na(n_provincias_con_casos), 0, n_provincias_con_casos)
  ) %>%
  arrange(ano, semana) %>%
  mutate(
    t = row_number(),
    periodo = paste0(ano, "-S", stringr::str_pad(semana, width = 2, pad = "0"))
  )

write.csv(serie_departamental, file.path(carpeta_salida, "01_base_departamental_semanal.csv"), row.names = FALSE)

###############################################################
# 5. ANALISIS DESCRIPTIVO DE LA SERIE
###############################################################

descriptivos_casos <- serie_departamental %>%
  summarise(
    n_semanas = n(),
    media = mean(casos, na.rm = TRUE),
    mediana = median(casos, na.rm = TRUE),
    desviacion_estandar = sd(casos, na.rm = TRUE),
    minimo = min(casos, na.rm = TRUE),
    p25 = quantile(casos, 0.25, na.rm = TRUE),
    p75 = quantile(casos, 0.75, na.rm = TRUE),
    p90 = quantile(casos, 0.90, na.rm = TRUE),
    maximo = max(casos, na.rm = TRUE),
    total_casos = sum(casos, na.rm = TRUE)
  )

casos_por_ano <- serie_departamental %>%
  group_by(ano) %>%
  summarise(
    semanas = n(),
    total_casos = sum(casos, na.rm = TRUE),
    media_semanal = mean(casos, na.rm = TRUE),
    mediana_semanal = median(casos, na.rm = TRUE),
    maximo_semanal = max(casos, na.rm = TRUE),
    .groups = "drop"
  )

estacionalidad_semana <- serie_departamental %>%
  group_by(semana) %>%
  summarise(
    media_casos = mean(casos, na.rm = TRUE),
    mediana_casos = median(casos, na.rm = TRUE),
    p90_casos = quantile(casos, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

top_semanas <- serie_departamental %>%
  arrange(desc(casos)) %>%
  slice(1:20)

write.csv(descriptivos_casos, file.path(carpeta_salida, "02_descriptivos_casos.csv"), row.names = FALSE)
write.csv(casos_por_ano, file.path(carpeta_salida, "03_casos_por_ano.csv"), row.names = FALSE)
write.csv(estacionalidad_semana, file.path(carpeta_salida, "04_estacionalidad_por_semana.csv"), row.names = FALSE)
write.csv(top_semanas, file.path(carpeta_salida, "05_top_20_semanas_mayor_incidente.csv"), row.names = FALSE)

# Grafico 1: serie temporal de casos
ggplot(serie_departamental, aes(x = t, y = casos)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "Casos semanales de dengue en Loreto, 2017-2025",
    x = "Orden temporal semanal",
    y = "Casos semanales"
  ) +
  theme_minimal()
ggsave(file.path(carpeta_salida, "grafico_01_serie_temporal_casos.png"), width = 10, height = 5, dpi = 300)

# Grafico 2: total anual
ggplot(casos_por_ano, aes(x = factor(ano), y = total_casos)) +
  geom_col() +
  labs(
    title = "Total anual de casos de dengue en Loreto",
    x = "Anio",
    y = "Total de casos"
  ) +
  theme_minimal()
ggsave(file.path(carpeta_salida, "grafico_02_total_anual_casos.png"), width = 8, height = 5, dpi = 300)

# Grafico 3: estacionalidad semanal promedio
ggplot(estacionalidad_semana, aes(x = semana, y = media_casos)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Estacionalidad promedio semanal del dengue en Loreto",
    x = "Semana epidemiologica",
    y = "Media de casos"
  ) +
  theme_minimal()
ggsave(file.path(carpeta_salida, "grafico_03_estacionalidad_semanal.png"), width = 9, height = 5, dpi = 300)

# ACF y PACF
png(file.path(carpeta_salida, "grafico_04_acf_casos.png"), width = 1100, height = 700)
acf(serie_departamental$casos, lag.max = 104, main = "ACF de casos semanales de dengue")
dev.off()

png(file.path(carpeta_salida, "grafico_05_pacf_casos.png"), width = 1100, height = 700)
pacf(serie_departamental$casos, lag.max = 104, main = "PACF de casos semanales de dengue")
dev.off()

###############################################################
# 6. CONSTRUCCION DE VARIABLES DERIVADAS Y REZAGOS
###############################################################

serie_modelado <- serie_departamental %>%
  arrange(ano, semana) %>%
  mutate(
    tmean_anom = tmean - tmean_clima,
    tmax_anom  = tmax  - tmax_clima,
    tmin_anom  = tmin  - tmin_clima,
    humr_anom  = humr  - humr_clima,
    ptot_anom  = ptot  - ptot_clima,
    sin_semana = sin(2 * pi * semana / 52),
    cos_semana = cos(2 * pi * semana / 52)
  )

# Rezagos de casos: corto plazo, mediano plazo y referencia anual.
lags_casos <- c(1:12, 26, 52)
for (k in lags_casos) {
  serie_modelado[[paste0("casos_lag", k)]] <- dplyr::lag(serie_modelado$casos, k)
}

# Medias moviles calculadas solo con semanas anteriores.
serie_modelado <- serie_modelado %>%
  mutate(
    casos_ma4  = slider::slide_dbl(dplyr::lag(casos, 1), mean, .before = 3,  .complete = TRUE, na.rm = TRUE),
    casos_ma8  = slider::slide_dbl(dplyr::lag(casos, 1), mean, .before = 7,  .complete = TRUE, na.rm = TRUE),
    casos_ma12 = slider::slide_dbl(dplyr::lag(casos, 1), mean, .before = 11, .complete = TRUE, na.rm = TRUE)
  )

# Rezagos climaticos candidatos. Se excluye el clima contemporaneo para evitar fuga temporal.
variables_clima_modelo <- c("tmean", "tmax", "tmin", "humr", "ptot",
                            "tmean_anom", "tmax_anom", "tmin_anom", "humr_anom", "ptot_anom")
lags_clima <- 1:12

for (v in variables_clima_modelo) {
  for (k in lags_clima) {
    serie_modelado[[paste0(v, "_lag", k)]] <- dplyr::lag(serie_modelado[[v]], k)
  }
}

# Base final de modelado: se excluyen filas sin rezagos completos.
base_modelado <- serie_modelado %>%
  filter(!is.na(casos_lag52)) %>%
  arrange(ano, semana)

write.csv(base_modelado, file.path(carpeta_salida, "06_base_modelado_con_rezagos.csv"), row.names = FALSE)

###############################################################
# 7. CORRELACION CRUZADA Y REDUNDANCIA ENTRE PREDICTORES
###############################################################

ccf_rezagos <- tidyr::expand_grid(
  variable = variables_clima_modelo,
  lag = lags_clima
) %>%
  mutate(
    predictor = paste0(variable, "_lag", lag),
    correlacion_con_casos = purrr::map_dbl(predictor, ~ cor(base_modelado$casos, base_modelado[[.x]], use = "complete.obs"))
  ) %>%
  arrange(desc(abs(correlacion_con_casos)))

write.csv(ccf_rezagos, file.path(carpeta_salida, "07_correlacion_cruzada_clima_casos.csv"), row.names = FALSE)

predictores_temporales <- c("sin_semana", "cos_semana")
predictores_casos <- c(paste0("casos_lag", lags_casos), "casos_ma4", "casos_ma8", "casos_ma12")
predictores_clima <- paste0(rep(variables_clima_modelo, each = length(lags_clima)), "_lag", rep(lags_clima, times = length(variables_clima_modelo)))

predictores_climatico_raw <- c(predictores_temporales, predictores_clima)
predictores_autoreg_raw <- c(predictores_temporales, predictores_casos)
predictores_combinado_raw <- c(predictores_temporales, predictores_casos, predictores_clima)

# Particiones temporales segun el plan de tesis.
base_train_cv <- base_modelado %>% filter(ano >= 2017, ano <= 2023)
base_validacion_2024 <- base_modelado %>% filter(ano == 2024)
base_prueba_2025 <- base_modelado %>% filter(ano == 2025)
base_train_2017_2024 <- base_modelado %>% filter(ano >= 2017, ano <= 2024)

if (nrow(base_train_cv) == 0 || nrow(base_validacion_2024) == 0 || nrow(base_prueba_2025) == 0) {
  stop("No hay observaciones suficientes en entrenamiento 2017-2023, validacion 2024 o prueba 2025. Revise la base.")
}

# Remocion de redundancia usando solo entrenamiento 2017-2023.
predictores_climatico <- remover_predictores_redundantes(base_train_cv, predictores_climatico_raw, corte = 0.95)
predictores_autoreg <- remover_predictores_redundantes(base_train_cv, predictores_autoreg_raw, corte = 0.95)
predictores_combinado <- remover_predictores_redundantes(base_train_cv, predictores_combinado_raw, corte = 0.95)

lista_predictores <- list(
  climatico = predictores_climatico,
  autoregresivo = predictores_autoreg,
  combinado = predictores_combinado
)

resumen_predictores <- tibble::tibble(
  modelo = names(lista_predictores),
  n_predictores = purrr::map_int(lista_predictores, length),
  predictores = purrr::map_chr(lista_predictores, ~ paste(.x, collapse = "; "))
)

write.csv(resumen_predictores, file.path(carpeta_salida, "08_predictores_finales_por_modelo.csv"), row.names = FALSE)

# Pares de alta correlacion del modelo combinado antes de depurar.
matriz_cor <- suppressWarnings(cor(base_train_cv[predictores_combinado_raw], use = "pairwise.complete.obs"))
matriz_cor[is.na(matriz_cor)] <- 0
pares_alta_cor <- which(abs(matriz_cor) >= 0.95 & upper.tri(matriz_cor), arr.ind = TRUE)
pares_alta_cor <- tibble::tibble(
  variable_1 = rownames(matriz_cor)[pares_alta_cor[, 1]],
  variable_2 = colnames(matriz_cor)[pares_alta_cor[, 2]],
  correlacion = matriz_cor[pares_alta_cor]
) %>%
  arrange(desc(abs(correlacion)))

write.csv(pares_alta_cor, file.path(carpeta_salida, "09_pares_predictores_alta_correlacion.csv"), row.names = FALSE)

###############################################################
# 8. VALIDACION TEMPORAL ROLLING ORIGIN Y TUNING DE XGBOOST
###############################################################

ventanas <- crear_ventanas_rolling(base_train_cv, initial = min(208, max(120, floor(nrow(base_train_cv) * 0.60))), assess = 26, skip = 13)

if (length(ventanas) < 2) {
  stop("No se pudieron crear suficientes ventanas rolling origin. Revise el numero de semanas disponibles.")
}

hiperparametros <- list(
  climatico = tunear_xgb(base_train_cv, ventanas, predictores_climatico, "XGBoost climatico"),
  autoregresivo = tunear_xgb(base_train_cv, ventanas, predictores_autoreg, "XGBoost autoregresivo"),
  combinado = tunear_xgb(base_train_cv, ventanas, predictores_combinado, "XGBoost combinado")
)

hiperparametros_tabla <- purrr::imap_dfr(hiperparametros, function(x, nm) {
  tibble::tibble(
    modelo = nm,
    max_depth = x$max_depth,
    eta = x$eta,
    subsample = x$subsample,
    colsample_bytree = x$colsample_bytree,
    min_child_weight = x$min_child_weight,
    lambda = x$lambda,
    alpha = x$alpha,
    nrounds = x$nrounds,
    metodo_tuning = x$metodo_tuning,
    mae_cv = x$mae_cv
  )
})

write.csv(hiperparametros_tabla, file.path(carpeta_salida, "10_hiperparametros_xgboost.csv"), row.names = FALSE)

###############################################################
# 9. VALIDACION EXTERNA INTERMEDIA: ANIO 2024
###############################################################

pred_validacion <- dplyr::bind_rows(
  predecir_modelo_final(base_train_cv, base_validacion_2024, predictores_climatico, hiperparametros$climatico, "XGBoost climatico"),
  predecir_modelo_final(base_train_cv, base_validacion_2024, predictores_autoreg, hiperparametros$autoregresivo, "XGBoost autoregresivo"),
  predecir_modelo_final(base_train_cv, base_validacion_2024, predictores_combinado, hiperparametros$combinado, "XGBoost combinado"),
  tibble::tibble(
    modelo = "Ingenuo estacional lag52",
    ano = base_validacion_2024$ano,
    semana = base_validacion_2024$semana,
    t = base_validacion_2024$t,
    observado = base_validacion_2024$casos,
    predicho = pmax(0, base_validacion_2024$casos_lag52),
    residuo = observado - predicho,
    abs_error = abs(residuo)
  )
)

metricas_validacion <- pred_validacion %>%
  group_by(modelo) %>%
  group_modify(~ calcular_metricas(.x$observado, .x$predicho, y_entrenamiento = base_train_cv$casos, m = 52)) %>%
  ungroup() %>%
  arrange(MAE)

write.csv(pred_validacion, file.path(carpeta_salida, "11_predicciones_validacion_2024.csv"), row.names = FALSE)
write.csv(metricas_validacion, file.path(carpeta_salida, "12_metricas_validacion_2024.csv"), row.names = FALSE)

# Seleccion del mejor modelo por MAE en 2024. Puede cambiarse a RMSE o MASE si el asesor lo exige.
mejor_modelo_nombre <- metricas_validacion %>%
  filter(modelo != "Ingenuo estacional lag52") %>%
  arrange(MAE) %>%
  slice(1) %>%
  pull(modelo)

message("\nMejor modelo en validacion 2024: ", mejor_modelo_nombre)

###############################################################
# 10. PRUEBA FINAL: ANIO 2025
###############################################################

pred_prueba <- dplyr::bind_rows(
  predecir_modelo_final(base_train_2017_2024, base_prueba_2025, predictores_climatico, hiperparametros$climatico, "XGBoost climatico"),
  predecir_modelo_final(base_train_2017_2024, base_prueba_2025, predictores_autoreg, hiperparametros$autoregresivo, "XGBoost autoregresivo"),
  predecir_modelo_final(base_train_2017_2024, base_prueba_2025, predictores_combinado, hiperparametros$combinado, "XGBoost combinado"),
  tibble::tibble(
    modelo = "Ingenuo estacional lag52",
    ano = base_prueba_2025$ano,
    semana = base_prueba_2025$semana,
    t = base_prueba_2025$t,
    observado = base_prueba_2025$casos,
    predicho = pmax(0, base_prueba_2025$casos_lag52),
    residuo = observado - predicho,
    abs_error = abs(residuo)
  )
)

metricas_prueba <- pred_prueba %>%
  group_by(modelo) %>%
  group_modify(~ calcular_metricas(.x$observado, .x$predicho, y_entrenamiento = base_train_2017_2024$casos, m = 52)) %>%
  ungroup() %>%
  arrange(MAE)

write.csv(pred_prueba, file.path(carpeta_salida, "13_predicciones_prueba_2025.csv"), row.names = FALSE)
write.csv(metricas_prueba, file.path(carpeta_salida, "14_metricas_prueba_2025.csv"), row.names = FALSE)

###############################################################
# 11. IMPORTANCIA DE VARIABLES DEL MODELO COMBINADO
###############################################################

# Se entrena el modelo combinado con 2017-2024 para reportar importancia final.
ajuste_combinado_final <- entrenar_predecir_xgb(
  train_df = base_train_2017_2024,
  new_df = base_prueba_2025,
  predictores = predictores_combinado,
  parametros = hiperparametros$combinado
)

importancia_combinado <- xgboost::xgb.importance(
  feature_names = predictores_combinado,
  model = ajuste_combinado_final$modelo
) %>%
  as_tibble()

write.csv(importancia_combinado, file.path(carpeta_salida, "15_importancia_variables_xgboost_combinado.csv"), row.names = FALSE)

if (nrow(importancia_combinado) > 0) {
  graf_importancia <- importancia_combinado %>%
    slice_max(Gain, n = min(20, nrow(importancia_combinado))) %>%
    mutate(Feature = reorder(Feature, Gain)) %>%
    ggplot(aes(x = Feature, y = Gain)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Importancia de variables - XGBoost combinado",
      x = "Variable",
      y = "Gain"
    ) +
    theme_minimal()

  ggsave(file.path(carpeta_salida, "grafico_06_importancia_variables_xgboost_combinado.png"),
         graf_importancia, width = 9, height = 7, dpi = 300)
}

###############################################################
# 12. ANALISIS DE RESIDUOS Y SEMANAS DE ALTA INCIDENCIA
###############################################################

umbral_p90 <- quantile(base_modelado$casos, 0.90, na.rm = TRUE)

pred_prueba_residuos <- pred_prueba %>%
  mutate(
    alta_incidencia = observado >= umbral_p90,
    error_porcentual = ifelse(observado == 0, NA_real_, 100 * residuo / observado)
  )

residuos_alta_incidencia <- pred_prueba_residuos %>%
  group_by(modelo, alta_incidencia) %>%
  summarise(
    n = n(),
    observado_promedio = mean(observado, na.rm = TRUE),
    predicho_promedio = mean(predicho, na.rm = TRUE),
    MAE = mean(abs_error, na.rm = TRUE),
    sesgo_promedio = mean(residuo, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(pred_prueba_residuos, file.path(carpeta_salida, "16_residuos_prueba_2025.csv"), row.names = FALSE)
write.csv(residuos_alta_incidencia, file.path(carpeta_salida, "17_resumen_residuos_alta_incidencia.csv"), row.names = FALSE)

# Grafico observado vs predicho en prueba final
ggplot(pred_prueba, aes(x = t)) +
  geom_line(aes(y = observado), linewidth = 0.8) +
  geom_line(aes(y = predicho), linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ modelo, scales = "free_y") +
  labs(
    title = "Casos observados vs predichos en prueba final 2025",
    x = "Orden temporal semanal",
    y = "Casos"
  ) +
  theme_minimal()
ggsave(file.path(carpeta_salida, "grafico_07_observado_vs_predicho_2025.png"), width = 12, height = 8, dpi = 300)

# Grafico de residuos del modelo combinado
pred_combinado_test <- pred_prueba_residuos %>% filter(modelo == "XGBoost combinado")

ggplot(pred_combinado_test, aes(x = t, y = residuo)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = alta_incidencia), size = 2) +
  labs(
    title = "Residuos del modelo XGBoost combinado en prueba 2025",
    x = "Orden temporal semanal",
    y = "Residuo: observado - predicho",
    shape = "Alta incidencia"
  ) +
  theme_minimal()
ggsave(file.path(carpeta_salida, "grafico_08_residuos_xgboost_combinado_2025.png"), width = 10, height = 5, dpi = 300)

# ACF de residuos del modelo combinado
png(file.path(carpeta_salida, "grafico_09_acf_residuos_xgboost_combinado_2025.png"), width = 1100, height = 700)
acf(pred_combinado_test$residuo, lag.max = min(30, nrow(pred_combinado_test) - 1),
    main = "ACF de residuos - XGBoost combinado 2025")
dev.off()

###############################################################
# 13. TABLA FINAL DE DECISION DE HIPOTESIS PREDICTIVA
###############################################################

metricas_comparativas <- metricas_prueba %>%
  mutate(
    supera_ingenuo_segun_MASE = ifelse(MASE < 1, "Si", "No"),
    interpretacion_MASE = dplyr::case_when(
      is.na(MASE) ~ "No evaluable",
      MASE < 1 ~ "Mejor que el modelo ingenuo estacional",
      MASE == 1 ~ "Equivalente al modelo ingenuo estacional",
      MASE > 1 ~ "Peor que el modelo ingenuo estacional"
    )
  )

write.csv(metricas_comparativas, file.path(carpeta_salida, "18_metricas_comparativas_decision_hipotesis.csv"), row.names = FALSE)

###############################################################
# 13.1 EXPORTAR TODOS LOS CSV EN UN SOLO ARCHIVO EXCEL
###############################################################

# Instalar y cargar openxlsx si no está disponible
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", dependencies = TRUE)
}
suppressPackageStartupMessages(library(openxlsx))

# Archivo Excel consolidado de salida
archivo_excel_salida <- file.path(
  carpeta_salida,
  "resultados_tesis_dengue_loreto_consolidado.xlsx"
)

# Buscar todos los CSV generados en la carpeta de salida
archivos_csv <- list.files(
  path = carpeta_salida,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(archivos_csv) == 0) {
  warning("No se encontraron archivos CSV en la carpeta de salida.")
} else {
  
  # Crear libro Excel
  wb <- openxlsx::createWorkbook()
  
  # Estilos
  estilo_titulo <- openxlsx::createStyle(
    textDecoration = "bold",
    fontSize = 12,
    fgFill = "#D9EAF7",
    border = "Bottom"
  )
  
  estilo_encabezado <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#B7DEE8",
    border = "TopBottomLeftRight",
    halign = "center",
    valign = "center"
  )
  
  estilo_normal <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    valign = "top"
  )
  
  # Función para crear nombres válidos de hoja Excel
  crear_nombre_hoja <- function(ruta_archivo, usados) {
    nombre <- tools::file_path_sans_ext(basename(ruta_archivo))
    
    # Quitar caracteres no permitidos por Excel
    nombre <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", nombre)
    
    # Máximo 31 caracteres
    nombre <- substr(nombre, 1, 31)
    
    # Evitar nombres repetidos
    nombre_base <- nombre
    contador <- 1
    
    while (nombre %in% usados) {
      sufijo <- paste0("_", contador)
      nombre <- substr(nombre_base, 1, 31 - nchar(sufijo))
      nombre <- paste0(nombre, sufijo)
      contador <- contador + 1
    }
    
    nombre
  }
  
  hojas_usadas <- character(0)
  
  # Crear hoja índice
  openxlsx::addWorksheet(wb, "INDICE")
  
  indice <- tibble::tibble(
    numero = seq_along(archivos_csv),
    archivo_csv = basename(archivos_csv),
    hoja_excel = NA_character_
  )
  
  # Exportar cada CSV como una hoja
  for (i in seq_along(archivos_csv)) {
    
    archivo <- archivos_csv[i]
    
    nombre_hoja <- crear_nombre_hoja(archivo, hojas_usadas)
    hojas_usadas <- c(hojas_usadas, nombre_hoja)
    indice$hoja_excel[i] <- nombre_hoja
    
    datos <- read.csv(
      archivo,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8"
    )
    
    openxlsx::addWorksheet(wb, nombre_hoja)
    
    # Título
    openxlsx::writeData(
      wb,
      sheet = nombre_hoja,
      x = paste0("Archivo fuente: ", basename(archivo)),
      startRow = 1,
      startCol = 1
    )
    
    openxlsx::addStyle(
      wb,
      sheet = nombre_hoja,
      style = estilo_titulo,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    
    # Tabla
    openxlsx::writeDataTable(
      wb,
      sheet = nombre_hoja,
      x = datos,
      startRow = 3,
      startCol = 1,
      tableStyle = "TableStyleMedium2",
      withFilter = TRUE
    )
    
    # Formato
    if (ncol(datos) > 0) {
      openxlsx::addStyle(
        wb,
        sheet = nombre_hoja,
        style = estilo_encabezado,
        rows = 3,
        cols = 1:ncol(datos),
        gridExpand = TRUE
      )
      
      openxlsx::setColWidths(
        wb,
        sheet = nombre_hoja,
        cols = 1:ncol(datos),
        widths = "auto"
      )
      
      openxlsx::freezePane(
        wb,
        sheet = nombre_hoja,
        firstActiveRow = 4,
        firstActiveCol = 1
      )
    }
  }
  
  # Escribir índice
  openxlsx::writeData(
    wb,
    sheet = "INDICE",
    x = "RESULTADOS CONSOLIDADOS - TESIS DENGUE LORETO",
    startRow = 1,
    startCol = 1
  )
  
  openxlsx::addStyle(
    wb,
    sheet = "INDICE",
    style = estilo_titulo,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  
  openxlsx::writeDataTable(
    wb,
    sheet = "INDICE",
    x = indice,
    startRow = 3,
    startCol = 1,
    tableStyle = "TableStyleMedium9",
    withFilter = TRUE
  )
  
  openxlsx::setColWidths(
    wb,
    sheet = "INDICE",
    cols = 1:ncol(indice),
    widths = "auto"
  )
  
  openxlsx::freezePane(
    wb,
    sheet = "INDICE",
    firstActiveRow = 4,
    firstActiveCol = 1
  )
  
  # Guardar Excel consolidado
  openxlsx::saveWorkbook(
    wb,
    archivo_excel_salida,
    overwrite = TRUE
  )
  
  cat("\nExcel consolidado generado correctamente:\n")
  cat(archivo_excel_salida, "\n")
}



###############################################################
# 14. RESUMEN AUTOMATICO EN CONSOLA
###############################################################

cat("\n============================================================\n")
cat("ANALISIS FINALIZADO\n")
cat("============================================================\n")
cat("Archivo analizado: ", archivo_excel, "\n", sep = "")
cat("Carpeta de salida: ", carpeta_salida, "\n", sep = "")
cat("Semanas consolidadas: ", nrow(serie_departamental), "\n", sep = "")
cat("Semanas usadas en modelado: ", nrow(base_modelado), "\n", sep = "")
cat("Entrenamiento/CV: 2017-2023 = ", nrow(base_train_cv), " semanas\n", sep = "")
cat("Validacion: 2024 = ", nrow(base_validacion_2024), " semanas\n", sep = "")
cat("Prueba final: 2025 = ", nrow(base_prueba_2025), " semanas\n", sep = "")
cat("\nMetricas de validacion 2024:\n")
print(metricas_validacion)
cat("\nMetricas de prueba final 2025:\n")
print(metricas_prueba)
cat("\nUmbral operativo de alta incidencia P90: ", round(umbral_p90, 2), " casos semanales\n", sep = "")
cat("\nRevise los archivos CSV y PNG generados en la carpeta de salida.\n")
cat("============================================================\n")
