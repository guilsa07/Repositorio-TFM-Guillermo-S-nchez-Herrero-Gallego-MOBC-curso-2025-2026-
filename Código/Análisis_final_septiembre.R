# ==============================================================================
# Análisis FINAL del TFM
# Guillermo Sánchez-Herrero Gallego
# ==============================================================================



# 0. Antes que nada ------------------------------------------------------------

cat("\014")    # Limpiar la consola
rm(list=ls())
# setwd(dirección/de/tu/directorio)
getwd()



# 1. Paquetes ------------------------------------------------------------------

# Importación, manipulación y visualización
library(tidyverse)     # Lectura, manipulación de datos (`dplyr`, `tidyr`) y gráficos (`ggplot2`)
library(openxlsx)      # Lectura y escritura de archivos Excel (`read.xlsx`, `write.xlsx`)
library(patchwork)     # Combinación y composición de múltiples figuras (`+`, `/`)

# Diversidad taxonómica y funcional
library(vegan)         # Diversidad taxonómica (Hill), disimilitud (Bray-Curtis), PERMANOVA (`adonis2`) y ord. (`metaMDS`)
library(FD)            # Matriz de distancias de Gower (`gowdis`) y medias ponderadas de rasgos (`functcomp`)
library(fundiversity)  # Cálculo rápido de índices funcionales (`fd_fric`, `fd_feve`, `fd_fdis`, `fd_raoq`)

# Análisis espacial y datos ornitológicos
library(sf)            # Procesamiento y manipulación de datos espaciales vectoriales (`st_read`, `st_transform`)
library(terra)         # Análisis y manejo de datos espaciales raster y vectoriales (`rast`, `extract`)
library(tmap)          # Visualización cartográfica y mapas temáticos (`tm_shape`, `tm_polygons`)
library(ebirdst)       # Descarga y análisis de modelos espacio-temporales de eBird (`ebirdst_download_status`)
library(exactextractr) # Extracción rápida de estadísticas de raster mediante polígonos (`exact_extract`)

# Estadística general y modelado
library(DescTools)     # Estadística descriptiva y pruebas complementarias (`PostHocTest`, `AUC`)
library(rstatix)       # Pruebas estadísticas integradas en tuberías/pipes (`shapiro_test`, `wilcox_test`, `anova_test`)
library(mgcv)          # Modelos Aditivos Generalizados (GAM/GAMM) y suavizados (`gam`, `s`)
library(emmeans)       # Cálculo de medias marginales estimadas y contrastes post-hoc (`emmeans`, `pairs`)

# Configuración global
theme_set(theme_bw(base_size = 12))
tmap_mode("plot")



# 2. Carga de los objetos, comprobaciones y preparativos -----------------------

# Variable explicativa -> Trayectoria
Trayectoria <- readRDS("Datos_codigo/Trayectorias.rds")

# Variable explicativa -> Uso actual
Uso_actual <- readRDS("Datos_codigo/Uso_actual.rds")

# Variable explicativa -> Estabilidad
Estabilidad <- readRDS("Datos_codigo/Estabilidad.rds")

# Variable explicativa -> Intensidad del cambio
Intensidad_cambio <- readRDS("Datos_codigo/Porcentaje_cambio.rds")

# Variable respuesta -> Rasgos funcionales de las especies
Rasgos <- read.csv ("Datos_codigo/Rasgos.csv", 
                    fileEncoding = "Latin1", 
                    sep = ",", 
                    stringsAsFactors = FALSE)

# Variable respuesta -> Abundancia bruta de las especies
Abundancia <- readRDS("Datos_codigo/Abundancia_bruta.rds")

# Referencias espaciales
grid_vect <- readRDS("Datos_codigo/grid_vect.rds")
grid_ebird <- readRDS("Datos_codigo/grid_ebird.rds")
recorte_SC <- readRDS("Datos_codigo/recorte_SC.rds")

# Comprobación (mapas de abundancia por especie)
Abundancia_sf <- grid_ebird %>%
  left_join(Abundancia, by = "ID")

tm_shape(Abundancia_sf) +
  tm_fill(
    fill = "Pica pica",
    fill.scale = tm_scale(values = "viridis"),
    fill.legend = tm_legend(title = "Abundancia")
  ) +
  tm_borders() +
  tm_scalebar(position = c("left", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_title("Abundancia reproductora de Pica pica")


# Conversión de ID de celda a formato character
Trayectoria$ID <- as.character(Trayectoria$ID)
Estabilidad$ID <- as.character(Estabilidad$ID)
Intensidad_cambio$ID <- as.character(Intensidad_cambio$ID)
Uso_actual$ID <- as.character(Uso_actual$ID)
grid_ebird$ID <- as.character(grid_ebird$ID)
Abundancia$ID <- as.character(Abundancia$ID)



# 3. Pruebas de independencia de las variables explicativas --------------------


# Trayectoria ~ Estabilidad 

# Base de datos
dependencia <- Trayectoria %>%
  st_drop_geometry() %>%
  select(ID, Tipo) %>%
  left_join(
    Estabilidad %>%
      st_drop_geometry() %>%
      select(ID, anio_estable),
    by = "ID"
  )

# Tabla de contingencia
tabla <- table(
  dependencia$Tipo,
  dependencia$anio_estable
)

# Prueba χ² (simulación Monte Carlo)
chi <- chisq.test(
  tabla,
  simulate.p.value = TRUE,
  B = 10000
)

chi

# Tamaño del efecto (V de Cramér)
DescTools::CramerV(tabla)

# Distribución porcentual por trayectoria
round(
  prop.table(tabla, margin = 1) * 100,
  1
)


# Porcentaje de cambio ~ Trayectoria 

# Base de datos
datos_kw <- Trayectoria %>%
  st_drop_geometry() %>%
  select(ID, Tipo) %>%
  left_join(
    Intensidad_cambio %>%
      st_drop_geometry() %>%
      select(ID, Porcentaje_cambio),
    by = "ID"
  )

# Prueba de Kruskal-Wallis
kruskal.test(
  Porcentaje_cambio ~ Tipo,
  data = datos_kw
)

# Tamaño del efecto (ε²)
datos_kw %>%
  kruskal_effsize(
    Porcentaje_cambio ~ Tipo
  )

# Estadísticos descriptivos
datos_kw %>%
  group_by(Tipo) %>%
  summarise(
    Mediana = median(Porcentaje_cambio),
    IQR = IQR(Porcentaje_cambio),
    n = n()
  ) %>%
  arrange(desc(Mediana))

# Representación gráfica
ggplot(
  datos_kw,
  aes(Tipo, Porcentaje_cambio)
) +
  geom_boxplot() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# Porcentaje de cambio ~ Estabilidad 

# Base de datos
datos_kw2 <- Estabilidad %>%
  st_drop_geometry() %>%
  select(ID, anio_estable) %>%
  left_join(
    Intensidad_cambio %>%
      st_drop_geometry() %>%
      select(ID, Porcentaje_cambio),
    by = "ID"
  )

# Prueba de Kruskal-Wallis
kruskal.test(
  Porcentaje_cambio ~ anio_estable,
  data = datos_kw2
)

# Tamaño del efecto (ε²)
datos_kw2 %>%
  kruskal_effsize(
    Porcentaje_cambio ~ anio_estable
  )

# Estadísticos descriptivos
datos_kw2 %>%
  group_by(anio_estable) %>%
  summarise(
    Mediana = median(Porcentaje_cambio),
    IQR = IQR(Porcentaje_cambio),
    n = n()
  )

# Representación gráfica
ggplot(
  datos_kw2,
  aes(anio_estable, Porcentaje_cambio)
) +
  geom_boxplot()



# 4. Limpieza y alineamiento de las tablas -------------------------------------

# Reordenar la tabla de rasgos para que coincida con la matriz de abundancia
rasgos <- Rasgos %>%
  arrange(match(Especie, names(Abundancia)[-1])) %>%
  column_to_rownames("Especie")

# Matriz de abundancia
abund <- Abundancia %>%
  column_to_rownames("ID") %>%
  as.matrix()

# Comprobación crítica
stopifnot(identical(colnames(abund), rownames(rasgos)))

cat(
  "Alineado:",
  ncol(abund), "especies x",
  nrow(abund), "celdas.\n"
)

# Creación de una submuestra de celdas
n_sub <- NULL      # <- Pon NULL para utilizar todas las celdas

set.seed(2026)

if (!is.null(n_sub) && n_sub < nrow(abund)) {
  
  filas <- sample(nrow(abund), n_sub)
  
  abund <- abund[filas, , drop = FALSE]
  
}

# Eliminación de celdas sin especies
abund <- abund[
  rowSums(abund, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

cat(
  "Matriz final:",
  nrow(abund), "celdas x",
  ncol(abund), "especies.\n"
)

# Determinación de los tipos de trayectorias según su endpoint
Endpoints <- list(
  Bosque = c(
    "Bosque permanente",
    "Recuperación forestal",
    "Sucesión de matorral a bosque",
    "Sucesión de espacio abierto a bosque",
    "Abandono agrícola (sucesión a bosque)"
  ),
  Matorral = c(
    "Matorral permanente",
    "Degradación del bosque (a matorral)",
    "Abandono agrícola (sucesión a matorral)"
  ),
  Agricola = c(
    "Agrícola permanente",
    "Expansión agrícola (origen matorral)",
    "Expansión agrícola (origen bosque)"
  )
)

# Número de celdas por trayectoria
trayectorias_n <- Trayectoria %>%
  st_drop_geometry() %>%
  count(Tipo)

# Filtrar trayectorias con ≥ 20 celdas manteniendo el orden original
Endpoints <- lapply(
  Endpoints,
  function(x) {
    x[x %in% trayectorias_n$Tipo[trayectorias_n$n >= 20]]
  }
)

# Ver las trayectorias finales de cada endpoint:
trayectorias_n %>% filter(Tipo %in% Endpoints$Bosque)
trayectorias_n %>% filter(Tipo %in% Endpoints$Matorral)
trayectorias_n %>% filter(Tipo %in% Endpoints$Agricola)



# 5. Diversidad taxonómica ------------------------------------------------------

# Riqueza (Hill q0)
R0 <- specnumber(abund)

# Diversidad de Shannon (Hill q1)
Q1 <- exp(
  diversity(
    abund,
    index = "shannon"
  )
)

# Diversidad de Simpson (Hill q2)
Q2 <- diversity(
  abund,
  index = "invsimpson"
)

# Coordenadas del centroide de cada celda
Coordenadas <- grid_ebird %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(
    ID = grid_ebird$ID
  )

# Base de datos para los análisis de diversidad taxonómica
Diversidad_taxonomica <- tibble(
  ID = rownames(abund),
  R0 = R0,
  Q1 = Q1,
  Q2 = Q2
) %>%
  left_join(
    Trayectoria %>%
      st_drop_geometry() %>%
      select(ID, Tipo),
    by = "ID"
  ) %>%
  left_join(
    Estabilidad %>%
      st_drop_geometry() %>%
      select(ID, anio_estable),
    by = "ID"
  ) %>%
  left_join(
    Intensidad_cambio %>%
      st_drop_geometry() %>%
      select(ID, Porcentaje_cambio),
    by = "ID"
  ) %>%
  left_join(
    Coordenadas,
    by = "ID"
  ) %>%
  rename(
    x = X,
    y = Y
  )

# Bases de datos por endpoint
Bosque <- Diversidad_taxonomica %>%
  filter(
    Tipo %in% Endpoints$Bosque
  ) %>%
  mutate(
    Tipo = factor(
      Tipo,
      levels = Endpoints$Bosque
    )
  )

Matorral <- Diversidad_taxonomica %>%
  filter(
    Tipo %in% Endpoints$Matorral
  ) %>%
  mutate(
    Tipo = factor(
      Tipo,
      levels = Endpoints$Matorral
    )
  )

Agricola <- Diversidad_taxonomica %>%
  filter(
    Tipo %in% Endpoints$Agricola
  ) %>%
  mutate(
    Tipo = factor(
      Tipo,
      levels = Endpoints$Agricola
    )
  )



# 6. Diversidad funcional basada en rasgos ------------------------------------

# Selección de rasgos funcionales
rasgos_fd <- data.frame(
  
  log_masa      = log10(rasgos$Masa),
  HWI           = rasgos$HWI,
  long_tarso    = rasgos$Longitud_tarso,
  long_ala      = rasgos$Longitud_ala,
  nivel_trofico = factor(rasgos$Nivel_trofico),
  nicho_trofico = factor(rasgos$Nicho_trofico),
  estilo_vida   = factor(rasgos$Estilo_vida),
  migracion     = factor(rasgos$Migracion),
  
  row.names = rownames(rasgos)
)

# Comprobación crítica: mismas especies y mismo orden que la matriz de abundancia
stopifnot(
  identical(
    rownames(rasgos_fd),
    colnames(abund)
  )
)

cat(
  nrow(rasgos_fd),
  "especies con rasgos funcionales correctamente alineadas.\n"
)

# Distancia funcional de Gower
d_gower <- gowdis(
  rasgos_fd
)

# PCoA corregido del espacio funcional
pco_fd <- cmdscale(
  d_gower,
  k = 4,
  eig = TRUE,
  add = TRUE
)

coords_fd <- pco_fd$points

colnames(coords_fd) <- paste0(
  "PC",
  1:ncol(coords_fd)
)

# Calidad del espacio funcional
eig_pos <- pco_fd$eig[
  pco_fd$eig > 0
]

calidad_fd <- sum(
  eig_pos[1:4]
) /
  sum(eig_pos)

cat(
  "Calidad del espacio funcional:",
  round(calidad_fd, 3),
  "\n"
)

# Cálculo de índices funcionales por celda
FRic <- fd_fric(
  coords_fd,
  abund
)

FDis <- fd_fdis(
  coords_fd,
  abund
)

# Unión de índices funcionales
fd_tab <- reduce(
  list(
    FRic,
    FDis
  ),
  full_join,
  by = "site"
) %>%
  rename(
    ID = site
  )

# Base de datos de diversidad funcional
Diversidad_funcional <- tibble(
  ID = rownames(abund),
  FRic = fd_tab$FRic,
  FDis = fd_tab$FDis
) %>%
  left_join(
    Trayectoria %>%
      st_drop_geometry() %>%
      select(
        ID,
        Tipo
      ),
    by = "ID"
  ) %>%
  left_join(
    Estabilidad %>%
      st_drop_geometry() %>%
      select(
        ID,
        anio_estable
      ),
    by = "ID"
  ) %>%
  left_join(
    Intensidad_cambio %>%
      st_drop_geometry() %>%
      select(
        ID,
        Porcentaje_cambio
      ),
    by = "ID"
  ) %>%
  left_join(
    Coordenadas,
    by = "ID"
  ) %>%
  rename(
    x = X,
    y = Y
  )

# Comprobación final
glimpse(
  Diversidad_funcional
)

summary(
  Diversidad_funcional %>%
    select(
      FRic,
      FDis,
      anio_estable,
      Porcentaje_cambio
    )
)



# 7. Pruebas de diversidad taxonómica para endpoint -> bosque ------------------

# Distribución de las variables respuesta según trayectoria
Bosque_long <- Bosque %>%
  select(
    ID,
    Tipo,
    R0,
    Q1,
    Q2
  ) %>%
  pivot_longer(
    cols = c(R0, Q1, Q2),
    names_to = "Indice",
    values_to = "Valor"
  )

ggplot(
  Bosque_long,
  aes(
    Tipo,
    Valor
  )
) +
  geom_boxplot(
    outlier.alpha = 0.2
  ) +
  facet_wrap(
    ~Indice,
    scales = "free_y"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# RIQUEZA (R0) =================================================================

# Modelos individuales
BOS_T_R0 <- bam(
  R0 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

BOS_E_R0 <- bam(
  R0 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

BOS_I_R0 <- bam(
  R0 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

# Modelos combinados
BOS_TE_R0 <- bam(
  R0 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

BOS_TI_R0 <- bam(
  R0 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

BOS_EI_R0 <- bam(
  R0 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

BOS_TEI_R0 <- bam(
  R0 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

# Comparación de modelos
Resumen_R0_Bosque <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(BOS_T_R0)$r.sq,
    summary(BOS_E_R0)$r.sq,
    summary(BOS_I_R0)$r.sq,
    summary(BOS_TE_R0)$r.sq,
    summary(BOS_TI_R0)$r.sq,
    summary(BOS_EI_R0)$r.sq,
    summary(BOS_TEI_R0)$r.sq
  ),
  Dev_exp = c(
    summary(BOS_T_R0)$dev.expl,
    summary(BOS_E_R0)$dev.expl,
    summary(BOS_I_R0)$dev.expl,
    summary(BOS_TE_R0)$dev.expl,
    summary(BOS_TI_R0)$dev.expl,
    summary(BOS_EI_R0)$dev.expl,
    summary(BOS_TEI_R0)$dev.expl
  ),
  AIC = c(
    AIC(BOS_T_R0),
    AIC(BOS_E_R0),
    AIC(BOS_I_R0),
    AIC(BOS_TE_R0),
    AIC(BOS_TI_R0),
    AIC(BOS_EI_R0),
    AIC(BOS_TEI_R0)
  )
)

Resumen_R0_Bosque

# Concurvity 
concurvity(BOS_T_R0, full = TRUE)
concurvity(BOS_E_R0, full = TRUE)
concurvity(BOS_I_R0, full = TRUE)
concurvity(BOS_TE_R0, full = TRUE)
concurvity(BOS_TI_R0, full = TRUE)
concurvity(BOS_EI_R0, full = TRUE)
concurvity(BOS_TEI_R0, full = TRUE)

# Comparaciones entre trayectorias
EMM_R0_Bosque <- emmeans(BOS_T_R0, "Tipo")
pairs(EMM_R0_Bosque, adjust = "tukey")


# SHANNON (Q1) =================================================================

# Modelos individuales
BOS_T_Q1 <- bam(
  Q1 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_E_Q1 <- bam(
  Q1 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_I_Q1 <- bam(
  Q1 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
BOS_TE_Q1 <- bam(
  Q1 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_TI_Q1 <- bam(
  Q1 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_EI_Q1 <- bam(
  Q1 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_TEI_Q1 <- bam(
  Q1 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos
Resumen_Q1_Bosque <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(BOS_T_Q1)$r.sq,
    summary(BOS_E_Q1)$r.sq,
    summary(BOS_I_Q1)$r.sq,
    summary(BOS_TE_Q1)$r.sq,
    summary(BOS_TI_Q1)$r.sq,
    summary(BOS_EI_Q1)$r.sq,
    summary(BOS_TEI_Q1)$r.sq
  ),
  Dev_exp = c(
    summary(BOS_T_Q1)$dev.expl,
    summary(BOS_E_Q1)$dev.expl,
    summary(BOS_I_Q1)$dev.expl,
    summary(BOS_TE_Q1)$dev.expl,
    summary(BOS_TI_Q1)$dev.expl,
    summary(BOS_EI_Q1)$dev.expl,
    summary(BOS_TEI_Q1)$dev.expl
  ),
  AIC = c(
    AIC(BOS_T_Q1),
    AIC(BOS_E_Q1),
    AIC(BOS_I_Q1),
    AIC(BOS_TE_Q1),
    AIC(BOS_TI_Q1),
    AIC(BOS_EI_Q1),
    AIC(BOS_TEI_Q1)
  )
)

Resumen_Q1_Bosque

# Concurvity 
concurvity(BOS_T_Q1, full = TRUE)
concurvity(BOS_E_Q1, full = TRUE)
concurvity(BOS_I_Q1, full = TRUE)
concurvity(BOS_TE_Q1, full = TRUE)
concurvity(BOS_TI_Q1, full = TRUE)
concurvity(BOS_EI_Q1, full = TRUE)
concurvity(BOS_TEI_Q1, full = TRUE)

# Comparaciones entre trayectorias
EMM_Q1_Bosque <- emmeans(BOS_T_Q1, "Tipo")
pairs(EMM_Q1_Bosque, adjust = "tukey")


# SIMPSON (Q2) =================================================================

# Modelos individuales
BOS_T_Q2 <- bam(
  Q2 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_E_Q2 <- bam(
  Q2 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_I_Q2 <- bam(
  Q2 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
BOS_TE_Q2 <- bam(
  Q2 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_TI_Q2 <- bam(
  Q2 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_EI_Q2 <- bam(
  Q2 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

BOS_TEI_Q2 <- bam(
  Q2 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos
Resumen_Q2_Bosque <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(BOS_T_Q2)$r.sq,
    summary(BOS_E_Q2)$r.sq,
    summary(BOS_I_Q2)$r.sq,
    summary(BOS_TE_Q2)$r.sq,
    summary(BOS_TI_Q2)$r.sq,
    summary(BOS_EI_Q2)$r.sq,
    summary(BOS_TEI_Q2)$r.sq
  ),
  Dev_exp = c(
    summary(BOS_T_Q2)$dev.expl,
    summary(BOS_E_Q2)$dev.expl,
    summary(BOS_I_Q2)$dev.expl,
    summary(BOS_TE_Q2)$dev.expl,
    summary(BOS_TI_Q2)$dev.expl,
    summary(BOS_EI_Q2)$dev.expl,
    summary(BOS_TEI_Q2)$dev.expl
  ),
  AIC = c(
    AIC(BOS_T_Q2),
    AIC(BOS_E_Q2),
    AIC(BOS_I_Q2),
    AIC(BOS_TE_Q2),
    AIC(BOS_TI_Q2),
    AIC(BOS_EI_Q2),
    AIC(BOS_TEI_Q2)
  )
)

Resumen_Q2_Bosque

# Concurvity 
concurvity(BOS_T_Q2, full = TRUE)
concurvity(BOS_E_Q2, full = TRUE)
concurvity(BOS_I_Q2, full = TRUE)
concurvity(BOS_TE_Q2, full = TRUE)
concurvity(BOS_TI_Q2, full = TRUE)
concurvity(BOS_EI_Q2, full = TRUE)
concurvity(BOS_TEI_Q2, full = TRUE)

# Comparaciones entre trayectorias
EMM_Q2_Bosque <- emmeans(BOS_T_Q2, "Tipo")
pairs(EMM_Q2_Bosque, adjust = "tukey")



# 8. Diversidad funcional: endpoint bosque -------------------------------------

# Base de datos del endpoint bosque
Bosque_FD <- Diversidad_funcional %>%
  filter(
    Tipo %in% Endpoints$Bosque
  ) %>%
  mutate(
    Tipo = factor(
      Tipo,
      levels = Endpoints$Bosque
    )
  )

# Comprobación de trayectorias incluidas
Bosque_FD %>%
  count(Tipo)


# FRic ========================================================================

# Modelos individuales
GAM_T_FRic_Bosque <- bam(
  FRic ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_E_FRic_Bosque <- bam(
  FRic ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_I_FRic_Bosque <- bam(
  FRic ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
GAM_TE_FRic_Bosque <- bam(
  FRic ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TI_FRic_Bosque <- bam(
  FRic ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_EI_FRic_Bosque <- bam(
  FRic ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TEI_FRic_Bosque <- bam(
  FRic ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos FRic
Resumen_FRic_Bosque <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(GAM_T_FRic_Bosque)$r.sq,
    summary(GAM_E_FRic_Bosque)$r.sq,
    summary(GAM_I_FRic_Bosque)$r.sq,
    summary(GAM_TE_FRic_Bosque)$r.sq,
    summary(GAM_TI_FRic_Bosque)$r.sq,
    summary(GAM_EI_FRic_Bosque)$r.sq,
    summary(GAM_TEI_FRic_Bosque)$r.sq
  ),
  Dev_exp = c(
    summary(GAM_T_FRic_Bosque)$dev.expl,
    summary(GAM_E_FRic_Bosque)$dev.expl,
    summary(GAM_I_FRic_Bosque)$dev.expl,
    summary(GAM_TE_FRic_Bosque)$dev.expl,
    summary(GAM_TI_FRic_Bosque)$dev.expl,
    summary(GAM_EI_FRic_Bosque)$dev.expl,
    summary(GAM_TEI_FRic_Bosque)$dev.expl
  ),
  AIC = c(
    AIC(GAM_T_FRic_Bosque),
    AIC(GAM_E_FRic_Bosque),
    AIC(GAM_I_FRic_Bosque),
    AIC(GAM_TE_FRic_Bosque),
    AIC(GAM_TI_FRic_Bosque),
    AIC(GAM_EI_FRic_Bosque),
    AIC(GAM_TEI_FRic_Bosque)
  )
)
Resumen_FRic_Bosque

# Concurvity FRic
concurvity(GAM_T_FRic_Bosque, full = TRUE)
concurvity(GAM_E_FRic_Bosque, full = TRUE)
concurvity(GAM_I_FRic_Bosque, full = TRUE)
concurvity(GAM_TE_FRic_Bosque, full = TRUE)
concurvity(GAM_EI_FRic_Bosque, full = TRUE)
concurvity(GAM_TI_FRic_Bosque, full = TRUE)
concurvity(GAM_TEI_FRic_Bosque, full = TRUE)

# Comparaciones entre trayectorias FRic
EMM_FRic_Bosque <- emmeans(
  GAM_T_FRic_Bosque,
  "Tipo"
)
pairs(EMM_FRic_Bosque, adjust = "tukey")


# FDis ========================================================================

# Modelos individuales
GAM_T_FDis_Bosque <- bam(
  FDis ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_E_FDis_Bosque <- bam(
  FDis ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_I_FDis_Bosque <- bam(
  FDis ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
GAM_TE_FDis_Bosque <- bam(
  FDis ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TI_FDis_Bosque <- bam(
  FDis ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_EI_FDis_Bosque <- bam(
  FDis ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TEI_FDis_Bosque <- bam(
  FDis ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Bosque_FD,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos FDis
Resumen_FDis_Bosque <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(GAM_T_FDis_Bosque)$r.sq,
    summary(GAM_E_FDis_Bosque)$r.sq,
    summary(GAM_I_FDis_Bosque)$r.sq,
    summary(GAM_TE_FDis_Bosque)$r.sq,
    summary(GAM_TI_FDis_Bosque)$r.sq,
    summary(GAM_EI_FDis_Bosque)$r.sq,
    summary(GAM_TEI_FDis_Bosque)$r.sq
  ),
  Dev_exp = c(
    summary(GAM_T_FDis_Bosque)$dev.expl,
    summary(GAM_E_FDis_Bosque)$dev.expl,
    summary(GAM_I_FDis_Bosque)$dev.expl,
    summary(GAM_TE_FDis_Bosque)$dev.expl,
    summary(GAM_TI_FDis_Bosque)$dev.expl,
    summary(GAM_EI_FDis_Bosque)$dev.expl,
    summary(GAM_TEI_FDis_Bosque)$dev.expl
  ),
  AIC = c(
    AIC(GAM_T_FDis_Bosque),
    AIC(GAM_E_FDis_Bosque),
    AIC(GAM_I_FDis_Bosque),
    AIC(GAM_TE_FDis_Bosque),
    AIC(GAM_TI_FDis_Bosque),
    AIC(GAM_EI_FDis_Bosque),
    AIC(GAM_TEI_FDis_Bosque)
  )
)
Resumen_FDis_Bosque

# Concurvity FDis
concurvity(GAM_T_FDis_Bosque, full = TRUE)
concurvity(GAM_E_FDis_Bosque, full = TRUE)
concurvity(GAM_I_FDis_Bosque, full = TRUE)
concurvity(GAM_TE_FDis_Bosque, full = TRUE)
concurvity(GAM_EI_FDis_Bosque, full = TRUE)
concurvity(GAM_TI_FDis_Bosque, full = TRUE)
concurvity(GAM_TEI_FDis_Bosque, full = TRUE)

# Comparaciones entre trayectorias FDis
EMM_FDis_Bosque <- emmeans(
  GAM_T_FDis_Bosque,
  "Tipo"
)
pairs(EMM_FDis_Bosque,adjust = "tukey")


# FIGURA — Medias marginales estimadas (EMM) ± IC95% -------------------------

# Obtención de EMM de FRic
EMM_FRic_plot <- as.data.frame(
  EMM_FRic_Bosque
) %>%
  mutate(
    Indice = "FRic"
  )

# 2. Obtención de EMM de FDis
EMM_FDis_plot <- as.data.frame(
  EMM_FDis_Bosque
) %>%
  mutate(
    Indice = "FDis"
  )

# Unión de ambos índices
EMM_FD_plot <- bind_rows(
  EMM_FRic_plot,
  EMM_FDis_plot
) %>%
  mutate(
    Indice = factor(
      Indice,
      levels = c("FRic", "FDis"),
      labels = c(
        "Riqueza funcional (FRic)",
        "Dispersión funcional (FDis)"
      )
    ),
    Tipo = factor(
      Tipo,
      levels = Endpoints$Bosque
    )
  )

# Gráfica
Figura_EMM <- ggplot(
  EMM_FD_plot,
  aes(
    x = Tipo,
    y = emmean
  )
) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(
    ~ Indice,
    scales = "free_y"
  ) +
  
  labs(
    x = "Tipo de trayectoria",
    y = "Media marginal estimada (IC95%)"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 10
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

Figura_EMM



# 9. Pruebas de diversidad taxonómica para endpoint -> matorral ----------------

# Distribución de las variables respuesta según trayectoria
Matorral_long <- Matorral %>%
  select(
    ID,
    Tipo,
    R0,
    Q1,
    Q2
  ) %>%
  pivot_longer(
    cols = c(R0, Q1, Q2),
    names_to = "Indice",
    values_to = "Valor"
  )

ggplot(
  Matorral_long,
  aes(
    Tipo,
    Valor
  )
) +
  geom_boxplot(
    outlier.alpha = 0.2
  ) +
  facet_wrap(
    ~Indice,
    scales = "free_y"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# RIQUEZA (R0) =================================================================

# Modelos individuales
MAT_T_R0 <- bam(
  R0 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

MAT_E_R0 <- bam(
  R0 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

MAT_I_R0 <- bam(
  R0 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

# Modelos combinados
MAT_TE_R0 <- bam(
  R0 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

MAT_TI_R0 <- bam(
  R0 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

MAT_EI_R0 <- bam(
  R0 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

MAT_TEI_R0 <- bam(
  R0 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

# Comparación de modelos
Resumen_R0_Matorral <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(MAT_T_R0)$r.sq,
    summary(MAT_E_R0)$r.sq,
    summary(MAT_I_R0)$r.sq,
    summary(MAT_TE_R0)$r.sq,
    summary(MAT_TI_R0)$r.sq,
    summary(MAT_EI_R0)$r.sq,
    summary(MAT_TEI_R0)$r.sq
  ),
  Dev_exp = c(
    summary(MAT_T_R0)$dev.expl,
    summary(MAT_E_R0)$dev.expl,
    summary(MAT_I_R0)$dev.expl,
    summary(MAT_TE_R0)$dev.expl,
    summary(MAT_TI_R0)$dev.expl,
    summary(MAT_EI_R0)$dev.expl,
    summary(MAT_TEI_R0)$dev.expl
  ),
  AIC = c(
    AIC(MAT_T_R0),
    AIC(MAT_E_R0),
    AIC(MAT_I_R0),
    AIC(MAT_TE_R0),
    AIC(MAT_TI_R0),
    AIC(MAT_EI_R0),
    AIC(MAT_TEI_R0)
  )
)

Resumen_R0_Matorral

# Concurvity 
concurvity(MAT_T_R0, full = TRUE)
concurvity(MAT_E_R0, full = TRUE)
concurvity(MAT_I_R0, full = TRUE)
concurvity(MAT_TE_R0, full = TRUE)
concurvity(MAT_TI_R0, full = TRUE)
concurvity(MAT_EI_R0, full = TRUE)
concurvity(MAT_TEI_R0, full = TRUE)

# Comparaciones entre trayectorias
EMM_R0_Matorral <- emmeans(MAT_T_R0, "Tipo")
pairs(EMM_R0_Matorral, adjust = "tukey")


# SHANNON (Q1) =================================================================

# Modelos individuales
MAT_T_Q1 <- bam(
  Q1 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_E_Q1 <- bam(
  Q1 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_I_Q1 <- bam(
  Q1 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
MAT_TE_Q1 <- bam(
  Q1 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_TI_Q1 <- bam(
  Q1 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_EI_Q1 <- bam(
  Q1 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_TEI_Q1 <- bam(
  Q1 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos
Resumen_Q1_Matorral <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(MAT_T_Q1)$r.sq,
    summary(MAT_E_Q1)$r.sq,
    summary(MAT_I_Q1)$r.sq,
    summary(MAT_TE_Q1)$r.sq,
    summary(MAT_TI_Q1)$r.sq,
    summary(MAT_EI_Q1)$r.sq,
    summary(MAT_TEI_Q1)$r.sq
  ),
  Dev_exp = c(
    summary(MAT_T_Q1)$dev.expl,
    summary(MAT_E_Q1)$dev.expl,
    summary(MAT_I_Q1)$dev.expl,
    summary(MAT_TE_Q1)$dev.expl,
    summary(MAT_TI_Q1)$dev.expl,
    summary(MAT_EI_Q1)$dev.expl,
    summary(MAT_TEI_Q1)$dev.expl
  ),
  AIC = c(
    AIC(MAT_T_Q1),
    AIC(MAT_E_Q1),
    AIC(MAT_I_Q1),
    AIC(MAT_TE_Q1),
    AIC(MAT_TI_Q1),
    AIC(MAT_EI_Q1),
    AIC(MAT_TEI_Q1)
  )
)

Resumen_Q1_Matorral

# Concurvity 
concurvity(MAT_T_Q1, full = TRUE)
concurvity(MAT_E_Q1, full = TRUE)
concurvity(MAT_I_Q1, full = TRUE)
concurvity(MAT_TE_Q1, full = TRUE)
concurvity(MAT_TI_Q1, full = TRUE)
concurvity(MAT_EI_Q1, full = TRUE)
concurvity(MAT_TEI_Q1, full = TRUE)

# Comparaciones entre trayectorias
EMM_Q1_Matorral <- emmeans(MAT_T_Q1, "Tipo")
pairs(EMM_Q1_Matorral, adjust = "tukey")


# SIMPSON (Q2) =================================================================

# Modelos individuales
MAT_T_Q2 <- bam(
  Q2 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_E_Q2 <- bam(
  Q2 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_I_Q2 <- bam(
  Q2 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
MAT_TE_Q2 <- bam(
  Q2 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_TI_Q2 <- bam(
  Q2 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_EI_Q2 <- bam(
  Q2 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

MAT_TEI_Q2 <- bam(
  Q2 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos
Resumen_Q2_Matorral <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(MAT_T_Q2)$r.sq,
    summary(MAT_E_Q2)$r.sq,
    summary(MAT_I_Q2)$r.sq,
    summary(MAT_TE_Q2)$r.sq,
    summary(MAT_TI_Q2)$r.sq,
    summary(MAT_EI_Q2)$r.sq,
    summary(MAT_TEI_Q2)$r.sq
  ),
  Dev_exp = c(
    summary(MAT_T_Q2)$dev.expl,
    summary(MAT_E_Q2)$dev.expl,
    summary(MAT_I_Q2)$dev.expl,
    summary(MAT_TE_Q2)$dev.expl,
    summary(MAT_TI_Q2)$dev.expl,
    summary(MAT_EI_Q2)$dev.expl,
    summary(MAT_TEI_Q2)$dev.expl
  ),
  AIC = c(
    AIC(MAT_T_Q2),
    AIC(MAT_E_Q2),
    AIC(MAT_I_Q2),
    AIC(MAT_TE_Q2),
    AIC(MAT_TI_Q2),
    AIC(MAT_EI_Q2),
    AIC(MAT_TEI_Q2)
  )
)

Resumen_Q2_Matorral

# Concurvity 
concurvity(MAT_T_Q2, full = TRUE)
concurvity(MAT_E_Q2, full = TRUE)
concurvity(MAT_I_Q2, full = TRUE)
concurvity(MAT_TE_Q2, full = TRUE)
concurvity(MAT_TI_Q2, full = TRUE)
concurvity(MAT_EI_Q2, full = TRUE)
concurvity(MAT_TEI_Q2, full = TRUE)

# Comparaciones entre trayectorias
EMM_Q2_Matorral <- emmeans(MAT_T_Q2, "Tipo")
pairs(EMM_Q2_Matorral, adjust = "tukey")



# 10. Diversidad funcional: endpoint Matorral -------------------------------------

# Base de datos del endpoint Matorral
Matorral_FD <- Diversidad_funcional %>%
  filter(
    Tipo %in% Endpoints$Matorral
  ) %>%
  mutate(
    Tipo = factor(
      Tipo,
      levels = Endpoints$Matorral
    )
  )

# Comprobación de trayectorias incluidas
Matorral_FD %>%
  count(Tipo)


# FRic ========================================================================

# Modelos individuales
GAM_T_FRic_Matorral <- bam(
  FRic ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_E_FRic_Matorral <- bam(
  FRic ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_I_FRic_Matorral <- bam(
  FRic ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
GAM_TE_FRic_Matorral <- bam(
  FRic ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TI_FRic_Matorral <- bam(
  FRic ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_EI_FRic_Matorral <- bam(
  FRic ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TEI_FRic_Matorral <- bam(
  FRic ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos FRic
Resumen_FRic_Matorral <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(GAM_T_FRic_Matorral)$r.sq,
    summary(GAM_E_FRic_Matorral)$r.sq,
    summary(GAM_I_FRic_Matorral)$r.sq,
    summary(GAM_TE_FRic_Matorral)$r.sq,
    summary(GAM_TI_FRic_Matorral)$r.sq,
    summary(GAM_EI_FRic_Matorral)$r.sq,
    summary(GAM_TEI_FRic_Matorral)$r.sq
  ),
  Dev_exp = c(
    summary(GAM_T_FRic_Matorral)$dev.expl,
    summary(GAM_E_FRic_Matorral)$dev.expl,
    summary(GAM_I_FRic_Matorral)$dev.expl,
    summary(GAM_TE_FRic_Matorral)$dev.expl,
    summary(GAM_TI_FRic_Matorral)$dev.expl,
    summary(GAM_EI_FRic_Matorral)$dev.expl,
    summary(GAM_TEI_FRic_Matorral)$dev.expl
  ),
  AIC = c(
    AIC(GAM_T_FRic_Matorral),
    AIC(GAM_E_FRic_Matorral),
    AIC(GAM_I_FRic_Matorral),
    AIC(GAM_TE_FRic_Matorral),
    AIC(GAM_TI_FRic_Matorral),
    AIC(GAM_EI_FRic_Matorral),
    AIC(GAM_TEI_FRic_Matorral)
  )
)
Resumen_FRic_Matorral

# Concurvity FRic
concurvity(GAM_T_FRic_Matorral, full = TRUE)
concurvity(GAM_E_FRic_Matorral, full = TRUE)
concurvity(GAM_I_FRic_Matorral, full = TRUE)
concurvity(GAM_TE_FRic_Matorral, full = TRUE)
concurvity(GAM_EI_FRic_Matorral, full = TRUE)
concurvity(GAM_TI_FRic_Matorral, full = TRUE)
concurvity(GAM_TEI_FRic_Matorral, full = TRUE)

# Comparaciones entre trayectorias FRic
EMM_FRic_Matorral <- emmeans(
  GAM_T_FRic_Matorral,
  "Tipo"
)
pairs(EMM_FRic_Matorral, adjust = "tukey")


# FDis ========================================================================

# Modelos individuales
GAM_T_FDis_Matorral <- bam(
  FDis ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_E_FDis_Matorral <- bam(
  FDis ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_I_FDis_Matorral <- bam(
  FDis ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
GAM_TE_FDis_Matorral <- bam(
  FDis ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TI_FDis_Matorral <- bam(
  FDis ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_EI_FDis_Matorral <- bam(
  FDis ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TEI_FDis_Matorral <- bam(
  FDis ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Matorral_FD,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos FDis
Resumen_FDis_Matorral <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(GAM_T_FDis_Matorral)$r.sq,
    summary(GAM_E_FDis_Matorral)$r.sq,
    summary(GAM_I_FDis_Matorral)$r.sq,
    summary(GAM_TE_FDis_Matorral)$r.sq,
    summary(GAM_TI_FDis_Matorral)$r.sq,
    summary(GAM_EI_FDis_Matorral)$r.sq,
    summary(GAM_TEI_FDis_Matorral)$r.sq
  ),
  Dev_exp = c(
    summary(GAM_T_FDis_Matorral)$dev.expl,
    summary(GAM_E_FDis_Matorral)$dev.expl,
    summary(GAM_I_FDis_Matorral)$dev.expl,
    summary(GAM_TE_FDis_Matorral)$dev.expl,
    summary(GAM_TI_FDis_Matorral)$dev.expl,
    summary(GAM_EI_FDis_Matorral)$dev.expl,
    summary(GAM_TEI_FDis_Matorral)$dev.expl
  ),
  AIC = c(
    AIC(GAM_T_FDis_Matorral),
    AIC(GAM_E_FDis_Matorral),
    AIC(GAM_I_FDis_Matorral),
    AIC(GAM_TE_FDis_Matorral),
    AIC(GAM_TI_FDis_Matorral),
    AIC(GAM_EI_FDis_Matorral),
    AIC(GAM_TEI_FDis_Matorral)
  )
)
Resumen_FDis_Matorral

# Concurvity FDis
concurvity(GAM_T_FDis_Matorral, full = TRUE)
concurvity(GAM_E_FDis_Matorral, full = TRUE)
concurvity(GAM_I_FDis_Matorral, full = TRUE)
concurvity(GAM_TE_FDis_Matorral, full = TRUE)
concurvity(GAM_EI_FDis_Matorral, full = TRUE)
concurvity(GAM_TI_FDis_Matorral, full = TRUE)
concurvity(GAM_TEI_FDis_Matorral, full = TRUE)

# Comparaciones entre trayectorias FDis
EMM_FDis_Matorral <- emmeans(
  GAM_T_FDis_Matorral,
  "Tipo"
)
pairs(EMM_FDis_Matorral,adjust = "tukey")


# FIGURA — Medias marginales estimadas (EMM) ± IC95% -------------------------

# Obtención de EMM de FRic
EMM_FRic_plot <- as.data.frame(
  EMM_FRic_Matorral
) %>%
  mutate(
    Indice = "FRic"
  )

# 2. Obtención de EMM de FDis
EMM_FDis_plot <- as.data.frame(
  EMM_FDis_Matorral
) %>%
  mutate(
    Indice = "FDis"
  )

# Unión de ambos índices
EMM_FD_plot <- bind_rows(
  EMM_FRic_plot,
  EMM_FDis_plot
) %>%
  mutate(
    Indice = factor(
      Indice,
      levels = c("FRic", "FDis"),
      labels = c(
        "Riqueza funcional (FRic)",
        "Dispersión funcional (FDis)"
      )
    ),
    Tipo = factor(
      Tipo,
      levels = Endpoints$Matorral
    )
  )

# Gráfica
Figura_EMM <- ggplot(
  EMM_FD_plot,
  aes(
    x = Tipo,
    y = emmean
  )
) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(
    ~ Indice,
    scales = "free_y"
  ) +
  
  labs(
    x = "Tipo de trayectoria",
    y = "Media marginal estimada (IC95%)"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 10
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

Figura_EMM



# 11. Pruebas de diversidad taxonómica para endpoint -> agricola ----------------

# Distribución de las variables respuesta según trayectoria
Agricola_long <- Agricola %>%
  select(
    ID,
    Tipo,
    R0,
    Q1,
    Q2
  ) %>%
  pivot_longer(
    cols = c(R0, Q1, Q2),
    names_to = "Indice",
    values_to = "Valor"
  )

ggplot(
  Agricola_long,
  aes(
    Tipo,
    Valor
  )
) +
  geom_boxplot(
    outlier.alpha = 0.2
  ) +
  facet_wrap(
    ~Indice,
    scales = "free_y"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# RIQUEZA (R0) =================================================================

# Modelos individuales
AGR_T_R0 <- bam(
  R0 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

AGR_E_R0 <- bam(
  R0 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

AGR_I_R0 <- bam(
  R0 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

# Modelos combinados
AGR_TE_R0 <- bam(
  R0 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

AGR_TI_R0 <- bam(
  R0 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

AGR_EI_R0 <- bam(
  R0 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

AGR_TEI_R0 <- bam(
  R0 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  family = nb(),
  discrete = TRUE
)

# Comparación de modelos
Resumen_R0_Agricola <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(AGR_T_R0)$r.sq,
    summary(AGR_E_R0)$r.sq,
    summary(AGR_I_R0)$r.sq,
    summary(AGR_TE_R0)$r.sq,
    summary(AGR_TI_R0)$r.sq,
    summary(AGR_EI_R0)$r.sq,
    summary(AGR_TEI_R0)$r.sq
  ),
  Dev_exp = c(
    summary(AGR_T_R0)$dev.expl,
    summary(AGR_E_R0)$dev.expl,
    summary(AGR_I_R0)$dev.expl,
    summary(AGR_TE_R0)$dev.expl,
    summary(AGR_TI_R0)$dev.expl,
    summary(AGR_EI_R0)$dev.expl,
    summary(AGR_TEI_R0)$dev.expl
  ),
  AIC = c(
    AIC(AGR_T_R0),
    AIC(AGR_E_R0),
    AIC(AGR_I_R0),
    AIC(AGR_TE_R0),
    AIC(AGR_TI_R0),
    AIC(AGR_EI_R0),
    AIC(AGR_TEI_R0)
  )
)

Resumen_R0_Agricola

# Concurvity 
concurvity(AGR_T_R0, full = TRUE)
concurvity(AGR_E_R0, full = TRUE)
concurvity(AGR_I_R0, full = TRUE)
concurvity(AGR_TE_R0, full = TRUE)
concurvity(AGR_TI_R0, full = TRUE)
concurvity(AGR_EI_R0, full = TRUE)
concurvity(AGR_TEI_R0, full = TRUE)

# Comparaciones entre trayectorias
EMM_R0_Agricola <- emmeans(AGR_T_R0, "Tipo")
pairs(EMM_R0_Agricola, adjust = "tukey")


# SHANNON (Q1) =================================================================

# Modelos individuales
AGR_T_Q1 <- bam(
  Q1 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_E_Q1 <- bam(
  Q1 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_I_Q1 <- bam(
  Q1 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
AGR_TE_Q1 <- bam(
  Q1 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_TI_Q1 <- bam(
  Q1 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_EI_Q1 <- bam(
  Q1 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_TEI_Q1 <- bam(
  Q1 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos
Resumen_Q1_Agricola <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(AGR_T_Q1)$r.sq,
    summary(AGR_E_Q1)$r.sq,
    summary(AGR_I_Q1)$r.sq,
    summary(AGR_TE_Q1)$r.sq,
    summary(AGR_TI_Q1)$r.sq,
    summary(AGR_EI_Q1)$r.sq,
    summary(AGR_TEI_Q1)$r.sq
  ),
  Dev_exp = c(
    summary(AGR_T_Q1)$dev.expl,
    summary(AGR_E_Q1)$dev.expl,
    summary(AGR_I_Q1)$dev.expl,
    summary(AGR_TE_Q1)$dev.expl,
    summary(AGR_TI_Q1)$dev.expl,
    summary(AGR_EI_Q1)$dev.expl,
    summary(AGR_TEI_Q1)$dev.expl
  ),
  AIC = c(
    AIC(AGR_T_Q1),
    AIC(AGR_E_Q1),
    AIC(AGR_I_Q1),
    AIC(AGR_TE_Q1),
    AIC(AGR_TI_Q1),
    AIC(AGR_EI_Q1),
    AIC(AGR_TEI_Q1)
  )
)

Resumen_Q1_Agricola

# Concurvity 
concurvity(AGR_T_Q1, full = TRUE)
concurvity(AGR_E_Q1, full = TRUE)
concurvity(AGR_I_Q1, full = TRUE)
concurvity(AGR_TE_Q1, full = TRUE)
concurvity(AGR_TI_Q1, full = TRUE)
concurvity(AGR_EI_Q1, full = TRUE)
concurvity(AGR_TEI_Q1, full = TRUE)

# Comparaciones entre trayectorias
EMM_Q1_Agricola <- emmeans(AGR_T_Q1, "Tipo")
pairs(EMM_Q1_Agricola, adjust = "tukey")


# SIMPSON (Q2) =================================================================

# Modelos individuales
AGR_T_Q2 <- bam(
  Q2 ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_E_Q2 <- bam(
  Q2 ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_I_Q2 <- bam(
  Q2 ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
AGR_TE_Q2 <- bam(
  Q2 ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_TI_Q2 <- bam(
  Q2 ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_EI_Q2 <- bam(
  Q2 ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

AGR_TEI_Q2 <- bam(
  Q2 ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos
Resumen_Q2_Agricola <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(AGR_T_Q2)$r.sq,
    summary(AGR_E_Q2)$r.sq,
    summary(AGR_I_Q2)$r.sq,
    summary(AGR_TE_Q2)$r.sq,
    summary(AGR_TI_Q2)$r.sq,
    summary(AGR_EI_Q2)$r.sq,
    summary(AGR_TEI_Q2)$r.sq
  ),
  Dev_exp = c(
    summary(AGR_T_Q2)$dev.expl,
    summary(AGR_E_Q2)$dev.expl,
    summary(AGR_I_Q2)$dev.expl,
    summary(AGR_TE_Q2)$dev.expl,
    summary(AGR_TI_Q2)$dev.expl,
    summary(AGR_EI_Q2)$dev.expl,
    summary(AGR_TEI_Q2)$dev.expl
  ),
  AIC = c(
    AIC(AGR_T_Q2),
    AIC(AGR_E_Q2),
    AIC(AGR_I_Q2),
    AIC(AGR_TE_Q2),
    AIC(AGR_TI_Q2),
    AIC(AGR_EI_Q2),
    AIC(AGR_TEI_Q2)
  )
)

Resumen_Q2_Agricola

# Concurvity 
concurvity(AGR_T_Q2, full = TRUE)
concurvity(AGR_E_Q2, full = TRUE)
concurvity(AGR_I_Q2, full = TRUE)
concurvity(AGR_TE_Q2, full = TRUE)
concurvity(AGR_TI_Q2, full = TRUE)
concurvity(AGR_EI_Q2, full = TRUE)
concurvity(AGR_TEI_Q2, full = TRUE)

# Comparaciones entre trayectorias
EMM_Q2_Agricola <- emmeans(AGR_T_Q2, "Tipo")
pairs(EMM_Q2_Agricola, adjust = "tukey")



# 12. Diversidad funcional: endpoint Agricola -------------------------------------

# Base de datos del endpoint Agricola
Agricola_FD <- Diversidad_funcional %>%
  filter(
    Tipo %in% Endpoints$Agricola
  ) %>%
  mutate(
    Tipo = factor(
      Tipo,
      levels = Endpoints$Agricola
    )
  )

# Comprobación de trayectorias incluidas
Agricola_FD %>%
  count(Tipo)


# FRic ========================================================================

# Modelos individuales
GAM_T_FRic_Agricola <- bam(
  FRic ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_E_FRic_Agricola <- bam(
  FRic ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_I_FRic_Agricola <- bam(
  FRic ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
GAM_TE_FRic_Agricola <- bam(
  FRic ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TI_FRic_Agricola <- bam(
  FRic ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_EI_FRic_Agricola <- bam(
  FRic ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TEI_FRic_Agricola <- bam(
  FRic ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos FRic
Resumen_FRic_Agricola <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(GAM_T_FRic_Agricola)$r.sq,
    summary(GAM_E_FRic_Agricola)$r.sq,
    summary(GAM_I_FRic_Agricola)$r.sq,
    summary(GAM_TE_FRic_Agricola)$r.sq,
    summary(GAM_TI_FRic_Agricola)$r.sq,
    summary(GAM_EI_FRic_Agricola)$r.sq,
    summary(GAM_TEI_FRic_Agricola)$r.sq
  ),
  Dev_exp = c(
    summary(GAM_T_FRic_Agricola)$dev.expl,
    summary(GAM_E_FRic_Agricola)$dev.expl,
    summary(GAM_I_FRic_Agricola)$dev.expl,
    summary(GAM_TE_FRic_Agricola)$dev.expl,
    summary(GAM_TI_FRic_Agricola)$dev.expl,
    summary(GAM_EI_FRic_Agricola)$dev.expl,
    summary(GAM_TEI_FRic_Agricola)$dev.expl
  ),
  AIC = c(
    AIC(GAM_T_FRic_Agricola),
    AIC(GAM_E_FRic_Agricola),
    AIC(GAM_I_FRic_Agricola),
    AIC(GAM_TE_FRic_Agricola),
    AIC(GAM_TI_FRic_Agricola),
    AIC(GAM_EI_FRic_Agricola),
    AIC(GAM_TEI_FRic_Agricola)
  )
)
Resumen_FRic_Agricola

# Concurvity FRic
concurvity(GAM_T_FRic_Agricola, full = TRUE)
concurvity(GAM_E_FRic_Agricola, full = TRUE)
concurvity(GAM_I_FRic_Agricola, full = TRUE)
concurvity(GAM_TE_FRic_Agricola, full = TRUE)
concurvity(GAM_EI_FRic_Agricola, full = TRUE)
concurvity(GAM_TI_FRic_Agricola, full = TRUE)
concurvity(GAM_TEI_FRic_Agricola, full = TRUE)

# Comparaciones entre trayectorias FRic
EMM_FRic_Agricola <- emmeans(
  GAM_T_FRic_Agricola,
  "Tipo"
)
pairs(EMM_FRic_Agricola, adjust = "tukey")


# FDis ========================================================================

# Modelos individuales
GAM_T_FDis_Agricola <- bam(
  FDis ~
    Tipo +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_E_FDis_Agricola <- bam(
  FDis ~
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_I_FDis_Agricola <- bam(
  FDis ~
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

# Modelos combinados
GAM_TE_FDis_Agricola <- bam(
  FDis ~
    Tipo +
    anio_estable +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TI_FDis_Agricola <- bam(
  FDis ~
    Tipo +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_EI_FDis_Agricola <- bam(
  FDis ~
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

GAM_TEI_FDis_Agricola <- bam(
  FDis ~
    Tipo +
    anio_estable +
    Porcentaje_cambio +
    s(x, y, bs = "tp", k = 50),
  data = Agricola_FD,
  method = "fREML",
  discrete = TRUE
)

# Comparación de modelos FDis
Resumen_FDis_Agricola <- tibble(
  Modelo = c(
    "T",
    "E",
    "I",
    "TE",
    "TI",
    "EI",
    "TEI"
  ),
  R2 = c(
    summary(GAM_T_FDis_Agricola)$r.sq,
    summary(GAM_E_FDis_Agricola)$r.sq,
    summary(GAM_I_FDis_Agricola)$r.sq,
    summary(GAM_TE_FDis_Agricola)$r.sq,
    summary(GAM_TI_FDis_Agricola)$r.sq,
    summary(GAM_EI_FDis_Agricola)$r.sq,
    summary(GAM_TEI_FDis_Agricola)$r.sq
  ),
  Dev_exp = c(
    summary(GAM_T_FDis_Agricola)$dev.expl,
    summary(GAM_E_FDis_Agricola)$dev.expl,
    summary(GAM_I_FDis_Agricola)$dev.expl,
    summary(GAM_TE_FDis_Agricola)$dev.expl,
    summary(GAM_TI_FDis_Agricola)$dev.expl,
    summary(GAM_EI_FDis_Agricola)$dev.expl,
    summary(GAM_TEI_FDis_Agricola)$dev.expl
  ),
  AIC = c(
    AIC(GAM_T_FDis_Agricola),
    AIC(GAM_E_FDis_Agricola),
    AIC(GAM_I_FDis_Agricola),
    AIC(GAM_TE_FDis_Agricola),
    AIC(GAM_TI_FDis_Agricola),
    AIC(GAM_EI_FDis_Agricola),
    AIC(GAM_TEI_FDis_Agricola)
  )
)
Resumen_FDis_Agricola

# Concurvity FDis
concurvity(GAM_T_FDis_Agricola, full = TRUE)
concurvity(GAM_E_FDis_Agricola, full = TRUE)
concurvity(GAM_I_FDis_Agricola, full = TRUE)
concurvity(GAM_TE_FDis_Agricola, full = TRUE)
concurvity(GAM_EI_FDis_Agricola, full = TRUE)
concurvity(GAM_TI_FDis_Agricola, full = TRUE)
concurvity(GAM_TEI_FDis_Agricola, full = TRUE)

# Comparaciones entre trayectorias FDis
EMM_FDis_Agricola <- emmeans(
  GAM_T_FDis_Agricola,
  "Tipo"
)
pairs(EMM_FDis_Agricola,adjust = "tukey")


# FIGURA — Medias marginales estimadas (EMM) ± IC95% -------------------------

# Obtención de EMM de FRic
EMM_FRic_plot <- as.data.frame(
  EMM_FRic_Agricola
) %>%
  mutate(
    Indice = "FRic"
  )

# 2. Obtención de EMM de FDis
EMM_FDis_plot <- as.data.frame(
  EMM_FDis_Agricola
) %>%
  mutate(
    Indice = "FDis"
  )

# Unión de ambos índices
EMM_FD_plot <- bind_rows(
  EMM_FRic_plot,
  EMM_FDis_plot
) %>%
  mutate(
    Indice = factor(
      Indice,
      levels = c("FRic", "FDis"),
      labels = c(
        "Riqueza funcional (FRic)",
        "Dispersión funcional (FDis)"
      )
    ),
    Tipo = factor(
      Tipo,
      levels = Endpoints$Agricola
    )
  )

# Gráfica
Figura_EMM <- ggplot(
  EMM_FD_plot,
  aes(
    x = Tipo,
    y = emmean
  )
) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(
    ~ Indice,
    scales = "free_y"
  ) +
  
  labs(
    x = "Tipo de trayectoria",
    y = "Media marginal estimada (IC95%)"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 10
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

Figura_EMM