# ==============================================================================
# Creación de los dataframes necesarios para el análisis
# ==============================================================================



# 0. Antes que nada ------------------------------------------------------------

cat("\014")    # Limpiar la consola
rm(list=ls())
# setwd(dirección/de/tu/directorio)
getwd()



# 1. Paquetes ------------------------------------------------------------------

library(tidyverse)     # Lectura, manipulación, ggplot2
library(vegan)         # Diversidad taxonómica (Hill), Bray-Curtis, PERMANOVA
library(FD)            # Gowdis() (distancia de Gower) y functcomp() (CWM)
library(fundiversity)  # Índices funcionales rápidos (FRic, FEve, FDis, Rao Q)
library(patchwork)     # Combinar figuras
library(sf)            # Manejo de datos espaciales vectoriales
library(terra)         # Análisis de datos espaciales raster y vectoriales
library(tmap)          # Creación y visualización de mapas temáticos
library(openxlsx)      # Lectura y escritura de archivos Excel
library(ebirdst)       # Descarga y análisis de modelos espacio-temporales de eBird
library(exactextractr)



# ------------------------------------------------------------------------------
# Variables respuesta
# ------------------------------------------------------------------------------



# 2. Preparación datos de abundancia de Ebird Status & Trends ------------------

# Carga del listado de especies utilizado para el análisis (cuya estructura es
# especie | codigo_ebirdst | gremio).
especies <- read.csv("Datos_codigo/especies_SC.csv",
                     fileEncoding = "Latin1",
                     sep = ";", 
                     stringsAsFactors = FALSE)

# Comprobaciones (si te interesa)
# head(especies)
# str (especies)
# View(especies)

# Descarga de los datos de EBird Status & Trends (solo una vez) 
# for(i in 1:nrow(especies)){
#   
#   codigo <- especies$codigo_ebirdst[i]
#   
#   cat("\nDescargando:", codigo, "\n")
#   
#   try(
#     
#     ebirdst_download_status(
#       species = codigo,
#       download_abundance = TRUE,
#       force = TRUE
#     ),
#     
#     silent = TRUE
#   )
# }



# 3. Creación del dataframe especies - rasgos (AVONET) -------------------------

# Carga de los datos de rasgos
rasgos <- readWorkbook(
  "Datos_codigo/ELEData/ELEData/TraitData/AVONET2_eBird.xlsx", 
  sheet = "AVONET2_eBird")

# Comprobaciones (si te interesa)
# head(rasgos)
# str(rasgos)
# View(rasgos)

# Creación del dataframe
especies_rasgos <- especies %>%
  left_join(
    rasgos %>%
      select(
        Species2,
        Mass,
        `Hand-Wing.Index`,
        Tarsus.Length,
        Wing.Length,
        Trophic.Level,
        Trophic.Niche,
        Primary.Lifestyle,
        Migration
      ) %>%
      rename(
        especie = Species2,
        Masa = Mass,
        HWI = `Hand-Wing.Index`,
        Longitud_tarso = Tarsus.Length,
        Longitud_ala = Wing.Length,
        Nivel_trofico = Trophic.Level,
        Nicho_trofico = Trophic.Niche,
        Estilo_vida = Primary.Lifestyle,
        Migracion = Migration
      ) %>%
      mutate(
        Nivel_trofico = recode(
          Nivel_trofico,
          "Carnivore" = "Carnívoro",
          "Omnivore" = "Omnívoro",
          "Herbivore" = "Herbívoro",
          "Scavenger" = "Carroñero"
        ),
        Nicho_trofico = recode(
          Nicho_trofico,
          "Vertivore" = "Vertebradófago",
          "Invertivore" = "Invertebradófago",
          "Scavenger" = "Carroñero",
          "Omnivore" = "Omnívoro",
          "Aquatic predator" = "Depredador acuático",
          "Frugivore" = "Frugívoro",
          "Herbivore aquatic" = "Herbívoro acuático",
          "Herbivore terrestrial" = "Herbívoro terrestre",
          "Nectarivore" = "Nectarívoro",
          "Granivore" = "Granívoro"
        ),
        Estilo_vida = recode(
          Estilo_vida,
          "Insessorial" = "Posador",
          "Generalist" = "Generalista",
          "Aerial" = "Aéreo",
          "Terrestrial" = "Terrestre",
          "Aquatic" = "Acuático"
        ),
        Migracion = recode(
          as.character(Migracion),
          "1" = "Residente",
          "2" = "Migradora parcial",
          "3" = "Migradora"
        )
      ),
    by = "especie"
  ) %>%
  filter(!is.na(Masa)) %>%
  rename(
    Especie = especie,
    Codigo_ebirdst = codigo_ebirdst,
    Gremio = gremio
  )

# Comprobaciones (si te interesa)
# head(especies_rasgos)
# str(especies_rasgos)
# View(especies_rasgos)

# Guardado del dataframe
write_csv(especies_rasgos,"Datos_codigo/Rasgos.csv")



# 4.Creación del dataframe celdas - especies -----------------------------------

# Carga del recorte del área de estudio y de la gradilla de 3 km x 3 km
grid_vect <- readRDS("Datos_codigo/grid_vect.rds")
grid_ebird <- readRDS ("Datos_codigo/grid_ebird.rds")
recorte_SC <- readRDS("Datos_codigo/recorte_SC.rds")

# Filtrado de especies sin rasgos disponibles
n_inicial <- nrow(especies)

especies <- especies %>%
  semi_join(
    especies_rasgos,
    by = c("codigo_ebirdst" = "Codigo_ebirdst")
  )

cat("Se eliminaron", n_inicial - nrow(especies), "especies por falta de rasgos.\n")

# Creación del dataframe
Abundancia_bruta <- grid_ebird %>%
  st_drop_geometry() %>%
  select(ID)

# Recorrido por todas las especies
for(i in 1:nrow(especies)){
  
  # Información de la especie
  nombre_cientifico <- especies$especie[i]
  codigo_ebird <- especies$codigo_ebirdst[i]
  
  cat("\n")
  cat("Procesando:", nombre_cientifico, "\n")
  
  # Carga del raster estacional de abundancia
  raster_sp <- load_raster(
    codigo_ebird,
    product = "abundance",
    period = "seasonal",
    metric = "mean",
    resolution = "3km"
  )
  
  # Selección de la capa de abundancia reproductora.
  # Para especies residentes se utiliza la capa "resident" y, en caso
  # contrario, la capa "breeding".
  capas <- names(raster_sp)
  
  if("resident" %in% capas){
    
    abundancia <- raster_sp[["resident"]]
    
  } else if("breeding" %in% capas){
    
    abundancia <- raster_sp[["breeding"]]
    
  } else{
    
    warning(paste("No se encontró capa reproductora para", nombre_cientifico))
    next
    
  }
  
  # Recorte al área de estudio
  abundancia_sc <- crop(
    abundancia,
    vect(recorte_SC)
  )
  
  # Extracción de abundancia media por celda
  abund_ext <- terra::extract(
    abundancia_sc,
    vect(grid_ebird),
    fun = mean,
    na.rm = TRUE
  )
  
  # Incorporación a la base de datos
  Abundancia_bruta[[nombre_cientifico]] <- abund_ext[,2]
  
}

# Ordenación de columnas
Abundancia_bruta <- Abundancia_bruta %>%
  select(ID, everything())

# Guardado RDS
saveRDS(
  Abundancia_bruta,
  "Datos_codigo/Abundancia_bruta.rds"
)

# Guardado CSV
write_csv(
  Abundancia_bruta,
  "Datos_codigo/Abundancia_bruta.csv"
)



# ------------------------------------------------------------------------------
# Variables explicativas
# ------------------------------------------------------------------------------



# 6. Preparación de los datos de Corine Land Cover -----------------------------

# Rasteres CORINE previamente reproyectados a WGS84 / Equal Earth Greenwich
CLC <- list(
  "1990" = rast("Datos_codigo/CLC_1990_wgs84_EEG.tif"),
  "2000" = rast("Datos_codigo/CLC_2000_wgs84_EEG.tif"),
  "2006" = rast("Datos_codigo/CLC_2006_wgs84_EEG.tif"),
  "2012" = rast("Datos_codigo/CLC_2012_wgs84_EEG.tif"),
  "2018" = rast("Datos_codigo/CLC_2018_wgs84_EEG.tif")
)

# Comprobación rápida
crs(CLC$`2018`)
levels(CLC$`2018`)

# Simplificación de categorías CORINE
CLC_original <- rast(
  "CLC/CLC_2018/DATA/U2018_CLC2018_V2020_20u1.tif"
)

niveles <- cats(CLC_original)[[1]]

niveles$Cod_3 <- niveles$CODE_18

niveles$Cod_Adaptado <- ifelse(
  substr(niveles$Cod_3, 1, 1) == "3",
  substr(niveles$Cod_3, 1, 2),
  paste0(substr(niveles$Cod_3, 1, 1), "0")
)

reclasificacion <- niveles[, c("Value", "Cod_Adaptado")]

# Aplicación a todos los años
CLC_simple <- lapply(
  CLC,
  classify,
  rcl = reclasificacion
)

# Etiquetas simplificadas
nombres_corine <- data.frame(
  Cod_Adaptado = c(10, 20, 31, 32, 33, 40, 50, 60),
  name = c(
    "Superficies artificiales",
    "Zonas agrícolas",
    "Bosques",
    "Matorrales y/o vegetación herbácea",
    "Espacios abiertos",
    "Zonas húmedas",
    "Superficies de agua",
    "SIN DATO"
  )
)

for(i in seq_along(CLC_simple)){
  levels(CLC_simple[[i]]) <- nombres_corine
}

# Colores para mapas
colores_corine <- c(
  "Superficies artificiales" = "#e31a1c",
  "Zonas agrícolas" = "#ff7f00",
  "Bosques" = "#006400",
  "Matorrales y/o vegetación herbácea" = "#66a61e",
  "Espacios abiertos" = "#b2df8a",
  "Zonas húmedas" = "#1f78b4",
  "Superficies de agua" = "#a6cee3",
  "SIN DATO" = "#d9d9d9"
)

# Área de estudio
recorte_SC <- st_as_sfc(
  st_bbox(
    c(
      xmin = -6.5,
      xmax = -2.5,
      ymin = 39.5,
      ymax = 41.5
    ),
    crs = 4326
  )
)

recorte_SC <- st_transform(
  recorte_SC,
  crs(CLC_simple[[1]])
)

# Recorte de todos los CORINE
CLC_SC <- lapply(
  CLC_simple,
  crop,
  y = recorte_SC
)

# Visualización rápida
tm_shape(CLC_SC$`2018`) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = colores_corine
    ),
    col.legend = tm_legend(
      title = "Uso del suelo"
    )
  ) +
  tm_scalebar(
    position = c("left", "bottom")
  ) +
  tm_compass(
    position = c("right", "top")
  ) +
  tm_layout(
    main.title = "Uso del suelo (CORINE 2018)",
    legend.outside = TRUE
  )



# 7. Variable explicativa: trayectoria dominante -------------------------------

# Construcción de una variable categórica que resume la trayectoria de cambio
# del uso del suelo seguida por cada celda entre 1990 y 2018.
#
# La trayectoria se calcula primero para cada píxel CORINE (100 m) y
# posteriormente se obtiene la trayectoria modal dentro de cada celda de la
# malla eBird.

# Conversión de la malla eBird a SpatVector
grid_vect <- vect(grid_ebird)

# Equivalencia entre categorías CORINE y códigos internos
codigos <- c(
  "10" = 1,
  "20" = 2,
  "31" = 3,
  "32" = 4,
  "33" = 5,
  "40" = 6,
  "50" = 7
)

equivalencias <- c(
  10, 20, 31, 32, 33, 40, 50
)

# Función para codificar la trayectoria de cada píxel
codificar_trayectoria <- function(x){
  
  if(any(is.na(x))){
    return(NA_integer_)
  }
  
  as.integer(
    paste0(
      codigos[as.character(x)],
      collapse = ""
    )
  )
  
}

# Cálculo de la trayectoria de cada píxel
trayectoria_raster <- app(
  c(
    CLC_SC[["1990"]],
    CLC_SC[["2000"]],
    CLC_SC[["2006"]],
    CLC_SC[["2012"]],
    CLC_SC[["2018"]]
  ),
  fun = codificar_trayectoria,
  cores = 1,
  filename = tempfile(fileext = ".tif"),
  overwrite = TRUE
)

# Raster con el identificador de cada celda eBird
grid_ID <- rasterize(
  grid_vect,
  trayectoria_raster,
  field = "ID"
)

# Obtención de la trayectoria modal por celda
trayectorias_modales <- terra::zonal(
  trayectoria_raster,
  grid_ID,
  fun = "modal",
  na.rm = TRUE
)

colnames(trayectorias_modales) <- c("ID", "Codigo")

# Función para recuperar las cinco clases CORINE
decodificar_trayectoria <- function(codigo){
  
  if(is.na(codigo)){
    return(rep(NA, 5))
  }
  
  codigo <- sprintf("%05d", codigo)
  
  equivalencias[
    as.integer(
      strsplit(codigo, "")[[1]]
    )
  ]
  
}

# Creación del dataframe de usos del suelo
datos_corine <- trayectorias_modales %>%
  mutate(
    Clases = lapply(Codigo, decodificar_trayectoria)
  ) %>%
  unnest_wider(
    Clases,
    names_sep = "_"
  ) %>%
  rename(
    CLC_1990 = Clases_1,
    CLC_2000 = Clases_2,
    CLC_2006 = Clases_3,
    CLC_2012 = Clases_4,
    CLC_2018 = Clases_5
  ) %>%
  select(
    ID,
    CLC_1990,
    CLC_2000,
    CLC_2006,
    CLC_2012,
    CLC_2018
  )

# Incorporación de la información a la malla eBird
grid_ebird <- grid_ebird %>%
  left_join(
    datos_corine,
    by = "ID"
  )

# Eliminación de celdas sin información
grid_ebird <- grid_ebird %>%
  filter(
    !is.na(CLC_1990),
    !is.na(CLC_2000),
    !is.na(CLC_2006),
    !is.na(CLC_2012),
    !is.na(CLC_2018)
  )

# Actualización de la versión SpatVector
grid_vect <- vect(grid_ebird)

# Creación de la trayectoria temporal
grid_ebird <- grid_ebird %>%
  rowwise() %>%
  mutate(
    Trayectoria = paste(
      CLC_1990,
      CLC_2000,
      CLC_2006,
      CLC_2012,
      CLC_2018,
      sep = "-"
    ),
    n_clases = n_distinct(
      c(
        CLC_1990,
        CLC_2000,
        CLC_2006,
        CLC_2012,
        CLC_2018
      )
    )
  ) %>%
  ungroup()

# Clasificación de las trayectorias de cambio
grid_ebird <- grid_ebird %>%
  mutate(
    Tipo = case_when(
      
      # Permanencias
      CLC_1990 == 10 & CLC_2000 == 10 & CLC_2006 == 10 &
        CLC_2012 == 10 & CLC_2018 == 10 ~
        "Artificial permanente",
      
      CLC_1990 == 20 & CLC_2000 == 20 & CLC_2006 == 20 &
        CLC_2012 == 20 & CLC_2018 == 20 ~
        "Agrícola permanente",
      
      CLC_1990 == 31 & CLC_2000 == 31 & CLC_2006 == 31 &
        CLC_2012 == 31 & CLC_2018 == 31 ~
        "Bosque permanente",
      
      CLC_1990 == 32 & CLC_2000 == 32 & CLC_2006 == 32 &
        CLC_2012 == 32 & CLC_2018 == 32 ~
        "Matorral permanente",
      
      CLC_1990 == 33 & CLC_2000 == 33 & CLC_2006 == 33 &
        CLC_2012 == 33 & CLC_2018 == 33 ~
        "Esp. abierto permanente",
      
      CLC_1990 == 40 & CLC_2000 == 40 & CLC_2006 == 40 &
        CLC_2012 == 40 & CLC_2018 == 40 ~
        "Zona húmeda permanente",
      
      CLC_1990 == 50 & CLC_2000 == 50 & CLC_2006 == 50 &
        CLC_2012 == 50 & CLC_2018 == 50 ~
        "Sup. de agua permanente",
      
      # Artificialización
      CLC_1990 == 20 & CLC_2018 == 10 & n_clases <= 2 ~
        "Artificialización (origen agrícola)",
      
      CLC_1990 == 31 & CLC_2018 == 10 & n_clases <= 2 ~
        "Artificialización (origen bosque)",
      
      CLC_1990 == 32 & CLC_2018 == 10 & n_clases <= 2 ~
        "Artificialización (origen matorral)",
      
      CLC_1990 == 33 & CLC_2018 == 10 & n_clases <= 2 ~
        "Artificialización (origen abierto)",
      
      CLC_2018 == 10 & n_clases >= 3 ~
        "Artificialización (trayectoria compleja)",
      
      CLC_2018 == 10 ~
        "Artificialización (otro)",
      
      # Recuperación forestal
      CLC_1990 == 31 &
        CLC_2018 == 31 &
        n_clases <= 3 ~
        "Recuperación forestal",
      
      # Degradación del bosque
      CLC_1990 == 31 &
        CLC_2018 == 32 &
        n_clases <= 3 ~
        "Degradación del bosque (a matorral)",
      
      CLC_1990 == 31 &
        CLC_2018 == 33 &
        n_clases <= 3 ~
        "Degradación del bosque (a espacio abierto)",
      
      # Sucesión natural
      CLC_1990 == 32 &
        CLC_2018 == 31 &
        n_clases <= 3 ~
        "Sucesión de matorral a bosque",
      
      CLC_1990 == 33 &
        CLC_2018 == 32 &
        n_clases <= 3 ~
        "Sucesión de espacio abierto a matorral",
      
      CLC_1990 == 33 &
        CLC_2018 == 31 &
        n_clases <= 3 ~
        "Sucesión de espacio abierto a bosque",
      
      # Abandono agrícola
      CLC_1990 == 20 &
        CLC_2018 == 31 &
        n_clases <= 3 ~
        "Abandono agrícola (sucesión a bosque)",
      
      CLC_1990 == 20 &
        CLC_2018 == 32 &
        n_clases <= 2 ~
        "Abandono agrícola (sucesión a matorral)",
      
      CLC_1990 == 20 &
        CLC_2018 == 33 &
        n_clases <= 2 ~
        "Abandono agrícola (sucesión a espacio abierto)",
      
      # Expansión agrícola
      CLC_1990 == 31 &
        CLC_2018 == 20 &
        n_clases <= 2 ~
        "Expansión agrícola (origen bosque)",
      
      CLC_1990 == 32 &
        CLC_2018 == 20 &
        n_clases <= 2 ~
        "Expansión agrícola (origen matorral)",
      
      CLC_1990 == 33 &
        CLC_2018 == 20 &
        n_clases <= 2 ~
        "Expansión agrícola (origen abierto)",
      
      CLC_2018 == 20 &
        n_clases >= 3 ~
        "Expansión agrícola (trayectoria compleja)",
      
      CLC_2018 == 20 ~
        "Expansión agrícola (otro)",
      
      # Resto de trayectorias
      TRUE ~ "Otros"
    )
  )

# Comprobaciones (si te interesa)
table(grid_ebird$Tipo)

Resumen_grid <- grid_ebird %>%
  st_drop_geometry() %>%
  count(Tipo, CLC_2018) %>%
  mutate(
    Estado_final = recode(
      as.character(CLC_2018),
      "10" = "Superficies artificiales",
      "20" = "Zonas agrícolas",
      "31" = "Bosques",
      "32" = "Matorrales y/o vegetación herbácea",
      "33" = "Espacios abiertos",
      "40" = "Zonas húmedas",
      "50" = "Superficies de agua"
    ),
    Porcentaje = n / sum(n) * 100
  ) %>%
  filter(n > 20)

# Variable específica para el color de la gráfica
Resumen_grid <- Resumen_grid %>%
  mutate(
    Estado_final_grafica = if_else(
      Tipo == "Otros",
      "Otros",
      Estado_final
    )
  )

# Gráfica
ggplot(
  Resumen_grid,
  aes(
    x = reorder(Tipo, -Porcentaje),
    y = Porcentaje,
    fill = Estado_final_grafica
  )
) +
  geom_col() +
  
  scale_fill_manual(
    values = c(
      colores_corine,
      "Otros" = "#808080"
    ),
    name = "Estado final"
  ) +
  
  labs(
    title = "Tipos de cambio de uso del suelo",
    subtitle = "Solo trayectorias con más de 20 celdas",
    x = "Tipo de trayectoria",
    y = "% de celdas"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(face = "bold"),
    
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Guardado de resultados
Trayectorias <- grid_ebird %>%
  select(
    ID,
    Trayectoria,
    Tipo,
    geometry
  )

saveRDS(
  Trayectorias,
  "Datos_codigo/Trayectorias.rds"
)



# 8. Variable explicativa: Uso actual ------------------------------------------

# Creación de una variable categórica que representa el uso del suelo dominante
# en cada celda durante el año 2018, a partir de la clasificación simplificada
# de CORINE Land Cover. Esta variable describe las condiciones actuales del
# paisaje, independientemente de la trayectoria seguida hasta alcanzarlas.

grid_ebird$Uso_actual <- factor(
  grid_ebird$CLC_2018,
  levels = c(10,20,31,32,33,40,50,60),
  labels = c(
    "Superficies artificiales",
    "Zonas agrícolas",
    "Bosques",
    "Matorrales y/o vegetación herbácea",
    "Espacios abiertos",
    "Zonas húmedas",
    "Superficies de agua",
    "SIN DATO"
  )
)

# Comprobaciones (si te interesa)
table(grid_ebird$Uso_actual)

# Guardado de resultados
Uso_actual <- grid_ebird %>%
  select(
    ID,
    Uso_actual,
    geometry
  )

saveRDS(
  Uso_actual,
  "Datos_codigo/Uso_actual.rds"
)



# 9. Variable explicativa: Estabilidad temporal -------------------------------

# Cálculo del año desde el que el uso del suelo permanece estable en cada
# píxel CORINE hasta 2018. Posteriormente se obtiene el año modal dentro de
# cada celda de la malla eBird.

# Función para calcular la estabilidad de un píxel
calcular_estabilidad <- function(x){
  
  if(any(is.na(x))){
    return(NA_integer_)
  }
  
  if(all(x == x[5])){
    return(1990)
  }
  
  if(all(x[2:5] == x[5])){
    return(2000)
  }
  
  if(all(x[3:5] == x[5])){
    return(2006)
  }
  
  if(all(x[4:5] == x[5])){
    return(2012)
  }
  
  2018
  
}

# Cálculo de la estabilidad para cada píxel CORINE
estabilidad_raster <- app(
  c(
    CLC_SC[["1990"]],
    CLC_SC[["2000"]],
    CLC_SC[["2006"]],
    CLC_SC[["2012"]],
    CLC_SC[["2018"]]
  ),
  fun = calcular_estabilidad,
  cores = 1,
  filename = tempfile(fileext = ".tif"),
  overwrite = TRUE
)

# Obtención del año modal por celda eBird
estabilidad_modal <- terra::zonal(
  estabilidad_raster,
  grid_ID,
  fun = "modal",
  na.rm = TRUE
)

colnames(estabilidad_modal) <- c("ID", "anio_estable")

# Incorporación a la malla eBird
grid_ebird <- grid_ebird %>%
  left_join(
    estabilidad_modal,
    by = "ID"
  )

# Eliminación de celdas sin información
grid_ebird <- grid_ebird %>%
  filter(
    !is.na(anio_estable)
  )

# Conversión a factor ordenado
grid_ebird$anio_estable <- factor(
  grid_ebird$anio_estable,
  levels = c(1990, 2000, 2006, 2012, 2018)
)

# Comprobaciones (si te interesa)
table(grid_ebird$anio_estable)

# Resumen de datos: Cálculo del porcentaje de celdas por año de estabilidad
resumen_estabilidad <- grid_ebird %>%
  count(anio_estable) %>%
  mutate(Porcentaje = (n / sum(n)) * 100)

# Generación de la gráfica
p1 <- ggplot(
  resumen_estabilidad,
  aes(
    x = anio_estable,
    y = Porcentaje
  )
) +
  geom_col() +
  labs(
    title = "Estabilidad temporal del uso del suelo",
    subtitle = "Año modal a partir del cual el uso permanece estable por celda",
    x = "Año de inicio de estabilidad",
    y = "% de celdas"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Guardado de resultados
Estabilidad <- grid_ebird %>%
  select(
    ID,
    anio_estable,
    geometry
  )

saveRDS(
  Estabilidad,
  "Datos_codigo/Estabilidad.rds"
)



# 10. Variable explicativa: Porcentaje de cambio -------------------------------

# Cálculo del porcentaje de superficie transformada dentro de cada celda de la
# malla eBird. Para ello, cada píxel CORINE se clasifica como estable si
# mantiene el mismo uso del suelo en los cinco años disponibles o como cambiado
# si presenta alguna modificación. Posteriormente se calcula el porcentaje de
# píxeles cambiados y estables en cada celda. 

# Creación de un stack con los cinco años CORINE
CLC_stack <- c(
  CLC_SC[["1990"]],
  CLC_SC[["2000"]],
  CLC_SC[["2006"]],
  CLC_SC[["2012"]],
  CLC_SC[["2018"]]
)

names(CLC_stack) <- c(
  "y1990",
  "y2000",
  "y2006",
  "y2012",
  "y2018"
)

# Conversión del stack a tabla de píxeles
CLC_pixels <- as.data.frame(
  CLC_stack,
  xy = TRUE,
  na.rm = TRUE
)

# Clasificación de los píxeles según si han cambiado de uso del suelo
CLC_pixels <- CLC_pixels %>%
  mutate(
    estable = ifelse(
      y1990 == y2018 &
        y2000 == y2018 &
        y2006 == y2018 &
        y2012 == y2018,
      1,
      0
    )
  )

# Conversión a SpatVector
CLC_points <- vect(
  CLC_pixels,
  geom = c("x", "y"),
  crs = crs(CLC_stack)
)

# Asignación de cada píxel CORINE a una celda eBird
CLC_points$ID <- terra::extract(
  grid_vect,
  CLC_points
)[,"ID"]

CLC_points <- CLC_points[
  !is.na(CLC_points$ID),
]

# Cálculo del porcentaje de cambio por celda
porcentaje_cambio_grid <- CLC_points %>%
  as.data.frame() %>%
  group_by(ID) %>%
  summarise(
    Porcentaje_cambio = mean(estable == 0) * 100,
    Porcentaje_estabilidad = mean(estable == 1) * 100,
    .groups = "drop"
  )

# Incorporación de las variables a la malla eBird
grid_ebird <- grid_ebird %>%
  left_join(
    porcentaje_cambio_grid,
    by = "ID"
  )

# Comprobaciones (si te interesa)
summary(grid_ebird$Porcentaje_cambio)
summary(grid_ebird$Porcentaje_estabilidad)

# Categorización en rangos y cálculo de porcentaje de celdas
resumen_rangos_cambio <- grid_ebird %>%
  mutate(
    Rango_cambio = cut(
      Porcentaje_cambio,
      breaks = c(-Inf, 0, 20, 40, 60, 80, 100),
      labels = c("Sin cambio (0%)", "1 - 20%", "21 - 40%", "41 - 60%", 
                 "61 - 80%", "81 - 100%"),
      include.lowest = TRUE
    )
  ) %>%
  count(Rango_cambio) %>%
  mutate(Porcentaje = (n / sum(n)) * 100)

# Generación de la gráfica
p2 <-ggplot(
  resumen_rangos_cambio,
  aes(
    x = Rango_cambio,
    y = Porcentaje
  )
) +
  geom_col() +
  labs(
    title = "Distribución del porcentaje de cambio por celda",
    subtitle = "Proporción de celdas según la intensidad de cambio",
    x = "Rango de % de píxeles cambiados",
    y = "% de celdas"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# EXTRA: Gráfica de estabilidad + gráfica de intensidad de cambio
p1 + p2

# Guardado de resultados
Porcentaje_cambio <- grid_ebird %>%
  select(
    ID,
    Porcentaje_cambio,
    Porcentaje_estabilidad,
    geometry
  )

saveRDS(
  Porcentaje_cambio,
  "Datos_codigo/Porcentaje_cambio.rds"
)

saveRDS(
  CLC_pixels,
  "Datos_codigo/Cambio_pixel.rds"
)



# 11. Mapas de las variables explicativas --------------------------------------

# Mapa de % de cambio de uso del suelo
mapa_cambio <- tm_shape(grid_ebird) +
  tm_polygons(
    fill = "Porcentaje_cambio",
    fill.scale = tm_scale_continuous(values = "viridis"),
    fill.legend = tm_legend(title = "Porcentaje de cambio"),
    col = NA
  ) +
  tm_scalebar(position = c("left", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_layout(
    main.title = "Porcentaje de píxeles cambiados",
    main.title.position = "center",
    main.title.fontface = "bold",
    legend.outside = TRUE,
    legend.outside.position = tm_pos_out("right", "center"),
    legend.title.fontface = "bold"
  )

# Mapa de estabilidad en el uso del suelo (modal)
mapa_estabilidad <- tm_shape(grid_ebird) +
  tm_polygons(
    fill = "anio_estable",
    fill.scale = tm_scale_categorical(values = "plasma"),
    fill.legend = tm_legend(title = "Último cambio"),
    col = NA
  ) +
  tm_scalebar(position = c("left", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_layout(
    main.title = "Año desde el que el uso actual permanece estable",
    main.title.position = "center",
    main.title.fontface = "bold",
    legend.outside = TRUE,
    legend.outside.position = tm_pos_out("right", "center"),
    legend.title.fontface = "bold"
  )

# Mapa de uso del suelo actual (modal)
mapa_uso_actual <- tm_shape(grid_ebird) +
  tm_polygons(
    fill = "Uso_actual",
    fill.scale = tm_scale_categorical(values = colores_corine),
    fill.legend = tm_legend(title = "Uso actual"),
    col = NA
  ) +
  tm_scalebar(position = c("left", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_layout(
    main.title = "Uso del suelo dominante (2018)",
    main.title.position = "center",
    main.title.fontface = "bold",
    legend.outside = TRUE,
    legend.outside.position = tm_pos_out("right", "center"),
    legend.title.fontface = "bold"
  )


colores_corine <- c(
  "Superficies artificiales" = "#e31a1c",
  "Zonas agrícolas" = "#ff7f00",
  "Bosques" = "#006400",
  "Matorrales y/o vegetación herbácea" = "#66a61e",
  "Espacios abiertos" = "#b2df8a",
  "Zonas húmedas" = "#1f78b4",
  "Superficies de agua" = "#a6cee3",
  "SIN DATO" = "#d9d9d9"
)


# Colores para las trayectorias de cambio
colores_trayectorias <- c(
  
  # Permanencias
  "Artificial permanente" = "#e31a1c",
  "Agrícola permanente" = "#ff7f00",
  "Bosque permanente" = "#006400",
  "Matorral permanente" = "#41ab5d",
  "Esp. abierto permanente" = "#a1d99b",
  "Zona húmeda permanente" = "#3182bd",
  "Sup. de agua permanente" = "#9ecae1",
  
  # Artificialización (morados)
  "Artificialización (origen agrícola)" = "#6a3d9a",
  "Artificialización (origen bosque)" = "#7b3294",
  "Artificialización (origen matorral)" = "#8e44ad",
  "Artificialización (origen abierto)" = "#9b59b6",
  "Artificialización (trayectoria compleja)" = "#c994c7",
  "Artificialización (otro)" = "#d4b9da",
  
  # Perturbación y recuperación forestal (verde oscuro)
  "Recuperación forestal" = "#006d2c",
  
  # Degradación forestal (marrones)
  "Degradación del bosque (a matorral)" = "#8c510a",
  "Degradación del bosque (a espacio abierto)" = "#bf812d",
  
  # Sucesión natural (verde-azulado)
  "Sucesión de matorral a bosque" = "#1b9e77",
  "Sucesión de espacio abierto a matorral" = "#4eb3a5",
  "Sucesión de espacio abierto a bosque" = "#66c2a4",
  
  # Abandono agrícola (verdes oliva)
  "Abandono agrícola (sucesión a bosque)" = "#66a61e",
  "Abandono agrícola (sucesión a matorral)" = "#8daa3a",
  "Abandono agrícola (sucesión a espacio abierto)" = "#b8c96f",
  
  # Expansión agrícola (amarillos)
  "Expansión agrícola (origen bosque)" = "#fdbf11",
  "Expansión agrícola (origen matorral)" = "#ffd92f",
  "Expansión agrícola (origen abierto)" = "#ffe680",
  "Expansión agrícola (trayectoria compleja)" = "#f6e8c3",
  "Expansión agrícola (otro)" = "#fff7bc",
  
  # Residual
  "Otros" = "grey60"
)

# Mapa de trayectorias dominantes
mapa_trayectorias <- tm_shape(grid_ebird) +
  tm_polygons(
    fill = "Tipo",
    fill.scale = tm_scale_categorical(values = colores_trayectorias),
    col = NA
  ) +
  tm_scalebar(position = c("left", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_layout(
    main.title = "Trayectoria dominante del uso del suelo (1990-2018)",
    main.title.position = "center",
    main.title.fontface = "bold",
    legend.outside = TRUE,
    legend.outside.position = tm_pos_out("right", "center"),
    legend.title.fontface = "bold"
  )

# Visualización
mapa_uso_actual
mapa_cambio
mapa_estabilidad
mapa_trayectorias
