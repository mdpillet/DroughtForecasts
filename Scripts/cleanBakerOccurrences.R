library(stringr)
# Set as needed
setwd("~/Downloads")
##
MABDB_2023_working <- read.csv("MABDB_2023_working.csv")
coordenadas <- MABDB_2023_working$LatLon
head(coordenadas)

# Crear una lista para almacenar los resultados
resultados <- list()

convert_to_decimal <- function(coordinate) {
  # Normalizar la entrada reemplazando diferentes tipos de comillas y espacios
  coordinate <- gsub("[`'´\"]", "'", coordinate)
  coordinate <- gsub("\\s+", " ", trimws(coordinate))
  
  # Definir una función auxiliar para convertir partes a decimal
  to_decimal <- function(deg, min = 0, sec = 0, dir = "N") {
    value <- as.numeric(deg) + as.numeric(min) / 60 + as.numeric(sec) / 3600
    if (dir %in% c("S", "W")) value <- -value
    return(value)
  }
  
  # Procesar cada formato
  if(grepl("^[NS]\\d+", coordinate)) { # Formatos con N/S al inicio
    if(grepl("'", coordinate)) { # Grados y minutos, posiblemente segundos
      parts <- str_extract_all(coordinate, "\\d+\\.\\d+|\\d+|[NSWE]")[[1]]
      if(length(parts) == 6) { # Grados y minutos decimales
        lat <- to_decimal(parts[2], parts[3], 0, parts[1])
        lon <- to_decimal(parts[5], parts[6], 0, parts[4])
      } else { # Grados, minutos y segundos
        lat <- to_decimal(parts[2], parts[3], parts[4], parts[1])
        lon <- to_decimal(parts[6], parts[7], parts[8], parts[5])
      }
    } else { # Solo grados decimales
      parts <- str_extract_all(coordinate, "\\d+\\.\\d+|[NSWE]")[[1]]
      lat <- to_decimal(parts[2], 0, 0, parts[1])
      lon <- to_decimal(parts[4], 0, 0, parts[3])
    }
  } else if(grepl("^-?\\d+\\.\\d+", coordinate)) { # Formato decimal sin dirección explícita
    parts <- str_extract_all(coordinate, "-?\\d+\\.\\d+")[[1]]
    lat <- as.numeric(parts[1])
    lon <- as.numeric(parts[2])
  } else {
    return("Formato de coordenada no reconocido")
  }
  
  return(c(lat = lat, lon = lon))
}


# Iterar a través de las coordenadas y aplicar la función convert_to_decimal
for (i in 1:length(coordenadas)) {
  resultado <- convert_to_decimal(coordenadas[i])
  resultados[[i]] <- resultado
}

# Convertir la lista de resultados en un data frame
resultados_df <- do.call(rbind, resultados)

# Asignar nombres a las columnas del data frame de resultados
colnames(resultados_df) <- c("lat", "lon")

# Combinar los resultados con el DataFrame original
final_df <- cbind(MABDB_2023_working, resultados_df)

# Mostrar las primeras filas del DataFrame final
head(final_df)

# Set as needed
write.csv(final_df, "exp1.csv")