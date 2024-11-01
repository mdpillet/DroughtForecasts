# Define data sources
relDir <- "D:/Research/DroughtForecasts/Data/"
iNat <- "iNaturalist/Processed/INAT.csv"
MABDB <- "MarcBaker/Processed/MABDB_final.csv"
MABHERB <- "MarcBaker/Processed/MABHERB_final.csv"
PGCOP <- "PabloGuerrero/Processed/PGCOP.csv"
PGCACT <- "PabloGuerrero/Processed/PGCACT.csv"
PGHERB <- "PabloGuerrero/Processed/PGHERB.csv"
LJM <- "LucasMajure/Processed/LJM.csv"
BIEN <- "BIEN/Processed/final.csv"
outPath <- "Occurrences/combined.csv"

# Read data sources and standardize column names
data1 <- read.csv(paste0(relDir, iNat), header = T)
data1 <- data1[, c("AccNoGlobal", "scientific_name", "latitude", "longitude")]
names(data1) <- c("AccNo", "Taxon", "lat", "lon")

data2 <- read.csv(paste0(relDir, MABDB), header = T)
data2$Taxon <- paste0(data2$Taxon1, " ", data2$Taxon2)
data2 <- data2[, c("AccNoGlobal", "Taxon", "lat", "lon")]
names(data2) <- c("AccNo", "Taxon", "lat", "lon")

data3 <- read.csv(paste0(relDir, MABHERB), header = T)
data3 <- data3[, c("AccNoGlobal", "scientificName", "decimalLatitude", "decimalLongitude")]
names(data3) <- c("AccNo", "Taxon", "lat", "lon")

data4 <- read.csv(paste0(relDir, PGCOP), header = T)
data4 <- data4[, c("CODIGO", "TAXON", "LATITUD", "LONGITUD")]
names(data4) <- c("AccNo", "Taxon", "lat", "lon")

data5 <- read.csv(paste0(relDir, PGCACT), header = T)
for (i in 1:nrow(data5)) {
  data5[i, "Taxon"] <- paste0(data5[i, "Genero"], " ", data5[i, "Especie"])
  if (data5[i, "subespecie"] != "") data5[i, "Taxon"] <- paste0(data5[i, "Taxon"], " ssp. ", data5[i, "subespecie"])
  if (data5[i, "variedad"] != "") data5[i, "Taxon"] <- paste0(data5[i, "Taxon"], " var. ", data5[i, "variedad"])
}
data5 <- data5[, c("AccNoGlobal", "Taxon", "latitudDecimal", "longitudDecimal")]
names(data5) <- c("AccNo", "Taxon", "lat", "lon")

data6 <- read.csv(paste0(relDir, PGHERB), header = T)
for (i in 1:nrow(data6)) {
  data6[i, "Taxon"] <- paste0(data6[i, "GENERO"], " ", data6[i, "SP"])
  if (data6[i, "Subespecie"] != "") data6[i, "Taxon"] <- paste0(data6[i, "Taxon"], " ssp. ", data6[i, "Subespecie"])
  if (data6[i, "Variedad"] != "") data6[i, "Taxon"] <- paste0(data6[i, "Taxon"], " var. ", data6[i, "Variedad"])
}
data6 <- data6[, c("AccNoGlobal", "Taxon", "LAT", "LONG")]
names(data6) <- c("AccNo", "Taxon", "lat", "lon")

data7 <- read.csv(paste0(relDir, BIEN), header = T)
data7 <- data7[, c("AccNoGlobal", "name_matched", "latitude", "longitude")]
names(data7) <- c("AccNo", "Taxon", "lat", "lon")

data8 <- read.csv(paste0(relDir, LJM), header = T)
data8 <- data8[, c("AccNo", "Taxon", "lat", "lon")]

# Combine data sources
occ <- rbind(data1, data2, data3, data4, data5, data6, data7, data8)

# Remove incomplete rows
occ <- occ[complete.cases(occ),]

# Export occurrences
write.csv(occ, paste0(relDir, outPath), row.names = F)