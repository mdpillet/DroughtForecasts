library(RCurl)

options(timeout = 36000)

#
# Get current data
#

# Set directory structure
relPath <- "D:/Research/Chapter3/"
outPath <- "Data/CHELSA/Raw/Current/"

# Get file names
fileNames <- read.table(paste0(relPath, "Data/CHELSA/current.txt"))

# Download files
for (i in 1:nrow(fileNames)) {
  print(i)
  outName <- strsplit(fileNames[i, 1], "/", fixed = T)[[1]][11]
  download.file(fileNames[i, 1], paste0(relPath, outPath, outName), mode = "wb")
}

#
# Get future data
#

# Set directory structure
outPath <- "Data/CHELSA/Raw/Future/"

# Get file names
fileNames <- read.table(paste0(relPath, "Data/CHELSA/future.txt"))

# Download files
for (i in 1:nrow(fileNames)) {
  print(i)
  outName <- strsplit(fileNames[i, 1], "/", fixed = T)[[1]][13]
  download.file(fileNames[i, 1], paste0(relPath, outPath, outName), mode = "wb")
}