# Set folder structure
outputPath <- "D:/Research/DroughtForecasts/Outputs/"

tmp <- list.dirs(path = outputPath, recursive = F)

remCount <- 0
for (i in 1:length(tmp)) {
  if (file.exists(paste0(tmp[i], "/ModelSetComparison/modelSetComparison.csv"))) {
    res <- read.csv(paste0(tmp[i], "/ModelSetComparison/modelSetComparison.csv"), header = T)
    if (sum(res$AverageValidationAUC <= 0.5) > 0) {
      remCount <- remCount + 1
      print(tmp[i])
    }
  }
}
print(remCount)
