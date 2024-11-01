library(todoBIEN)

# Set directory structure
relPath <- "F:/Chapter3/"
outPath <- "Data/BIEN/Raw/"
outFile <- "occurrences_10-30-2023.csv"

# Fetch occurrences
start_time <- Sys.time()
occurrences <- BIEN_occurrence_family(family = "Cactaceae",
                                      cultivated = F,
                                      new.world = T,
                                      observation.type = T,
                                      all.taxonomy = T,
                                      native.status = T,
                                      natives.only = T,
                                      political.boundaries = T,
                                      collection.info = T,
                                      only.geovalid = T,
                                      user = "", # Username and password removed
                                      password = "", # Username and password removed
                                      schema = "analytical_db")
end_time <- Sys.time()
time_diff <- end_time - start_time

# Export occurrence data
# write.csv(occurrences, paste0(relPath, outPath, outFile), row.names = F)