# DroughtForecasts
Drought forecasts for cacti
 
Environmental data preparation
* getCHELSAData.R: downloads current and future CHELSA data
* getSoilGridsData.R: downloads SoilGrids data.
* averageSoilGridsData.R: crops and creates weighted average of SoilGrids depth slices.
* resampleCHELSAData.R: resamples CHELSA data to match origin and resolution of SoilGrids data.
* resampleBIO2FIData.R: resamples BIO2FI data to match origin and resolution of SoilGrids data, creates historical and future averages of BIO2FI data to match periods of CHELSA data.
* stackPredictorData.R: stacks predictor layers.
* maskPredictorData.R: standardizes all predictor files to have identical masking.
* subsetPredictorData.R: excludes severe drought variables.
* renamePredictorData.R: updates variable names.

Occurrence data preparation
* getBIENOccurences.R: downloads occurrence data from BIEN.
* combineiNaturalistData.R: combines iNaturalist public and obscured data.
* cleaniNaturalistData.R: cleans iNaturalist data.
* cleanBakerOccurrences.R: standardizes coordinates for Marc Baker’s data.
* Data deriving from Lucas Majure, Pablo Guerrero and Marc Baker underwent manual cleaning, including coordinate standardization, removal of records not identified to species, coordinate correction (e.g. lat./lon. switched), removal of records with incomplete or unresolvable coordinates, and correction of certain characters not accepted by TNRS.
* combineOccurrenceData.R: combines occurrence data from different data sources.
* cleanOccurrenceData.R: cleans of coordinate data and subsets New World data.
* standardizeTaxonomy.R: standardizes taxonomy with Caryophyllales.org.
* See StandardizationNotes.txt for more information on manual taxonomic cleaning.
* checkInvasives.R: extracts country names for occurrences, standardizes country names with GNRS, and runs these against NSR, then, only for species with 10 or more records, removes records from countries where species are thought to be non-native, as well as species with considerable uncertainty about their native range.
* See InvasiveNotes.txt for more information on manual curation of invasive status.
* splitSpecies.R: splits subsetted occurrence records by species.
* occTest.R: subsets occurrences based on automated tests using occTest package, and separates out species with 10 or more final occurrences.

Create models and maps
* thinning.R: performs spatial thinning of filtered occurrence data.
* functions.R: contains auxiliary functions.
* workflow_sensitivityAnalysis.R: performs modeling workflow for a random subset of species.
* workflow.R: performs modeling workflow.

Species-level analyses
* sensitivityAnalysis.R: performs sensitivity analysis of models and maps to parameter choices.
* modelQuality.R: analyzes model quality.
* variableSelection.R: summarizes variable selection.
* calcRangeChanges.R: calculates range changes.
* analyzeRangeChanges.R: analyzes range changes.
* uncertaintyAnalysis: variance decomposition of modeling choices.

Richness-level analyses
* stackMaps.R: creates richness maps and corresponding summary maps.
* expertMap.R: creates richness map based on IUCN SSC CSSG Global Cactus Assessment expert maps and calculates correlations with model-based maps.
* bioregions.R: summarizes richness changes by bioregion following Calvente et al. (2023), Appendix S11.

iNaturalist data
* iNaturalist.zip: archive including all raw iNaturalist occurrences used for this study. Secured data are not included. See https://doi.org/10.15468/dd.6jxctw.