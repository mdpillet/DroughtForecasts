# DroughtForecasts
Drought forecasts for cacti
 
Environmental data preparation
* getCHELSAData.R: download current and future CHELSA data
* getSoilGridsData.R: downloads SoilGrids data.
* averageSoilGridsData.R: crop, create weighted average of SoilGrids depth slices.
* resampleCHELSAData.R: resample CHELSA data to match origin and resolution of SoilGrids data.
* resampleBIO2FIData.R: resample BIO2FI data to match origin and resolution of SoilGrids data, create historical and future averages of BIO2FI data to match periods of CHELSA data.
* stackPredictorData.R: stack predictor layers.
* maskPredictorData.R: standardize all predictor files to have identical masking.

Occurrence data preparation
* getBIENOccurences.R: download occurrence data from BIEN.
* combineiNaturalistData.R: combine iNaturalist public and obscured data.
* cleaniNaturalistData.R: clean iNaturalist data.
* cleanBakerOccurrences.R: code used to standardize coordinates for Marc Baker’s data.
* Data deriving from Lucas Majure, Pablo Guerrero and Marc Baker underwent manual cleaning, including coordinate standardization, removal of records not identified to species, coordinate correction (e.g. lat./lon. switched), removal of records with incomplete or unresolvable coordinates, and correction of certain characters not accepted by TNRS.
* combineOccurrenceData.R: combine occurrence data from different data sources.
* cleanOccurrenceData.R: automated cleaning of coordinate data, subset New World data.
* standardizeTaxonomy.R: standardize taxonomy with Caryophyllales.org.
* See StandardizationNotes.txt for more information on manual taxonomic cleaning.
* checkInvasives.R: extracts country names for occurrences, standardizes country names with GNRS, and runs these against NSR, then, only for species with 10 or more records, removes records from countries where species are thought to be non-native, as well as species with considerable uncertainty about their native range.
* See InvasiveNotes.txt for more information on manual curation of invasive status.
* splitSpecies.R: split subsetted occurrence records by species.
* occTest.R: subsets occurrences based on automated tests using occTest package, and separates out species with 10 or more final occurrences.

Create models and maps
* functions.R: auxiliary functions.
* workflow.R: performs modeling workflow.

Species-level analyses
* rangeChanges.R: analyzes range changes.
* modelComparison.R: compares model quality and summarizes variable inclusion rates.

Richness-level analyses
* stackMaps.R: creates richness maps and corresponding summary maps.
* expertMap.R: creates richness map based on IUCN SSC CSSG Global Cactus Assessment expert maps and calculates correlations with model-based maps.
* bioregions.R: summarizes richness changes by bioregion following Calvente et al. (2023), Appendix S11.

iNaturalist data
* iNaturalist.zip: archive including all raw iNaturalist occurrences used for this study. Secured data are not included. See https://doi.org/10.15468/dd.6jxctw.