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
* NOTE: For a further 24 species, sample size requirements were not met.
* NOTE: For 23 species, models failed to properly converge. These were excluded from analyses. Species: Cipocereus crassisepalus, Brachycereus nesioticus, Cereus mirabella, Discocactus bahiensis, Frailea pumila, Gymnocalycium schroederianum, Harrisia tortuosa, Maihueniopsis conoidea, Melocactus inconcinnus, M. levitestatus, M. pachyacanthus, Opuntia aurantiaca, O. bonaerensis, O. nemoralis, O. undulata, Parodia carambeiensis, P. oxycostata, Pelecyphora zilziana, Pterocactus fischeri, Rhipsalis grandiflora, R. shaferi, Selenicereus pteranthus, Strophocactus wittii.
* NOTE: For 10 species, models could not be properly compared due to the inability to calculate AICc. These were excluded from analyses. Species: Discocactus zehntneri, Espostoa frutescens, Gymnocalycium castellanosii, Mammillaria carmenae, M. moelleriana, M. peninsularis, Opuntia repens, Rapicactus subterraneus, Rhipsalis crispata, R. olivifera.
* removePoorModels.R: checks for species with AUC below or equal to 0.5, and removes those.
* NOTE: For 35 species, AUC was below or equal to 0.5. These were excluded from analyses. Species: Austrocactus bertinii, Cereus fricii, C. phatnospermus, Cleistocactus baumannii, Cylindropuntia alcahes, C. waltoniorum, C. x viridiflora, Epithelantha bokei, Eriosyce iquiquensis, Espostoopsis dybowskii, Ferocactus uncinatus, Grusonia schottii, Gymnocalycium hyptiacanthum, G. quehlianum, G. schickendantzii, Lepismium waringianum, Maihueniopsis hickenii, Melocactus bahiensis, Opuntia anacantha, O. aureispina, O. cespitosa, O. feracantha, O. rioplatensis, O. robinsonii, O. x debreczyi, Parodia langsdorfii, P. linkii, P. microsperma, Pilosocerereus flavipulvinatus, P. piauhyensis, Praecereus saxicola, Rhipsalis dissimilis, R. hileiabaiana, Sclerocactus papyracanthus, Tephrocactus articulatus.

Species-level analyses
* rankTopModels.R: ranks performance by AICc of variable sets (A, AS, AE, ASE) with or without GCMs considered.
* compareModelPerformance.R: analyzes model performance.
* rangeChanges.R: analyzes range changes.

Richness-level analyses
* stackMaps.R: creates richness maps and corresponding summary maps.
* expertMap.R: creates richness map based on IUCN SSC CSSG Global Cactus Assessment expert maps and calculates correlations with model-based maps.
* bioregions.R: summarizes richness changes by bioregion following Calvente et al. (2023), Appendix S11.