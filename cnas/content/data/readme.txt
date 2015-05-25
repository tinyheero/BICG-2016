OncoSNP Configuration files

Standard:
 
-paramsfile hyperparameters.dat (standard settings)
-trainingstatesfile trainingStates.dat (do not modify)
-tumourstatesfile tumourStates.dat (contain copy number states up to CN6)
-levelsfile levels-610.dat (for InfiniumHD arrays)

Options:

-tumourstatesfile tumourStates8.dat (contain copy number states up to CN8)
-tumourstatesfile tumourStatesCLL.dat (restricted copy number set for CLL analysis)
-tumourstatesfile tumourStatesNormal.dat (restricted copy number set for germline analysis)

-levelsfile levels.dat (for Infinium I/II arrays)

-levelsfile levels-tQN.dat (for data preprocessed using tQN -- experimental)

-levelsfile levels-affy.dat (settings for Affymetrix array analysis after preprocessing using PennCNV-Affy toolkit -- experimental)

-paramsfile hyperparameters-affy.dat (settings for Affymetrix array analysis -- experimental)

-paramsfile hyperparameters-affy.dat (settings for intra-tumour heterogeneity mode only)

