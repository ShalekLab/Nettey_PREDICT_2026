# This function will filter variable genes from a Seurat object or a
# user-provided list of genes. 

require(Seurat)

FeaturesFilter <- function(seurat = NULL, features = NULL, 
                           removeGenes = "", removeIG = F, removeTR = F, 
                           removeMito = F, removeRibo = F, removeHLA = F){
  library(dplyr)
  
  # get variable genes
  if(is.null(features) & !is.null(seurat)){
    message("Defaulting to Seurat object variable features for filtering.")
    vargenes <- VariableFeatures(seurat)
  } else {
    vargenes <- features
  }
  
  
  # remove specified genes (only if genes are specified and found in the list)
  genesToRemove <- which(vargenes %in% removeGenes)
  
  if (!identical(genesToRemove,integer(0))) {
    vargenes <- vargenes[-which(vargenes %in% removeGenes)]
  }
  
  # remove non-constant immunoglobulin genes
  if (removeIG) {
    IGgenes <- vargenes[grepl("^IG[HKL][VJ]",vargenes)]
    
    # # add additional, novel genes
    # AC244452.3 = IGHVII-44-2D
    # AC135068.8 = IGHV1OR15-1
    # AC135068.9 = IGHV1OR15-2
    # AC135068.2 = IGHV4-4
    # AC134981.1 = IGHV1OR15-4
    # AC233755.1 = IGHV4-38-2
    # AC233755.2 = IGHV3-43D
    # AC141272.1 = Unknown
    # AC140658.1 = Unknown
    # AC136428.1 = Unknown
    # AC136616.2 = Unknown

    IGgenes <- c(IGgenes, 
                 "AC244452.3", "AC135068.8", "AC135068.9", 
                 "AC135068.2", "AC134981.1", "AC233755.1", 
                 "AC233755.2", "AC141272.1", "AC140658.1",
                 "AC136428.1", "AC136616.2")
    
    vargenes <- vargenes[!(vargenes %in% IGgenes)]
  }
  
  # remove non-constant T cell receptor genes
  if (removeTR) {
    TRgenes <- vargenes[grepl("^TR[ABDG][VJ]",vargenes)]
    
    vargenes <- vargenes[!(vargenes %in% TRgenes)]
  }
  
  # remove mitochondrial genes
  if (removeMito) {
    mitogenes <- vargenes[grep("^MT-",vargenes)]
    
    vargenes <- vargenes[!(vargenes %in% mitogenes)]
  }
  
  # remove ribosomal protein genes
  if (removeRibo) {
    ribogenes <- vargenes[grepl("^RP[SL]", vargenes)]
    
    vargenes <- vargenes[!(vargenes %in% ribogenes)]
  }
  
  # remove HLA genes
  if (removeHLA) {
    HLAgenes <- vargenes[grepl("^HLA",vargenes)]
    
    vargenes <- vargenes[!(vargenes %in% HLAgenes)]
  }
  
  return(vargenes)
}
