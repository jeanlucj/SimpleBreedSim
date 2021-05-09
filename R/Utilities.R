#' initializeProgram function
#'
#' Function to initialize the simulation
#' Read the founder genetic architecture and population genetic parameters
#' Create simulation data structures
#'
#' @param founderFile String name of the text file with parameters governing
#' the founders and the genetic architecture
#' @param schemeFile String name of the text file with parameters governing
#' the variety development pipeline and selection functions
#' @return A named list containing the `AlphaSimR` SP file, bsd a list of
#' breeding scheme data, breedingPop the breeding population, varietyCand
#' an empty pop for variety candidates, phenoRecords an empty tibble for
#' phenotype records, and inventory an empty tibble for variety candidates
#'
#' @details Call this function at the beginning of the simulation
#'
#' @examples
#' bsd <- initializeProgram("FounderCtrlFile.txt")
#'
#' @export
initializeProgram <- function(founderFile, schemeFile){
  # Read parameters to create founders
  parmNames <- c("nChr", "nFounders", "effPopSize", "quickHaplo", "segSites", "nQTL", "nSNP", "genVar", "gxeVar", "gxyVar", "gxlVar", "gxyxlVar", "meanDD", "varDD", "relAA")
  bsd <- readControlFile(founderFile, parmNames)

  # Create haplotypes for founder population of outbred individuals
  if (bsd$quickHaplo){
    founderHap <- quickHaplo(nInd=bsd$nFounders, nChr=bsd$nChr, segSites=bsd$segSites)
  } else{
    founderHap <- runMacs2(nInd=bsd$nFounders, nChr=bsd$nChr, segSites=bsd$segSites, Ne=bsd$effPopSize)
  }

  # Global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderHap)
  SP$restrSegSites(minQtlPerChr=1, minSnpPerChr=10, overlap=FALSE)
  # Additive, dominance, and epistatic trait architecture
  SP$addTraitADE(nQtlPerChr=bsd$nQTL, var=bsd$genVar, meanDD=bsd$meanDD, varDD=bsd$varDD, relAA=bsd$relAA, useVarA=FALSE)
  # Observed SNPs per chromosome
  SP$addSnpChip(bsd$nSNP)

  # Create the founders
  breedingPop <- newPop(founderHap, simParam=SP)

  # Read parameters about trial types
  parmNames <- c("nStages", "stageNames", "nReps", "nLocs", "errVars", "minInventory", "inventoryOut")
  bsdNew <- readControlFile(schemeFile, parmNames)
  bsd <- c(bsd, bsdNew)

  # Add miscelaneous data structures
  bsd$year <- 0
  bsd$nextTrialID <- 1

  return(bsd)
}

#' readControlFile function
#'
#' Function to read a text control file
#'
#' The text file should be organized as follows
#' 1. Any text after a comment symbol # will be ignored
#' 2. Control parameter names should be on their own line
#' 3. Parameter values should be on the following line. If multiple parameter values are needed they should be separated by white space but on the same line
#' @param fileName The name of the text file to be read. Must include the path to the file
#' @param parmNames A string vector with the names of the control parameters that will be searched in the text file
#' @return A named list of the parameter values read from the control file
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' params <- readControlFile("./inputDir/ctrlFile.txt", c("nStages", "nParents", "nCrosses"))
#'
#' @export
readControlFile <- function(fileName, parmNames){
  ctrlLines <- readLines(fileName)
  ctrlLines <- sapply(ctrlLines, function(st) strsplit(st, "#", fixed=T)[[1]][1])
  getParm <- function(parmName){
    parmRow <- grep(parmName, ctrlLines)+1
    parms <- unlist(strsplit(ctrlLines[parmRow], "[[:space:]]"))
    names(parms) <- NULL
    parmsNum <- suppressWarnings(as.numeric(parms))
    if (!any(is.na(parmsNum))) parms <- parmsNum
    return(parms)
  }
  parms <- lapply(parmNames, getParm)
  names(parms) <- parmNames
  return(parms)
}
