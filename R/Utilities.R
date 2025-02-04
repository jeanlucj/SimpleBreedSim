#' initializeProgram function
#'
#' Function to initialize the simulation.
#' Read the founder genetic architecture and population genetic parameters.
#' Create simulation data structures.
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
  parmNames <- c("nChr", "nFounders", "effPopSize", "quickHaplo", "segSites", "nQTL", "nSNP", "genVar", "gxeVar", "gxyVar", "gxlVar", "gxyxlVar", "meanDD", "varDD", "relAA", "quickHaplo")
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
  bsd$breedingPop <- AlphaSimR::newPop(founderHap, simParam=SP)

  # Read parameters about trial types
  parmNames <- c("nStages", "stageNames", "nReps", "nLocs", "errVars",
                 "seedNeeded", "seedProduced", "optiContEffPop",
                 "nBreedingPopProg", "nStartVarietyCand")
  bsdNew <- readControlFile(schemeFile, parmNames)
  bsd <- c(bsd, bsdNew)

  # Add miscelaneous data structures
  bsd$year <- 0
  bsd$nextTrialID <- 1
  bsd <- calcDerivedParms(bsd)
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

#' calcDerivedParms function
#'
#' Once you have read in parameters from a control file, or set them yourself, there are still a few derived parameters that are needed.  This function calculates them.
#'
#' @param bsd List of breeding program data
#' @return A list bsd that extends the input with a few derived parameters
#'
#' @details This function is only called internally by other functions used to specify the pipeline
#'
calcDerivedParms <- function(bsd){
  # Function to check if a parameter has no value
  nv <- function(parm){
    is.null(parm) | length(parm) == 0
  }

  # Some parms have to be logical
  makeLogical <- function(parm){
    if (nv(parm)) parm <- FALSE else parm <- as.logical(parm)
    if (is.na(parm)) parm <- FALSE
    return(parm)
  }

  # Prevent errors having to do with inconsistent parameters
  if (bsd$nSNP + bsd$nQTL >= bsd$segSites){
    print("The number of segregating sites (segSites) has to be greater than the number of SNPs (nSNP) and the number of QTL (nQTL). segSites set 10% bigger than nSNP + nQTL")
    bsd$segSites <- ceiling((bsd$nSNP + bsd$nQTL) * 1.1)
  }

  bsd$quickHaplo <- makeLogical(bsd$quickHaplo)

  # Genetic architecture defaults
  if (nv(bsd$meanDD)) bsd$meanDD <- 0
  if (nv(bsd$varDD)) bsd$varDD <- 0
  if (nv(bsd$relAA)) bsd$relAA <- 0

  # Check that these vectors are of the right length
  rightLength <- function(objName) length(get(objName, bsd)) == bsd$nStages
  vecNames <- c("stageNames", "nReps", "nLocs", "errVars",
                "seedNeeded", "seedProduced")
  rl <- sapply(vecNames, rightLength)
  if (any(!rl)){
    stop(paste("These vectors do not have the right length:",
               paste(vecNames[!rl], collapse=" ")))
  }

  # Defaults for GxE variance
  if (any(nv(bsd$gxyVar), nv(bsd$gxlVar), nv(bsd$gxyxlVar))){
    if (!nv(bsd$gxeVar)){
      if (nv(bsd$gxyVar)) bsd$gxyVar <- bsd$gxeVar / 3
      if (nv(bsd$gxlVar)) bsd$gxlVar <- bsd$gxeVar / 3
      if (nv(bsd$gxyxlVar)) bsd$gxyxlVar <- bsd$gxeVar / 3
    } else{
      if (nv(bsd$gxyVar)) bsd$gxyVar <- 0
      if (nv(bsd$gxlVar)) bsd$gxlVar <- 0
      if (nv(bsd$gxyxlVar)) bsd$gxyxlVar <- 0
    }
  }

  # Make sure everything has names
  names(bsd$nReps) <- names(bsd$nLocs) <- names(bsd$errVars) <-
    names(bsd$seedNeeded) <- names(bsd$seedProduced)<- bsd$stageNames

  return(bsd)
}
