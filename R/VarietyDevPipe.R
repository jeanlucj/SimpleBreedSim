
#' runVDPtrial function
#'
#' Function to run a variety development pipeline stage
#'
#' @param trialType String name of the kind of trial to run. Must link back to
#' a trial type with parameters in bsd
#' @param entries A vector of entry IDs
#' @param bsd The breeding scheme data list
#'
#' @return Updated bsd
#'
#' @details Add records to phenoRecords and update the inventory
#'
#' @examples
#' params <- runVDPtrial("CET", c("1", "2", "3", "4"), bsd)
#'
#' @export
runVDPtrial <- function(trialType, entries, bsd){
  # Calculate the error variance resulting from the number of locations and reps
  # Note that the number of years for a trial is always 1
  nL <- bsd$nLoc[trialType]; nR <- bsd$nRep[trialType]
  eV <- bsd$errVars[trialType]
  impliedVarE <- bsd$gxlVar/nL + gxyVar + gxyxlVar/nL + eV/nL/nR
  # Phenotypic evaluation of experimental lines
  newPheno <- AlphaSimR::setPheno(bsd$varietyCandidates[entries],
                                  varE=impliedVarE, simParam=SP, onlyPheno=T)
  # Add the information to phenoRecords
  bsd <- addPhenoRecords(newPheno, trialType, bsd)
  # Manage the inventory
  bsd$inventory[entries] <- bsd$inventory[entries] + bsd$seedProduced[trialType]
  # Manage the trialID number
  bsd$nextTrialID <- bsd$nextTrialID + 1

  return(bsd)
}

#' choseTrialEntries function
#'
#' Function to select which varieties to advance to the next stage assuming
#' they are IID (i.e., estimating genotypic values), by truncation selection
#'
#' @param trialType String the trial type the entries are going to (inventory)
#' @param nToSelect The number of IDs to return that have the highest BLUPs
#' @param bsd The breeding scheme data list
#' @return List with entries vector of variety IDs selected and updated bsd
#' @details Accesses all data in phenoRecords to pick the highest among
#' candidates with enough seed inventory.
#'
#' @examples
#' entries <- choseTrialEntries(phenoRecords, trialType, nToSelect, bsd)]
#' bsd <- entries$bsd
#' entries <- entries$entries
#'
#' @export
choseTrialEntries <- function(trialType, nToSelect, bsd){
  nCandidates <- bsd$phenoRecords$id %>% unique %>% length
  if (nrow(bsd$phenoRecords) > nCandidates){ # There is some replication
    crit <- iidPhenoEval(bsd$phenoRecords)
  } else{
    crit <- bsd$phenoRecords$pheno
    names(crit) <- bsd$phenoRecords$id
  }
  minToEnter <- bsd$seedNeeded[trialType]
  hasInventory <- bsd$inventory[bsd$inventory >= minToEnter] %>% names
  if (length(hasInventory) < nToSelect){
    stop("There are too few variety candidates with enough inventory")
  }
  crit <- crit[hasInventory]
  entries <- crit[crit %>% order(decreasing=T)[1:nToSelect]] %>% names
  # Manage the inventory
  bsd$inventory[entries] <- bsd$inventory[entries] - bsd$seedNeeded[trialType]
  return(list(entries=entries, bsd=bsd))
}

#' makeVarietyCandidates function
#'
#' @param bsd List of breeding program data
#'
#' @return Updated bsd with new variety candidates in bsd$varietyCandidates
#' @details Here, creates DHs evenly distributed among the breeding progeny
#' in the last generation
#' @examples
#' bsd <- makeVarietyCandidates(bsd)]
#'
#' @export
makeVarietyCandidates <- function(bsd){
  nBP <- AlphaSimR::nInd(bsd$breedingPop)
  lastGen <- nBP - (bsd$nBreedingPopProg-1):0
  lastGenBreedProg <- bsd$breedingPop[lastGen]
  nDH <- bsd$nStartVarietyCand %/% bsd$nBreedingPopProg
  newVarCand <- makeDH(lastGenBreedProg, nDH=nDH, simParam=bsd$SP)
  if (bsd$nStartVarietyCand %% bsd$nBreedingPopProg > 0){
    nExtra <- bsd$nStartVarietyCand %% bsd$nBreedingPopProg
    whichPar <- sample(bsd$nBreedingPopProg, nExtra)
    newVarCand <- c(newVarCand,
                    makeDH(lastGenBreedProg[whichPar], nDH=1, simParam=bsd$SP)
    )
  }
  bsd$varietyCandidates <- c(bsd$varietyCandidates, newVarCand)
  return(bsd)
}

#' iidPhenoEval function
#'
#' Function to estimate variety candidate BLUPs assuming they are IID
#'
#' @param phenoRecords Tibble of phenotypic observations on variety candidates.
#' @return Named real vector of the BLUPs of all individuals in phenoRecords,
#' (names are the individual ids), with appropriate weights by error variance
#' @details Given all the phenotypic records calculate the best prediction of
#' the genotypic value for each individual using all its records
#'
#' @examples
#' iidBLUPs <- iidPhenoEval(phenoRecords)
#'
#' @export
iidPhenoEval <- function(phenoRecords){
  require(lme4)
  # Make error variances into weights
  phenoRecords <- phenoRecords %>% dplyr::mutate(wgt=1/errVar)
  fm <- lmer(pheno ~ (1|id), weights=wgt, data=phenoRecords)
  blup <- as.matrix(ranef(fm)[[1]])[,1] # Make into matrix to get names
  names(blup) <- (names(blup) %>% strsplit(":", fixed=T) %>% unlist %>%
                    matrix(nrow=2))[1,]
  # Ensure output has variation: needed for optimal contributions
  if (sd(blup) == 0){
    namesBlup <- names(blup)
    blup <- tapply(phenoRecords$pheno, phenoRecords$id, mean)
    names(blup) <- namesBlup
  }
  return(blup)
}

#' addPhenoRecords function
#'
#' Function to add new phenotypes from a VDP trial to phenoRecords
#'
#' @param pheno Named vector of phenotypic values
#' @param trialType String name of the trial type to get the error variance
#' @param bsd The breeding scheme data list
#' @return Updates bsd with new phenotypes in phenoRecords
#'
#' @details Call this function from the runVDPstage
#'
#' @examples
#' newPhenoRec <- addPhenoRecords(newPheno, "CET", bsd)
#'
#' @export
addPhenoRecords <- function(pheno, trialType, bsd){
  newRec <- tibble(year=bsd$year, trialID=bsd$nextTrialID, trialType=trialType,
                   id=rownames(pheno), pheno=pheno, error=bsd$errVars[trialType])
  bsd$phenoRecords <- bsd$phenoRecords %>% bind_rows(newRec)
  return(bsd)
}
