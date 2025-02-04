% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarietyDevPipe.R
\name{choseTrialEntries}
\alias{choseTrialEntries}
\title{choseTrialEntries function}
\usage{
choseTrialEntries(trialType, nToSelect, bsd)
}
\arguments{
\item{trialType}{String the trial type the entries are going to (inventory)}

\item{nToSelect}{The number of IDs to return that have the highest BLUPs}

\item{bsd}{The breeding scheme data list}
}
\value{
List with entries vector of variety IDs selected and updated bsd
}
\description{
Function to select which varieties to advance to the next stage assuming
they are IID (i.e., estimating genotypic values), by truncation selection
}
\details{
Accesses all data in phenoRecords to pick the highest among
candidates with enough seed inventory.
}
\examples{
entries <- choseTrialEntries(phenoRecords, trialType, nToSelect, bsd)]
bsd <- entries$bsd
entries <- entries$entries

}
