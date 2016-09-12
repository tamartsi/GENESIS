\name{admixMapMixedModel}
\alias{admixMapMixedModel}

\title{admixMapMixedModel}
\description{
Run admixture analyses
}
\usage{
admixMapMixedModel(admixDataList, cholSigmaInv, outcome, covar.vec = NULL, scan.exclude = NULL, snpStart = NULL, snpEnd = NULL, block.size = 5000, verbose = TRUE)
}

\arguments{
  \item{admixDataList}{named list of \code{\link{GenotypeData}} objects}
  \item{cholSigmaInv}{cholSigmaInv matrix, output from \code{\link{estVarComp}}}
  \item{outcome}{outcome variable name}
  \item{covar.vec}{vector of covariate variable names}
  \item{scan.exclude}{scanIDs to exclude}
  \item{snpStart}{index of first SNP to run}
  \item{snpEnd}{index of last SNP to run}
  \item{block.size}{integer number of SNPs to loop over}
  \item{verbose}{logical indicating whether to report progress information}
}


\value{data frame with admixture mapping results}

\author{Matt Conomos}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{estVarComp}}}
