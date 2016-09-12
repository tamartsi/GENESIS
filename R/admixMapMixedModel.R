admixMapMixedModel <- function(admixDataList,
                               cholSigmaInv,
                               outcome,
                               covar.vec = NULL,
                               scan.exclude = NULL,
                               snpStart = NULL,
                               snpEnd = NULL,
                               #impute.local = TRUE,
                               block.size = 5000,
                               verbose = TRUE){

  # if admixDataList is GenodtypeData (one file), convert to a list
  if(class(admixDataList) == "GenotypeData"){
    admixDataList <- list(admixDataList)
  }

  # how many ancestry components to be tested
  v <- length(admixDataList)

  # if admixDataList doesn't have names, assign them
  if(is.null(names(admixDataList))){
    names(admixDataList) <- paste("Anc",1:v,sep="")
  }

  # set snpStart and snpEnd
  if(is.null(snpStart)){
    snpStart <- 1
  }
  if(is.null(snpEnd)){
    snpEnd <- nsnp(admixDataList[[1]])
  }
  
  # get scanIDs
  scanID <- getScanID(admixDataList[[1]])

  # get chromosome information
  chr <- getChromosome(admixDataList[[1]], index=snpStart:snpEnd)

  # check that items match in all elements of admixDataList
  if(v > 1){
    for(i in 2:v){
      # scanIDs
      if(!all(scanID == getScanID(admixDataList[[i]]) )){
        stop("scanIDs do not match for all elements of admixDataList")
      }
      # chromosomes
      if(!all(chr == getChromosome(admixDataList[[i]], index=snpStart:snpEnd) )){
        stop("Chromosomes do not match for all elements of admixDataList")
      }
    }
  }

  # set which samples to keep
  keep <- rep(TRUE, nscan(admixDataList[[1]]))
  # samples excluded from entire analysis
  if(!is.null(scan.exclude)){
    keep <- keep & !(scanID %in% scan.exclude)
  }
  
  # Sex chromosome checks
  Xcode <- rep(NA,v); Ycode <- rep(NA,v)
  for(i in 1:v){
    Xcode[i] <- XchromCode(admixDataList[[i]])
    Ycode[i] <- YchromCode(admixDataList[[i]])
  }
  if(any(Xcode %in% chr)){
    # check for matching
    if(length(unique(Xcode)) != 1){
      stop("X chromosome codes do not match for all elements of admixDataList")
    }
    # check for sex variable
    for(i in 1:v){
      if(!hasSex(admixDataList[[i]])){
        stop("Sex values for the samples are required for X chromosome SNPs")
      }
    }
  }
  if(any(Ycode %in% chr)){
    # check for matching
    if(length(unique(Ycode)) != 1){
      stop("Y chromosome codes do not match for all elements of admixDataList")
    }
    if(!all(chr == Ycode[1])){
      stop("Y chromosome must be analyzed separately")
    }
    # check for sex variable
    for(i in 1:v){
      if(!hasSex(admixDataList[[i]])){
        stop("Sex values for the samples are required for Y chromosome SNPs")
      }
      # check for matching
      if(!all( getSex(admixDataList[[1]]) == getSex(admixDataList[[i]]) )){
        stop("Sex values for the samples do not match for all elements of admixDataList (required for Y chr)")
      }
    }    
    # only keep males
    keep <- keep & (getSex(admixDataList[[1]]) == "M")
  }  
    
  
  # read in outcome and covariate data
  if(verbose) message("Reading in Phenotype and Covariate Data...")
  if(!is.null(covar.vec)){
    cvnames <- unique(unlist(strsplit(covar.vec,"[*:]")))
    # read in data
    dat <- as.data.frame(getScanVariable(admixDataList[[1]], c(outcome,cvnames)))
    # check for matching
    if(v > 1){
      for(i in 2:v){
        if(!identical( dat, as.data.frame(getScanVariable(admixDataList[[i]], c(outcome,cvnames))) )){
          stop("Outcome or covariate data does not match for all elements of admixDataList")
        }
      }
    }
    # identify samples with any missing data
    keep <- keep & apply(dat,1,function(x){ all(!is.na(x)) })
    # subset samples
    dat <- as.data.frame(dat[keep,])
    # outcome
    Y <- dat[,outcome]
    # create design matrix
    model.formula <- as.formula(paste(paste(outcome,"~"), paste(covar.vec,collapse="+")))
    W <- model.matrix(model.formula, data=dat)    
    
  }else{
    # read in data
    dat <- getScanVariable(admixDataList[[1]],outcome)
    if(v > 1){
      for(i in 2:v){
        if(!identical( dat, getScanVariable(admixDataList[[i]], outcome) )){
          stop("Outcome data does not match for all elements of admixDataList")
        }
      }
    }
    # identify samples with any missing data
    keep <- keep & !is.na(dat)    
    # outcome
    Y <- dat[keep]
    # design matrix
    W <- matrix(1,nrow=length(Y),ncol=1)
  }
  k <- ncol(W)
  scanID <- scanID[keep]

  # which samples to remove from cholSigmaInv
  if(!all(scanID %in% colnames(cholSigmaInv))){
    stop("All of the included Samples must be in the cholSigmaInv matrix")
  }
  chol.idx <- which(!(colnames(cholSigmaInv) %in% scanID))
  cholSigmaInv <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)
  

  # sample size, assuming no missing genotypes
  n <- length(scanID)
  if(verbose) message("Running analysis with ", n, " Samples")
  
  # number of SNPs in the segment
  nsnp.seg <- snpEnd - snpStart + 1
  # determine number of SNP blocks
  nblocks <- ceiling(nsnp.seg/block.size)

  # set up results matrix
  # add in frequency of each ancestry at the SNP
  nv <- c("snpID","chr","n")
  if(v > 1){
    for(i in 1:v){
      nv <- append(nv, paste(names(admixDataList)[i],c(".freq",".Est",".SE"), sep=""))
    }
    nv <- append(nv, c("Joint.Stat", "Joint.pval"))
  }else{
    nv <- append(nv, c("freq", "Est", "SE", "Stat", "pval"))
  }
  res <- matrix(NA, nrow=nsnp.seg, ncol=length(nv), dimnames=list(NULL, nv))
    
  
  # chromosome
  res[,"chr"] <- chr
  
  # since we haven't done any loops here yet:
  keep.previous <- rep(TRUE, n)
  
  if(verbose) message("Beginning Calculations...")
  # loop through blocks
  for(b in 1:nblocks){
    
    # keep track of time for rate reporting
    startTime <- Sys.time()
    
    snp.start.pos <- snpStart + (b-1)*block.size
    nsnp.block <- block.size
    if(snp.start.pos + nsnp.block > snpEnd){
      nsnp.block <- snpEnd - snp.start.pos + 1
    }
    snp.end.pos <- snp.start.pos + nsnp.block - 1
    
    bidx <- ((b-1)*block.size+1):((b-1)*block.size+nsnp.block)
    
    # get local ancestry for the block
    local <- array(NA, dim=c(n,nsnp.block,v)) # indices: scan, snp, ancestry
    for(i in 1:v){
      tmplocal <- getGenotype(admixDataList[[i]], snp=c(snp.start.pos, nsnp.block), scan=c(1,-1), transpose=TRUE)
      # subset
      tmplocal <- as.matrix(tmplocal)[keep, , drop=F]
      # put into the array
      local[,,i] <- tmplocal; rm(tmplocal)
    }

    # ancestral frequency
    freq <- 0.5*colMeans(local, dims=1, na.rm=T) # matrix:  rows are SNPs, columns are ancestries
    # for X chr
    if(XchromCode(admixDataList[[1]]) %in% chr[bidx]){
      # which are on X chr
      Xidx <- chr[bidx] == XchromCode(admixDataList[[1]])
      # males
      m <- (getSex(admixDataList[[1]]) == "M")[keep]
      f <- (getSex(admixDataList[[1]]) == "F")[keep]
      # calculate allele freq for X
      freq[Xidx,] <- (0.5*colSums(local[m,Xidx,], dims=1, na.rm=T) + colSums(local[f,Xidx,], dims=1, na.rm=T)) / (colSums(!is.na(local[m,Xidx,]), dims=1) + 2*colSums(!is.na(local[f,Xidx,]), dims=1) )
    }
    
    # add to output
    if(v > 1){
      for(i in 1:v){
        res[bidx,paste(names(admixDataList)[i],".freq", sep="")] <- freq[,i]
      }
    }else{
      freq <- as.vector(freq)
      res[bidx,"freq"] <- freq
    }
    
    # impute missing local ancestry values
    #if(impute.local){      
    #  miss.idx <- which(is.na(local))
    #  if(length(miss.idx) > 0){
    #    freq.idx <- ceiling(miss.idx/n)
    #    local[miss.idx] <- 2*freq[freq.idx]  # double check that this indexes in the array correctly (pretty sure it does)
    #  }
    #}
    # this is imputing as population average.  would it be better to impute to individual's global ancestry?
    ### no missing data currently - worry about this later ###
    
    # check for missingness
    check <- rowSums(is.na(local), dims=1)    
    if (!all(check %in% c(0, v*nsnp.block))) {
      stop("sporadic missing in block size > 1")
    }
    
    keep.geno <- check == 0
    
    # get rid of missing for this block
    # can probably put this in an if statement
    local <- as.array(local)[keep.geno, , , drop=F]
    
    # sample size
    n <- sum(keep.geno)
    res[bidx, "n"] <- n
    
    if (b == 1 | !all(keep.previous == keep.geno)){
      # calculate matrices
      
      if (b > 1) warning("recalculating matrices!")
      #if(verbose) message("Pre-Computing some Matrices for Analysis...")      
      
      # subsetting
      W.block <- W[keep.geno, , drop=F]
      Y.block <- Y[keep.geno]
      
      # here we have to subset the matrix, not the inverse
      # this is a fancy way of getting the inverse of the subset without having to get the original matrix
      chol.idx <- which(!(colnames(cholSigmaInv) %in% scanID[keep.geno]))
      C.block <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)
      
      # W: covariate matrix
      # C: cholesky decomposition of sigma inverse
      # sigma inverse: inverse phenotype covariance matrix
      CW <- crossprod(C.block, W.block)
      Mt <- C.block - tcrossprod(tcrossprod(C.block,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW) # this is a matrix used to adjust phenotype and local ancestry for fixed effect covariates AND "decorrelating" the phenotype and genotype -- ie adjusting out covariance structure given by sigma matrix.
      Ytilde <- crossprod(Mt,Y.block) # so ytilde is the phenotype adjusted for the covariates/correlation structure
      sY2 <- sum(Ytilde^2)
    }
    
        
    # perform regressions
    if(v == 1){
      local <- local[,,1]
      Xtilde <- crossprod(Mt,local) # adjust local for correlation structure and fixed effects
      XtX <- colSums(Xtilde^2) # vector of X^T SigmaInv X (for each SNP)
      # filter monomorphic SNPs
      XtX[which(freq == 0 || freq == 1)] <- NA
      beta <- as.vector(crossprod(Xtilde,Ytilde)/XtX)
      Vbeta <- (sY2/XtX - beta^2)/(n - k - 1) # RSS/XtX
      Stat <- beta^2/Vbeta

      res[bidx,"Est"] <- beta
      res[bidx,"SE"] <- sqrt(Vbeta)
      res[bidx,"Stat"] <- Stat
      res[bidx,"pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)      

    }else{
      Joint.Stat <- rep(NA, ncol(local))
      Est <- matrix(NA, nrow=ncol(local), ncol=v)
      SE <- matrix(NA, nrow=ncol(local), ncol=v)
      for(g in 1:ncol(local)){
        if(identical(local[,g,], local[,(g-1),])){
          Joint.Stat[g] <- Joint.Stat[g-1]
          Est[g,] <- Est[g-1,]
          SE[g,] <- SE[g-1,]                     
          next
        }
        # filter monomorphic or missing SNPs
        if(any(freq[g,]==1) || sum(freq[g,]==0)){ next }
        Xtilde <- crossprod(Mt,local[,g,])
        XtX <- crossprod(Xtilde)
        XtXinv <- tryCatch( chol2inv(chol(XtX)), error=function(e){TRUE})
        # check that the error function above hasn't been called (which returns TRUE instead of the inverse matrix)
        if(is.logical(XtXinv)){ next }
        XtY <- crossprod(Xtilde, Ytilde)
        betas <- crossprod(XtXinv, XtY)
        RSS <- as.numeric((sY2 - crossprod(XtY,betas))/(n - k - v))
        Vbetas <- XtXinv*RSS

        Est[g,] <- betas
        SE[g,] <- sqrt(diag(Vbetas))
        Joint.Stat[g] <- tryCatch( crossprod(betas,crossprod(XtX,betas))/RSS, error=function(e){NA})
      } # g loop

      # collect results
      for(i in 1:v){
        res[bidx,paste(names(admixDataList)[i],".Est", sep="")] <- Est[,i]
        res[bidx,paste(names(admixDataList)[i],".SE", sep="")] <- SE[,i]
      }
      res[bidx,"Joint.Stat"] <- Joint.Stat
      res[bidx,"Joint.pval"] <- pchisq(Joint.Stat, df=v, lower.tail=FALSE)
    } # else

    endTime <- Sys.time()
    rate <- format(endTime - startTime, digits=4)

    keep.previous <- keep.geno

    if(verbose) message(paste("Block", b, "of", nblocks, "Completed -", rate))
  } # end block loop
  
  # results data frame
  res <- as.data.frame(res)
  
  # add in snpID
  res$snpID <- getSnpID(admixDataList[[1]], snpStart:snpEnd)
  
  return(res)
}


.subsetCholSigmaInv <- function(cholSigmaInv, chol.idx) {
  if(length(chol.idx) > 0){
    # subset cholSigmaInv
    SigmaInv <- tcrossprod(cholSigmaInv)
    for(i in sort(chol.idx, decreasing=TRUE)){
      SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
    }
    cholSigmaInv <- t(chol(SigmaInv))
  }
  
  cholSigmaInv
}
