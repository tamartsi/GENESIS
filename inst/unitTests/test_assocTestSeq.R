library(SeqVarTools)
library(Biobase)
library(SNPRelate)
library(logistf)
library(gdsfmt)

.testObject <- function(){
    showfile.gds(closeall=TRUE, verbose=FALSE)
    
    gdsfile <- seqExampleFileName("gds")
    gds.seq <- seqOpen(gdsfile)
    
    data(pedigree)
    pedigree <- pedigree[match(seqGetData(gds.seq, "sample.id"), pedigree$sample.id),]
    pedigree$outcome <- rnorm(nrow(pedigree))
    pedigree$status <- rbinom(nrow(pedigree), 1, 0.4)
    
    SeqVarData(gds.seq, sampleData=AnnotatedDataFrame(pedigree))
}

.testKing <- function(gds){
    ibd <- snpgdsIBDKING(gds, verbose=FALSE)
    kc <- ibd$kinship
    rownames(kc) <- ibd$sample.id
    colnames(kc) <- ibd$sample.id
    kc
}

.testGRM <- function(seqData, ...){
    kinship <- .testKing(seqData)
    mypcair <- pcair(seqData, kinMat=kinship, divMat=kinship, verbose=FALSE, ...)
    mypcrel <- pcrelate(seqData, pcMat=mypcair$vectors[,1:2], training.set=mypcair$unrels, verbose=FALSE, ...)
    pcrelateMakeGRM(mypcrel)
}

.testNullModel <- function(seqData, type="mixed", binary=FALSE, ...){
    scanData <- sampleData(seqData)
    if (binary) {
        outcome <- "status"
        family <- binomial
    } else {
        outcome <- "outcome"
        family <- gaussian
    }

    if (type == "mixed") {
        grm <- .testGRM(seqData, ...)
        fitNullMM(scanData, outcome=outcome, covars="sex", covMatList=grm, family=family, verbose=FALSE, ...)
    } else {
        fitNullReg(scanData, outcome=outcome, covars="sex", family=family, verbose=FALSE, ...)
    }
}

.testBurden <- function(seqData) {
    rowSums(altDosage(seqData), na.rm=TRUE)
}


test_projection_mm <- function(){
    seqData <- .testObject()

    nullmod <- .testNullModel(seqData, type="mixed")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Score")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Wald")
    proj <- GENESIS:::.calculateProjection(nullmod, test="SKAT", burden.test=NULL)

    seqClose(seqData)
}

test_projection_reg <- function(){
    seqData <- .testObject()

    nullmod <- .testNullModel(seqData, type="reg")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Score")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Wald")
    proj <- GENESIS:::.calculateProjection(nullmod, test="SKAT", burden.test="")

    seqClose(seqData)
}

test_projection_logistic <- function(){
    seqData <- .testObject()

    nullmod <- .testNullModel(seqData, type="reg", binary=TRUE)
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Score")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Wald")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Firth")
    proj <- GENESIS:::.calculateProjection(nullmod, test="SKAT", burden.test="")

    seqClose(seqData)
}

# fails unpredictably
## test_projection_gmmat <- function(){
##     seqData <- .testObject()

##     nullmod <- .testNullModel(seqData, type="mixed", binary=TRUE)
##     proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test="Score")
##     proj <- GENESIS:::.calculateProjection(nullmod, test="SKAT", burden.test="")

##     seqClose(seqData)
## }

test_burden_score <- function(){
    seqData <- .testObject()
    burden <- .testBurden(seqData)

    burden.test <- "Score"
    nullmod <- .testNullModel(seqData, type="reg")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test=burden.test)
    bt <- GENESIS:::.runBurdenTest(burden, proj, burden.test=burden.test)

    seqClose(seqData)
}

test_burden_wald <- function(){
    seqData <- .testObject()
    burden <- .testBurden(seqData)

    burden.test <- "Wald"
    nullmod <- .testNullModel(seqData, type="reg")
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test=burden.test)
    bt <- GENESIS:::.runBurdenTest(burden, proj, burden.test=burden.test)

    dat <- cbind(pData(sampleData(seqData)), burden)
    mod <- lm("outcome ~ sex + burden", data=dat)
    chk <- coef(summary(mod))["burden", c("Estimate", "Std. Error")]
    checkEquals(chk, bt[c("Est", "SE")], checkNames=FALSE)

    seqClose(seqData)
}

test_burden_wald_bin <- function(){
    seqData <- .testObject()
    burden <- .testBurden(seqData)

    burden.test <- "Wald"
    nullmod <- .testNullModel(seqData, type="reg", binary=TRUE)
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test=burden.test)
    bt <- GENESIS:::.runBurdenTest(burden, proj, burden.test=burden.test)

    dat <- cbind(pData(sampleData(seqData)), burden)
    mod <- glm("status ~ sex + burden", data=dat, family=binomial)
    chk <- coef(summary(mod))["burden", c("Estimate", "Std. Error")]
    checkEquals(chk, bt[c("Est", "SE")], checkNames=FALSE)

    seqClose(seqData)
}

test_burden_firth <- function(){
    seqData <- .testObject()
    burden <- .testBurden(seqData)

    burden.test <- "Firth"
    nullmod <- .testNullModel(seqData, type="reg", binary=TRUE)
    proj <- GENESIS:::.calculateProjection(nullmod, test="Burden", burden.test=burden.test)
    bt <- GENESIS:::.runBurdenTest(burden, proj, burden.test=burden.test)

    dat <- cbind(pData(sampleData(seqData)), burden)
    mod <- logistf(as.formula("status ~ sex + burden"), data=dat)
    checkEquals(coef(mod)["burden"], bt["Est"], checkNames=FALSE)
    
    seqClose(seqData)
}


test_monomorphs <- function(){
    seqData <- .testObject()
    nullmod <- fitNullReg(sampleData(seqData), outcome="outcome")
    
    ## check that function works when list of supplied variants includes monomorphs
    af <- seqAlleleFreq(seqData)
    mono <- which(af == 1)
    agg <- list(data.frame(variant.id=(mono[1]-2):(mono[1]+2), allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg, verbose=FALSE)

    ## now check if monomorph is the only variant in a block
    agg <- list(data.frame(variant.id=1:10, allele.index=1),
                data.frame(variant.id=mono, allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg, verbose=FALSE)
    
    seqClose(seqData)
}


test_sexchrom <- function(){
    # make a test file with chr="X"
    gdsfile <- seqExampleFileName("gds")
    tmpfile <- tempfile()
    file.copy(gdsfile, tmpfile)
    gds <- openfn.gds(tmpfile, readonly=FALSE)
    node <- index.gdsn(gds, "chromosome")
    compression.gdsn(node, "")
    chr <- read.gdsn(node)
    write.gdsn(node, rep("X", length(chr)))
    closefn.gds(gds)
    gds.seq <- seqOpen(tmpfile)
    
    data(pedigree)
    pedigree <- pedigree[match(seqGetData(gds.seq, "sample.id"), pedigree$sample.id),]
    pedigree$outcome <- rnorm(nrow(pedigree))
    
    seqData <- SeqVarData(gds.seq, sampleData=AnnotatedDataFrame(pedigree))
    nullmod <- fitNullReg(sampleData(seqData), outcome="outcome", verbose=FALSE)
    
    agg <- list(data.frame(variant.id=1:10, allele.index=1),
                data.frame(variant.id=11:20, allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg, verbose=FALSE)
    checkTrue(all(do.call(rbind, assoc$variantInfo)[,"chr"] == "X"))
    
    assoc <- assocTestSeqWindow(seqData, nullmod, variant.include=1:50, verbose=FALSE)
    checkTrue(all(assoc$results[,"chr"] == "X"))
    checkTrue(all(assoc$variantInfo[,"chr"] == "X"))
    
    seqClose(seqData)
    unlink(tmpfile)
}


test_indexNotValue <- function(){
    seqData <- .testObject()
    nullmod <- .testNullModel(seqData, type="reg")
    assoc <- assocTestSeqWindow(seqData, nullmod, window.size=1000, window.shift=500, verbose=FALSE)

    # make new object with different variant.id
    annot <- sampleData(seqData)
    seqClose(seqData)
    gdsfile <- seqExampleFileName("gds")
    tmpfile <- tempfile()
    file.copy(gdsfile, tmpfile)
    gds <- openfn.gds(tmpfile, readonly=FALSE)
    node <- index.gdsn(gds, "variant.id")
    compression.gdsn(node, "")
    var <- read.gdsn(node)
    write.gdsn(node, var+1000)
    closefn.gds(gds)
    gds.seq <- seqOpen(tmpfile)
    seqData <- SeqVarData(gds.seq, sampleData=annot)
    assoc2 <- assocTestSeqWindow(seqData, nullmod, window.size=1000, window.shift=500, verbose=FALSE)
    checkEquals(assoc$results, assoc2$results)
    checkEquals(assoc$variantInfo[,-1], assoc2$variantInfo[,-1])
    
    seqClose(seqData)
    unlink(tmpfile)
}


test_variantInclude <- function() {
    seqData <- .testObject()
    var.id <- seqGetData(seqData, "variant.id")
    var.incl <- sample(var.id, 1000, replace=TRUE)
    res <- GENESIS:::getVariantInclude(seqData, var.incl, chromosome=NULL)
    freq <- SeqVarTools::alleleFrequency(seqData, use.names=TRUE)
    mono <- as.integer(names(freq[freq == 1]))
    checkEquals(res$value, sort(setdiff(var.incl, mono)))
    seqClose(seqData)
}


test_aggVarList <- function(){
    seqData <- .testObject()
    nullmod <- .testNullModel(seqData, type="reg")
    
    agg <- list(data.frame(variant.id=11:25, allele.index=1),
                data.frame(variant.id=1:15, allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg, verbose=FALSE)
    assoc1 <- assocTestSeq(seqData, nullmod, agg[1], verbose=FALSE)
    assoc2 <- assocTestSeq(seqData, nullmod, agg[2], verbose=FALSE)
    checkEquals(assoc$results[1,], assoc1$results)
    checkEquals(unlist(assoc$results[2,]), unlist(assoc2$results)) # avoid comparing row names
    
    seqClose(seqData)
}


test_skat <- function(){
    seqData <- .testObject()
    nullmod <- .testNullModel(seqData, type="reg")
    assoc <- assocTestSeqWindow(seqData, nullmod, window.size=1000, window.shift=500, chromosome=1, test="SKAT", rho=c(0,0.5,1), verbose=FALSE)
    seqClose(seqData)
}
