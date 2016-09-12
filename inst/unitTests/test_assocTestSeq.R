library(SeqVarTools)
library(Biobase)
library(SNPRelate)
library(logistf)
library(gdsfmt)
library(coxme)
library(bigQF)

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
    pcrelateMakeGRM(mypcrel, scaleKin=1)
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
    nullmod <- .testNullModel(seqData, type="reg")
    
    ## check that function works when list of supplied variants includes monomorphs
    af <- seqAlleleFreq(seqData)
    mono <- which(af == 1)
    agg <- list(data.frame(variant.id=(mono[1]-2):(mono[1]+2), allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg)

    ## now check if monomorph is the only variant in a block
    agg <- list(data.frame(variant.id=1:10, allele.index=1),
                data.frame(variant.id=mono, allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg)
    
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
    nullmod <- .testNullModel(seqData, type="reg")
    
    agg <- list(data.frame(variant.id=1:10, allele.index=1),
                data.frame(variant.id=11:20, allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg)
    checkTrue(all(do.call(rbind, assoc$variantInfo)[,"chr"] == "X"))
    
    assoc <- assocTestSeqWindow(seqData, nullmod, variant.include=1:50)
    checkTrue(all(assoc$results[,"chr"] == "X"))
    checkTrue(all(assoc$variantInfo[,"chr"] == "X"))
    
    seqClose(seqData)
    unlink(tmpfile)
}


test_indexNotValue <- function(){
    seqData <- .testObject()
    nullmod <- .testNullModel(seqData, type="reg")
    assoc <- assocTestSeqWindow(seqData, nullmod, window.size=1000, window.shift=500)

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
    assoc2 <- assocTestSeqWindow(seqData, nullmod, window.size=1000, window.shift=500)
    checkEquals(assoc$results, assoc2$results)
    checkEquals(assoc$variantInfo[,-1], assoc2$variantInfo[,-1])
    
    seqClose(seqData)
    unlink(tmpfile)
}


test_famSKAT_null <- function(){
    seqData <- .testObject()
    grm <- .testGRM(seqData)
    dat <- pData(sampleData(seqData))
    dat$sex <- unname(c(M=1, F=2)[dat$sex])
    names(dat)[names(dat) == "sample.id"] <- "id" # famSKAT only works if random effect is called 'id'
    mod <- coxme::lmekin(outcome ~ sex + (1|id), data=dat, varlist=grm*2,
                         x=TRUE, y=TRUE, method="REML")

    nullmod <- fitNullfamSKAT(sampleData(seqData), outcome="outcome", covars="sex", covMatList=grm)

    checkEquals(abs(mod$coefficients$fixed["sex"]), abs(nullmod$lmekin$coefficients$fixed["sexM"]), checkNames=FALSE)
    checkEquals(mod$coefficients$random, nullmod$lmekin$coefficients$random, checkNames=FALSE)
    
    seqClose(seqData)
}

test_famSKAT <- function(){
    seqData <- .testObject()
    grm <- .testGRM(seqData)
    nullmod <- fitNullfamSKAT(sampleData(seqData), outcome="outcome", covars="sex", covMatList=grm)
    
    g <- altDosage(seqData)
    g <- g[,colSums(is.na(g)) == 0]

    ## without using update.famSKAT
    f <- bigQF::famSKAT(g, nullmod$lmekin, nullmod$covMatList, weights= GENESIS:::.weightFunction(c(1,25)))
    Q <- f$Q()

    f <-  GENESIS:::.makeFastSKATTest(nullmod, weight.beta=c(1,25))
    res <-  GENESIS:::.runFastSKATTest(f, g, pval.method="kuonen", neig=10)
    checkEquals(Q, res["Q_0"], checkNames=FALSE)
    
    ## compare to regular SKAT
    nullmod <- .testNullModel(seqData)
    proj <- GENESIS:::.calculateProjection(nullmod, "SKAT", NULL)
    scores <- as.vector(crossprod(g, proj$resid))
    freq <- colMeans(g/2)
    weights <- GENESIS:::.weightFunction(c(1,25))(freq)
    res0 <- GENESIS:::.runSKATTest(scores, g, weights, rho=0, pval.method="kuonen", optimal=FALSE)
    
    seqClose(seqData)
}
    
test_fastSKAT <- function(){
    seqData <- .testObject()
    grm <- .testGRM(seqData)
    nullmod <- fitNullfamSKAT(sampleData(seqData), outcome="outcome", covars="sex", covMatList=grm)
    agg <- list(data.frame(variant.id=seqGetData(seqData, "variant.id"), allele.index=1))
    assoc <- assocTestSeq(seqData, nullmod, agg, test="fastSKAT", fastSKAT.neig=10)
    
    ## compare to regular SKAT
    nullmod0 <- .testNullModel(seqData)
    assoc0 <- assocTestSeq(seqData, nullmod0, agg, test="SKAT")

    ## window
    assoc <- assocTestSeqWindow(seqData, nullmod, test="fastSKAT", window.size=100000, window.shift=50000, fastSKAT.neig=3)
    assoc0 <- assocTestSeqWindow(seqData, nullmod0, test="SKAT", window.size=100000, window.shift=50000)
       
    seqClose(seqData)
}
