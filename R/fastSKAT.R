fitNullfamSKAT <- function(scanData,						
                           outcome, 
                           covars = NULL,
                           covMatList,
                           scan.include = NULL) {
    
    if(class(scanData) == "AnnotatedDataFrame"){
        scanData <- pData(scanData)
        names(scanData)[names(scanData) == "sample.id"] <- "scanID"
    }
    if(class(scanData) != "data.frame"){
        stop("scanData should either be a data.frame or an AnnotatedDataFrame")
    }

    # create design matrices and outcome vector
    dat <- createDesignMatrix(scanData = scanData, outcome = outcome, covars = covars, scan.include = scan.include)

    # if covMatList is a list, convert to a matrix
    if(is(covMatList, "list")){
        if (length(covMatList) > 1) warning("Using only the first element of covMatList")
        covMatList <- covMatList[[1]]
    }

    df <- data.frame(Y=dat$Y, dat$W[,-1,drop=FALSE], id=rownames(dat$W), stringsAsFactors=FALSE)
    model <- as.formula(paste("Y ~", paste(colnames(dat$W)[-1], collapse=" + "), "+ (1|id)"))
    mod <- coxme::lmekin(model, data=df, varlist=covMatList*2, x=TRUE, y=TRUE, method="REML")

    return(list(lmekin=mod, covMatList=covMatList, scanID=df$id))
}


.weightFunction <- function(weight.beta) {
    function(maf) dbeta(maf, weight.beta[1], weight.beta[2])
}


.runFastSKATTest <- function(model, geno.adj, weight.beta, pval.method, neig) {
    f <- bigQF::famSKAT(geno.adj, model$lmekin, model$covMatList, weights=.weightFunction(weight.beta))
    Q <- f$Q()
    meth <- if (pval.method == "kuonen") "saddlepoint" else "integration"
    neig <- min(c(neig, dim(geno.adj)))
    p <- bigQF::pQF(Q, f, neig=neig, convolution.method=meth)
    out <- c("Q_0"=Q, "pval_0"=p)
    return(out)
}
