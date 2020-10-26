#' Function to do cross-validation for zero-inflated Poisson classification.
#' @description Perform cross-validation for the function that
#' implements the "sparse zero-inflated Poisson linear discriminant
#' analysis classifier",which is similar to linear discriminant analysis
#' but assumes a zero-inflated Poisson model rather than a Gaussian model
#' for the data. The classifies soft-thresholds the estimated effect of
#' each feature in order to achieve sparsity. This cross-validation
#' function selects the proper value of the tuning parameter that controls
#' the level of soft-thresholding.
#' @usage ZIPDA.cv(x, y, rhos = NULL, beta = 1, nfolds = 5, prob0=NULL,
#' type=c("mle","deseq","quantile"),folds = NULL, transform=TRUE, alpha=NULL,
#' prior=NULL)
#' @param x A n-by-p training data matrix; n observations and p features.
#' @param y A numeric vector of class labels of length n: 1, 2, ...., K if
#' there are K classes.Each element of y corresponds to a row of x;
#' i.e. these are the class labels for the observationsin x.
#' @param rhos A vector of tuning parameters to try out in cross-validation.
#' Rho controls the level of shrinkage performed, i.e. the number of features
#' that are not involved in the classifier. When rho=0 then all features
#' are involved in the classifier, and when rho is very large no features
#' are involved. If rhos=NULL then a vector of rho values will be chosen
#' automatically.
#' @param beta A smoothing term. A Gamma(beta,beta) prior is used to fit the
#' zero-inflated Poisson model.Recommendation is to leave it at 1, the default
#' value.
#' @param nfolds The number of folds in the cross-validation; default is
#' 5-fold cross-validation.
#' @param prob0 The probability that the read is 0
#' @param type How should the observations be normalized within the
#' zero-inflated
#' Poisson model, i.e. how should the size factors be estimated?
#' Options are "quantile" or "deseq" (more robust) or "mle" (less robust).
#' In greater detail: "quantile" is quantile normalization approach
#' of Bullard et al 2010 BMC Bioinformatics, "deseq"  is median of the
#' ratio of an observation to a pseudoreference obtained by taking the
#' geometric mean, described in Anders and Huber 2010 Genome
#' Biology and implemented in Bioconductor package "DESeq", and "mle" is
#' the sum of counts for each sample; this is the maximum likelihood
#' estimate under a simple Poisson model.
#' @param prior Vector of length equal to the number of classes,
#' representing prior probabilities for each class. If NULL then
#' uniform priors are used (i.e.each class is equally likely).
#' @param transform Should data matrices x and xte first be power transformed 
#' so that it more  closely fits the zero-inflated Poisson model? TRUE or FALSE.
#' Power transformation is especially  useful if the data are overdispersed
#' relative to the zero-inflated Poisson model.
#' @param alpha If transform=TRUE, this determines the power to which the data
#' matrices x and xte are transformed. If alpha=NULL then the
#' transformation that makes the zero-inflated Poisson model best fit the data
#' matrix x is computed. (Note that alpha is computed based on x, not based
#' on xte). Or a value of alpha, 0<alpha<=1, can be  entered by the user.
#' @param folds Instead of specifying the number of folds in cross-validation,
#' one can explicitly specify the folds. To do this, input a list of length
#' r(to perform r-fold cross-validation). The rth element of the list
#' should be vector containing the indices of the test observations in the
#' rth fold.
#' @return list(.) A list of output, "errs" represents A matrix of dimension
#' (number of folds)-by-(length of rhos)."bestrho" represents The tuning
#' parameter value resulting in the lowest overall cross-validation error
#' rate for. "rhos" represent the vector of rho values used in
#' cross-validation."nnonzero" represents A matrix of dimension (number
#' of folds)-by-(length of rhos)."folds" represents Cross-validation folds
#' used. "alpha" represents Power transformation used (if transform=TRUE).
#' @examples
#' library(SummarizedExperiment)
#' dat <- newCountDataSet(n=40,p=500, K=4, param=10, sdsignal=0.1,drate=0.4)
#' x <- t(assay(dat$sim_train_data))
#' y <- as.numeric(colnames(dat$sim_train_data))
#' xte <- t(assay(dat$sim_test_data))
#' prob<-estimatep(x=x, y=y, xte=x, beta=1, type="mle", prior=NULL)
#' prob0<-estimatep(x=x, y=y, xte=xte, beta=1,type="mle", prior=NULL)
#' cv.out <- ZIPDA.cv(x=x, y=y, prob0=t(prob))
#' out <- ZIPLDA(x=x, y=y, xte=xte, rho=cv.out$bestrho, prob0=t(prob0))
#' @export


ZIPDA.cv <-
    function(x,y,rhos = NULL, beta = 1, nfolds = 5, prob0 = NULL,
             type = c("mle","deseq","quantile"), folds = NULL, transform = TRUE,
             alpha = NULL, prior = NULL){
        type <- match.arg(type)
        if (!transform && !is.null(alpha)) stop("You have asked for NO
                            transformation but have entered alpha.")
        if (transform && is.null(alpha)) 
            alpha <- PoiClaClu::FindBestTransform(x)
        if (transform) {
            if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
            x <- x^alpha
        }
        if (is.null(rhos)) {
            ns <- PoiClaClu::NullModel(x,type = type)$n
            uniq <- sort(unique(y))
            maxrho <- rep(NA, length(uniq))
            for (k in seq_len(length(uniq))) {
                a <- colSums(x[y == uniq[k],]) + beta
                b <- colSums(ns[y == uniq[k],]) + beta
                maxrho[k] <- max(abs(a/b - 1)*sqrt(b),na.rm = TRUE)
            }
            rhos <- seq(0, max(maxrho,na.rm = TRUE)*(2/3), len = 30)
        }
        if (is.null(folds)) folds <- balanced.folds(y,nfolds = nfolds)
        nfolds <- length(folds)
        errs <- nnonzero <- matrix(NA, nrow = nfolds, ncol = length(rhos))
        for (i in seq_len(nfolds)) {
            cat(i, fill = FALSE)
            tr <- -folds[[i]]
            te <- folds[[i]]
            out <- ZIPLDA(x[tr,], y[tr], x[te,], rhos = rhos, beta = beta,
                          prob0 = prob0, type = "mle", prior = prior, 
                          transform = FALSE)
            for (j in seq_len(length(rhos))) {
                errs[i,j] <- sum(out[[j]]$ytehat != y[te])
                nnonzero[i,j] <- sum(colSums(out[[j]]$ds != 1) != 0)
            }
        }
        cat(fill = TRUE)
        save <- list(errs = errs, bestrho = rhos
                     [max(which(colMeans(errs) == min(colMeans(errs))))],
                     rhos = rhos, nnonzero = nnonzero, folds = folds,
                     alpha = alpha, type = type)
        return(save)
    }
