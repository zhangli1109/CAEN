#' Estimate the probability that the read is 0 in a Zero-inflated
#' Poisson model.
#' @description Estimate the probability that the read is 0 in
#' a Zero-inflated Poisson model.
#' @usage estimatep(x,y,xte=NULL,beta=1,
#' type=c("mle","deseq","quantile"), prior=NULL)
#' @param x MUST be a n times p matrix - i.e. observations on the rows and
#' features on the columns
#' @param y A numeric vector of class labels of length n: 1, 2, ...., K
#' if there are K classes.Each element of y corresponds to a row of x;
#' i.e. these are the class labels for the observations in x.
#' @param xte A m-by-p data matrix: m test observations and p features.
#' The classifier fit on the training data set x will be tested on this
#' data set. If NULL, then testing will be performed on the training set.
#' @param beta A standardized parameter
#' @param type the method of normality
#' @param prior vector of length equal to the number of classes, representing
#' prior probabilities for each class.If NULL then uniform priors are used
#' (i.e. each class is equally likely)
#' @return p the probability that the read is 0 in a
#'  Zero-inflated Poisson model
#' @examples
#' library(SummarizedExperiment)
#' dat <- newCountDataSet(n=40,p=500, K=4, param=10, sdsignal=0.1,drate=0.4)
#' x <- t(assay(dat$sim_train_data))
#' y <- as.numeric(colnames(dat$sim_train_data))
#' xte <- t(assay(dat$sim_test_data))
#' prob <- estimatep(x=x, y=y, xte=x, beta=1, type="mle", prior=NULL)
#' @export

estimatep <- function(x, y, xte=NULL, beta=1,
                    type=c("mle","deseq","quantile"), prior=NULL){
    
    if (is.null(xte)) {
    xte <- x
    warning("Since no xte was provided, testing was performed
        on training data set.")
    }
    type <- match.arg(type)
    if (is.null(prior)) prior <- rep(1/length(unique(y)),
                                    length(unique(y)))
    null.out <- PoiClaClu::NullModel(x, type = type)
    ns <- null.out$n
    nste <- PoiClaClu::NullModelTest(null.out,x,xte,type = type)$nste
    ds <-  GetDn(ns,x,y,beta)
    mu <- matrix(NA, nrow = nrow(x), ncol = length(x[1,]))
    for (i in seq_len(nrow(x))) {
    dstar = ds[y[i],]
    mu[i,] <- nste[i,]*dstar
    }
    G <- length(x[1,])
    lib <- rowSums(x)
    x1 <- t(x)
    mu1 <- as.vector(t(mu))
    librep <- rep(lib,rep(G,length(lib)))/(lib[1])
    x2 <- as.vector(x1)
    y <- x2
    y[y != 0] <- 1
    xreg <- cbind(y,librep,mu1)
    glm.out <- glm(y ~ librep + mu1,family = binomial("logit"),
                    data = data.frame(xreg))
    coef <- as.matrix(glm.out$coefficients)
    inter <- rep(1,G)
    muu <- t(mu)
    xte1 <- t(xte)
    p <- xte1
    for (i in seq_len(length(xte1[1,]))) {
        libsize <- rep(sum(xte1[,i]), G)/lib[1]
        estx1 <- cbind(inter, libsize, muu[,1])
        dd <- estx1 %*% coef
        dd[dd > 50] <- 50
        dd[dd < (-50)] <- -50
        p1 <- exp(dd)/(1 + exp(dd))
        p[,i] <- 1 - p1
    }
    return(p)
}
