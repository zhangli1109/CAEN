#' Generate a simulated sequencing data set using a negative binomial model
#' @description Generate two nxp data sets: a training set and a test set,
#' as well as outcome vectors y and yte of length n indicating the class
#' labels of the training and test observations.
#' @usage newCountDataSet(n, p, K, param, sdsignal,drate)
#' @param n Number of observations desired.
#' @param p Number of features desired. Note that drate of the features
#'  will differ between classes, though some of those differences may be small.
#' @param K Number of classes desired. Note that the function requires that n be
#'  at least equal to 4K.i.e. there must be at least 4 observations per class on
#'  average.
#' @param param The dispersion parameter for the negative binomial distribution.
#' The negative binomial distribution is parameterized using "mu" and "size" in
#' the R function "rnbinom". That is, Y ~ NB(mu, param) means that E(Y)=mu and
#' Var(Y) = mu+mu^2/param.So when param is very large this is essentially a
#' Poisson distribution, and when param is smaller then there is a lot of
#' overdispersion relative to the Poisson distribution.
#' @param sdsignal The extent to which the classes are different. If this equals
#' zero then there are no class differences and if this is large then the
#' classes are very different.
#' @param drate The proportion of differentially expressed genes
#' @return list(.) A list of output, "x" represents training data of nxq
#' data matrix.
#' May have q<p because features with 0 total counts are removed."y" represents
#' class labels for the n observations in x."xte" represents nxq data
#' matrix of test observations; the q features are those with >0 total counts
#' in x. So q<=p."yte"represnets class labels for the n observation in xte.
#'"truesf" denotes size factors for training observations."isDE" represnts
#' the differential gene label.
#' @examples
#' library(SummarizedExperiment)
#' sim_data<-newCountDataSet(n=20,p=100,K=4,param=10,sdsignal=2,drate=0.4)
#' @export
#' @import SummarizedExperiment
newCountDataSet<- function (n, p, K, param, sdsignal, drate)
{
    if (n < 4 * K)
    stop("We require n to be at least 4*K.")
    q0 <- rexp(p, rate = 1/20)#lamuda
    isDE <-  runif(p)< drate    #runif(p)
    classk <- matrix(NA, nrow = K, ncol = p)
        for (k in seq(from = 1, to = K)) {
            lfc <- rnorm(p, sd = sdsignal)#ln d
            classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)
            }
    #truesf <- runif(n) * 5 + 50#s
    #truesfte <- runif(n) * 5 + 50
    truesf <- runif(n) * 5 + 50#s
    truesfte <- runif(n) * 5+ 50
    conds <- sample(c(rep(seq(from = 1, to = K), 4),
            sample(seq(from = 1, to = K),n - 4 * K, replace = TRUE)))
    condste <- sample(c(rep(seq(from = 1, to = K), 4),
            sample(seq(from = 1, to = K), n - 4 * K, replace = TRUE)))
    x <- xte <- matrix(NA, nrow = n, ncol = p)
    #prob<-rep(0,length=n)
    prob<- runif(n,0.05,0.5)
    for (i in seq(from = 1, to = n)) {
        for (k in seq(from = 1, to = K)) {
            if (conds[i] == k)
            x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, ],size = param)
            ind0<-which(rank(x[i, ])<round(p*prob[i]))
            #ind0<-sample(1:p, round(p*prob[i]), replace = FALSE)
            x[i,ind0] <-0
            if (condste[i] == k)
                xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, ],
                        size = param)
            ind0<-which(rank(xte[i, ])<round(p*prob[i]))
            #ind0<-sample(1:p, round(p*prob[i]), replace = FALSE)
            xte[i,ind0] <-0
    }
}
    rm <- apply(x, 2, sum) == 0
    return(list(x = x[,!rm], xte = xte[,!rm], y = conds, yte = condste,
        truesf = truesf, truesfte = truesfte, isDE=isDE[!rm]))
}




