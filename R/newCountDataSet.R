#' Generate a simulated sequencing data set using a negative binomial model
#' @description Generate two nxp data sets: a training set and a test set,
#' as well as outcome vectors y and yte of length n indicating the class
#' labels of the training and test observations.
#' @usage newCountDataSet(n, p, K, param, sdsignal,drate)
#' @param n Number of observations desired.
#' @param p Number of features desired. Note that drate of the features
#'  will differ between classes, though some of those differences may be small.
#' @param K Number of classes desired. Note that the function requires that n 
#' be at least equal to 4K.i.e. there must be at least 4 observations per class 
#' on average.
#' @param param The dispersion parameter for the negative binomial distribution.
#' The negative binomial distribution is parameterized using "mu" and "size" in
#' the R function "rnbinom". That is, Y ~ NB(mu, param) means that E(Y)=mu and
#' Var(Y) = mu+mu^2/param.So when param is very large this is essentially a
#' Poisson distribution, and when param is smaller then there is a lot of
#' overdispersion relative to the Poisson distribution.
#' @param sdsignal The extent to which the classes are different. If this 
#' equals zero then there are no class differences and if this is large then the
#' classes are very different.
#' @param drate The proportion of differentially expressed genes
#' @return list(.) A list of output, 
#' "sim_train_data" represents training data of q*n
#' data matrix. 
#' "sim_test_data" represents test data of q*n data matrix.
#' The colnames of this two matrix are class labels for the n observations
#' May have q<p because features with 0 total counts are removed.
#' The q features are those with >0 total counts in dataset. So q <= p.
#' "truesf" denotes size factors for training observations."isDE" represnts
#' the differential gene label.
#' @examples
#' dat <- newCountDataSet(n=40,p=500, K=4, param=10, sdsignal=0.1,drate=0.4)
#' @export


newCountDataSet <- function(n, p, K, param, sdsignal, drate)
{
    if (n < 4 * K)
    stop("We require n to be at least 4*K.")
    q0 <- rexp(p, rate = 1/20)
    isDE <- runif(p) < drate    
    classk <- matrix(NA, nrow = K, ncol = p)
        for (k in seq_len(K)) {
            lfc <- rnorm(p, sd = sdsignal)
            classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)
        }
    truesf <- runif(n) * 5 + 50
    truesfte <- runif(n) * 5 + 50
    conds <- sample(c(rep(seq_len(K), 4),
            sample(seq_len(K),n - 4 * K, replace = TRUE)))
    condste <- sample(c(rep(seq_len(K), 4),
            sample(seq_len(K), n - 4 * K, replace = TRUE)))
    x <- xte <- matrix(NA, nrow = n, ncol = p)
    prob <- runif(n,0.05,0.5)
    for (i in seq_len(n)) {
        for (k in seq_len(K)) {
            if (conds[i] == k)
            x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, ],size = param)
            ind0 <- which(rank(x[i, ]) < round(p*prob[i]))
            x[i,ind0] <- 0
            if (condste[i] == k)
                xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, ],
                        size = param)
            ind0 <- which(rank(xte[i, ]) < round(p*prob[i]))
            xte[i,ind0] <- 0
    }
}
    rm <- apply(x, 2, sum) == 0
    sim_train_data <- 
      SummarizedExperiment::SummarizedExperiment(assays 
                                                 = list(as.matrix(t(x[,!rm]))))
    rownames(sim_train_data) <- seq_len(nrow(t(x[,!rm])))
    colnames(sim_train_data) <- conds
    sim_test_data <- 
      SummarizedExperiment::SummarizedExperiment(assays 
                                                 = list(as.matrix(t(xte[,!rm]))))
    rownames(sim_test_data) <- seq_len(nrow(t(xte[,!rm])))
    colnames(sim_test_data) <- condste
    return(list(sim_train_data = sim_train_data, sim_test_data = sim_test_data,
        truesf = truesf, truesfte = truesfte, isDE = isDE[!rm]))
}




