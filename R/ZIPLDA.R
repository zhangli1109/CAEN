#' Classify observations using a Zero-inflated Poisson model.
#' @description Classify observations using a Zero-inflated Poisson model.
#' @usage ZIPLDA(x, y, xte=NULL, rho = 0, beta = 1, rhos = NULL, prob0=NULL,
#' type=c("mle","deseq","quantile"),prior = NULL, transform=TRUE, alpha=NULL)
#' @param x A n-by-p training data matrix; n observations and p features.
#' Used to train the classifier.
#' @param y A numeric vector of class labels of length n: 1, 2, ...., K
#' if there are K classes.Each element of y corresponds to a row of x;
#' i.e. these are the class labels for the observations in x.
#' @param xte A m-by-p data matrix: m test observations and p features.
#' The classifier fit on the training data set x will be tested on this
#' data set. If NULL, then testing will be performed on the training set.
#' @param rho Tuning parameter controlling the amount of soft thresholding
#'  performed, i.e. the level of sparsity, i.e. number of nonzero features
#'  in classifier. Rho=0 means that there is no soft-thresolding,
#'  i.e. all  features used in classifier. Larger rho means that fewer
#'  features will be used.
#' @param beta A smoothing term. A Gamma(beta,beta) prior is used to fit
#'  the Zero-inflated Poisson model.Recommendation is to just leave it
#'  at 1, the default value.
#' @param rhos A vector of tuning parameters that control the amount of
#' soft thresholding performed. If "rhos" is provided then a number of models
#'  will be fit (one for each element of "rhos"), and a number of predicted
#'  class labels will be output (one for each element of "rhos").
#' @param prob0  The probability that the read is 0
#' @param type How should the observations be normalized within the
#' Zero-inflated Poisson model, i.e. how should the size factors be estimated?
#'  Options are "quantile" or "deseq" (more robust) or "mle" (less robust).
#' In greater detail: "quantile" is quantile normalization approach of
#' Bullard et al 2010 BMC Bioinformatics,"deseq" is median of the ratio of
#' an observation to a pseudoreference obtained by taking the geometric mean,
#' described in Anders and Huber 2010 Genome Biology and implemented in
#' Bioconductor package "DESeq", and "mle" is the sum of counts for each
#' sample; this is the maximum likelihood estimate under a simple Zero-inflated
#' Poisson model.
#' @param prior vector of length equal to the number of classes, representing
#'  prior probabilities for each class.If NULL then uniform priors are used
#' (i.e. each class is equally likely)
#' @param transform should data matrices x and xte first be power transformed
#' so that it more closely fits the Zero-inflated Poisson model? TRUE or FALSE.
#' Power transformation is especially useful if the data are overdispersed
#' relative to the Zero-inflated Poisson model.
#' @param alpha if transform=TRUE, this determines the power to which the
#' data matrices x and xte are transformed.If alpha=NULL then the transformation
#' that makes the Zero-inflated Poisson model best fit the data matrix x is
#' computed.(Note that alpha is computed based on x, not based on xte).
#' Or a value of alpha, 0<alpha<=1, can be entered by the user.
#' @return list(.) A list of output, "ytehat" represents The predicted class
#' labels for each of the test observations(rows of xte)."discriminant"
#' represents A m-by-K matrix, where K is the number of classes. The (i,k)
#' element is large if the ith element of xte belongs to class k."ds" A
#' K-by-p matrix indicating the extent to which each feature is under-or
#' over-expressed in each class. The (k,j) element is >1 if feature j is
#' over-expressed in class k, and is <1 if feature j is under-expressed in
#' class k. When rho is large then many of the elemtns of this matrix
#' are shrunken towards 1(no over- or under-expression)."alpha" represents
#' Power transformation used (if transform=TRUE).
#'@examples
#'dat <- newCountDataSet(n=40,p=500,sdsignal=0.1,K=4,param=10,drate=0.4)
#'cv.out <- ZIPDA.cv(dat$x,dat$y)
#'out <- ZIPLDA(dat$x,dat$y,dat$xte,rho=cv.out$bestrho)
#' @export

ZIPLDA<-
    function(x,y,xte=NULL,rho=0,beta=1,rhos=NULL,prob0=NULL,
        type=c("mle","deseq","quantile"),
        prior=NULL, transform=TRUE, alpha=NULL){
    if(is.null(xte)){
        xte <- x
        warning("Since no xte was provided, testing was
            performed on training data set.")
    }
    if(!is.null(rho) && length(rho)>1)
        stop("Can only enter 1 value of rho. If you would like
            to enter multiple values, use rhos argument.")
    type <- match.arg(type)
    if(!transform && !is.null(alpha)) stop("You have asked for
                NO transformation but have entered alpha.")
    if(transform && is.null(alpha)) alpha <- FindBestTransform(x)
    if(transform){
        if(alpha<=0 || alpha>1) stop("alpha must be between 0 and 1")
        x <- x^alpha
        xte <- xte^alpha
    }
    if(is.null(prior)) prior <- rep(1/length(unique(y)), length(unique(y)))
    if(is.null(rho)&&is.null(rhos)) stop("Must enter rho or rhos.")
    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    signx<-sign(xte==0)
    if(is.null(rhos)){
        ds <- GetD(ns,x,y,rho,beta)
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in seq(from = 1, to = length(uniq))){
        for(i in seq(from = 1, to = nrow(xte))){
            dstar = ds[k,]
            part2=nste[i,]*dstar
            part1=prob0[i,]+(1-prob0[i,])*exp(-part2)
            part1[part1==0]=1
            discriminant[i,k] <-sum(signx[i,]*log(part1))+
                sum(xte[i,]*(1-signx[i,])*log(dstar))-sum((1-signx[i,])*part2)+
                log(prior[k])

        }
    }
    save <- list(ns=ns,nste=nste,ds=ds,discriminant=discriminant,
                ytehat=uniq[apply(discriminant,1,which.max)],
                alpha=alpha,rho=rho,x=x,y=y,xte=xte,type=type)
    return(save)
    } else {
        save <- list()
        ds.list <- GetD(ns,x,y,rho=NULL, rhos=rhos,beta)
        for(rho in rhos){
        ds <- ds.list[[which(rhos==rho)]]
        discriminant <- matrix(NA, nrow=nrow(xte), ncol=length(uniq))
        for(k in seq(from = 1, to = length(uniq))){
            for(i in seq(from = 1, to = nrow(xte))){
            dstar = ds[k,]
            part2=nste[i,]*dstar
            part1=prob0[i]+(1-prob0[i])*exp(-part2)
            part1[part1==0]=1
            discriminant[i,k] <-sum(signx[i,]*log(part1))+sum(xte[i,]*
                    (1-signx[i,])*log(dstar))-sum((1-signx[i,])*part2)
                    +log(prior[k])
        }
        }
        save[[which(rhos==rho)]] <- list(ns=ns,nste=nste,ds=ds,
                discriminant=discriminant,
                ytehat=uniq[apply(discriminant,1,which.max)],
                alpha=alpha, rho=rho,x=x,y=y,xte=xte,type=type)
        }
        return(save)
        }
}


