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
#'@examples
#'dat <- newCountDataSet(n=40,p=500,sdsignal=0.1,K=4,param=10,drate=0.4)
#'prob<-estimatep(dat$x,dat$y,dat$x,beta=1,
#'type=c("mle","deseq","quantile"),prior=NULL)
#'prob0<-estimatep(dat$x,dat$y,dat$xte,beta=1,
#'type=c("mle","deseq","quantile"),prior=NULL)
#'cv.out <- ZIPDA.cv(dat$x,dat$y,prob0=t(prob))
#'out <- ZIPLDA(dat$x,dat$y,dat$xte,rho=cv.out$bestrho,prob0=t(prob0))
#' @export

estimatep<-function(x,y,xte=NULL,beta=1,
                    type=c("mle","deseq","quantile"), prior=NULL){
    if(is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed
        on training data set.")
    }
    type <- match.arg(type)
    if(is.null(prior)) prior <- rep(1/length(unique(y)),
                                    length(unique(y)))

    null.out <- NullModel(x, type=type)
    ns <- null.out$n
    nste <- NullModelTest(null.out,x,xte,type=type)$nste
    uniq <- sort(unique(y))
    signx<-sign(xte==0)
    ds <-  GetDn(ns,x,y,beta)

    mu <- matrix(NA, nrow=nrow(x), ncol=length(x[1,]))
    # for(k in 1:length(y)){
    for(i in seq(from = 1, to = nrow(x))){
    dstar = ds[y[i],]
    mu[i,]<-nste[i,]*dstar
    }
    #}
    #mu1<-t(mu)


    G<-length(x[1,])
    lib<-rowSums(x)
    x1<-t(x)
    #stx1<-t(testx)

    mu1<-as.vector(t(mu))
    librep<-rep(lib,rep(G,length(lib)))/(lib[1])
    x2<-as.vector(x1)

    y<-x2
    y[y!=0]<-1
    xreg<-cbind(y,librep,mu1)
    glm.out<-glm(y~ librep+ mu1,family=binomial("logit"),
                    data=data.frame(xreg))
    summary(glm.out)

    coef<-as.matrix(glm.out$coefficients)
    inter<-rep(1,G)
    #libsize<-rep(sum(x1[,2]),G)/lib[1]
    muu<-t(mu)
    #estx1<-cbind(inter,libsize,muu[,1])
    #p<-exp(estx1%*% coef)/(1+exp(estx1%*% coef))


    xte1=t(xte)
    p<-xte1
    for(i in seq(from = 1, to = length(xte1[1,]))){
        libsize<-rep(sum(xte1[,i]),G)/lib[1]
        estx1<-cbind(inter,libsize,muu[,1])
        dd<-estx1%*% coef
        dd[dd>50]<-50
        dd[dd<(-50)]<--50
        p1<-exp(dd)/(1+exp(dd))
        p[,i]<-1-p1
    }
    return(p)
}
