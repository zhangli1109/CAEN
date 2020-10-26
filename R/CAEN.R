#' Compute the correlation coefficient of gene with category number to identify
#' differentially expressed genes
#' @description To Compute the correlation coefficient of gene with category
#' number to identify  differentially expressed genes.
#' @usage CAEN(dataTable, y, K, gene_no_list)
#' @param dataTable Matrix or data.frame containing read counts and the row 
#' denotes the gene
#' @param y the category for each sample
#' @param K the number of class
#' @param gene_no_list the number of differentially expressed genes you want
#'  to select
#' @return list(.) A list of computed correlation coefficient and the first some
#' differentially expressed genes , where "r" represents correlation coefficient
#'  between gene and category number, and "np" represents the top differential
#'  feature label.
#' @examples
#' library(SummarizedExperiment)
#' dat <- newCountDataSet(n=40,p=500,sdsignal=0.1,K=4,param=10,drate=0.4)
#' x <- t(assay(dat$sim_train_data))                  
#' y <- as.numeric(colnames(dat$sim_train_data))      
#' xte <- t(assay(dat$sim_test_data))                 
#' prob<-estimatep(x=x,y=y,xte=x,beta=1,
#' type=c("mle","deseq","quantile"),prior=NULL)      
#' prob0<-estimatep(x=x,y=y,xte=xte,beta=1,
#' type=c("mle","deseq","quantile"),prior=NULL)   
#' myscore<-CAEN(dataTable=assay(dat$sim_train_data),
#' y=as.numeric(colnames(dat$sim_train_data)),
#' K=4,gene_no_list=100)
#' @export
#' @import SummarizedExperiment

CAEN <- function(dataTable, y, K, gene_no_list)
{      
    if (K < 2)
    stop("We require K to be at least 2.")
    c <- rep(0, length = K)
    for (ii in seq_len(K)) {
        c[ii] <- sum(y == ii)
    }
    newdata <- dataTable
    p <- nrow(newdata)
    n <- ncol(newdata)
    newdata <- rbind(newdata, y)
    newdata <- newdata[, order(newdata[p + 1, ])]
    trainx <- newdata[-(p + 1), ]
    smyscore <- rep(0, length = nrow(trainx))
    ntrainx <- matrix(0, nrow = nrow(trainx), ncol = n)
    newylabel <- matrix(0, nrow = nrow(trainx), ncol = n)
    classmean <- seq_len(K)
    for (i in seq_len(nrow(trainx))) {
        for (iii in seq_len(K)) {
            classmean[iii] <- mean(dataTable[i, y == iii])
    }
    if (length(which(classmean == 0)) >= 2)
        next
    crank <-  rank(classmean)
    c1 <- rep(crank[1], c[1])
    for (ii in seq(from = 2, to = K)) {
        c1 <- c(c1,rep(crank[ii], c[ii]))
    }
    newylabel[i, ] <- c1
    b <- rbind(trainx[i, ], newylabel[i, ])
    ntrainx[i, ] <- b[1, order(b[2, ])]
    label <- c()
    nn <- c()
    nn[1] <- 0
    for (j in seq_len(K)) {
        bj <- which(newylabel[i, ] == j)
        nn[j + 1] <- nn[j] + length(bj)
        label <- c(label, rank(trainx[i, bj]) + nn[j])
    }
    rank <- rank(ntrainx[i, ])
    rankclass <- rank(label)
    sdIX <- sd(rank)
    sdIY <- sd(rankclass)
    sums <- 0
    for (jj in seq_len(n)) {
        sums <- sums + (rank[jj] - rankclass[jj])^2
    }
    fenzi <- (sdIX^2 + sdIY^2 - sums/n)/2
    smyscore[i] <- fenzi/(sdIX*sdIY)
    ssorttrainx <- sort.list(smyscore, decreasing = TRUE)
    ssorttrainx <- ssorttrainx[seq_len(gene_no_list)]
    }
    return(list(r = smyscore, np = ssorttrainx))

}

