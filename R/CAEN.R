#' Compute the correlation coefficient of gene with category number to identify
#' differentially expressed genes
#' @description To Compute the correlation coefficient of gene with category
#' number to identify  differentially expressed genes.
#' @usage CAEN(data,y,K,gene_no_list)
#' @param data Matrix or data.frame containing read counts and the row denotes
#' the gene
#' @param y the category for each sample
#' @param K the number of class
#' @param gene_no_list the number of differentially expressed genes you want
#'  to select
#' @return list(.) A list of computed correlation coefficient and the first some
#' differentially expressed genes , "r"represents correlation coefficient
#'  between gene and category number, "np" represents the top differential
#'  feature label.
#' @examples
#' library(SummarizedExperiment)
#' data(sim_data)
#' CAEN(data=assay(sim_data),y=y,K=4,gene_no_list=100)
#' @export
#' @import SummarizedExperiment

CAEN <- function (data,y,K,gene_no_list)
{


    if (K < 2)
    stop("We require K to be at least 2.")
    c<-rep(0,length=K)
    for(ii in seq(from = 1, to = K)){
        c[ii]=sum(y==ii)
    }

    #dat$x
    ##ZIPLDA
    #prob11<-estp(dat$x)
    #prob01<-estp(dat$xte)
    newdata<-data
    p=nrow(newdata)
    n=ncol(newdata)
    newdata <- rbind(newdata,y)
    #dim(newdata)
    newdata=newdata[,order(newdata[p+1,])]
    trainx=newdata[-(p+1),]
    #[,trlabel]
    smyscore <- rep(0,length=nrow(trainx))
    ntrainx<-matrix(0,nrow=nrow(trainx),ncol = n)
    newylabel<-matrix(0,nrow=nrow(trainx),ncol = n)
    classmean<-seq(from = 1, to = K)
    for (i in seq(from = 1, to = nrow(trainx))){
        for(iii in seq(from = 1, to = K)){
            classmean[iii]=mean(data[i,y==iii])
    }
    if(length(which(classmean==0))>=2)
        next
    crank= rank(classmean)
    c1=rep(crank[1],c[1])
    for(ii in seq(from = 2, to = K)){
        c1<- c(c1,rep(crank[ii],c[ii]))
    }
    newylabel[i,]=c1
    b = rbind(trainx[i,],newylabel[i,])
    ntrainx[i,]=b[1,order(b[2,])]#new x
    label<-c()
    nn<-c()
    nn[1]<-0
    for(j in seq(from = 1, to = K)){
        bj=which(newylabel[i,]==j)
        nn[j+1]=nn[j]+length(bj)
        label<-c(label,rank(trainx[i,bj])+nn[j])
    }
    rank<-rank(ntrainx[i,])
    rankclass<-rank(label)
    sdIX=sd(rank)
    sdIY=sd(rankclass)
    sum<-0
    for(jj in seq(from = 1, to = n)){
        sum=sum+(rank[jj]-rankclass[jj])^2
    }
    fenzi<-(sdIX^2+sdIY^2-sum/n)/2
    smyscore[i]=fenzi/(sdIX*sdIY)
    ssorttrainx <- sort.list(smyscore, decreasing=TRUE)
    ssorttrainx <- ssorttrainx[seq(from = 1, to = gene_no_list)]

    }
    return(list(r=smyscore,np=ssorttrainx))

}

