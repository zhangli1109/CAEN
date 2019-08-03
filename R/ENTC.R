#' Compute the correlation coefficient of gene with category number to identify
#' feature genes
#' @description To Compute the correlation coefficient of gene with category
#' number to identify  differentially expressed genes.
#' @usage ENTC(data,y,K,gene_no_list)
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
#' dat <- newCountDataSet(n=40,p=500,sdsignal=0.1,K=4,param=10,drate=0.4)
#' ENTC(t(dat$x),dat$y,K=4,gene_no_list=100)
#' @export

ENTC <- function (data,y,K,gene_no_list)
{


    if (K < 2)
    s
