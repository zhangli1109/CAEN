# 02/07/11 -- Helper functions.

balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
    totals <- table(y)
    fmax <- max(totals)
    nfolds <- min(nfolds, fmax)
    yids <- split(seq(y), y)
    bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
    for (i in seq(totals)) {
        bigmat[seq(totals[i]), i] <- sample(yids[[i]])
    }
    smallmat <- matrix(bigmat, nrow = nfolds) 
    smallmat <- permute.rows(t(smallmat)) 
    x <- apply(smallmat, 2, function(x) x[!is.na(x)])
    if (is.matrix(x)) {
        xlist <- list()
        for (i in seq(from = 1, to = ncol(x))) {
            xlist[[i]] <- x[,i]
        }
        return(xlist)
    }
    return(x)
}

permute.rows <- function(x){
    dd <- dim(x)
    n <- dd[1]
    p <- dd[2]
    mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
    matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}


Soft <- function(x,a){
    return(sign(x)*pmax(abs(x) - a,0))
}

GetD <- function(ns, x, y, rho,beta,rhos=NULL){
    if (!is.null(rho) && !is.null(rhos))
        stop("do you want to use rho or rhos in GetD function???")
    if (is.null(rhos)) {
        uniq <- sort(unique(y))
        ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
        for (k in seq(from = 1, to = length(uniq))) {
            a <- colSums(x[y == uniq[k],]) + beta
            b <- colSums(ns[y == uniq[k],]) + beta
            ds[k,] <- 1 + Soft(a/b - 1,rho/sqrt(b))
        }
        return(ds)
    } else {
        uniq <- sort(unique(y))
        ds.list <- list()
        for (rho in rhos) {
            ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
            for (k in seq(from = 1, to = length(uniq))) {
                a <- colSums(x[y == uniq[k],]) + beta
                b <- colSums(ns[y == uniq[k],]) + beta
                ds[k,] <- 1 + Soft(a/b - 1,rho/sqrt(b))
            }
        ds.list[[which(rhos == rho)]] <- ds
        }
    return(ds.list)
    }
}


GetDn <- function(ns, x, y, beta){
    uniq <- sort(unique(y))
    ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
    for (k in seq(from = 1, to = length(uniq))) {
        a <- colSums(x[y == uniq[k],]) + beta
        b <- colSums(ns[y == uniq[k],]) + beta
        ds[k,] <- a/b
    }
    return(ds)
}


