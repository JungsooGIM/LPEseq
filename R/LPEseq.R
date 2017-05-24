#' A function testing differential expression between two conditions
#' @param expr_x    a numeric vector or matrix in first condition
#' @param expr_y    a numeric vector or matrix in second condition
#' @param n.bin     the number of quantile bins (default 100)
#' @param df        degrees of freedom in smooth.spline() function (default 10)
#' @param d         expression difference in LPEseq.outlier() function (defalut 0.5; use 1.2 for biological replicates)
#' @details See the reference paper
#' @return \item{dataframe}{dataframe containing mean of each condition, local pooled variance, statistics, p-value, and BH adjusted p-value}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
LPEseq.test <- function(expr_x, expr_y, n.bin=100, df=10, d=1.2) {
  n.bin = n.bin
  df = df
  n.x = ncol(as.data.frame(expr_x))
  n.y = ncol(as.data.frame(expr_y))
  n.gene = nrow(as.data.frame(expr_x))
  if(n.x == 1 & n.y == 1){
    mu.x <- expr_x
    mu.y <- expr_y
    tmp_dat <- cbind(expr_x, expr_y)
    tmp.var <- LPEseq.var(tmp_dat, n.bin=n.bin, d=d, df=df)
    var.x <- LPEseq.predict.var(mu.x, tmp.var)
    var.y <- LPEseq.predict.var(mu.y, tmp.var)
    std.dev <- sqrt(var.x + var.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x == 1 & n.y == 2){
    mu.x <- expr_x
    mu.y <- apply(expr_y, 1, mean)
    tmp.dat <- cbind(mu.x, mu.y)
    tmp.var <- LPEseq.var(tmp.dat, n.bin=n.bin, df=df, d=d)
    var.x <- LPEseq.predict.var(mu.x, tmp.var)
    basevar.y <- lpe.var(expr_y, n.bin=n.bin, df=df)
    var.y <- basevar.y[,2]
    std.dev <- sqrt(var.x + var.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
      data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x == 1 & n.y > 2){
    mu.x <- expr_x
    mu.y <- apply(expr_y, 1, mean)
    tmp.dat <- cbind(mu.x, mu.y)
    tmp.var <- LPEseq.var(tmp.dat, n.bin=n.bin, d=d, df=df)
    var.x <- LPEseq.predict.var( mu.x, tmp.var)
    basevar.y <- lpe.var(expr_y, n.bin=n.bin, df=df)
    var.y <- basevar.y[,2]
    std.dev <- sqrt(var.x + (pi/2)*var.y/n.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x == 2 & n.y == 1){
    basevar.x <- lpe.var(expr_x, n.bin=n.bin, df=df)
    var.x <- basevar.x[,2]
    mu.x <- basevar.x[,1]
    mu.y <- expr_y
    tmp.dat <- cbind(mu.x, mu.y)
    tmp.var <- LPEseq.var(tmp.dat, n.bin=n.bin, d=d, df=df)
    var.y <- LPEseq.predict.var( mu.y, tmp.var)
    std.dev <- sqrt(var.x + var.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x == 2 & n.y == 2){
    basevar.x <- lpe.var(expr_x, n.bin=n.bin, df=df)
    basevar.y <- lpe.var(expr_y, n.bin=n.bin, df=df)
    var.x <- basevar.x[,2]
    var.y <- basevar.y[,2]
    mu.x <- basevar.x[,1]
    mu.y <- basevar.y[,1]
    std.dev <- sqrt(var.x+var.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x == 2 & n.y > 2){
    basevar.x <- lpe.var(expr_x, n.bin=n.bin, df=df)
    var.x <- basevar.x[,2]
    mu.x <- basevar.x[,1]
    basevar.y <- lpe.var(expr_y, n.bin=n.bin, df=df)
    mu.y <- basevar.y[,1]
    var.y <- basevar.y[,2]
    std.dev <- sqrt(var.x + (pi/2)*var.y/n.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x > 2 & n.y == 1){
    basevar.x <- lpe.var(expr_x, n.bin=n.bin, df=df)
    var.x <- basevar.x[,2]
    mu.x <- basevar.x[,1]
    mu.y <- expr_y
    tmp.dat <- cbind(mu.x, mu.y)
    tmp.var <- LPEseq.var(tmp.dat, n.bin=n.bin, d=d, df=df)
    var.y <- LPEseq.predict.var(mu.y, tmp.var)
    std.dev <- sqrt(var.x + var.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x > 2 & n.y == 2){
    basevar.x <- lpe.var(expr_x, n.bin=n.bin, df=df)
    var.x <- basevar.x[,2]
    mu.x <- basevar.x[,1]
    basevar.y <- lpe.var(expr_y, n.bin=n.bin, df=df)
    mu.y <- basevar.y[,1]
    var.y <- basevar.y[,2]
    std.dev <- sqrt((pi/2)*var.x/n.x + var.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
  if(n.x > 2 & n.y > 2){
    basevar.x <- lpe.var(expr_x, n.bin=n.bin, df=df)
    basevar.y <- lpe.var(expr_y, n.bin=n.bin, df=df)
    var.x <- basevar.x[,2]
    var.y <- basevar.y[,2]
    mu.x <- basevar.x[,1]
    mu.y <- basevar.y[,1]
    std.dev <- sqrt((pi/2)*var.x/n.x + (pi/2)*var.y/n.y)
    z.stats <- (mu.y-mu.x)/std.dev
    p.val <- as.numeric(2*(1-pnorm(abs(z.stats))))
    ix.na <- which(mu.x == 0 & mu.y == 0)
    p.val[ix.na] <- NA
    adj.p.fdr <- p.adjust(p.val, method="fdr")
    data.out <- data.frame(
      mu.x = mu.x,
      mu.y = mu.y,
      pooled.std.dev = std.dev,
      z.stats = z.stats,
      p.value = p.val,
      q.value = adj.p.fdr
    )
    return(data.out)
  }
}


#' A function Normalising input data according to their total count values
#' @param expr_dat  count value matrix (should not be count values; FPKM or RPKM can be applied)
#' @param method  one of c("mean", "median") for summary of the column sums (default "mean")
#' @details See the reference paper
#' @return \item{matrix}{log2 transformed normalized matrix}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
LPEseq.normalise <- function(expr_dat, method="mean"){
  colSum <- apply(expr_dat, 2, sum)
  if(method=="mean"){
    meanVec <- colSum/mean(colSum)
  }else if(method=="median"){
    meanVec <- colSum/median(colSum)
  }else{
      stop("method should be one of \"mean\" or \"median\"")
  }
  normData <- expr_dat
  for(i in 1:ncol(expr_dat)){
    normData[,i] <- expr_dat[,i]/meanVec[i]
  }
  return(log2(normData+1))
}



#' A function performing MA transformation of data without replicate
#' @param dat  expression data with non-replicate per each condition
#' @details See the reference paper
#' @return \item{matrix}{MA matrix}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
LPEseq.matrans <- function(dat){
  if(ncol(dat)!=2){
    stop("The number of column should be two")
  }
  MA <- dat
  colnames(MA) = c("A", "M")
  MA[,1] <- apply(dat, 1, mean)
  MA[,2] <- dat[,1] - dat[,2]
  return(MA)
}


#' A function performing MA transformation of data with replicates
#' @param y  expression data with replicates per each condition
#' @details See the reference paper
#' @return \item{matrix}{MA matrix}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
am.trans <- function (y){
  if (ncol(y) > 5)
    y <- y[, sample(1:ncol(y), 5)]
  n <- ncol(y)
  if (n < 2) {
    stop("There are no replicated arrays!")
  }
  A <- c()
  M <- c()
  cc <- permute(1:n)
  for (i in 1:(n - 1)) {
    A <- c(A, c((y + y[, cc[i, ]])/2), recursive = TRUE)
    M <- c(M, c(y - y[, cc[i, ]]), recursive = TRUE)
  }
  return(cbind(A, M))
}


#' A function performing permutation used in am.trans function
#' @param a  numeric value indicating the number of samples
#' @details See the reference paper
#' @return \item{numeric}{permuted number}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
permute <- function(a){
  aa <- matrix(NA, length(a) - 1, length(a))
  for (i in 1:(length(a) - 1)) {
    aa[i, ] <- a[c((i + 1):length(a), 1:i)]
  }
  return(aa)
}


#' A function removing outliers
#' @param expr_dat  expression data with non-replicate per each condition
#' @param d  expression difference between condition (log2 scaled)
#' @details See the reference paper
#' @return \item{matrix}{intput matrix outlier-removed}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
LPEseq.outlier <- function(expr_dat, d=d){
  ix.outlier <- which(abs(apply(expr_dat, 1, diff)) >= d)
  expr_dat_sub <- expr_dat[-ix.outlier,]
  return(expr_dat_sub)
}


#' A function evaluating local pooled error variance with non-replicated data
#' @param y  expression data with non-replicate per each condition
#' @param n.bin  the number of quantile bins
#' @param df  degrees of freedom of smooth.spline() function
#' @param d  expression difference of LPEseq.outlier() function
#' @details See the reference paper
#' @return \item{list}{output of smooth.spline() function}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
LPEseq.var <- function(y, n.bin=n.bin, df=df, d=d){
  y <- LPEseq.outlier(y, d=d)
  qnt <- 1/n.bin
  AM <- am.trans(y)
  A <- AM[,1]
  M <- AM[,2]
  mu.x <- y[,1]
  mu.y <- y[,2]
  quantile.A <- quantile(A, probs = seq(0,1,qnt), na.rm=TRUE)
  quan.n <- length(quantile.A)-1
  var.M <- rep(0, length=quan.n-1)
  medianAs <- rep(0, length=quan.n-1)
  if(sum(A==min(A)) > (qnt*length(A))){
    tmpA <- A[!(A==min(A))]
    quantile.A <- c(min(A), quantile(tmpA, probs=seq(qnt, 1, qnt), na.rm=TRUE))
  }
  for(i in 2:(quan.n+1)){
    n.i <- length(!is.na(M[A>=quantile.A[i-1]&A<quantile.A[i]]))
    if(n.i > 1){
      mult.factor <- 0.5*((n.i-0.5)/(n.i-1))
      tmp.M <- M[A>=quantile.A[i-1] & A<quantile.A[i]]
      var.M[i-1] <- mult.factor * var(tmp.M)
      medianAs[i-1] <- median(A[A>=quantile.A[i-1] & A<quantile.A[i]], na.rm=TRUE)
    }
  }
  if(any(is.na(var.M))){
    for(i in (quan.n-1):1){
      if(is.na(var.M[i])){
        var.M[i] <- ifelse(!is.na(var.M[i+1]), var.M[i+1])
      }
    }
  }
  var.M[1:which(var.M==max(var.M))] <- max(var.M)
  base.var <- cbind(A=medianAs, var.M=var.M)
  if(IQR(base.var[,1])<=0){
    tol=0.00001
    sm.spline <- smooth.spline(base.var[,1], base.var[,2], df=df, tol=tol)
  }
  else{
    sm.spline <- smooth.spline(base.var[,1], base.var[,2], df=df)
  }
  return(sm.spline)
}


#' A function predicting per-gene variance using LPEseq.var() output
#' @param gene_expr  gene expression whose variance to be estimated
#' @param var.spline  output of LPEseq.var() function
#' @details See the reference paper
#' @return \item{numeric}{estimated variance}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
LPEseq.predict.var <- function(gene_expr, var.spline){
  tmp.predict <- fixbounds.predict.smooth.spline(var.spline, gene_expr)$y
  min.var <- min(var.spline$y)
  if(any(tmp.predict < min.var)){
    tmp.predict[tmp.predict < min.var] <- min.var
  }
  return(tmp.predict)
}



#' A function evaluating LPE variance of data with replicates
#' @param y  gene expressino whose variance to be estimated
#' @param n.bin  the number of bins
#' @param df  degree of freedom in smoothing (arg of smooth.spline())
#' @details See the reference paper
#' @return \item{matrix}{evaluated LPE}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
lpe.var <- function(y, n.bin=n.bin, df=df){
  qnt <- 1/n.bin
  AM <- am.trans(y)
  A <- AM[,1]
  M <- AM[,2]
  median.y <- apply(y, 1, median)
  quantile.A <- quantile(A, probs = seq(0,1,qnt), na.rm=TRUE)
  quan.n <- length(quantile.A)-1
  var.M <- rep(0, length=quan.n-1)
  medianAs <- rep(0, length=quan.n-1)
  if(sum(A==min(A)) > (qnt*length(A))){
    tmpA <- A[!(A==min(A))]
    quantile.A <- c(min(A), quantile(tmpA, probs=seq(qnt, 1, qnt), na.rm=TRUE))
  }
  for(i in 2:(quan.n+1)){
    n.i <- length(!is.na(M[A>=quantile.A[i-1]&A<quantile.A[i]]))
    if(n.i > 1){
      mult.factor <- 0.5*((n.i-0.5)/(n.i-1))
      var.M[i-1] <- mult.factor * var(M[A>=quantile.A[i-1] & A<quantile.A[i]], na.rm=TRUE)
      medianAs[i-1] <- median(A[A>=quantile.A[i-1] & A<quantile.A[i]], na.rm=TRUE)
    }
  }
  if(any(is.na(var.M))){
    for(i in (quan.n-1):1){
      if(is.na(var.M[i])){
        var.M[i] <- ifelse(!is.na(var.M[i-1]), mean(var.M[i+1], var.M[i-1]), var.M[i+1])
      }
    }
  }
  var.M[1:which(var.M==max(var.M))] <- max(var.M)
  base.var <- cbind(A=medianAs, var.M=var.M)
  sm.spline <- smooth.spline(base.var[,1], base.var[,2], df=df)
  min.Var <- min(base.var[,2])
  var.genes <- fixbounds.predict.smooth.spline(sm.spline, median.y)$y
  if(any(var.genes < min.Var))
    var.genes[var.genes < min.Var] <- min.Var
  basevar.step1 <- cbind(A=median.y, var.M=var.genes)
  return(basevar.step1)
}



#' A function predicting LPE variance using the output of lpe.var()
#' @param object  output of lpe.var()
#' @param x       a numeric value (or a vector) whose variance be estimated
#' @param deriv   for missing input handling (default = 0)
#' @details See the reference paper
#' @return \item{numeric}{estimated variance}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
fixbounds.predict.smooth.spline <- function (object, x, deriv = 0){
  if (missing(x)) {
    if (deriv == 0) {
      return(object[c("x", "y")])
    }
    else {
      x <- object$x
    }
  }
  if (is.null(object)) {
    stop("not a valid smooth.spline object")
  }
  else {
    out <- predict(object, x, deriv)
    maxpredY <- object$y[object$x == max(object$x)]
    out$y[out$x > max(object$x)] <- maxpredY
    minpredY <- object$y[object$x == min(object$x)]
    out$y[out$x < min(object$x)] <- minpredY
    invisible(out)
  }
}



#' A function generating simulation data
#' @param n.gene  the number of genes (default 20000)
#' @param n.cond  the number of conditions (default 2)
#' @param n.deg  the number of differentially expressed genes (DEGs) (default 0)
#' @param eff  count difference of DEGs between conditions  (default 1000)
#' @param n.rep  the number of replicates in each condition (default 3)
#' @param disp  dispersion parameter for NB distribution (default 0.25)
#' @details See the reference paper
#' @return \item{matrix}{simulated datasets}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
generateData <- function(n.gene=20000, n.cond=2, n.deg=0, eff=1000, n.rep=3, disp=0.25){
  tmp.size <- 1/disp
  tmp.mu <- NBparameter[sample(1:nrow(NBparameter), n.gene, replace=TRUE),1]
  tmp.dat <- matrix(0, nrow=n.gene, ncol=n.cond*n.rep+1)
  for(i in 1:n.gene){
    tmp.dat[i,1:(n.cond*n.rep)] <- rnbinom(n.cond*n.rep, size=tmp.size, mu=tmp.mu[i])
  }
  if(n.deg!=0){
    deg.ix <- sample(1:n.gene, n.deg)
    deg.ix.u <- sample(deg.ix, round(n.deg/2))
    deg.ix.l <- deg.ix[!deg.ix %in% deg.ix.u]
    for(i in 1:length(deg.ix.u)){
      tmp.dat[deg.ix.u[i], (n.rep+1):(n.cond*n.rep)] <- rnbinom(n.rep, size=tmp.size, mu=tmp.mu[deg.ix.u[i]]+eff)
    }
    for(i in 1:length(deg.ix.l)){
      tmp.dat[deg.ix.l[i], 1:n.rep] <- rnbinom(n.rep, size=tmp.size, mu=tmp.mu[deg.ix.l[i]]+eff)
    }
    tmp.dat[deg.ix, n.cond*n.rep+1] <- 1
  }
  rowVec <- paste("gene", 1:n.gene, sep="_")
  colVec <- c(paste("condition1", 1:n.rep, sep="."), paste("condition2", 1:n.rep, sep="."), "DEG")
  rownames(tmp.dat) <- rowVec
  colnames(tmp.dat) <- colVec
  return(tmp.dat)
}

#' A function plotting average intensity versus variance
#' @param data  expression data matrix
#' @param avg  measure for average value across the columns of the data matrix (default "mean")
#' @param w.value  parameter for lowess weight value (default 0.1)
#' @param logged  logical variable indicating log tranformation (default True)
#' @details See the reference paper
#' @return \item{plot}{average-variance plot}
#' @export
#' @author Jungsoo Gim
#' @references J. Gim and T. Park, 2014 (http://bibs.snu.ac.kr/software/LPEseq)
AVplot <- function(data, avg="mean", w.value=0.1, logged=TRUE){
  if(logged==1){
    tmp.avg <- apply(data, 1, avg)
    tmp.avg <- cbind(tmp.avg, apply(data, 1, var))
  }else{
    tmp.avg <- apply(log2(data+1), 1, avg)
    tmp.avg <- cbind(tmp.avg, apply(log2(data+1), 1, var))
  }
  plot(tmp.avg, xlab="Average Intensity", ylab="Variance", pch="*", cex=1, col="grey50", main="Average-Variance plot")
  lines(lowess(tmp.avg, f=w.value), col="red", lwd=3)
}

