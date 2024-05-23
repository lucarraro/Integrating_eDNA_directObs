eval_wilcox_paired <- function(df, subset, groupby=NULL, val=NULL){

  if (is.null(groupby)) groupby <- "frac_e"
  if (is.null(val)) val <- c(0.25,0.5,0.75,1)

  wt <- t(combn(val,2))
  wt <- cbind(wt, numeric(nrow(wt)))

  for (i in 1:nrow(wt)){
    tmp <- wilcox.test(df$MAE_p[subset & df[[groupby]]==wt[i,1]],
                       df$MAE_p[subset & df[[groupby]]==wt[i,2]],
                       paired=TRUE)
    wt[i,3] <- tmp$p.value
  }

  invisible(wt)
}

draw_significance <- function(wt, add, coeff=NULL, height=NULL, offset=NULL){

  if (is.null(coeff)) coeff=4
  if (is.null(height)) height=10
  if (is.null(offset)) offset=0.3
  significant.pair <- which(wt[,3]<0.001)
  X1 <- add + coeff*wt[significant.pair, 1] #+ 0.1
  X2 <- add + coeff*wt[significant.pair, 2] #- 0.1

 if (length(X1)>0){
  for (i in 1:length(X1)){lines(c(X1[i],X2[i]),
                                (height+offset-offset*i)*c(1,1))}}

}

best_median <- function(df, subset, add, groupby=NULL, val=NULL,
                        split.median=FALSE){

  if (is.null(groupby)) groupby <- "frac_e"
  if (is.null(val)) val <- c(0.25,0.5,0.75,1)
  if (split.median==FALSE){
    mm <- numeric(length(val))
    for (i in 1:length(val)) mm[i] <- median(df$MAE_p[subset & df[[groupby]]==val[i]])
  } else {
    mm <- numeric(2)
    mm[1] <- median(df$MAE_p[subset & df[[groupby]] <= median(df[[groupby]])])
    mm[2] <- median(df$MAE_p[subset & df[[groupby]] > median(df[[groupby]])])
  }
  x <- which.min(mm)
  y <- mm[x]
  points(x+add,y,bg="yellow",pch=23,cex=0.7)
}
