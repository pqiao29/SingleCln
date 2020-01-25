indptSTM3_select <- function(dat1, dat2, dat3, n){
 
  K <- max(dat1[, 4])
  
  dat1_tilled <- tilling(dat1, n)
  dat2_tilled <- tilling(dat2, n)
  dat3_tilled <- tilling(dat3, n)
  
  all_models <- t(as.matrix(expand.grid(rep(list(1:0), K))))
  all_models <- rbind(rep(1, ncol(all_models)), all_models)
  return(Select3(dat1_tilled, dat2_tilled, dat3_tilled, n, all_models))
}