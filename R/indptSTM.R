tilling <- function(data, n){
  cordi_x <- data[, 2]
  cordi_y <- data[, 3]
  
  Mx <- max(cordi_x) + 1.e-8
  mx <- min(cordi_x) - 1.e-8
  My <- max(cordi_y) + 1.e-8
  my <- min(cordi_y) - 1.e-8
  l.x <- (Mx - mx)/n
  l.y <- (My - my)/n
  
  loc_x <- findInterval(cordi_x, seq(mx, Mx, l.x))
  loc_y <- findInterval(cordi_y, seq(my, My, l.y))
  tmp_tile  <- (loc_y - 1) * n + loc_x
  
  cbind(data[, c(1, 4)], tmp_tile)
}

indptSTM3 <- function(data1, data2, data3, n){
  dat1 <- tilling(data1, n)
  dat2 <- tilling(data2, n)
  dat3 <- tilling(data3, n)
  return(indptSTM_cpp3(dat1, dat2, dat3, n))
}