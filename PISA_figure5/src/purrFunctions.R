tTest <- function(data){
  
  test_01 <- t.test(x = data$`00um`, y = data$`01um`)
  test_20 <- t.test(x = data$`00um`, y = data$`20um`)
  pval_01 <- test_01$p.value
  pval_20 <- test_20$p.value
  delta_01 <- mean(data$`01um`) - mean(data$`00um`)
  delta_20 <- mean(data$`20um`) - mean(data$`00um`)
  sumSD_01 <- sd(data$`01um`) + sd(data$`00um`)
  sumSD_20 <- sd(data$`20um`) + sd(data$`00um`)
  
  return(data.frame(
    delta_01,
    delta_20,
    pval_01,
    pval_20,
    sumSD_01,
    sumSD_20))
}