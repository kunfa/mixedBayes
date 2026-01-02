make_Zblock <- function(k, n1,slope=TRUE) {
  if(slope){
    k_1 = c(1:k)
    k_1 = k_1 - mean(k_1)
    z = cbind(rep(1,k),k_1)
  } else {
    z = as.matrix(rep(1,k))
  }
  c <- ncol(z)
  Zblock <- matrix(0, nrow = n1 * k, ncol = n1 * c)
  for (i in 1:n1) {
    rows <- ((i - 1) * k + 1):(i * k)
    cols <- ((i - 1) * c + 1):(i * c)
    Zblock[rows, cols] <- z
  }
  Zblock
}
