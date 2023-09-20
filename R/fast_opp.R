#' @export
f.grad.rcpp <- function(X, W, H, beta, alpha, y, delta, theta, j){
  f_grad_rcpp(X, W, H, beta, alpha, y, delta, theta, j)
}