grad_hess <- function(func, x, method.args = list(), ...) {
    out <- numDeriv::genD(func, x, method.args = method.args, ...)
    p <- length(x)
    gr <- out$D[1:p]
    he <- matrix(0, p, p)
    ut <- upper.tri(he, diag = TRUE)
    he[ut] <- out$D[-seq_len(p)]
    lt <- lower.tri(he, diag = TRUE)
    he[lt] <- t(he)[lt]
    list(grad = gr, hess = he)
}
