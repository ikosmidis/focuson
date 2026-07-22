.focus_core <- function(theta,
                        V,
                        on,
                        correction,
                        components_fun = NULL,
                        on_gradient = NULL,
                        on_hessian = NULL,
                        ...) {
    correction <- match.arg(correction, c("no", "median", "mean"))
    d1_psi <- if (is.null(on_gradient)) {
        numDeriv::grad(on, theta, ...)
    } else {
        on_gradient(theta, ...)
    }
    stopifnot(is.numeric(d1_psi),
              length(d1_psi) == length(theta),
              !anyNA(d1_psi))
    d1_psi <- as.numeric(d1_psi)
    hat_psi <- on(theta, ...)
    muffin <- drop(V %*% d1_psi)
    var_psi <- sum(d1_psi * muffin)
    out <- list(estimate = unname(hat_psi),
                se = sqrt(var_psi))
    if (identical(correction, "no")) {
        return(out)
    }
    d2_psi <- if (is.null(on_hessian)) {
        numDeriv::hessian(on, theta, ...)
    } else {
        on_hessian(theta, ...)
    }
    stopifnot(is.numeric(d2_psi),
              identical(dim(d2_psi), c(length(theta), length(theta))),
              !anyNA(d2_psi))
    constant_on <- all(d1_psi == 0) && all(d2_psi == 0)
    if (constant_on) {
        return(out)
    }
    if (!is.function(components_fun)) {
        stop("`components_fun` must be a function for mean or median correction.")
    }
    correction_components <- components_fun()
    bias <- correction_components$bias
    if (is.null(bias)) {
        stop("`bias` is required for mean or median correction.")
    }
    mean_b <- sum(d1_psi * bias) + 0.5 * sum(d2_psi * V)
    if (identical(correction, "mean")) {
        out$estimate <- unname(hat_psi - mean_b)
        return(out)
    }
    P <- correction_components$P
    Q <- correction_components$Q
    if (is.null(P) || is.null(Q)) {
        stop("`P` and `Q` are required for median correction.")
    }
    cheese <- lapply(seq_along(P), function(k) {
        muffin[k] * (P[[k]] / 3 + Q[[k]] / 2)
    })
    cheese <- -Reduce("+", cheese) + 0.5 * d2_psi
    skew <- sum((cheese %*% muffin) * muffin) / var_psi
    out$estimate <- unname(hat_psi - mean_b + skew)
    out
}
