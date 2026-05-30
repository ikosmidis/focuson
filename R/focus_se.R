#' @export
focus_on_all <- function(object, ...) {
    UseMethod("focus_on_all")
}

#' @export
focus_se <- function(object, ...) {
    UseMethod("focus_se")
}


#' @export
focus_on_all.focus_list_glm <- function(object, ...) {
    fit <- object$object
    theta <- coef(fit, model = "full")
    if (fit$family$family %in% c("poisson", "binomial")) {
        theta <- theta[names(coef(fit, model = "mean"))]
    }
    on_coordinate <- function(theta, j) {
        theta[j]
    }
    on_coordinate_gradient <- function(theta, j) {
        out <- numeric(length(theta))
        out[j] <- 1
        out
    }
    on_coordinate_hessian <- function(theta, j) {
        matrix(0, nrow = length(theta), ncol = length(theta))
    }
    out <- vapply(seq_along(theta), function(j) {
        res <- focus(
            fit,
            on = on_coordinate,
            correction = object$correction,
            on_gradient = on_coordinate_gradient,
            on_hessian = on_coordinate_hessian,
            j = j
        )
        c(estimate = coef(res), se = res$se)
    }, numeric(2))
    colnames(out) <- names(theta)
    out
}


#' @export
focus_se.focus_list_glm <- function(object, tol_deriv = 1e-10, tol_opt = 1e-06, ...) {
    all_coefs <- focus_on_all(object, ...)
    estimates <- all_coefs["estimate", ]
    ses <- all_coefs["se", ]
    on_estimate <- coef(object)
    on_fun <- object$on$on
    dots <- object$dots
    ders <- do.call(numDeriv::grad,
                    c(list(func = on_fun, x = estimates), dots))
    id <- which.max(abs(ders))
    if (!is.finite(ders[id]) || abs(ders[id]) <= tol_deriv) {
        stop("Could not identify an active parameter to replace in the reference parameterization. ",
             "All derivatives of the focus function at the corrected estimates are zero, ",
             "non-finite, or below `tol_deriv`.")
    }
    fun <- function(co, id) {
        theta <- estimates
        theta[id] <- co
        (do.call(on_fun, c(list(theta), dots)) - on_estimate)^2
    }
    scale <- max(abs(ses[id]), abs(estimates[id]), 1e-02, na.rm = TRUE)
    lims <- estimates[id] + c(-10, 10) * scale
    opt <- optimize(fun, lims, id = id)
    if (!is.finite(opt$objective) || sqrt(opt$objective) > tol_opt) {
        stop("Could not reconstruct the reference parameter vector by replacing parameter `",
            names(estimates)[id], "` over interval [", lims[1], ", ", lims[2], "].")
    }
    estimates[id] <- opt$minimum
    estimates
}
