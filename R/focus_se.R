#' Focus on all model parameters
#'
#' Compute focused estimates, and associated delta-method standard errors, for
#' each coordinate of the model parameter vector used by a [`focus()`] result.
#'
#' @param object An object returned by [`focus()`].
#' @param ... Additional arguments passed to methods.
#'
#' @return
#' For the currently implemented methods, a numeric matrix with one column per
#' parameter. The first row contains the focused estimates and the second row
#' contains their standard errors.
#'
#' @details
#' It applies the same focus correction recorded in `object` to each coordinate
#' of the model parameter vector. The `focus_list_glm` and `focus_engine_list`
#' methods use analytic coordinate gradients and Hessians for those coordinate
#' functions.
#'
#' @seealso [`focus()`], [`focus_se()`]
#'
#' @export
focus_on_all <- function(object, ...) {
    UseMethod("focus_on_all")
}

#' Compute standard errors for focus objects
#'
#' Compute a delta-method standard error for a focus estimate after
#' reconstructing the model parameter vector at which the standard error is
#' evaluated.
#'
#' @param object An object returned by [`focus()`].
#' @param control A list of control parameters used by methods; see
#'   [`focus_se_control()`].
#' @param ... Additional arguments passed to methods.
#'
#' @return
#' For the `focus_list_glm` method, a list with components:
#' \describe{
#'   \item{`se`}{Numeric scalar, the delta-method standard error.}
#'   \item{`theta`}{Numeric vector, the reconstructed model parameter vector.}
#'   \item{`V`}{Covariance matrix evaluated at `theta`.}
#'   \item{`gradient`}{Gradient of the focus function evaluated at `theta`.}
#'   \item{`replace`}{Integer index of the parameter coordinate replaced when
#'     reconstructing `theta`.}
#' }
#'
#' @details
#' For `focus_list_glm` and `focus_engine_list` objects, the method first
#' obtains focused estimates for all model parameters using [`focus_on_all()`].
#' It then selects an active coordinate of the original focus function,
#' replaces that coordinate so that the focus function matches the reported
#' focus estimate, evaluates the model information at the reconstructed
#' parameter vector, and applies the delta method. The `focus_engine_list`
#' method requires a `V_function` function.
#'
#' @seealso [`focus()`], [`focus_on_all()`]
#'
#' @export
focus_se <- function(object, control = list(), ...) {
    UseMethod("focus_se")
}


#' @rdname focus_on_all
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
            se_at = "naive",
            j = j
        )
        c(estimate = coef(res), se = res$se)
    }, numeric(2))
    colnames(out) <- names(theta)
    out
}

#' @rdname focus_on_all
#' @export
focus_on_all.focus_engine_list <- function(object, ...) {
    theta <- object$theta
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
        res <- focus_engine(
            theta = theta,
            components = object$components,
            on = on_coordinate,
            correction = object$correction,
            estimator = object$estimator,
            on_gradient = on_coordinate_gradient,
            on_hessian = on_coordinate_hessian,
            se_at = "naive",
            j = j
        )
        c(estimate = coef(res), se = res$se)
    }, numeric(2))
    colnames(out) <- names(theta)
    out
}


#' @rdname focus_se
#' @details
#' The current `focus_list_glm` and `focus_engine_list` methods recognize the
#' following `control` entries:
#' \describe{
#'   \item{`tol_deriv`}{Numeric scalar. Tolerance used to decide whether a
#'     derivative is large enough for the corresponding parameter to be treated
#'     as active.}
#'   \item{`tol_opt`}{Numeric scalar. Tolerance for the one-dimensional
#'     reconstruction of the replaced parameter coordinate. The optimization
#'     objective is the squared reconstruction error divided by
#'     `max(abs(estimate), 1)^2`, so `tol_opt` is applied to the corresponding
#'     unitless absolute error.}
#' }
#' @export
focus_se.focus_list_glm <- function(object, control = list(), ...) {
    if (!is.list(control)) {
        stop("`control` must be a list.")
    }
    control <- do.call(focus_se_control, control)
    tol_deriv <- control$tol_deriv
    tol_opt <- control$tol_opt
    all_coefs <- focus_on_all(object, ...)
    estimates <- all_coefs["estimate", ]
    ses <- all_coefs["se", ]
    on_estimate <- coef(object)
    on_fun <- object$on$on
    dots <- object$dots
    objective_scale <- max(abs(on_estimate), 1)
    ders <- do.call(numDeriv::grad,
                    c(list(func = on_fun, x = estimates), dots))
    id <- which.max(abs(ders))
    if (!is.finite(ders[id]) || abs(ders[id]) <= tol_deriv) {
        stop("Could not identify an active parameter to replace. ",
             "All derivatives of the focus function at the corrected estimates are zero, ",
             "non-finite, or below `tol_deriv`.")
    }
    fun <- function(co, id) {
        theta <- estimates
        theta[id] <- co
        ((do.call(on_fun, c(list(theta), dots)) - on_estimate) /
             objective_scale)^2
    }
    scale <- max(abs(ses[id]), abs(estimates[id]), 1e-02, na.rm = TRUE)
    lims <- estimates[id] + c(-10, 10) * scale
    opt <- optimize(fun, lims, id = id)
    if (!is.finite(opt$objective) || sqrt(opt$objective) > tol_opt) {
        stop("Could not reconstruct the model parameter vector by replacing parameter `",
            names(estimates)[id], "` over interval [", lims[1], ", ", lims[2], "].")
    }
    estimates[id] <- opt$minimum
    afuns <- enrichwith::get_auxiliary_functions(object$object)
    p_mean <- length(coef(object$object, model = "mean"))
    if (length(estimates) == p_mean) {
        info <- afuns$information(coefficients = estimates)
    } else {
        info <- afuns$information(coefficients = estimates[seq_len(p_mean)],
                                  dispersion = estimates[-seq_len(p_mean)])
    }
    V <- solve(info)
    d1_psi <- do.call(numDeriv::grad,
                      c(list(func = on_fun, x = estimates), dots))
    se <- sqrt(sum(d1_psi * drop(V %*% d1_psi)))
    list(
        se = se,
        theta = estimates,
        V = V,
        gradient = d1_psi,
        replace = id
    )
}

#' @rdname focus_se
#' @param V_function Function returning the covariance matrix at a supplied
#'   model parameter vector. Required for `focus_engine_list` objects.
#' @export
focus_se.focus_engine_list <- function(object, V_function, control = list(), ...) {
    if (missing(V_function) || !is.function(V_function)) {
        stop("`V_function` must be a function.")
    }
    if (!is.list(control)) {
        stop("`control` must be a list.")
    }
    control <- do.call(focus_se_control, control)
    tol_deriv <- control$tol_deriv
    tol_opt <- control$tol_opt
    all_coefs <- focus_on_all(object)
    estimates <- all_coefs["estimate", ]
    ses <- all_coefs["se", ]
    on_estimate <- coef(object)
    on_fun <- object$on$on
    dots <- object$dots
    objective_scale <- max(abs(on_estimate), 1)
    ders <- do.call(numDeriv::grad,
                    c(list(func = on_fun, x = estimates), dots))
    id <- which.max(abs(ders))
    if (!is.finite(ders[id]) || abs(ders[id]) <= tol_deriv) {
        stop("Could not identify an active parameter to replace. ",
             "All derivatives of the focus function at the corrected estimates are zero, ",
             "non-finite, or below `tol_deriv`.")
    }
    fun <- function(co, id) {
        theta <- estimates
        theta[id] <- co
        ((do.call(on_fun, c(list(theta), dots)) - on_estimate) /
             objective_scale)^2
    }
    scale <- max(abs(ses[id]), abs(estimates[id]), 1e-02, na.rm = TRUE)
    lims <- estimates[id] + c(-10, 10) * scale
    opt <- optimize(fun, lims, id = id)
    if (!is.finite(opt$objective) || sqrt(opt$objective) > tol_opt) {
        stop("Could not reconstruct the model parameter vector by replacing parameter `",
            names(estimates)[id], "` over interval [", lims[1], ", ", lims[2], "].")
    }
    estimates[id] <- opt$minimum
    V <- V_function(estimates, ...)
    d1_psi <- do.call(numDeriv::grad,
                      c(list(func = on_fun, x = estimates), dots))
    se <- sqrt(sum(d1_psi * drop(V %*% d1_psi)))
    list(
        se = se,
        theta = estimates,
        V = V,
        gradient = d1_psi,
        replace = id
    )
}

#' Control parameters for corrected focus standard errors
#'
#' Construct the control list used by [`focus_se()`] when reconstructing the
#' model parameter vector for corrected standard errors.
#'
#' @param tol_deriv Numeric scalar. Tolerance used to decide whether a
#'   derivative is large enough for the corresponding parameter to be treated as
#'   active.
#' @param tol_opt Numeric scalar. Tolerance for the one-dimensional
#'   reconstruction of the replaced parameter coordinate. The optimization
#'   objective is the squared reconstruction error divided by
#'   `max(abs(estimate), 1)^2`, so `tol_opt` is applied to the corresponding
#'   unitless absolute error.
#'
#' @return A list with components `tol_deriv` and `tol_opt`.
#'
#' @seealso [`focus_se()`]
#'
#' @export
focus_se_control <- function(tol_deriv = 1e-10, tol_opt = 1e-04) {
    list(tol_deriv = tol_deriv, tol_opt = tol_opt)
}
