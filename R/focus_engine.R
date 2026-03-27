#' Low-level engine for focus estimation and inference
#'
#' Estimate and perform inference on a scalar function of a supplied
#' parameter vector, assumed to be the maximum likelihood estimate of
#' the model parameters, using the inverse of the expected information
#' matrix at that estimaate, and, when needed, auxiliary quantities
#' for bias correction.
#'
#' @param theta Numeric parameter vector, assumed to be the maximum
#'     likelihood estimate of the model parameters.
#' @param V the inverse of the expected information at `theta`.
#' @param on A function specifying the scalar parameter of interest.
#'     It must take the parameter vector as first argument and return
#'     a single numeric value. The default, `function(theta, ...)
#'     theta[1]`, focuses on the first component of the parameter
#'     vector.
#' @param correction Character string specifying the bias correction
#'     method.  One of `"no"`, `"median"`, or `"mean"`.
#' @param alpha Numeric scalar in `(0, 1)`. Target miscoverage level
#'     for the Wald confidence interval.
#' @param P Optional list of matrices \eqn{P_t}, one matrix per
#'     parameter, as defined in Kosmidis (2014, expression
#'     (5)). Required when `correction = "median"` or `correction =
#'     "mean"`.
#' @param Q Optional list of matrices \eqn{Q_t}, one matrix per
#'     parameter, as defined in Kosmidis (2014, expression
#'     (5)). Required when `correction = "median"` or `correction =
#'     "mean"`.
#' @param on_gradient Optional function returning the gradient of
#'     `on`; if `NULL` (default), the gradient is computed numerically using
#'     [numDeriv::grad()].
#' @param on_hessian Optional function returning the Hessian of `on`;
#'     if `NULL` (default), the Hessian is computed numerically using
#'     [numDeriv::hessian()].
#' @param ... Additional arguments passed to `on`, `on_gradient`, and
#'     `on_hessian`.
#'
#' @return
#' A list with the same components as [focus()] for Wald inference:
#' `estimate`, `se`, `ci_type`, and `confint`.
#'
#' @details
#' `focus_engine()` assumes that the supplied `theta` is the maximum likelihood
#' estimate of the model parameters, and that any requested bias correction is
#' to be based on this estimator and its supplied covariance matrix `V`.
#'
#' `focus_engine()` is a low-level helper intended for advanced use and
#' therefore performs only minimal validation of its inputs. In particular,
#' when `correction` is `"mean"` or `"median"`, `P` and `Q` are assumed to be
#' supplied in the format required by the computation.
#'
#' Confidence intervals returned by `focus_engine()` are Wald-type intervals
#' using the delta-method and the supplied covariance matrix `V`.
#'
#' @seealso [focus()], [ci_control()]
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196. \doi{10.1002/wics.1296}.
#'
#' @export
focus_engine <- function(theta,
                         V,
                         on = function(theta, ...) theta[1],
                         correction = "median",
                         alpha = 0.05,
                         P = NULL,
                         Q = NULL,
                         on_gradient = NULL,
                         on_hessian = NULL, ...) {
    theta <- as.numeric(theta)
    stopifnot(length(alpha) == 1L, is.numeric(alpha), !is.na(alpha), alpha > 0, alpha < 1)
    stopifnot(is.null(on_gradient) || is.function(on_gradient))
    stopifnot(is.null(on_hessian) || is.function(on_hessian))
    correction <- match.arg(correction, c("no", "median", "mean"))
    d1_psi <- if (is.null(on_gradient)) numDeriv::grad(on, theta, ...) else on_gradient(theta, ...)
    d1_psi <- as.numeric(d1_psi)
    hat_psi <- on(theta, ...)
    muffin <- drop(V %*% d1_psi)
    var_psi <- sum(d1_psi * muffin)
    if (identical(correction, "no")) {
        out <- hat_psi
    } else {
        stopifnot(!is.null(P), !is.null(Q))
        d2_psi <- if (is.null(on_hessian)) numDeriv::hessian(on, theta, ...) else on_hessian(theta, ...)
        constant_on <- all(d1_psi == 0) && all(d2_psi == 0)
        if (constant_on) {
            out <- hat_psi
        } else {
            A <- vapply(seq_along(P), function(t) 0.5 * sum(diag(V %*% (P[[t]] + Q[[t]]))), numeric(1))
            bias <- -drop(V %*% A)
            mean_b <- sum(d1_psi * bias) + 0.5 * sum(d2_psi * V)
            if (identical(correction, "mean")) {
                out <- hat_psi - mean_b
            }
            if (identical(correction, "median")) {
                cheese <- lapply(seq_along(P), function(k) muffin[k] * (P[[k]] / 3 + Q[[k]] / 2))
                cheese <- - Reduce("+", cheese) + 0.5 * d2_psi
                skew <- sum((cheese %*% muffin) * muffin) / var_psi
                out <- hat_psi - mean_b + skew
            }
        }
    }
    se <- sqrt(var_psi)
    confint <- out + c(-1, 1) * qnorm(1 - alpha / 2) * se
    names(confint) <- c("lower", "upper")
    attr(confint, "alpha") <- alpha
    out <- unname(out)
    attr(out, "correction") <- correction
    list(estimate = out, se = se, ci_type = "wald", confint = confint)
}
