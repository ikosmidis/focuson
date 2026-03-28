#' Low-level engine for focus estimation and inference
#'
#' Estimate and perform inference on a scalar function of a supplied
#' parameter vector, assumed to be either the maximum likelihood estimate
#' of the model parameters or a reduced mean-bias estimator, using the
#' inverse of the expected information matrix at that estimate, and, when
#' needed, auxiliary quantities for bias correction.
#'
#' @param theta Numeric parameter vector, assumed to be either the maximum
#'     likelihood estimate of the model parameters or a reduced mean-bias
#'     estimator of them.
#' @param V the inverse of the expected information at `theta`.
#' @param on A function specifying the scalar parameter of interest.
#'     It must take the parameter vector as first argument and return
#'     a single numeric value. The default, `function(theta, ...)
#'     theta[1]`, focuses on the first component of the parameter
#'     vector.
#' @param correction Character string specifying the bias correction
#'     method.  One of `"no"`, `"median"`, or `"mean"`.
#' @param estimator Character string specifying the estimator represented by
#'     the supplied `theta`. One of `"ML"` (maximum likelihood, default) or
#'     `"meanBR"` (reduced mean-bias estimator).
#' @param alpha Numeric scalar in `(0, 1)`. Target miscoverage level
#'     for the Wald confidence interval.
#' @param P Optional list of matrices \eqn{P_t}, one matrix per
#'     parameter, as defined in Kosmidis (2014, expression
#'     (5)). Required when `correction = "median"`, and also when
#'     `correction = "mean"` and `estimator = "ML"`.
#' @param Q Optional list of matrices \eqn{Q_t}, one matrix per
#'     parameter, as defined in Kosmidis (2014, expression
#'     (5)). Required when `correction = "median"`, and also when
#'     `correction = "mean"` and `estimator = "ML"`.
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
#' `focus_engine()` assumes that the supplied `theta` is either the maximum
#' likelihood estimator of the model parameters or a reduced mean-bias
#' estimator of them, as indicated by `estimator`.
#'
#' If `estimator = "ML"`, then the first-order bias term required for
#' `correction = "mean"` or `correction = "median"` is computed from
#' the supplied `P`, `Q`, and `V`. If `estimator = "meanBR"`, then
#' that first-order bias term is taken to be zero.
#'
#' If `estimator = "meanBR"` and `correction = "mean"`, then `P` and
#' `Q` are not needed, because the correction depends only on the
#' Hessian of `on` and the supplied covariance matrix `V`.
#'
#' `focus_engine()` is a low-level helper intended for advanced use and
#' therefore performs only minimal validation of its inputs. In particular,
#' whenever they are needed by the requested correction, `P` and `Q` are
#' assumed to be supplied in the format required by the computation.
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
                         estimator = "ML",
                         alpha = 0.05,
                         P = NULL,
                         Q = NULL,
                         on_gradient = NULL,
                         on_hessian = NULL, ...) {
    theta <- as.numeric(theta)
    correction <- match.arg(correction, c("no", "median", "mean"))
    estimator <- match.arg(estimator, c("ML", "meanBR"))
    d1_psi <- if (is.null(on_gradient)) numDeriv::grad(on, theta, ...) else on_gradient(theta, ...)
    d1_psi <- as.numeric(d1_psi)
    hat_psi <- on(theta, ...)
    muffin <- drop(V %*% d1_psi)
    var_psi <- sum(d1_psi * muffin)
    if (identical(correction, "no")) {
        out <- hat_psi
    } else {
        d2_psi <- if (is.null(on_hessian)) numDeriv::hessian(on, theta, ...) else on_hessian(theta, ...)
        constant_on <- all(d1_psi == 0) && all(d2_psi == 0)
        if (constant_on) {
            out <- hat_psi
        } else {
            bias <- if (identical(estimator, "ML")) {
                A <- vapply(seq_along(P), function(t) 0.5 * sum(diag(V %*% (P[[t]] + Q[[t]]))), numeric(1))
                -drop(V %*% A)
            } else {
                0
            }
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
