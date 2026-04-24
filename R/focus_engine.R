#' Estimate components for [focus_engine()]
#'
#' Estimate the model-side components used by [focus_engine()] through repeated
#' simulation at a supplied parameter vector.
#'
#' @param theta Numeric parameter vector at which the components are estimated.
#' @param loglik A function returning the log-likelihood evaluated at a supplied
#'   parameter vector and dataset. It is expected to have an interface of the
#'   form `loglik(theta, data, ...)` and to return a single numeric value.
#' @param score Optional function returning the gradient of the log-likelihood
#'   at a supplied parameter vector and dataset. It is expected to have an
#'   interface of the form `score(theta, data, ...)`. If `NULL` (default), the
#'   score is computed numerically from `loglik`.
#' @param information Optional function returning the observed information
#'   matrix, namely the negative Hessian of the log-likelihood, at a supplied
#'   parameter vector and dataset. It is expected to have an interface of the
#'   form `information(theta, data, ...)`. If `NULL` (default), the information
#'   matrix is computed numerically from `loglik`.
#' @param simulate A function that simulates one dataset at the supplied
#'   parameter vector `theta`. It is expected to have an interface of the form
#'   `simulate(theta, ...)` and to return one simulated dataset at each call.
#' @param nsim Integer. The number of simulated datasets used to estimate the
#'   components.
#' @param parallelize Logical. If `TRUE`, use
#'   [future.apply::future_lapply()] for the repeated derivative
#'   calculations. This requires the suggested package
#'   \pkg{future.apply}.
#' @param ... Additional arguments passed to `simulate`, `loglik`, and, when
#'   supplied, to `score` and `information`.
#'
#' @return
#' A named list suitable for the `components` argument of [focus_engine()].
#' It contains `V`, `P`, and `Q`.
#'
#' @details
#' The supplied `simulate` function is called repeatedly as
#' `simulate(theta, ...)`, once for each simulated dataset.
#'
#' The supplied `loglik` function is then used on each simulated dataset
#' through calls of the form `loglik(theta, data = data, ...)`.
#'
#' The returned `V` is obtained by averaging the observed information matrices
#' across the simulated datasets and inverting the result.
#'
#' The returned `P` and `Q` are estimated by averaging the corresponding
#' simulation-based quantities across the simulated datasets.
#'
#' @seealso [focus_engine()]
#'
#' @export
estimate_focus_components <- function(theta,
                                      loglik,
                                      score = NULL,
                                      information = NULL,
                                      simulate,
                                      nsim = 1000,
                                      parallelize = FALSE,
                                      ...) {
    theta <- as.numeric(theta)
    if (parallelize && !requireNamespace("future.apply", quietly = TRUE)) {
        stop("Package `future.apply` is required when `parallelize = TRUE`.")
    }
    if (is.null(score) && is.null(information)) {
        simu_one <- function(i) {
            data <- simulate(theta, ...)
            ders <- grad_hess(loglik, theta, data = data, ...)
            list(S = ders$grad,
                 I = -ders$hess)
        }
    } else {
        s_fun <- if (is.null(score)) {
                     function(x, data, ...) numDeriv::grad(loglik, x, data = data, ...)
                 } else {
                     score
             }
        i_fun <- if (is.null(information)) {
                     function(x, data, ...) -numDeriv::hessian(loglik, x, data = data, ...)
                 } else {
                     information
                 }
        simu_one <- function(i) {
            data <- simulate(theta, ...)
            list(S = s_fun(theta, data, ...),
                 I = i_fun(theta, data, ...))
        }
    }

    derivatives <- if (parallelize) {
        future.apply::future_lapply(
            seq_len(nsim),
            simu_one,
            future.seed = TRUE
        )
    } else {
        lapply(seq_len(nsim), simu_one)
    }

    Ihat <- Reduce("+", lapply(derivatives, `[[`, "I")) / nsim
    out <- list(V = solve(Ihat))
    sc <- mean(diag(out$V))
    if (!is.finite(sc) || sc < 1e-6 || sc > 1) sc <- 1
    out$P <- lapply(seq_along(theta), function(t) {
        Reduce(
            "+",
            lapply(derivatives, function(der) tcrossprod(der$S) * der$S[t] * sc)
        ) / (sc * nsim)
    })

    out$Q <- lapply(seq_along(theta), function(t) {
        Reduce(
            "+",
            lapply(derivatives, function(der) -der$I * der$S[t] * sc)
        ) / (sc * nsim)
    })

    out
}

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
#' @param components A named list containing the model-side quantities needed
#'     by the correction. At present, this can include `V`, `P`, and `Q`,
#'     where `V` is the inverse of the expected information at `theta`, and
#'     `P` and `Q` are lists of matrices \eqn{P_t} and \eqn{Q_t} in
#'     Kosmidis (2014, expression (5)); see below for when each is required.
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
#' A list similar to the output of [focus()], with components:
#' \describe{
#'   \item{`estimate`}{Numeric scalar, the estimate of the quantity defined by `on`.}
#'   \item{`se`}{Numeric scalar, the delta-method standard error.}
#'   \item{`correction`}{Character string recording the bias correction method used.}
#'   \item{`on`}{A list containing the supplied `on`, `on_gradient`, and
#'     `on_hessian` functions.}
#'   \item{`dots`}{A list with the additional arguments supplied through `...`.}
#'   \item{`call`}{The matched call to `focus_engine()`.}
#' }
#' The returned object has class `"focus_engine_list"` and inherits from
#' `"focus_list"`.
#'
#' @details
#' `focus_engine()` assumes that the supplied `theta` is either the maximum
#' likelihood estimator of the model parameters or a reduced mean-bias
#' estimator of them, as indicated by `estimator`.
#'
#' If `estimator = "ML"`, then the first-order bias term required for
#' `correction = "mean"` or `correction = "median"` is computed from
#' the supplied `components$P`, `components$Q`, and `components$V`. If
#' `estimator = "meanBR"`, then that first-order bias term is taken to be
#' zero.
#'
#' If `estimator = "meanBR"` and `correction = "mean"`, then `P` and
#' `Q` are not needed, because the correction depends only on the
#' Hessian of `on` and the supplied covariance matrix `components$V`.
#'
#' `focus_engine()` is a low-level helper intended for advanced use and
#' therefore performs only minimal validation of its inputs. In particular,
#' `components$V` is always assumed to be supplied, and whenever they are
#' needed by the requested correction, `components$P` and `components$Q` are
#' assumed to be supplied in the format required by the computation.
#'
#' Standard errors are computed using the delta method, with covariance
#' matrix and gradients evaluated at the supplied `theta`. They are therefore
#' not re-evaluated at the corrected estimates.
#'
#' Confidence intervals can be constructed from the returned `estimate`
#' and `se`, for example through the normal approximation.
#'
#' @seealso [focus()]
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196. \doi{10.1002/wics.1296}.
#'
#' @export
focus_engine <- function(theta,
                         components,
                         on = function(theta, ...) theta[1],
                         correction = "median",
                         estimator = "ML",
                         on_gradient = NULL,
                         on_hessian = NULL, ...) {
    cl <- match.call()
    dots <- list(...)
    theta <- as.numeric(theta)
    correction <- match.arg(correction, c("no", "median", "mean"))
    estimator <- match.arg(estimator, c("ML", "meanBR"))
    stopifnot(is.list(components), !is.null(components$V))
    V <- components$V
    P <- components$P
    Q <- components$Q
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
    out <- unname(out)
    out <- list(
        call = cl,
        on = list(on = on,
                  on_gradient = on_gradient,
                  on_hessian = on_hessian),
        dots = dots,
        correction = correction,
        estimate = out,
        se = se
    )
    class(out) <- c("focus_engine_list", "focus_list", class(out))
    out
}
