#' Estimate components for [focus_engine()]
#'
#' Estimate the model-side components used by [focus_engine()] through
#' repeated simulation at a supplied parameter vector.
#'
#' @param theta Numeric parameter vector at which the components are
#'     estimated.
#' @param loglik A function returning the log-likelihood evaluated at
#'     a supplied parameter vector and dataset. It is expected to have
#'     an interface of the form `loglik(theta, data, ...)` and to
#'     return a single numeric value.
#' @param score Optional function returning the gradient of the
#'     log-likelihood at a supplied parameter vector and dataset. It
#'     is expected to have an interface of the form `score(theta,
#'     data, ...)`. If `NULL` (default), the score is computed
#'     numerically from `loglik`.
#' @param information Optional function returning the observed
#'     information matrix, namely the negative Hessian of the
#'     log-likelihood, at a supplied parameter vector and dataset. It
#'     is expected to have an interface of the form
#'     `information(theta, data, ...)`. If `NULL` (default), the
#'     information matrix is computed numerically from `loglik`, or,
#'     if `score` is supplied, as the negative Jacobian of `score`.
#' @param simulate A function that simulates one dataset at the
#'     supplied parameter vector `theta`. It is expected to have an
#'     interface of the form `simulate(theta, ...)` and to return one
#'     simulated dataset at each call.
#' @param nsim Integer. The number of simulated datasets used to
#'     estimate the components.
#' @param parallelize Logical. If `TRUE`, use
#'     [future.apply::future_lapply()] for the repeated derivative
#'     calculations. This requires the suggested package
#'     \pkg{future.apply}.
#' @param ... Additional arguments passed to `simulate`, `loglik`,
#'     and, when supplied, to `score` and `information`.
#'
#' @return A named list suitable for the `components` argument of
#'     [focus_engine()], with elements `V`, `P`, and `Q`.
#'
#' @details
#' The supplied `simulate` function is called repeatedly as
#' `simulate(theta, ...)`, once for each simulated dataset.
#'
#' The supplied `loglik` function is then used on each simulated dataset
#' through calls of the form `loglik(theta, data = data, ...)`.
#'
#' If both `score` and `information` are `NULL`, then both quantities
#' are obtained numerically from `loglik`.
#'
#' If `score` is supplied but `information` is `NULL`, then the
#' simulated score contributions are computed using `score`, and the
#' simulated information matrices are computed as the negative
#' Jacobian of `score`.
#'
#' If `information` is supplied, it is used directly, regardless of
#' whether `score` is also supplied.
#'
#' The returned `V` is obtained by averaging the observed information
#' matrices across the simulated datasets and inverting the result.
#'
#' The returned `P` and `Q` are estimated by averaging the
#' corresponding simulation-based quantities across the simulated
#' datasets.
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
    no_score <- is.null(score)
    no_info <- is.null(information)
    if (no_score && no_info) {
        simu_one <- function(i) {
            data <- simulate(theta, ...)
            ders <- grad_hess(loglik, theta, data = data, ...)
            list(S = ders$grad,
                 I = -ders$hess)
        }
    } else {
        if (no_score) {
            s_fun <- function(x, data, ...) numDeriv::grad(loglik, x, data = data, ...)
        } else {
            s_fun <- score
        }
        if (no_info) {
            if (no_score) {
                i_fun <- function(x, data, ...) -numDeriv::hessian(loglik, x, data = data, ...)
            } else {
                i_fun <- function(x, data, ...) -numDeriv::jacobian(score, x, data = data, ...)
            }
        } else {
            i_fun <- information
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
