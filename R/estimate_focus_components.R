.diagnostics_matrix <- function(draws, mean_mat, nsim = length(draws)) {
    second_moment <- Reduce("+", lapply(draws, function(x) x^2)) / nsim
    var_mat <- if (nsim > 1L) {
        pmax((second_moment - mean_mat^2) * nsim / (nsim - 1), 0)
    } else {
        mean_mat * 0
    }
    mcse_mat <- sqrt(var_mat / nsim)
    frob_mean <- sqrt(sum(mean_mat^2))
    list(
        mcse_max = max(mcse_mat),
        mcse_frobenius = sqrt(sum(mcse_mat^2)),
        rel_mcse_frobenius = sqrt(sum(mcse_mat^2)) /
            max(frob_mean, .Machine$double.eps)
    )
}

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
#' @param diagnostics Logical. If `TRUE`, return simple Monte Carlo
#'     diagnostics for the estimated components.
#' @param ... Additional arguments passed to `simulate`, `loglik`,
#'     and, when supplied, to `score` and `information`.
#'
#' @return An object of class `"focus_components"`, suitable for the
#'     `components` argument of [focus_engine()], with elements `V`,
#'     `P`, and `Q`. If `diagnostics = TRUE`, an additional element
#'     `diagnostics` is included.
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
#' The returned `V` is obtained by averaging the outer products of the
#' simulated score vectors across the simulated datasets and inverting
#' the result.
#'
#' The returned `P` and `Q` are estimated by averaging the
#' corresponding simulation-based quantities across the simulated
#' datasets.
#'
#' If `diagnostics = TRUE`, simple Monte Carlo diagnostics are
#' returned for `I`, `P`, and `Q`, based on sample variability of the
#' simulated contributions.
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
                                      diagnostics = FALSE,
                                      ...) {
    cl <- match.call()
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
                 I = -ders$hess,
                 SS = tcrossprod(ders$grad))
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
            S <- s_fun(theta, data, ...)
            I <- i_fun(theta, data, ...)
            SS <- tcrossprod(S)
            list(S = S, I = I, SS = SS)
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

    Ihat <- Reduce("+", lapply(derivatives, function(der) der$SS)) / nsim
    out <- list(V = solve(Ihat))
    sc <- mean(diag(out$V))
    if (!is.finite(sc) || sc < 1e-6 || sc > 1) sc <- 1
    p_diag <- vector("list", length(theta))
    out$P <- lapply(seq_along(theta), function(t) {
        draws_t <- lapply(derivatives, function(der) der$SS * der$S[t] * sc)
        mean_t <- Reduce(
            "+",
            draws_t
        ) / (sc * nsim)
        if (diagnostics) {
            p_diag[[t]] <<- .diagnostics_matrix(draws_t, mean_t, nsim = nsim)
        }
        mean_t
    })

    q_diag <- vector("list", length(theta))
    out$Q <- lapply(seq_along(theta), function(t) {
        draws_t <- lapply(derivatives, function(der) -der$I * der$S[t] * sc)
        mean_t <- Reduce(
            "+",
            draws_t
        ) / (sc * nsim)
        if (diagnostics) {
            q_diag[[t]] <<- .diagnostics_matrix(draws_t, mean_t, nsim = nsim)
        }
        mean_t
    })

    if (diagnostics) {
        I_diag <- .diagnostics_matrix(lapply(derivatives, function(der) der$SS), Ihat, nsim = nsim)
        I_diag$kappa <- kappa(Ihat)
        I_diag$min_eigen <- min(Re(eigen((Ihat + t(Ihat)) / 2,
                                         only.values = TRUE)$values))
        out$diagnostics <- list(I = I_diag,
                                P = p_diag,
                                Q = q_diag)
    }

    out$call <- cl
    out$meta <- list(score_supplied = !no_score,
                     information_supplied = !no_info,
                     parallelize = parallelize,
                     diagnostics = diagnostics,
                     nsim = nsim)
    class(out) <- c("focus_components", class(out))
    out
}

#' Estimate components for [focus_engine()] in a full exponential family
#'
#' Estimate the model-side components used by [focus_engine()] in the
#' special case where the reference parameterization is that of a full
#' exponential family, so that `Q = 0` and `V` can be obtained directly
#' from the observed information evaluated at the supplied estimate.
#'
#' @param theta Numeric parameter vector at which the components are
#'     estimated.
#' @param data Observed dataset at which the observed information is
#'     evaluated.
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
#'     estimate `P`.
#' @param parallelize Logical. If `TRUE`, use
#'     [future.apply::future_lapply()] for the repeated score
#'     calculations. This requires the suggested package
#'     \pkg{future.apply}.
#' @param diagnostics Logical. If `TRUE`, return simple Monte Carlo
#'     diagnostics for the estimated components.
#' @param ... Additional arguments passed to `simulate`, `loglik`,
#'     and, when supplied, to `score` and `information`.
#'
#' @return An object of class `"focus_components"`, suitable for the
#'     `components` argument of [focus_engine()], with elements `V`,
#'     `P`, and `Q`. The returned `Q` consists of zero matrices. If
#'     `diagnostics = TRUE`, an additional element `diagnostics` is
#'     included.
#'
#' @details
#' This helper assumes that the supplied parameterization is that of a
#' full exponential family. In that case, `Q` is taken to be zero and
#' `V` is obtained by inverting the observed information evaluated at
#' the supplied `theta` and `data`. Only `P` is estimated by Monte
#' Carlo simulation.
#'
#' @seealso [focus_engine()], [estimate_focus_components()]
#'
#' @export
estimate_focus_components_fef <- function(theta,
                                          data,
                                          loglik,
                                          score = NULL,
                                          information = NULL,
                                          simulate,
                                          nsim = 1000,
                                          parallelize = FALSE,
                                          diagnostics = FALSE,
                                          ...) {
    cl <- match.call()
    theta <- as.numeric(theta)
    if (parallelize && !requireNamespace("future.apply", quietly = TRUE)) {
        stop("Package `future.apply` is required when `parallelize = TRUE`.")
    }
    no_score <- is.null(score)
    no_info <- is.null(information)

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

    Ihat <- i_fun(theta, data, ...)
    out <- list(V = solve(Ihat))
    sc <- mean(diag(out$V))
    if (!is.finite(sc) || sc < 1e-6 || sc > 1) sc <- 1

    simu_one <- function(i) {
        sim_data <- simulate(theta, ...)
        S <- s_fun(theta, sim_data, ...)
        SS <- tcrossprod(S)
        list(S = S, SS = SS)
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

    p_diag <- vector("list", length(theta))
    out$P <- lapply(seq_along(theta), function(t) {
        draws_t <- lapply(derivatives, function(der) der$SS * der$S[t] * sc)
        mean_t <- Reduce("+", draws_t) / (sc * nsim)
        if (diagnostics) {
            p_diag[[t]] <<- .diagnostics_matrix(draws_t, mean_t, nsim = nsim)
        }
        mean_t
    })

    zero_mat <- matrix(0, nrow = length(theta), ncol = length(theta))
    q_diag <- replicate(length(theta),
                        .diagnostics_matrix(list(zero_mat), zero_mat, nsim = 1L),
                        simplify = FALSE)
    out$Q <- replicate(length(theta), zero_mat, simplify = FALSE)

    if (diagnostics) {
        I_diag <- .diagnostics_matrix(list(Ihat), Ihat, nsim = 1L)
        I_diag$kappa <- kappa(Ihat)
        I_diag$min_eigen <- min(Re(eigen((Ihat + t(Ihat)) / 2,
                                         only.values = TRUE)$values))
        out$diagnostics <- list(I = I_diag,
                                P = p_diag,
                                Q = q_diag)
    }

    out$call <- cl
    out$meta <- list(score_supplied = !no_score,
                     information_supplied = !no_info,
                     parallelize = parallelize,
                     diagnostics = diagnostics,
                     nsim = nsim,
                     fef = TRUE)
    class(out) <- c("focus_components", class(out))
    out
}

#' Estimate components for [focus_engine()] under iid sampling
#'
#' Estimate the model-side components used by [focus_engine()] by
#' simulating single iid observations and exploiting the additive
#' structure of the score and information.
#'
#' @param theta Numeric parameter vector at which the components are
#'     estimated.
#' @param n Integer. The intended sample size of the iid dataset for
#'     which the components are required.
#' @param loglik A function returning the log-likelihood contribution
#'     of a single observation evaluated at a supplied parameter
#'     vector. It is expected to have an interface of the form
#'     `loglik(theta, data, ...)` and to return a single numeric
#'     value.
#' @param score Optional function returning the score contribution of a
#'     single observation. It is expected to have an interface of the
#'     form `score(theta, data, ...)`. If `NULL` (default), the score
#'     contribution is computed numerically from `loglik`.
#' @param information Optional function returning the observed
#'     information contribution of a single observation, namely the
#'     negative Hessian contribution of the log-likelihood. It is
#'     expected to have an interface of the form `information(theta,
#'     data, ...)`. If `NULL` (default), the information contribution
#'     is computed numerically from `loglik`, or, if `score` is
#'     supplied, as the negative Jacobian of `score`.
#' @param simulate A function that simulates one iid observation at the
#'     supplied parameter vector `theta`. It is expected to have an
#'     interface of the form `simulate(theta, ...)` and to return one
#'     simulated observation at each call.
#' @param nsim Integer. The number of simulated iid observations used
#'     to estimate the components.
#' @param parallelize Logical. If `TRUE`, use
#'     [future.apply::future_lapply()] for the repeated derivative
#'     calculations. This requires the suggested package
#'     \pkg{future.apply}.
#' @param diagnostics Logical. If `TRUE`, return simple Monte Carlo
#'     diagnostics for the estimated components.
#' @param ... Additional arguments passed to `simulate`, `loglik`, and,
#'     when supplied, to `score` and `information`.
#'
#' @return An object of class `"focus_components"`, suitable for the
#'     `components` argument of [focus_engine()], with elements `V`,
#'     `P`, and `Q`. If `diagnostics = TRUE`, an additional element
#'     `diagnostics` is included.
#'
#' @details
#' The supplied `simulate` function is called repeatedly as
#' `simulate(theta, ...)`, once for each simulated iid observation.
#'
#' The supplied `loglik` function is then used on each simulated
#' observation through calls of the form `loglik(theta, data = data,
#' ...)`.
#'
#' If both `score` and `information` are `NULL`, then both quantities
#' are obtained numerically from `loglik`.
#'
#' If `score` is supplied but `information` is `NULL`, then the
#' simulated score contributions are computed using `score`, and the
#' simulated information contributions are computed as the negative
#' Jacobian of `score`.
#'
#' If `information` is supplied, it is used directly, regardless of
#' whether `score` is also supplied.
#'
#' The returned `V`, `P`, and `Q` target the corresponding full-sample
#' components under iid sampling by multiplying the simulated
#' one-observation quantities by `n`.
#'
#' @seealso [focus_engine()], [estimate_focus_components()]
#'
#' @export
estimate_focus_components_iid <- function(theta,
                                          n,
                                          loglik,
                                          score = NULL,
                                          information = NULL,
                                          simulate,
                                          nsim = 1000,
                                          parallelize = FALSE,
                                          diagnostics = FALSE,
                                          ...) {
    cl <- match.call()
    theta <- as.numeric(theta)
    n <- as.integer(n)
    if (length(n) != 1L || is.na(n) || n < 1L) {
        stop("`n` must be a positive integer.")
    }
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
                 I = -ders$hess,
                 SS = tcrossprod(ders$grad))
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
            S <- s_fun(theta, data, ...)
            I <- i_fun(theta, data, ...)
            SS <- tcrossprod(S)
            list(S = S, I = I, SS = SS)
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

    Ihat <- n * Reduce("+", lapply(derivatives, function(der) der$SS)) / nsim
    out <- list(V = solve(Ihat))
    sc <- mean(diag(out$V))
    if (!is.finite(sc) || sc < 1e-6 || sc > 1) sc <- 1
    p_diag <- vector("list", length(theta))
    out$P <- lapply(seq_along(theta), function(t) {
        draws_t <- lapply(derivatives, function(der) n * der$SS * der$S[t] * sc)
        mean_t <- Reduce("+", draws_t) / (sc * nsim)
        if (diagnostics) {
            p_diag[[t]] <<- .diagnostics_matrix(draws_t, mean_t, nsim = nsim)
        }
        mean_t
    })

    q_diag <- vector("list", length(theta))
    out$Q <- lapply(seq_along(theta), function(t) {
        draws_t <- lapply(derivatives, function(der) -n * der$I * der$S[t] * sc)
        mean_t <- Reduce("+", draws_t) / (sc * nsim)
        if (diagnostics) {
            q_diag[[t]] <<- .diagnostics_matrix(draws_t, mean_t, nsim = nsim)
        }
        mean_t
    })

    if (diagnostics) {
        I_diag <- .diagnostics_matrix(lapply(derivatives, function(der) n * der$SS), Ihat, nsim = nsim)
        I_diag$kappa <- kappa(Ihat)
        I_diag$min_eigen <- min(Re(eigen((Ihat + t(Ihat)) / 2,
                                         only.values = TRUE)$values))
        out$diagnostics <- list(I = I_diag,
                                P = p_diag,
                                Q = q_diag)
    }

    out$call <- cl
    out$meta <- list(score_supplied = !no_score,
                     information_supplied = !no_info,
                     parallelize = parallelize,
                     diagnostics = diagnostics,
                     nsim = nsim,
                     iid = TRUE,
                     n = n)
    class(out) <- c("focus_components", class(out))
    out
}

#' @export
print.focus_components <- function(x,
                                   digits = max(3L, getOption("digits") - 2L),
                                   ...) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
    meta <- x$meta
    derivative_mode <-
        if (isTRUE(meta$score_supplied) && isTRUE(meta$information_supplied)) {
            "supplied score, supplied information"
        } else if (isTRUE(meta$score_supplied) && !isTRUE(meta$information_supplied)) {
            "supplied score, Jacobian-based information"
        } else if (!isTRUE(meta$score_supplied) && isTRUE(meta$information_supplied)) {
            "numerical score, supplied information"
        } else {
            "numerical score and information"
        }
    cat("Monte Carlo component estimates\n")
    cat("Parameters:", length(x$P), "\n")
    if (isTRUE(meta$iid)) {
        cat("Sampling:", "iid contributions", "\n")
        cat("Target sample size:", meta$n, "\n")
    }
    if (isTRUE(meta$fef)) {
        cat("Structure:", "full exponential family shortcut", "\n")
    }
    cat("Simulations:", meta$nsim, "\n")
    cat("Derivatives:", derivative_mode, "\n")
    cat("Diagnostics:", if (isTRUE(meta$diagnostics)) "yes" else "no", "\n")
    cat("Parallelized:", if (isTRUE(meta$parallelize)) "yes" else "no", "\n")
    if (isTRUE(meta$diagnostics) && !is.null(x$diagnostics)) {
        cat("\nDiagnostics\n")
        I_diag <- x$diagnostics$I
        cat("I: kappa = ", format(signif(I_diag$kappa, digits)),
            ", min eigen = ", format(signif(I_diag$min_eigen, digits)),
            ", MCSE(Frobenius) = ", format(signif(I_diag$mcse_frobenius, digits)),
            ", relative MCSE(Frobenius) = ", format(signif(I_diag$rel_mcse_frobenius, digits)),
            "\n", sep = "")
        p_rel <- vapply(x$diagnostics$P, `[[`, numeric(1), "rel_mcse_frobenius")
        q_rel <- vapply(x$diagnostics$Q, `[[`, numeric(1), "rel_mcse_frobenius")
        cat("P: max relative MCSE(Frobenius) = ",
            format(signif(max(p_rel), digits)), "\n")
        cat("Q: max relative MCSE(Frobenius) = ",
            format(signif(max(q_rel), digits)), "\n")
    }
    invisible(x)
}
