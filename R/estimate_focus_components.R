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

    diagnostics_matrix <- function(draws, mean_mat) {
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
            p_diag[[t]] <<- diagnostics_matrix(draws_t, mean_t)
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
            q_diag[[t]] <<- diagnostics_matrix(draws_t, mean_t)
        }
        mean_t
    })

    if (diagnostics) {
        I_diag <- diagnostics_matrix(lapply(derivatives, function(der) der$SS), Ihat)
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
