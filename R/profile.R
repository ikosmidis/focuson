#' Profile likelihood confidence intervals
#'
#' Compute a profile likelihood confidence interval for a scalar function of a
#' model parameter vector using the endpoint equations of Venzon and
#' Moolgavkar (1988).
#'
#' @param loglik Function returning the log-likelihood at a supplied
#'     parameter vector. It must take the parameter vector as its
#'     first argument.
#' @param score Optional function returning the score vector at a
#'     supplied parameter vector. It must take the parameter vector as
#'     its first argument.  If `NULL`, the score is computed
#'     numerically from `loglik`.
#' @param information Optional function returning the observed or
#'     expected information matrix at a supplied parameter vector. It
#'     must take the parameter vector as its first argument. If
#'     `NULL`, the information is computed numerically as minus the
#'     Hessian of `loglik`.
#' @param mle Numeric vector. The unrestricted maximum likelihood
#'     estimate of the full model parameter vector.
#' @param on Function specifying the scalar parameter of interest. It
#'     must take the parameter vector as its first argument and return
#'     a numeric scalar.
#' @param on_gradient Optional function returning the gradient of `on`
#'     with respect to the parameter vector. It must take the
#'     parameter vector as its first argument. If `NULL`, the gradient
#'     is computed numerically.
#' @param likelihood_args List of additional arguments passed to
#'     `loglik`, `score`, and `information`.
#' @param nleqslv_args List of additional arguments passed to
#'     [nleqslv::nleqslv()].
#' @param level Confidence level. Default is `0.95`.
#' @param ... Additional arguments passed to `on` and `on_gradient`.
#'
#' @return
#' A numeric vector of length 2 with names `"lower"` and `"upper"`. The result
#' has attributes:
#' \describe{
#'   \item{`"type"`}{The string `"profile"`.}
#'   \item{`"max|fvec|"`}{The maximum absolute endpoint equation residual for
#'     each endpoint.}
#'   \item{`"messages"`}{The convergence messages returned by
#'     [nleqslv::nleqslv()] for the lower and upper endpoint solves.}
#' }
#'
#' @details
#'
#' The endpoint equations are
#' \deqn{2\{\ell(\hat\theta) - \ell(\theta)\} - c = 0}
#' and
#' \deqn{\nabla \ell(\theta) - \lambda \nabla g(\theta) = 0,}
#' where \eqn{\ell(\theta)} is `loglik`, \eqn{\hat\theta} is `mle`,
#' \eqn{g} is the function supplied through `on`, and \eqn{c} is the
#' `level` quantile of a chi-squared distribution with one degree of
#' freedom.
#'
#' The function solves these equations twice, from Wald-type starting values
#' in opposite directions, and returns the corresponding values of `on`. This
#' is a low-level routine: endpoint convergence diagnostics are returned as
#' attributes, and callers can decide how strictly to enforce them.
#'
#' @references
#' Venzon D J, Moolgavkar S H (1988). A method for computing
#' profile-likelihood-based confidence intervals. *Journal of the Royal
#' Statistical Society: Series C (Applied Statistics)*, **37**, 87--94.
#'
#' @examples
#' y <- c(-1, 0, 1, 2, 3)
#' loglik <- function(theta, y) {
#'     sum(dnorm(y, mean = theta[1], sd = 1, log = TRUE))
#' }
#' score <- function(theta, y) {
#'     sum(y - theta[1])
#' }
#' information <- function(theta, y) {
#'     matrix(length(y), 1, 1)
#' }
#' profile_ci(loglik = loglik,
#'            score = score,
#'            information = information,
#'            mle = mean(y),
#'            likelihood_args = list(y = y))
#'
#' @export
profile_ci <- function(loglik,
                       score = NULL,
                       information = NULL,
                       mle,
                       on = function(theta) theta[1],
                       on_gradient = NULL,
                       likelihood_args = list(),
                       nleqslv_args = list(),
                       level = 0.95,
                       ...) {
    quant <- qchisq(level, 1)
    if (!is.function(loglik))
        stop("`loglik` must be a function.")
    if (!is.null(score) && !is.function(score))
        stop("`score` must be a function.")
    if (!is.null(information) && !is.function(information))
        stop("`information` must be a function.")
    if (!is.function(on))
        stop("`on` must be a function.")
    if (!is.null(on_gradient) && !is.function(on_gradient))
        stop("`on_gradient` must be a function.")
    if (!is.numeric(mle) || length(mle) == 0 || any(!is.finite(mle)))
        stop("`mle` must be a finite numeric vector.")
    if (!is.list(likelihood_args))
        stop("`likelihood_args` must be a list.")
    if (!is.list(nleqslv_args))
        stop("`nleqslv_args` must be a list.")
    ll <- function(theta) do.call(loglik, c(list(theta), likelihood_args))
    if (is.null(on_gradient))
        on_gradient  <- function(theta, ...) numDeriv::grad(on, theta, ...)
    if (is.null(score))
        sc <- function(theta) numDeriv::grad(ll, theta)
    else
        sc <- function(theta) do.call(score, c(list(theta), likelihood_args))
    info <- if (is.null(information)) {
        function(theta) -numDeriv::hessian(ll, theta)
    } else {
        function(theta) do.call(information, c(list(theta), likelihood_args))
    }
    maxloglik <- ll(mle)
    if (!is.numeric(maxloglik) || length(maxloglik) != 1 || !is.finite(maxloglik))
        stop("`loglik` must return a finite numeric scalar at `mle`.")
    score_mle <- sc(mle)
    if (!is.numeric(score_mle) || length(score_mle) != length(mle) || any(!is.finite(score_mle)))
        stop("`score` must return a finite numeric vector of length `length(mle)` at `mle`.")
    on_grad <- on_gradient(mle, ...)
    if (!is.numeric(on_grad) || length(on_grad) != length(mle) || any(!is.finite(on_grad)))
        stop("`on_gradient` must return a finite numeric vector of length `length(mle)` at `mle`.")
    if (max(abs(score_mle)) > 1e-02)
        warning("`mle` is not the maximizer of the log-likelihood;",
                "the maximum absolute value of the gradient of the log-likelihood at `mle` is ",
                round(max(abs(score_mle)), 2), ".")
    npars <- length(mle) + 1
    eq <- function(pars) {
        theta <- pars[1:(npars - 1)]
        l <- ll(theta)
        s <- sc(theta)
        c(2 * (maxloglik - l) - quant, s - pars[npars] * on_gradient(theta, ...))
    }
    info_mle <- info(mle)
    if (!is.numeric(info_mle) || !identical(dim(info_mle), c(length(mle), length(mle))) ||
        any(!is.finite(info_mle)))
        stop("`information` must return a finite numeric matrix with one row and one column per parameter.")
    iinfo <- try(solve(info_mle), silent = TRUE)
    if (inherits(iinfo, "try-error"))
        stop("Could not invert the information matrix at `mle`.")
    step <- drop(iinfo %*% on_grad)
    var_on <- sum(on_grad * step)
    if (!is.finite(var_on) || var_on <= 0)
        stop("Could not construct profile endpoint starting values; ",
             "the delta-method variance of `on` at `mle` is non-positive or non-finite.")
    lam <- sqrt(quant / var_on)
    trans <- lam * step
    endpoints <- lapply(c(-1, 1), function(s) {
        do.call(nleqslv, c(list(x = c(mle + s * trans, - s * lam), fn = eq), nleqslv_args))
    })
    ci <- sapply(endpoints, function(end) on(end$x[1:(npars - 1)], ...))
    names(ci) <- c("lower", "upper")
    attr(ci, "max|fvec|") <- c(lower = max(abs(endpoints[[1]]$fvec)),
                               upper = max(abs(endpoints[[2]]$fvec)))
    attr(ci, "messages") <- c(lower = endpoints[[1]]$message,
                              upper = endpoints[[2]]$message)
    attr(ci, "type") <- "profile"
    ci
}
