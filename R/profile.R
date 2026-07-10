#' df <- data.frame(ldose = rep(0:5, 2),
#'                  numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16),
#'                  sex = factor(rep(c("M", "F"), c(6, 6))))
#' df <- df |> transform(numalive = 20 - numdead)
#' budworm.lg <- glm(cbind(numalive, numdead) ~ sex*ldose, family = binomial, data = df)
#' pr1 <- profile(budworm.lg)
#' plot(pr1)
#' pairs(pr1)
#'
#' ll <- function(theta, object) {
#'     dfun <- get_dmodel_function(object)
#'     sum(dfun(coefficients = theta, log = TRUE))
#' }
#' sc <- function(theta, object) {
#'     sc <- get_score_function(object)
#'     sc(coefficients = theta)
#' }
#' info <- function(theta, object) {
#'     ii <- get_information_function(object)
#'     ii(coefficients = theta)
#' }
#' profile_ci(ll, score = sc,
#'            mle = coef(budworm.lg), on = function(theta) theta[2],
#'            likelihood_args = list(object = budworm.lg))
profile_ci <- function(loglik,
                       score = NULL,
                       mle,
                       on = function(theta) theta[1],
                       on_gradient = NULL,
                       likelihood_args = list(),
                       nleqslv_args = list(),
                       level = 0.95,
                       ...) {
    quant <- qchisq(level, 1)
    ll <- function(theta) do.call(loglik, c(list(theta), likelihood_args))
    if (is.null(on_gradient))
        on_gradient  <- function(theta, ...) numDeriv::grad(on, theta, ...)
    if (is.null(score))
        sc <- function(theta) numDeriv::grad(ll, theta)
    else
        sc <- function(theta) do.call(score, c(list(theta), likelihood_args))
    maxloglik <- ll(mle)
    score_mle <- sc(mle)
    on_grad <- on_gradient(mle, ...)
    if (max(abs(score_mle)) > 1e-02)
        warning("`mle` is not the maximizer of the log-likelihood; the maximum absolute value of the gradient of the log-likelihood at `mle` is ", round(max(abs(score_mle)), 2), ".")
    npars <- length(mle) + 1
    eq <- function(pars) {
        theta <- pars[1:(npars - 1)]
        l <- ll(theta)
        s <- sc(theta)
        c(2 * (maxloglik - l) - quant, s - pars[npars] * on_gradient(theta, ...))
    }
    iinfo <- solve(-numDeriv::hessian(ll, mle))
    step <- drop(iinfo %*% on_grad)
    lam <- sqrt(quant / sum(on_grad * step))
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
