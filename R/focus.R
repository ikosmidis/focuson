#' Focus on a scalar function of the model parameters
#'
#' Estimates a scalar function of the model parameters, from a fitted
#' object of class [`"glm"`][glm] or [`"brglmFit"`][brglmFit].
#'
#' The `focus()` method is intended for estimation of and inference on
#' a user-defined quantity of interest, represented by a scalar
#' function of the full parameter vector, and additional arguments.
#'
#' @param object A fitted model object of class [`"glm"`][glm] or
#'     [`"brglmFit"`][brglmFit].
#' @param on A function specifying the scalar parameter of
#'     interest. It must take the model parameter vector as first
#'     argument and return a single numeric value. The default,
#'     `function(theta, ...) theta[1]`, focuses on the first component
#'     of the parameter vector.
#' @param correction The type of correction to be used. Possible
#'     values are `"no"` (`on` is simply evaluated at the estimates
#'     from the supplied object), `"mean"` (explicit reduced-mean-bias
#'     estimation of the `on` parameter), and `"median"` (default; explicit
#'     reduced-median-bias estimation of the `on` parameter).
#' @param ... Additional arguments passed to `on`.
#'
#' @return
#' A numeric scalar: the estimate  of the quantity defined by `on`.
#'
#' @details
#' If object is glm then refit using brglmFit with type = "ML" starting at the glm coefs
#' If object is brglmFit and type is not ML, "AS_mean", or "correction" then refit using AS_mean (just to avoid potential separation issues and have bias = 0)
#'
#' @examples
#' \dontrun{
#' library(brglm2)
#'
#' data("endometrial", package = "brglm2")
#' fit <- glm(HG ~ NV + PI + EH,
#'            data = endometrial,
#'            family = binomial("logit"),
#'            method = "brglm_fit")
#'
#' ## Focus on the first regression parameter
#' focus(fit)
#'
#' ## Focus on the second coefficient
#' focus(fit, on = function(theta) theta[2])
#'
#' ## which are exactly the same as the estiamte from one step
#' ## quasi-Fisher scoring with `type = "AS_median"` starting from the
#' ## reduced-bias estimates
#' fit1FS <- update(fit, maxit = 1, max_step_factor = 1, type = "AS_median",
#'                 start = coef(fit))
#' coef(fit1FS)[1:2]
#'
#' ## Focus on an odds ratio
#' focus(fit, on = function(theta, ...) exp(theta["NV"]))
#'
#' ## Focus on a contrast
#' focus_fun <- function(theta, a = 1, b = 2) theta[a] - theta[b]
#' focus(fit, on = focus_fun, a = 2, b = 3)
#'
#' ## Or using mean bias correction
#' focus(fit, on = focus_fun, a = 2, b = 3, correction = "mean")
#'
#' ## The mean-bias reduced estimate is exactly the same as it is
#' ## equivariant to linear transformations
#' coef(fit)[2] - coef(fit)[3]
#'
#' }
#'
#' @export
focus <- function(object, on = function(theta, ...) theta[1], correction = "median", ...) {
    if (isTRUE(inherits(object, "glm", TRUE) == 1)) {
        object <- update(object, method = "brglmFit", type = "ML", start = coef(object))
    }
    if (inherits(object, "brglmFit")) {
        otype <- object$type
        if (!(otype %in% c("ML", "AS_mean", "correction"))) {
            object <- update(object, type = "AS_mean")
        }
    }
    correction <- match.arg(correction, c("no", "median", "mean"))
    V <- vcov(object, model = "full")
    theta <- coef(object, model = "full")
    if (object$family$family %in% c("poisson", "binomial")) {
        cnams <- names(coef(object, model = "mean"))
        V <- V[cnams, cnams]
        theta <- theta[cnams]
    }
    d1_psi <- numDeriv::grad(on, theta, ...)
    hat_psi <- on(theta, ...)
    muffin <-  drop(V %*% d1_psi)
    var_psi <- sum(d1_psi * muffin)
    if (identical(correction, "no")) {
        out <- hat_psi
    } else {
        d2_psi <- numDeriv::hessian(on, theta, ...)
        afuns <- enrichwith::get_auxiliary_functions(object)
        P <- afuns$Pmat()
        Q <- afuns$Qmat()
        bias <- if (object$type %in% c("ML")) afuns$bias() else 0
        mean_b <- sum(d1_psi * bias) + 0.5 * sum(d2_psi * V) # 1st term of the mean bias expansion
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
    ## If the model has only one parameter (and some other cases) it
    ## is posible to evaluate the SE at the corrected estimates but
    ## not more generally. So, we always use V(theta) and d1_psi with
    ## theta at the original estimates (i.e. ML or
    ## AS_mean/correction).
    attr(out, "se") <- sqrt(var_psi)
    out
}
