#' Focus on a scalar function of the model parameters
#'
#' Estimate and perform inference on a scalar function of the model
#' parameters from a fitted model object of class [`"glm"`][stats::glm]
#' or [`"brglmFit"`][brglm2::brglmFit].
#'
#' The function evaluates a user-specified scalar function of the full
#' parameter vector and optionally applies mean or median bias correction.
#'
#' @param object A fitted model object of class [`"glm"`][stats::glm] or
#'   [`"brglmFit"`][brglm2::brglmFit].
#' @param on A function specifying the scalar parameter of interest.
#'   It must take the model parameter vector as first argument and return
#'   a single numeric value. The default,
#'   `function(theta, ...) theta[1]`, focuses on the first component of
#'   the parameter vector.
#' @param correction Character string specifying the bias correction method.
#'   One of:
#'   \describe{
#'     \item{`"no"`}{No correction.}
#'     \item{`"median"`}{Median bias correction; default.}
#'     \item{`"mean"`}{Mean bias correction.}
#'   }
#' @param on_gradient Optional function returning the gradient of `on` with
#'   respect to the model parameter vector. It must take the parameter vector
#'   as first argument and return a numeric vector of the same length.
#'   If `NULL` (default), the gradient is computed numerically using
#'   [numDeriv::grad()]. Apart from basic checks on type and dimension,
#'   `focus()` does not verify that the supplied gradient is correct.
#' @param on_hessian Optional function returning the Hessian matrix of `on`
#'   with respect to the model parameter vector. It must take the parameter
#'   vector as first argument and return a numeric matrix with one row and one
#'   column per model parameter. If `NULL` (default), the Hessian is computed
#'   numerically using [numDeriv::hessian()]. Apart from basic checks on type
#'   and dimension, `focus()` does not verify that the supplied Hessian is
#'   correct.
#' @param ... Additional arguments passed to `on`, `on_gradient`, and
#'   `on_hessian`.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{`estimate`}{Numeric scalar, the estimate of the quantity defined by `on`.
#'     The returned value has a `"correction"` attribute recording the bias
#'     correction method used.}
#'   \item{`se`}{Numeric scalar, the delta-method standard error.}
#'   \item{`object`}{The fitted model object used internally by `focus()`,
#'     after any refitting described below.}
#'   \item{`on`}{A list containing the supplied `on`, `on_gradient`, and
#'     `on_hessian` functions.}
#'   \item{`dots`}{A list with the additional arguments supplied through `...`.}
#' }
#'
#' @details
#'
#' If the primary class of `object` is [`"glm"`][stats::glm], then
#' `focus()` refits the model using [brglm2::brglmFit()] with `type = "ML"`
#' and starting values `coef(object)`.
#'
#' If the primary class of `object` is [`"brglmFit"`][brglm2::brglmFit] and
#' `object$type` is not one of `"ML"`, `"AS_mean"`, or `"correction"`,
#' then `focus()` refits the model with `type = "AS_mean"`. This ensures that
#' the first-order term in the bias expansion of the maximum likelihood
#' estimator is removed, and helps avoid separation issues in categorical
#' response models.
#'
#' The current implementation assumes that `object` supports
#' [stats::coef()] and [stats::vcov()] for the full parameter vector
#' and the inverse of the expected information.
#'
#' Let \eqn{\psi(\theta)} denote the scalar function specified by `on`.
#' The plug-in estimator is `on(theta, ...)`, where `theta` is the estimated
#' parameter vector. Mean and median bias corrections are computed using
#' first- and second-order derivatives of \eqn{\psi(\theta)} together with
#' model-specific auxiliary quantities.
#'
#' If `on_gradient` or `on_hessian` are supplied, then `focus()` uses them in
#' place of numerical derivatives. Apart from basic checks on type and
#' dimension, their validity is not verified and is the user's
#' responsibility.
#'
#' Arguments in `...` are reserved for `on`, `on_gradient`, and
#' `on_hessian`. However, names that coincide with formal arguments of
#' `focus()` itself, such as `correction` and `object`, are matched
#' before `...` is formed and therefore cannot be passed through `...`.
#'
#' Standard errors are computed using the delta method, with
#' covariance matrix and gradients evaluated at the estimated
#' parameters from `object` or the refit version of it, as described
#' above. They are therefore not re-evaluated at the corrected
#' estimates.
#'
#' Confidence intervals can be obtained from the returned object using
#' [confint()].
#'
#' @seealso [confint.focus_list()], [hulc_ci()]
#'
#' @examples
#'
#' library("brglm2")
#'
#' data("endometrial", package = "brglm2")
#' endo <- glm(HG ~ NV + PI + EH,
#'             data = endometrial,
#'             family = binomial("logit"),
#'             method = "brglmFit")
#'
#' ## Focus on the first regression parameter
#' focus(endo)
#'
#' ## Focus on the second coefficient
#' focus(endo, on = function(theta) theta[2])
#'
#' ## These are exactly the same as the estimates from one-step
#' ## quasi-Fisher scoring with type = "AS_median" starting from the
#' ## reduced-bias estimates
#' endo1FS <- update(endo, maxit = 1, max_step_factor = 1, type = "AS_median",
#'                   start = coef(endo))
#' coef(endo1FS)[1:2]
#'
#' ## Focus on an odds ratio
#' focus(endo, on = function(theta) exp(theta["NV"]))
#'
#' ## Focus on a contrast
#' focus_fun <- function(theta, i = 1, j = 2) theta[i] - theta[j]
#' focus(endo, on = focus_fun, i = 2, j = 3)
#'
#' ## The same contrast with analytic derivatives
#' focus_grad <- function(theta, i = 1, j = 2) {
#'   out <- rep(0, length(theta))
#'   out[i] <- 1
#'   out[j] <- -1
#'   out
#' }
#' focus_hessian <- function(theta, i = 1, j = 2) {
#'   matrix(0, nrow = length(theta), ncol = length(theta))
#' }
#' focus(endo, on = focus_fun, on_gradient = focus_grad, on_hessian = focus_hessian,
#'       i = 2, j = 3)
#'
#' ## Mean bias correction
#' focus(endo, on = focus_fun, i = 2, j = 3, correction = "mean")
#'
#' ## The mean-bias reduced estimate is exactly the same as it is
#' ## equivariant to linear transformations
#' coef(endo)[2] - coef(endo)[3]
#'
#' ## Wald confidence interval
#' confint(focus(endo), level = 0.9)
#'
#' \dontrun{
#' ## HulC confidence interval
#' set.seed(678)
#' confint(focus(endo), method = "hulc", level = 0.9,
#'         Delta = 0, check_statistic = FALSE)
#' }
#'
#' @export
focus <- function(object, ...) {
    UseMethod("focus")
}

#' @export
focus.default <- function(object, ...) {
    stop("No `focus()` method is available for objects of class ",
         paste(class(object), collapse = "/"), ".")
}

#' @export
print.focus_list <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
    if (!is.null(x$call)) {
        cat("Call:\n")
        print(x$call)
        cat("\n")
    }
    out <- cbind(
        unname(x$estimate),
        unname(x$se)
    )
    rownames(out) <- "on"
    colnames(out) <- c("Estimate", "Std. Error")
    cat("Focus estimate\n")
    printCoefmat(out, digits = digits, P.values = FALSE, has.Pvalue = FALSE)
    cat("\n")
    cat("Type of correction:", attr(x$estimate, "correction"), "\n")
    invisible(x)
}


#' @export
focus.glm <- function(object,
                      on = function(theta, ...) theta[1],
                      correction = "median",
                      on_gradient = NULL,
                      on_hessian = NULL, ...) {
    cl <- match.call()
    dots <- list(...)
    stopifnot(is.null(on_gradient) || is.function(on_gradient))
    stopifnot(is.null(on_hessian) || is.function(on_hessian))
    is_glm <- identical(class(object)[1], "glm")
    is_brglmFit <- inherits(object, "brglmFit")
    if (!(is_glm || is_brglmFit)) {
        stop("`focus.glm()` supports only objects of primary class `glm` or `brglmFit`.")
    }
    if (is_glm) {
        object <- update(object, method = "brglmFit", type = "ML", start = coef(object))
    }
    if (is_brglmFit) {
        if (!(object$type %in% c("ML", "AS_mean", "correction"))) {
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
    d1_psi <- if (is.null(on_gradient)) numDeriv::grad(on, theta, ...) else on_gradient(theta, ...)
    stopifnot(is.numeric(d1_psi), length(d1_psi) == length(theta), !anyNA(d1_psi))
    d1_psi <- as.numeric(d1_psi)
    hat_psi <- on(theta, ...)
    muffin <-  drop(V %*% d1_psi)
    var_psi <- sum(d1_psi * muffin)
    if (identical(correction, "no")) {
        out <- hat_psi
    } else {
        d2_psi <- if (is.null(on_hessian)) {
            numDeriv::hessian(on, theta, ...)
        } else {
            on_hessian(theta, ...)
        }
        stopifnot(is.numeric(d2_psi), identical(dim(d2_psi), c(length(theta), length(theta))), !anyNA(d2_psi))
        constant_on <-
            all(d1_psi == 0) &&
            all(d2_psi == 0)
        if (constant_on) {
            out <- hat_psi
        } else {
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
    }
    ## If the model has only one parameter (and some other cases) it
    ## is posible to evaluate the SE at the corrected estimates but
    ## not more generally. So, we always use V(theta) and d1_psi with
    ## theta at the original estimates (i.e. ML or
    ## AS_mean/correction).
    se <- sqrt(var_psi)
    out <- unname(out)
    attr(out, "correction") <- correction
    out <- list(
        call = cl,
        object = object,
        on = list(on = on,
                  on_gradient = on_gradient,
                  on_hessian = on_hessian),
        dots = dots,
        estimate = out,
        se = se)
    class(out) <- c("focus_list", class(out))
    out
}


#' Focus on a scalar function of the model parameters for new data
#'
#' Evaluate a scalar function of the model parameters after refitting
#' a model object on a supplied dataset.
#'
#' This is a convenience wrapper around [focus()], intended for use in
#' resampling procedures such as the bootstrap, where the model is
#' repeatedly refit on different datasets.
#'
#' @param data A data frame containing the variables required to refit
#'   the model.
#' @param object A fitted model object of class [`"glm"`][stats::glm] or
#'   [`"brglmFit"`][brglm2::brglmFit].
#' @param on A function specifying the scalar parameter of interest.
#'   It must take the model parameter vector as first argument and return
#'   a single numeric value. The default,
#'   `function(theta, ...) theta[1]`, focuses on the first component of
#'   the parameter vector.
#' @param correction Character string specifying the bias correction method.
#'   One of `"no"`, `"median"`, or `"mean"`.
#' @param on_gradient Optional function returning the gradient of `on`; see
#'   [focus()].
#' @param on_hessian Optional function returning the Hessian of `on`; see
#'   [focus()].
#' @param ... Additional arguments passed to [focus()]. These can be used
#'   to supply further arguments to `on`, `on_gradient`, and
#'   `on_hessian`.
#'
#' @return
#' A numeric scalar: the estimate of the quantity defined by `on`.
#'
#' @details
#' The function refits `object` using `data` via [update()], and then
#' applies [focus()] to the refitted model. Only the `estimate`
#' component of the result from [focus()] is returned.
#'
#' This function is particularly useful in resampling settings, where a
#' statistic function returning a single numeric value is required.
#'
#' @examples
#' library("brglm2")
#'
#' data("endometrial", package = "brglm2")
#'
#' endo <- glm(HG ~ NV + PI + EH,
#'             data = endometrial,
#'             family = binomial("logit"),
#'             method = "brglmFit")
#'
#' ## Focus on the first coefficient using the original data
#' focus_statistic(endometrial, endo, correction = "no")
#'
#' ## This is the same point estimate as
#' focus(endo, correction = "no")$estimate
#'
#' \dontrun{
#' if (requireNamespace("boot", quietly = TRUE)) {
#'   boot_fun <- function(data, indices) {
#'     d <- data[indices, ]
#'     focus_statistic(d, endo, correction = "mean")
#'   }
#'   boot::boot(endometrial, boot_fun, R = 100)
#' }
#' }
#'
#' @seealso [focus()]
#'
#' @export
focus_statistic <- function(data, object,
                            on = function(theta, ...) theta[1],
                            correction = "median",
                            on_gradient = NULL,
                            on_hessian = NULL, ...) {
    object <- do.call(update, list(object = object, data = data))
    focus(object, on = on, correction = correction,
          on_gradient = on_gradient, on_hessian = on_hessian, ...)$estimate
}
