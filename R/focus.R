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
#' @param control A list, as returned by a call to
#'     [focus_control()]. See Details.
#' @param ... Additional arguments passed to `on`.
#'
#' @return
#' A numeric scalar: the estimate  of the quantity defined by `on`.
#'
#' @details
#'
#' If the primary class of the `object` is [`"glm"`][glm], then
#' `focus()` will refit using [brglm2::brglmFit()] with `type = "ML"` starting
#' at `coef(object)`.
#'
#' If the primary class of the `object` is [`"brglmFit"`][brglmFit],
#' and `object$type` is not one of "ML", "AS_mean", or "correction",
#' then `focus()` refits using `type = "AS_mean"`. This ensures that
#' the first-term in the bias expansion of the maximum likelihood
#' estimator is eliminated, and avoids potential separation issues for
#' categorical response models.
#'
#' @examples
#' \dontrun{
#' library("brglm2")
#'
#' data("endometrial", package = "brglm2")
#' endo <- glm(HG ~ NV + PI + EH,
#'            data = endometrial,
#'            family = binomial("logit"),
#'            method = "brglm_fit")
#'
#' ## Focus on the first regression parameter
#' focus(endo)
#'
#' ## Focus on the second coefficient
#' focus(endo, on = function(theta) theta[2])
#'
#' ## which are exactly the same as the estiamte from one step
#' ## quasi-Fisher scoring with `type = "AS_median"` starting from the
#' ## reduced-bias estimates
#' endo1FS <- update(endo, maxit = 1, max_step_factor = 1, type = "AS_median",
#'                 start = coef(endo))
#' coef(endo1FS)[1:2]
#'
#' ## Focus on an odds ratio
#' focus(endo, on = function(theta, ...) exp(theta["NV"]))
#'
#' ## Focus on a contrast
#' focus_fun <- function(theta, a = 1, b = 2) theta[a] - theta[b]
#' focus(endo, on = focus_fun, a = 2, b = 3)
#'
#' ## Or using mean bias correction
#' focus(endo, on = focus_fun, a = 2, b = 3, correction = "mean")
#'
#' ## The mean-bias reduced estimate is exactly the same as it is
#' ## equivariant to linear transformations
#' coef(endo)[2] - coef(endo)[3]
#'
#' }
#'
#' @export
focus <- function(object, on = function(theta, ...) theta[1], control = focus_control(), ...) {
    if (isTRUE(inherits(object, "glm", TRUE) == 1)) {
        object <- update(object, method = "brglmFit", type = "ML", start = coef(object))
    }
    if (inherits(object, "brglmFit")) {
        if (!(object$type %in% c("ML", "AS_mean", "correction"))) {
            object <- update(object, type = "AS_mean")
        }
    }
    correction <- control$correction
    ci <- control$ci
    alpha <- control$alpha
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
    se <- sqrt(var_psi)
    if (identical(ci, "wald")) {
        confint <- out + c(-1, 1) * qnorm(1 - alpha / 2) * se
        names(confint) <- c("lower", "upper")
    }
    list(estimate = unname(out), se = se, confint = confint, ci_type = ci)
}

#' Configure [focus()] inference options
#'
#' Create a control object specifying bias correction and confidence interval
#' construction options for [focus()] inference procedures.
#'
#' This function standardizes user inputs and returns a list of settings that
#' can be passed to downstream estimation or inference routines.
#'
#' @param correction Character string specifying the bias correction method.
#'   One of:
#'   \describe{
#'     \item{`"no"`}{`on` in [focus()] is simply evaluated at the estimates from the supplied object.}
#'     \item{`"mean"`}{explicit reduced-mean-bias estimation of the `on` parameter.}
#'     \item{`"median"`}{default; explicit reduced-median-bias estimation of the `on` parameter.}
#'   }
#' @param ci Character string specifying the confidence interval method.
#'   One of:
#'   \describe{
#'     \item{`"wald"`}{Wald-type confidence interval.}
#'     \item{`"hulc"`}{HulC confidence interval (see [hulc_ci()]).}
#'   }
#' @param alpha Numeric scalar in `(0, 1)`. Target miscoverage level.
#' @param Delta Numeric scalar in `[0, 1/2)`. Median bias used by the
#'   HulC method; see [hulc_ci()]. Ignored unless `ci = "hulc"`.
#' @param randomize Logical. Whether to use randomized HulC intervals; see [hulc_ci()].
#'   Ignored unless `ci = "hulc"`.
#' @param check_statistic Logical. Whether to validate the statistic on the
#'   smallest partition when using HulC; see [hulc_ci()].. Ignored unless `ci = "hulc"`.
#'
#' @return
#' A named list with components:
#' \describe{
#'   \item{`correction`}{Selected bias correction method.}
#'   \item{`ci`}{Selected confidence interval method.}
#'   \item{`alpha`}{Miscoverage level.}
#'   \item{`Delta`}{Median bias of the estimator of `on`.}
#'   \item{`randomize`}{Logical flag for randomized HulC.}
#'   \item{`check_statistic`}{Logical flag controlling statistic validation.}
#' }
#'
#' @details
#'
#' This function performs input validation and matches arguments using
#' `match.arg()`. The returned object is intended to be passed to the
#' `control` argument of the [focus()], and controls estimation and
#' inference behavior.
#'
#' Parameters `Delta`, `randomize`, and `check_statistic` are only
#' relevant when `ci = "hulc"`, but are included in the returned
#' object for convenience.
#'
#' @seealso [focus()] [hulc_ci()]
#'
#' @examples
#' focus_control()
#'
#' focus_control(
#'   correction = "mean",
#'   ci = "hulc",
#'   alpha = 0.1,
#'   Delta = 0.05
#' )
#'
#' @export
focus_control <- function(correction = "median", ci = "wald", alpha = 0.05,
                          Delta = 0, randomize = TRUE, check_statistic = TRUE) {
    stopifnot(length(alpha) == 1, is.numeric(alpha), alpha > 0, alpha < 1)
    stopifnot(length(Delta) == 1, is.numeric(Delta), Delta >= 0, Delta < 1/2)
    correction <- match.arg(correction, c("no", "median", "mean"))
    ci <- match.arg(ci, c("wald", "hulc"))
    list(correction = correction, ci = ci, alpha = alpha)
}

#' Focus on a scalar function of the model parameters for new data
#'
#' Evaluates a scalar function of the model parameters after refitting
#' a model object on a supplied dataset. This is a convenience wrapper
#' around [focus()], intended for use in resampling procedures (e.g.,
#' bootstrap) where the model needs to be refit on different datasets.
#'
#' @param data A data frame (or object coercible to a data frame)
#'     containing the variables required to refit the model.
#' @param object A fitted model object of class [`"glm"`][glm] or
#'     [`"brglmFit"`][brglmFit].
#' @param control A list, as returned by a call to
#'     [focus_control()]. See Details.
#' @param on A function specifying the scalar parameter of
#'     interest. It must take the model parameter vector as first
#'     argument and return a single numeric value. The default,
#'     `function(theta, ...) theta[1]`, focuses on the first component
#'     of the parameter vector.
#' @param control A list, as returned by a call to
#'     [focus_control()]. See Details.
#' @param ... Additional arguments passed to `on`.
#'
#' @return
#' A numeric scalar: the estimate of the quantity defined by `on`,
#' with an attribute `"se"` giving the corresponding standard error.
#'
#' @details
#' The function refits `object` using `data` via [update()], and then
#' applies [focus()] to the refitted model. This ensures that the
#' quantity of interest is evaluated using parameter estimates
#' corresponding to the supplied dataset.
#'
#' This function is particularly useful in resampling settings, such as
#' bootstrap or cross-validation, where `data` varies across iterations.
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
#' focus_statistic(endometrial, endo)
#'
#' ## which is the same as
#' focus(endo)
#'
#' \dontrun{
#'
#' ## Example in a bootstrap setting
#' if (requireNamespace("boot", quietly = TRUE)) {
#'   boot_fun <- function(data, indices) {
#'     d <- data[indices, ]
#'     focus_statistic(d, endo, correction = "mean")
#'   }
#'   boot::boot(endometrial, boot_fun, R = 100)
#' }
#'
#' }
#'
#' @seealso [focus()]
#'
#' @export
focus_statistic <- function(data, object, on = function(theta, ...) theta[1], control = focus_control(), ...) {
    object <- do.call(update, list(object = object, data = data))
    focus(object, on = on, control = focus_control(), ...)
}


