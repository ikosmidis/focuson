#' Focus on a scalar function of the model parameters
#'
#' Estimate and perform inference on a scalar function of the model
#' parameters from a fitted model object of class [`"glm"`][stats::glm]
#' or [`"brglmFit"`][brglm2::brglmFit].
#'
#' The function evaluates a user-specified scalar function of the full
#' parameter vector and optionally applies mean or median bias correction.
#' Confidence intervals can be constructed using either Wald or HULC methods.
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
#'     \item{`"no"`}{No bias correction.}
#'     \item{`"median"`}{Median bias correction.}
#'     \item{`"mean"`}{Mean bias correction.}
#'   }
#' @param ci Character string specifying the confidence interval method.
#'   One of:
#'   \describe{
#'     \item{`"wald"`}{Wald-type confidence interval based on the delta method.}
#'     \item{`"hulc"`}{HulC confidence interval; see [hulc_ci()].}
#'   }
#' @param alpha Numeric scalar in `(0, 1)`. Target miscoverage level for the
#'   confidence interval.
#' @param control_ci A list of control options for confidence interval
#'   construction as returned by [ci_control()].
#' @param ... Additional arguments passed to `on`.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{`estimate`}{Numeric scalar, the estimate of the quantity defined by `on`.
#'     The returned value has a `"correction"` attribute recording the bias
#'     correction method used.}
#'   \item{`se`}{Numeric scalar, the delta-method standard error.}
#'   \item{`ci_type`}{Character string indicating the confidence interval method.}
#'   \item{`confint`}{Numeric vector of length 2 with names `"lower"` and `"upper"` if
#'     `ci = "wald"`, or the result of [hulc_ci()] if `ci = "hulc"`.}
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
#' Let \eqn{\psi(\theta)} denote the scalar function specified by `on`.
#' The plug-in estimator is `on(theta, ...)`, where `theta` is the fitted
#' parameter vector. Mean and median bias corrections are computed using
#' first- and second-order derivatives of \eqn{\psi(\theta)} together with
#' model-specific auxiliary quantities.
#'
#' Standard errors are computed using the delta method,
#' with covariance matrix evaluated at the original fitted parameter vector.
#' They are therefore not generally re-evaluated at the bias-corrected estimate.
#'
#' When `ci = "hulc"`, the HULC interval is computed by applying [hulc_ci()]
#' to the model frame using `focus_statistic()` as the statistic evaluated on
#' each partition.
#'
#' @seealso [ci_control()], [hulc_ci()]
#'
#' @examples
#'
#' library("brglm2")
#'
#' data("endometrial", package = "brglm2")
#' endo <- glm(HG ~ NV + PI + EH,
#'             data = endometrial,
#'             family = binomial("logit"),
#'             method = "brglm_fit")
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
#' focus_fun <- function(theta, a = 1, b = 2) theta[a] - theta[b]
#' focus(endo, on = focus_fun, a = 2, b = 3)
#'
#' ## Mean bias correction
#' focus(endo, on = focus_fun, a = 2, b = 3, correction = "mean")
#'
#' ## The mean-bias reduced estimate is exactly the same as it is #
#' ## equivariant to linear transformations
#' coef(endo)[2] - coef(endo)[3]
#'
#' \dontrun{
#' ## HulC confidence interval
#' focus(endo, ci = "hulc", control_ci = ci_control(Delta = 0.05))
#' }
#'
#' @export
focus <- function(object,
                  on = function(theta, ...) theta[1],
                  correction = "median",
                  ci = "wald",
                  alpha = 0.05,
                  control_ci = ci_control(), ...) {
    if (!inherits(control_ci, "ci_control")) {
        stop("Please use `ci_control()` to construct the list for `control_ci`")
    }
    if (identical(class(object)[1], "glm")) {
        object <- update(object, method = "brglmFit", type = "ML", start = coef(object))
    }
    if (inherits(object, "brglmFit")) {
        if (!(object$type %in% c("ML", "AS_mean", "correction"))) {
            object <- update(object, type = "AS_mean")
        }
    }
    correction <- match.arg(correction, c("no", "median", "mean"))
    ci <- match.arg(ci, c("wald", "hulc"))
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
        attr(confint, "alpha") <- alpha
    }
    if (identical(ci, "hulc")) {
        confint <- hulc_ci(model.frame(object),
                           statistic = focus_statistic,
                           object = object,
                           on = on,
                           correction = correction,
                           alpha = alpha,
                           Delta = control_ci$Delta,
                           randomize = control_ci$randomize,
                           check_statistic = control_ci$check_statistic, ...)
    }
    out <- unname(out)
    attr(out, "correction") <- correction
    list(estimate = out, se = se, ci_type = ci, confint = confint)
}

#' Control options for confidence intervals returned by [focus()]
#'
#' Create a control object specifying options for confidence intervals
#' returned by [focus()]
#'
#' @param Delta Numeric scalar in `[0, 1/2)`. Median bias of the
#'     estimator of the `on` function for when `ci = "hulc"` in
#'     [focus()]; see [hulc_ci()] for details.
#'
#' @param randomize Logical. If `TRUE`, use the randomization to
#'     determine the number of partitions when `ci = "hulc"` in
#'     [focus()]; see [hulc_ci()] for details.
#'
#' @param check_statistic Logical. If `TRUE` and `ci = "hulc"` in
#'     [focus()], check whether the estimator of the `on` function can
#'     be evaluated on the smallest partition before computing the
#'     interval; see [hulc_ci()] for details.
#'
#' @param ... Currently unused. Included for future extensions.
#'
#' @return
#'
#' A named list with components `Delta`, `randomize`,
#' `check_statistic` with the values supplied.
#'
#' @details
#'
#' This function does not perform extensive validation. It is
#' primarily a convenience constructor to group options related to
#' confidence intervals.
#'
#' @seealso [hulc_ci()], [focus()]
#'
#' @examples
#' ci_control()
#'
#' ci_control(Delta = 0.05, randomize = FALSE)
#'
#'
#' @export
ci_control <- function(Delta = 0, randomize = TRUE, check_statistic = TRUE, ...) {
    out <- list(Delta = Delta, randomize = randomize, check_statistic = check_statistic)
    class(out) <- c("ci_control", class(out))
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
#' @param ... Additional arguments passed to `on`.
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
#' focus_statistic(endometrial, endo)
#'
#' ## This is the same point estimate as
#' focus(endo)$estimate
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
                            correction = "median", ...) {
    object <- do.call(update, list(object = object, data = data))
    focus(object, on = on, correction = correction, ...)$estimate
}
