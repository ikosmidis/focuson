#' Plot a focus profile log-likelihood
#'
#' Plot a focus profile log-likelihood either on the log-likelihood scale or the
#' signed likelihood-root scale.
#'
#' @param x An object of class `"profile_focus_list"`, as returned by
#'     [profile_focus()] or [profile.focus_list_glm()].
#' @param level Confidence level used to draw the horizontal cutoff
#'     and, if `ci = TRUE`, vertical confidence limits.
#' @param signed Logical. If `TRUE`, plot the signed likelihood root;
#'     otherwise plot the profile log-likelihood. Default is `FALSE`.
#' @param what Character. The quantity to plot. `"pl"` (default) selects the
#'     ordinary profile likelihood, `"mpl"` selects the modified profile
#'     likelihood, and `"rstar"` selects the modified signed
#'     likelihood-ratio statistic. The latter two require an object returned
#'     by [modified_profile_focus()]. Selecting `"rstar"` sets `signed = TRUE`
#'     internally.
#' @param interpolation Character. Interpolation method used between
#'     computed profile points. `"linear"` (default) uses [stats::approxfun()]
#'     and `"cubic"` uses [stats::splinefun()].
#' @param ci Logical. If `TRUE`, display the confidence limits for
#'     `level`. Default is `FALSE`. The corresponding cutoff is
#'     displayed irrespective of `ci`.
#' @param ... Additional graphical arguments passed to
#'     [graphics::plot.default()].
#'
#' @details
#'
#' When `ci = TRUE`, `level` must lie within the range covered on both
#' sides of the computed profile. If it does not, try recomputing the
#' profile with a larger `max_level` or a wider `focus_range`, as
#' appropriate, or inspecting the profile using `approach =
#' "focus_grid"` in the profile methods.  In particular, two-sided
#' confidence limits cannot be displayed from a focus-grid profile
#' lying entirely on one side of the MLE.  When `ci = FALSE`, the
#' corresponding cutoff is drawn but confidence limits are not
#' computed, so `level` need not be covered by the profile.
#' Confidence limits are extracted with [confint.profile_focus_list()]
#' using the selected interpolation method.
#'
#' If `signed = TRUE` and `what = "mpl"`, the plotted quantity is the signed
#' likelihood root based on the modified profile likelihood. With
#' `what = "rstar"`, the plotted quantity is the modified signed
#' likelihood-ratio statistic \eqn{r^*}. It has no corresponding
#' log-likelihood representation, so `signed` is treated as `TRUE`
#' irrespective of its supplied value.
#'
#' @return Called for its side effect of drawing a plot.
#'
#' @seealso [profile_focus()], [profile.focus_list_glm()],
#'     [modified_profile_focus()], [confint.profile_focus_list()],
#'     [profile_ci()]
#'
#' @export
plot.profile_focus_list <- function(x, level = 0.95, signed = FALSE,
                                    interpolation = c("linear", "cubic"),
                                    ci = FALSE,
                                    what = c("pl", "mpl", "rstar"), ...) {
    what <- match.arg(what)
    interpolation <- match.arg(interpolation)
    lin <- identical(interpolation, "linear")
    ci <- isTRUE(ci)
    if (what == "rstar")
        signed <- TRUE
    if (what != "pl" &&
        !inherits(x, "modified_profile_focus_list"))
        stop("`what = \"", what, "\"` requires an object returned by ",
             "`modified_profile_focus()`.")
    if (what == "pl") {
        profile_loglik <- x$loglik
        profile_maximum <- c(psi = unname(attr(x, "mle")),
                             loglik = unname(attr(x, "max_loglik")))
        signed_root <- .profile_signed(x)
    } else if (what == "mpl") {
        profile_loglik <- x$modified_loglik
        profile_maximum <- .modified_profile_maximum(x)
        signed_root <- sign(profile_maximum["psi"] - x$psi) *
            sqrt(2 * pmax(profile_maximum["loglik"] - profile_loglik, 0))
    } else {
        profile_loglik <- profile_maximum <- NULL
        signed_root <- x$rstar
    }
    qua <- qnorm(0.5 + level/2)
    keep <- is.finite(x$psi) & is.finite(signed_root)
    root_values <- signed_root[keep]
    focus_values <- x$psi[keep]
    fn <- if (lin) approxfun(x = root_values, y = focus_values) else
        splinefun(root_values, y = focus_values)
    r <- seq(min(root_values), max(root_values), length.out = 201)
    psi <- fn(r)
    if (ci)
        limits <- confint(x, level = level, method = what,
                          interpolation = interpolation)
    if (signed) {
        ylab <- switch(
            what,
            pl = "Signed likelihood root",
            mpl = "Signed likelihood root based on modified profile likelihood",
            rstar = "Modified signed likelihood-ratio statistic"
        )
        plot.default(x$psi, signed_root, pch = 21, bg = "lightgray",
                     col = "lightgray", cex = 0.8,
                     xlab = expression(psi), ylab = ylab, ...)
        points(psi, r, type = "l")
        abline(h = c(-qua, qua), lty = 3, col = "lightgray")
        if (what == "rstar") {
            if (min(root_values) <= 0 && max(root_values) >= 0)
                points(fn(0), 0, pch = 21, bg = "lightgray")
        } else {
            points(profile_maximum["psi"], 0, pch = 21, bg = "lightgray")
        }
        hei <- min(signed_root, na.rm = TRUE)
    } else {
        ylab <- if (what == "pl") "Profile log-likelihood" else
            "Modified profile log-likelihood"
        plot.default(x$psi, profile_loglik, pch = 21, bg = "lightgray",
                     col = "lightgray", cex = 0.8,
                     xlab = expression(psi), ylab = ylab, ...)
        points(psi, profile_maximum["loglik"] - r^2/2, type = "l")
        cutoff <- profile_maximum["loglik"] - qchisq(level, 1) / 2
        abline(h = cutoff, lty = 3, col = "lightgray")
        points(profile_maximum["psi"], profile_maximum["loglik"],
               pch = 21, bg = "lightgray")
        hei <- max(profile_loglik, na.rm = TRUE)
    }
    if (ci) {
        abline(v = limits, lty = 1, col = "lightgray")
        text(limits, hei, label = format(round(limits, 2), nsmall = 2, digits = 2))
        text(max(x$psi, na.rm = TRUE), hei, label = paste(round(100 * level, 2), "% CI", sep=""),
             adj = 1)
    }
}
