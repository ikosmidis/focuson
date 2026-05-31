library("brglm2")

data("coalition", package = "brglm2")
coalition_fit <- glm(duration ~ fract + numst2,
                     family = Gamma,
                     data = coalition,
                     method = "brglmFit",
                     type = "ML")

focus_square <- focus(coalition_fit,
                      on = function(theta, k) theta[k]^2,
                      correction = "median",
                      k = 2)
se_square <- focus_se(focus_square)
theta_square <- se_square$theta
se_square_control <- focus_se(focus_square, control = list(tol_opt = 1e-5))
se_square_partial_control <- focus_se(focus_square, control = list(tol_op = 1e-5))

expect_true(abs(unname(theta_square[2]^2) - unname(coef(focus_square))) < 1e-8)
expect_true(is.numeric(se_square$se))
expect_equal(dim(se_square$V), rep(length(theta_square), 2))
expect_equal(length(se_square$gradient), length(theta_square))
expect_equal(length(se_square$replace), 1L)
expect_true(is.numeric(se_square_control$se))
expect_true(is.numeric(se_square_partial_control$se))
expect_error(focus_se(focus_square, control = 1),
             pattern = "`control` must be a list")

all_square <- focus_on_all(focus_square)

expect_equal(rownames(all_square), c("estimate", "se"))
expect_equal(colnames(all_square), names(coef(focus_square$object, model = "full")))
expect_equal(length(theta_square), ncol(all_square))

focus_square_default <- focus(coalition_fit,
                              on = function(theta, k) theta[k]^2,
                              correction = "median",
                              k = 2)
ci_square_default <- confint(focus_square_default)
ci_square_compatible <- confint(focus_square_default,
                                se_at = "compatible",
                                se_control = list(tol_opt = 1e-5))

expect_equal(focus_square_default$estimate,
             focus_square$estimate,
             check.attributes = FALSE)
expect_equal(attr(ci_square_default, "se_at"), "supplied")
expect_equal(attr(ci_square_compatible, "se_at"), "compatible")
expect_true(is.null(attr(ci_square_default, "se_info")))
expect_true(is.list(attr(ci_square_compatible, "se_info")))
expect_equal(attr(ci_square_compatible, "se_info")$theta,
             se_square_control$theta,
             tolerance = 1e-8,
             check.attributes = FALSE)
expect_error(confint(focus_square_default,
                     se_at = "compatible",
                     se_control = 1),
             pattern = "`se_control` must be a list")

focus_exp <- focus(coalition_fit,
                   on = function(theta) exp(theta[2]),
                   correction = "median")
se_exp <- focus_se(focus_exp)

expect_true(abs(unname(exp(se_exp$theta[2]) - coef(focus_exp))) < 1e-4)
expect_true(is.numeric(se_exp$se))

expect_warning(
    ci_exp_fallback <- confint(focus_exp,
                               se_at = "compatible",
                               se_control = list(tol_opt = 1e-10))
)
expect_equal(attr(ci_exp_fallback, "se_at"), "supplied")
expect_true(is.null(attr(ci_exp_fallback, "se_info")))

warpbreaks_fit <- glm(breaks ~ wool + tension,
                      data = warpbreaks,
                      family = poisson,
                      method = "brglmFit",
                      type = "ML")

focus_warpbreaks <- focus(warpbreaks_fit,
                          on = function(theta, k) theta[k]^2,
                          correction = "median",
                          k = 2)
se_warpbreaks <- focus_se(focus_warpbreaks, control = list(tol_opt = 1e-5))

expect_true(is.numeric(se_warpbreaks$se))
expect_equal(dim(se_warpbreaks$V), rep(length(se_warpbreaks$theta), 2))
expect_equal(length(se_warpbreaks$theta), length(coef(focus_warpbreaks$object, model = "mean")))


data("endometrial", package = "brglm2")
endo <- glm(HG ~ NV + PI + EH,
            data = endometrial,
            family = binomial("probit"),
            method = "brglmFit")

focus_fun <- function(theta, i = 1, j = 2) theta[i] - theta[j]
fdiff0 <- focus(endo, on = focus_fun, i = 2, j = 3)

focus_all <- focus_on_all(fdiff0)

expect_warning(
    one_step_endo <- update(endo, type = "AS_median", start = coef(endo), maxit = 1, max_step_factor = 0)
)

expect_equal(focus_all["estimate", ], coef(one_step_endo), tol = 1e-06)

focus_all_correction <- sapply(1:4, function(j) {
    out <- focus(endo, on = function(theta) theta[j])
    c(estimate = coef(out), se = focus_se(out)$se)
    })

expect_equal(coef(summary(one_step_endo))[, 1:2], t(focus_all_correction),
             check.attributes = FALSE, tol = 1e-06)
