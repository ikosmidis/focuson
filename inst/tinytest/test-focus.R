library("brglm2")

data("endometrial", package = "brglm2")
endo <- glm(HG ~ NV + PI + EH,
            data = endometrial,
            family = binomial("logit"),
            method = "brglmFit")

set.seed(678)
out <- focus(
    endo,
    on = function(theta, randomize = 2) randomize * theta[1],
    ci = "hulc",
    randomize = 2,
    control_ci = ci_control(check_statistic = FALSE)
)

expect_identical(out$ci_type, "hulc")
expect_true(is.numeric(out$estimate))
expect_true(is.numeric(out$confint))
expect_identical(attr(out$confint, "type"), "hulc")


set.seed(678)
out <- focus(
    endo,
    on = function(theta, check_statistic = 2) check_statistic * theta[1],
    correction = "no",
    ci = "hulc",
    check_statistic = 2,
    control_ci = ci_control(check_statistic = FALSE)
)

expect_identical(out$ci_type, "hulc")
expect_true(is.numeric(out$estimate))
expect_identical(attr(out$confint, "type"), "hulc")
expect_equal(
    unname(out$estimate),
    2 * unname(focus(endo, on = function(theta) theta[1], correction = "no")$estimate)
)


out_wald <- focus(endo, correction = "mean")
expect_equal(confint(out_wald), out_wald$confint)
expect_equal(
    unname(confint(out_wald, level = 0.9)),
    unname(out_wald$estimate) + c(-1, 1) * qnorm(0.95) * out_wald$se,
    check.attributes = FALSE
)

on_index <- function(theta, i = 1) theta[i]
grad_index <- function(theta, i = 1) {
    out <- rep(0, length(theta))
    out[i] <- 1
    out
}
hess_index <- function(theta, i = 1) matrix(0, nrow = length(theta), ncol = length(theta))

out_wald_i2 <- focus(
    endo,
    on = on_index,
    on_gradient = grad_index,
    on_hessian = hess_index,
    i = 2,
    correction = "no"
)
set.seed(678)
ci_hulc_method <- confint(
    out_wald_i2,
    method = "hulc",
    randomize = FALSE,
    check_statistic = FALSE
)
set.seed(678)
ci_hulc_focus <- focus(
    endo,
    on = on_index,
    on_gradient = grad_index,
    on_hessian = hess_index,
    i = 2,
    correction = "no",
    ci = "hulc",
    control_ci = ci_control(randomize = FALSE, check_statistic = FALSE)
)$confint

expect_equal(ci_hulc_method, ci_hulc_focus, check.attributes = FALSE)
expect_identical(attr(ci_hulc_method, "type"), "hulc")


for (correction in c("no", "mean", "median")) {
    out <- focus(endo, on = function(theta) pi, correction = correction)

    expect_equal(unname(out$estimate), pi, check.attributes = FALSE)
    expect_equal(out$se, 0)
    expect_equal(unname(out$confint), c(pi, pi), check.attributes = FALSE)

    set.seed(678)
    out_hulc <- focus(
        endo,
        on = function(theta) pi,
        correction = correction,
        ci = "hulc",
        control_ci = ci_control(check_statistic = FALSE)
    )

    expect_equal(unname(out_hulc$estimate), pi, check.attributes = FALSE)
    expect_equal(out_hulc$se, 0)
    expect_equal(unname(out_hulc$confint), c(pi, pi), check.attributes = FALSE)
    expect_identical(attr(out_hulc$confint, "type"), "hulc")
}


expect_error(focus(endo, alpha = 0))
expect_error(focus(endo, alpha = 1))
expect_error(focus(endo, alpha = 2))


f_base <- focus(endo, on = function(theta) theta[1], correction = "median")$estimate
f_scaled <- focus(endo,  on = function(theta) 1e-9 * theta[1],  correction = "median")$estimate

expect_equal(
    unname(f_scaled),
    1e-9 * unname(f_base),
    tolerance = 1e-12,
    check.attributes = FALSE
)


grad_or <- function(theta) {
    out <- rep(0, length(theta))
    out[1] <- exp(theta[1])
    out
}

hessian_or <- function(theta) {
    out <- matrix(0, nrow = length(theta), ncol = length(theta))
    out[1, 1] <- exp(theta[1])
    out
}

on_or <- function(theta) exp(theta[1])

out_num <- focus(
    endo,
    on = on_or,
    correction = "median"
)
out_ana <- focus(
    endo,
    on = on_or,
    correction = "median",
    on_gradient = grad_or,
    on_hessian = hessian_or
)

expect_equal(out_ana$estimate, out_num$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(out_ana$se, out_num$se, tolerance = 1e-8)
expect_equal(out_ana$confint, out_num$confint, tolerance = 1e-8, check.attributes = FALSE)
expect_identical(out_ana$on$on, on_or)
expect_identical(out_ana$on$on_gradient, grad_or)
expect_identical(out_ana$on$on_hessian, hessian_or)

set.seed(678)
out_hulc <- focus(
    endo,
    on = on_or,
    correction = "median",
    ci = "hulc",
    on_gradient = grad_or,
    on_hessian = hessian_or,
    control_ci = ci_control(check_statistic = FALSE)
)

expect_identical(out_hulc$ci_type, "hulc")
expect_true(is.numeric(out_hulc$estimate))
expect_identical(attr(out_hulc$confint, "type"), "hulc")
