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


f_base <- focus(
    endo,
    on = function(theta) theta[1],
    correction = "median"
)$estimate
f_scaled <- focus(
    endo,
    on = function(theta) 1e-9 * theta[1],
    correction = "median"
)$estimate

expect_equal(
    unname(f_scaled),
    1e-9 * unname(f_base),
    tolerance = 1e-12,
    check.attributes = FALSE
)
