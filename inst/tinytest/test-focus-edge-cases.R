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
