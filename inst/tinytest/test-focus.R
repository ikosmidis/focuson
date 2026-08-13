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
    randomize = 2
)

expect_true(is.numeric(out$estimate))
expect_true(is.numeric(out$se))
set.seed(678)
ci_hulc <- confint(out, method = "hulc", check_statistic = FALSE)
expect_true(is.numeric(ci_hulc))
expect_identical(attr(ci_hulc, "type"), "hulc")


set.seed(678)
out <- focus(
    endo,
    on = function(theta, check_statistic = 2) check_statistic * theta[1],
    correction = "no",
    check_statistic = 2
)

set.seed(678)
ci_hulc <- confint(out, method = "hulc", check_statistic = FALSE)
expect_identical(attr(ci_hulc, "type"), "hulc")
expect_equal(
    unname(out$estimate),
    2 * unname(focus(endo, on = function(theta) theta[1], correction = "no")$estimate)
)


out_wald <- focus(endo, correction = "mean")
expect_identical(out_wald$correction, "mean")
expect_identical(unname(coef(out_wald)), unname(out_wald$estimate))
expect_equal(drop(vcov(out_wald)), out_wald$se^2)
expect_true(inherits(confint(out_wald), "focus_ci"))
expect_equal(
    as.numeric(confint(out_wald)),
    unname(out_wald$estimate) + c(-1, 1) * qnorm(0.975) * out_wald$se,
    check.attributes = FALSE
)
expect_equal(
    as.numeric(confint(out_wald, level = 0.9)),
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
ci_hulc_direct <- hulc_ci(
    data = focuson:::.refit_data(out_wald_i2$object),
    statistic = function(data) {
        focus_statistic(
            data = data,
            object = out_wald_i2$object,
            on = on_index,
            correction = "no",
            on_gradient = grad_index,
            on_hessian = hess_index,
            i = 2
        )
    },
    level = 0.95,
    randomize = FALSE,
    check_statistic = FALSE
)

expect_equal(ci_hulc_method, ci_hulc_direct, check.attributes = FALSE)
expect_identical(attr(ci_hulc_method, "type"), "hulc")


for (correction in c("no", "mean", "median")) {
    out <- focus(endo, on = function(theta) pi, correction = correction)

    expect_equal(unname(out$estimate), pi, check.attributes = FALSE)
    expect_equal(out$se, 0)
    expect_equal(as.numeric(confint(out)), c(pi, pi))

    set.seed(678)
    ci_hulc <- confint(
        out,
        method = "hulc",
        check_statistic = FALSE
    )

    expect_equal(as.numeric(ci_hulc), c(pi, pi))
    expect_identical(attr(ci_hulc, "type"), "hulc")
}


f_base <- focus(endo, on = function(theta) theta[1], correction = "median")$estimate
f_scaled <- focus(endo, on = function(theta) 1e-9 * theta[1], correction = "median")$estimate

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

expect_identical(out_ana$correction, "median")
expect_equal(out_ana$estimate, out_num$estimate, tolerance = 1e-8, check.attributes = FALSE)
expect_equal(out_ana$se, out_num$se, tolerance = 1e-8)
expect_equal(confint(out_ana), confint(out_num), tolerance = 1e-8, check.attributes = FALSE)
expect_identical(out_ana$on$on, on_or)
expect_identical(out_ana$on$on_gradient, grad_or)
expect_identical(out_ana$on$on_hessian, hessian_or)

set.seed(678)
ci_hulc <- confint(
    out_ana,
    method = "hulc",
    check_statistic = FALSE
)

expect_true(is.numeric(out_ana$estimate))
expect_identical(attr(ci_hulc, "type"), "hulc")


budworm_fit <- glm(cbind(numalive, numdead) ~ sex * ldose,
                   family = binomial,
                   data = transform(
                       data.frame(ldose = rep(0:5, 2),
                                  numdead = c(1, 4, 9, 13, 18, 20,
                                              0, 2, 6, 10, 12, 16),
                                  sex = factor(rep(c("M", "F"), c(6, 6)))),
                       numalive = 20 - numdead
                   ))
budworm_glm_ci <- confint(profile(budworm_fit), parm = 2, level = 0.95)
budworm_focus <- focus(budworm_fit,
                       on = on_index,
                       on_gradient = grad_index,
                       on_hessian = hess_index,
                       correction = "no",
                       i = 2)
budworm_profile_ci <- confint(budworm_focus, method = "pl", level = 0.95)
expect_equal(as.numeric(budworm_profile_ci), as.numeric(budworm_glm_ci),
             tolerance = 1e-04)
expect_identical(attr(budworm_profile_ci, "type"), "pl")

budworm_focus_median <- focus(budworm_fit,
                              on = on_index,
                              on_gradient = grad_index,
                              on_hessian = hess_index,
                              correction = "median",
                              i = 2)
budworm_profile_ci_median <- confint(budworm_focus_median,
                                     method = "pl",
                                     level = 0.95)
expect_equal(as.numeric(budworm_profile_ci_median),
             as.numeric(budworm_glm_ci), tolerance = 1e-04)
expect_identical(attr(budworm_profile_ci_median, "type"), "pl")

budworm_refit_data <- focuson:::.refit_data(budworm_focus$object)
expect_true(all(c("numalive", "numdead", "sex", "ldose") %in% names(budworm_refit_data)))
budworm_refit_stat <- focus_statistic(budworm_refit_data,
                                      budworm_focus$object,
                                      on = on_index,
                                      on_gradient = grad_index,
                                      on_hessian = hess_index,
                                      correction = "no",
                                      i = 2)
expect_equal(budworm_refit_stat, budworm_focus$estimate,
             tolerance = 1e-08, check.attributes = FALSE)
