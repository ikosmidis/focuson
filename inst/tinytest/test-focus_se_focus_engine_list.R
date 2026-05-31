library("brglm2")
library("enrichwith")

data("coalition", package = "brglm2")
coalition_fit <- glm(duration ~ fract + numst2,
                     family = Gamma,
                     data = coalition,
                     method = "brglmFit",
                     type = "ML")

afuns <- enrichwith::get_auxiliary_functions(coalition_fit)
theta <- coef(coalition_fit, model = "full")
components <- list(
    V = vcov(coalition_fit, model = "full"),
    P = afuns$Pmat(),
    Q = afuns$Qmat()
)

engine_out <- focus_engine(theta = theta,
                           components = components,
                           on = function(theta) exp(theta[1]),
                           correction = "mean")

engine_all <- focus_on_all(engine_out)


coalition_bc <- update(coalition_fit, type = "correction")


expect_equal(rownames(engine_all), c("estimate", "se"))
expect_equal(colnames(engine_all), names(theta))
expect_equal(ncol(engine_all), length(theta))
expect_equal(engine_all["estimate", ], coef(coalition_bc, model = "full"))

coordinate_one <- focus_engine(
    theta = theta,
    components = components,
    on = function(theta, j) theta[j],
    on_gradient = function(theta, j) {
        out <- numeric(length(theta))
        out[j] <- 1
        out
    },
    on_hessian = function(theta, j) {
        matrix(0, nrow = length(theta), ncol = length(theta))
    },
    correction = "mean",
    j = 1
)

expect_equal(engine_all["estimate", 1],
             coordinate_one$estimate,
             check.attributes = FALSE)
expect_equal(engine_all["se", 1],
             coordinate_one$se,
             check.attributes = FALSE)

p_mean <- length(coef(coalition_fit, model = "mean"))
V_function <- function(theta) {
    solve(afuns$information(coefficients = theta[seq_len(p_mean)],
                            dispersion = theta[-seq_len(p_mean)]))
}

engine_square <- focus_engine(theta = theta,
                              components = components,
                              on = function(theta, k) theta[k]^2,
                              correction = "median",
                              k = 2)
engine_square_se <- focus_se(engine_square, V_function = V_function)

expect_true(is.numeric(engine_square_se$se))
expect_equal(dim(engine_square_se$V), rep(length(engine_square_se$theta), 2))
expect_equal(length(engine_square_se$gradient), length(engine_square_se$theta))
expect_equal(length(engine_square_se$replace), 1L)

engine_square_corrected <- confint(engine_square,
                                   se_at = "corrected",
                                   V_function = V_function)

expect_equal(attr(engine_square_corrected, "se_at"), "corrected")
expect_true(is.list(attr(engine_square_corrected, "se_info")))
expect_equal(attr(engine_square_corrected, "se_info")$theta,
             engine_square_se$theta,
             tolerance = 1e-8,
             check.attributes = FALSE)

expect_error(focus_se(engine_square),
             pattern = "`V_function` must be a function")

expect_warning(
    engine_square_fallback <- confint(engine_square,
                                      se_at = "corrected")
)
expect_equal(attr(engine_square_fallback, "se_at"), "naive")
expect_true(is.null(attr(engine_square_fallback, "se_info")))

coalition_mean <- update(coalition_fit, type = "AS_mean")
theta_mean <- coef(coalition_mean, model = "full")
components_mean <- list(V = vcov(coalition_mean, model = "full"),
                        P = afuns$Pmat(coef(coalition_mean, "mean"),
                                       coef(coalition_mean, "dispersion")),
                        Q = afuns$Qmat(coef(coalition_mean, "mean"),
                                       coef(coalition_mean, "dispersion")))

engine_mean <- focus_engine(theta = theta_mean,
                            components = components_mean,
                            on = function(theta) theta[1] - theta[4],
                            correction = "median",
                            estimator = "meanBR")

engine_all_mean <- focus_on_all(engine_mean)

expect_equal(rownames(engine_all_mean), c("estimate", "se"))
expect_equal(colnames(engine_all_mean), names(theta_mean))
expect_equal(ncol(engine_all_mean), length(theta_mean))

coordinate_mean_one <- focus_engine(
    theta = theta_mean,
    components = components_mean,
    on = function(theta, j) theta[j],
    on_gradient = function(theta, j) {
        out <- numeric(length(theta))
        out[j] <- 1
        out
    },
    on_hessian = function(theta, j) {
        matrix(0, nrow = length(theta), ncol = length(theta))
    },
    correction = "median",
    estimator = "meanBR",
    j = 1
)

expect_equal(engine_all_mean["estimate", 1],
             coordinate_mean_one$estimate,
             check.attributes = FALSE)
expect_equal(engine_all_mean["se", 1],
             coordinate_mean_one$se,
             check.attributes = FALSE)

expect_warning(
    coalition_xmed <- update(coalition_mean,
                             type = "AS_median",
                             maxit = 1,
                             max_step_factor = 0,
                             start = theta_mean)
)

expect_equal(coef(coalition_xmed, model = "full"), engine_all_mean["estimate", ],
             tol = 1e-03)

data("endometrial", package = "brglm2")
endo <- glm(HG ~ NV + PI + EH,
            data = endometrial,
            family = binomial("probit"),
            method = "brglmFit", type = "AS_mean")

afuns <- enrichwith::get_auxiliary_functions(endo)
comps <- list(V = solve(afuns$information()),
              P = afuns$Pmat(),
              Q = afuns$Qmat())

focus_fun <- function(theta, i = 1, j = 2) theta[i] - theta[j]

fdiff0 <- focus(endo, on = focus_fun, i = 2, j = 3)
fdiff1 <- focus_engine(coef(endo), comps, focus_fun,
                       estimator = "meanBR", i = 2, j = 3)

expect_equal(coef(fdiff0), coef(fdiff1))
engine_ests <- focus_on_all(fdiff1)

expect_warning(
    one_step_endo <- update(endo, type = "AS_median", start = coef(endo), maxit = 1, max_step_factor = 0)
)

expect_equal(coef(one_step_endo), engine_ests["estimate", ], tol = 1e-06)

focus_all_correction <- sapply(1:4, function(j) {
    out <- focus_engine(coef(endo), comps, on = function(theta) theta[j], estimator = "meanBR")
    c(estimate = coef(out), se = focus_se(out, V_function = function(theta) solve(afuns$information(theta)))$se)
    })

expect_equal(coef(summary(one_step_endo))[, 1:2], t(focus_all_correction),
             check.attributes = FALSE, tol = 1e-06)
