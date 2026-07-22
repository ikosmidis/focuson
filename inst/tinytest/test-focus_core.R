theta <- log(2)
V <- matrix(0.25, 1, 1)
bias_ml <- 0.1
P <- list(matrix(0.6, 1, 1))
Q <- list(matrix(-0.2, 1, 1))

on <- function(theta, power) exp(theta)^power
on_gradient <- function(theta, power) power * exp(theta)^power
on_hessian <- function(theta, power) {
    matrix(power^2 * exp(theta)^power, 1, 1)
}

## At theta = log(2) and power = 1:
## hat_psi = 2, d1_psi = d2_psi = 2, and var_psi = 1.
## For ML, mean_b = 2 * 0.1 + 0.5 * 2 * 0.25 = 0.45.
## The median adjustment is 0.2375.
core_no <- focuson:::.focus_core(
    theta = theta,
    V = V,
    on = on,
    correction = "no",
    components_fun = function() stop("supplier was called"),
    power = 1
)
expect_equal(core_no$estimate, 2)
expect_equal(core_no$se, 1)

core_mean_ml <- focuson:::.focus_core(
    theta = theta,
    V = V,
    components_fun = function() list(bias = bias_ml),
    on = on,
    correction = "mean",
    on_gradient = on_gradient,
    on_hessian = on_hessian,
    power = 1
)
expect_equal(core_mean_ml$estimate, 1.55)
expect_equal(core_mean_ml$se, 1)

core_median_ml <- focuson:::.focus_core(
    theta = theta,
    V = V,
    components_fun = function() {
        list(bias = bias_ml, P = P, Q = Q)
    },
    on = on,
    correction = "median",
    on_gradient = on_gradient,
    on_hessian = on_hessian,
    power = 1
)
expect_equal(core_median_ml$estimate, 1.7875)
expect_equal(core_median_ml$se, 1)

zero_bias <- numeric(length(theta))
core_mean_br <- focuson:::.focus_core(
    theta = theta,
    V = V,
    components_fun = function() list(bias = zero_bias),
    on = on,
    correction = "mean",
    on_gradient = on_gradient,
    on_hessian = on_hessian,
    power = 1
)
expect_equal(core_mean_br$estimate, 1.75)

core_median_br <- focuson:::.focus_core(
    theta = theta,
    V = V,
    components_fun = function() {
        list(bias = zero_bias, P = P, Q = Q)
    },
    on = on,
    correction = "median",
    on_gradient = on_gradient,
    on_hessian = on_hessian,
    power = 1
)
expect_equal(core_median_br$estimate, 1.9875)

core_numeric <- focuson:::.focus_core(
    theta = theta,
    V = V,
    components_fun = function() {
        list(bias = bias_ml, P = P, Q = Q)
    },
    on = on,
    correction = "median",
    power = 1
)
expect_equal(core_numeric$estimate, core_median_ml$estimate,
             tolerance = 1e-06)
expect_equal(core_numeric$se, core_median_ml$se,
             tolerance = 1e-06)

constant <- focuson:::.focus_core(
    theta = theta,
    V = V,
    on = function(theta) 2,
    correction = "median",
    components_fun = function() stop("supplier was called"),
    on_gradient = function(theta) 0,
    on_hessian = function(theta) matrix(0, 1, 1)
)
expect_identical(constant$estimate, 2)
expect_identical(constant$se, 0)

expect_error(
    focuson:::.focus_core(
        theta = theta,
        V = V,
        on = on,
        correction = "mean",
        components_fun = function() list(),
        on_gradient = on_gradient,
        on_hessian = on_hessian,
        power = 1
    ),
    pattern = "`bias` is required"
)
expect_error(
    focuson:::.focus_core(
        theta = theta,
        V = V,
        components_fun = function() list(bias = zero_bias),
        on = on,
        correction = "median",
        on_gradient = on_gradient,
        on_hessian = on_hessian,
        power = 1
    ),
    pattern = "`P` and `Q` are required"
)
