library("betareg")

data("GasolineYield", package = "betareg")

parameter <- function(theta, j)
    theta[j]

parameter_gradient <- function(theta, j)
    replace(numeric(length(theta)), j, 1)

parameter_hessian <- function(theta, j)
    matrix(0, length(theta), length(theta))

fit <- betareg(yield ~ temp | temp, data = GasolineYield, type = "ML")
theta <- coef(fit, model = "full")
aux <- enrichwith::get_auxiliary_functions(fit)
response <- model.response(model.frame(fit))

loglik <- function(theta, data) {
    sum(aux$dmodel(response = data, coefficients = theta, log = TRUE))
}

score <- function(theta, data) {
    aux$score(coefficients = theta, response = data)
}

information <- function(theta, data) {
    aux$information(coefficients = theta,
                    type = "observed",
                    response = data)
}

simulate <- function(theta) {
    aux$simulate(coefficients = theta, nsim = 1L)[, 1L]
}

fitted_focus <- focus(
    fit,
    correction = "no",
    on = parameter,
    on_gradient = parameter_gradient,
    on_hessian = parameter_hessian,
    j = 2
)

set.seed(1)
modified_low <- modified_profile_focus(
    loglik = loglik,
    score = score,
    information = information,
    simulate = simulate,
    mle = theta,
    on = parameter,
    on_gradient = parameter_gradient,
    on_hessian = parameter_hessian,
    likelihood_args = list(data = response),
    nsim = 50,
    grid_size = 2,
    max_level = 0.8,
    j = 2
)

set.seed(1)
modified_fitted <- modified_profile(
    fitted_focus,
    nsim = 50,
    grid_size = 2,
    max_level = 0.8
)

expect_equal(modified_fitted, modified_low, tolerance = 1e-10)
expect_true(inherits(modified_fitted, "modified_profile_focus_list"))
expect_identical(attr(modified_fitted, "nsim"), 50)
noncentral <- modified_fitted$signed_root != 0
expect_true(all(is.finite(modified_fitted$modified_loglik)))
expect_true(all(is.finite(modified_fitted$rstar[noncentral])))
