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

fitted_focus <- focus(
    fit,
    on = parameter,
    on_gradient = parameter_gradient,
    on_hessian = parameter_hessian,
    j = 2
)

profile_low <- profile_focus(
    loglik = loglik,
    score = score,
    information = information,
    mle = theta,
    on = parameter,
    on_gradient = parameter_gradient,
    on_hessian = parameter_hessian,
    likelihood_args = list(data = response),
    grid_size = 3,
    max_level = 0.9,
    j = 2
)
profile_fitted <- profile(fitted_focus, grid_size = 3, max_level = 0.9)

expect_equal(profile_fitted, profile_low, tolerance = 1e-10)
expect_true(inherits(profile_fitted, "profile_focus_list"))
expect_equal(nrow(profile_fitted), 7L)
expect_equal(ncol(profile_fitted$theta), length(theta))
expect_identical(colnames(profile_fitted$theta), names(theta))

ci_low <- profile_ci(
    loglik = loglik,
    score = score,
    information = information,
    mle = theta,
    on = parameter,
    on_gradient = parameter_gradient,
    on_hessian = parameter_hessian,
    likelihood_args = list(data = response),
    level = 0.9,
    j = 2
)
ci_fitted <- confint(fitted_focus, method = "pl", level = 0.9)

expect_equal(ci_fitted, ci_low, tolerance = 1e-10)
expect_identical(attr(ci_fitted, "type"), "pl")

focus_range <- theta[2] + c(-2, 2) * sqrt(vcov(fit, model = "full")[2, 2])
profile_grid <- profile(
    fitted_focus,
    approach = "focus_grid",
    focus_range = focus_range,
    grid_size = 3
)

expect_identical(attr(profile_grid, "approach"), "focus_grid")
expect_equal(range(profile_grid$psi), sort(focus_range), tolerance = 1e-6)
expect_true(max(abs(as.numeric(profile_grid$theta[, 2]) -
                    profile_grid$psi)) < 1e-6)
expect_true(all(is.finite(profile_grid$loglik)))

for (type in c("BC", "BR")) {
    adjusted_fit <- betareg(yield ~ temp | temp,
                            data = GasolineYield,
                            type = type)
    adjusted_focus <- focus(
        adjusted_fit,
        on = parameter,
        on_gradient = parameter_gradient,
        on_hessian = parameter_hessian,
        j = 2
    )
    ml_fit <- update(adjusted_fit,
                     type = "ML",
                     start = coef(adjusted_fit, model = "full"))
    adjusted_profile <- profile(adjusted_focus,
                                grid_size = 2,
                                max_level = 0.8)

    expect_equal(attr(adjusted_profile, "theta_mle"),
                 coef(ml_fit, model = "full"),
                 tolerance = 1e-10)
    expect_equal(attr(adjusted_profile, "max_loglik"),
                 as.numeric(logLik(ml_fit)),
                 tolerance = 1e-10)
}
