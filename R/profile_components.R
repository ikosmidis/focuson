.profile_components_glm <- function(fitted) {
    fit <- fitted$object
    if (!identical(fit$type, "ML"))
        fit <- update(fit, type = "ML", start = coef(fit, model = "mean"))
    p_mean <- length(coef(fit, model = "mean"))
    theta <- coef(fit, model = "full")
    if (fit$family$family %in% c("poisson", "binomial"))
        theta <- theta[names(coef(fit, model = "mean"))]
    aux <- enrichwith::get_auxiliary_functions(fit)
    response <- stats::model.response(stats::model.frame(fit))
    m_inds <- seq_len(p_mean)
    d_ind <- p_mean + 1L
    get_dispersion <- function(theta) {
        if (length(theta) == p_mean)
            1
        else
            theta[d_ind]
    }
    loglik <- function(theta, data) {
        sum(aux$dmodel(
            response = data,
            coefficients = theta[m_inds],
            dispersion = get_dispersion(theta),
            log = TRUE
        ))
    }
    score <- function(theta, data) {
        aux$score(
            coefficients = theta[m_inds],
            dispersion = get_dispersion(theta),
            response = data
        )
    }
    information <- function(theta, data) {
        aux$information(
            coefficients = theta[m_inds],
            dispersion = get_dispersion(theta),
            type = "observed",
            response = data
        )
    }
    simulate <- function(theta) {
        aux$simulate(
            coefficients = theta[m_inds],
            dispersion = get_dispersion(theta),
            nsim = 1L
        )[[1L]]
    }
    list(loglik = loglik,
         score = score,
         information = information,
         simulate = simulate,
         mle = theta,
         likelihood_args = list(data = response))
}

.profile_components_betareg <- function(fitted) {
    fit <- fitted$object
    if (!identical(fit$type, "ML"))
        fit <- update(fit, type = "ML", start = coef(fit, model = "full"))
    theta <- coef(fit, model = "full")
    aux <- enrichwith::get_auxiliary_functions(fit)
    response <- stats::model.response(stats::model.frame(fit))
    loglik <- function(theta, data) {
        sum(aux$dmodel(response = data,
                       coefficients = theta,
                       log = TRUE))
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
    list(loglik = loglik,
         score = score,
         information = information,
         simulate = simulate,
         mle = theta,
         likelihood_args = list(data = response))
}
