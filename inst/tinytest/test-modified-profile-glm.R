budworm <- data.frame(
    ldose = rep(0:5, 2),
    numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16),
    sex = factor(rep(c("M", "F"), c(6, 6)))
)
budworm$numalive <- 20 - budworm$numdead
budworm_fit <- glm(
    cbind(numalive, numdead) ~ sex * ldose,
    family = binomial,
    data = budworm,
    method = brglm2::brglmFit,
    type = "ML"
)
budworm_focus <- focus(
    budworm_fit,
    correction = "no",
    on = function(theta) theta[1],
    on_gradient = function(theta)
        replace(numeric(length(theta)), 1L, 1),
    on_hessian = function(theta)
        matrix(0, length(theta), length(theta))
)

set.seed(1)
budworm_modified <- modified_profile(
    budworm_focus,
    nsim = 50,
    grid_size = 2,
    max_level = 0.8
)

expect_true(inherits(budworm_modified, "modified_profile_focus_list"))
expect_true(inherits(budworm_modified, "profile_focus_list"))
expect_equal(nrow(budworm_modified), 5)
expect_identical(attr(budworm_modified, "nsim"), 50)
expect_identical(ncol(budworm_modified$theta),
                 length(coef(budworm_fit)))
noncentral <- budworm_modified$signed_root != 0
expect_true(all(is.finite(budworm_modified$modified_loglik)))
expect_true(all(is.finite(budworm_modified$rstar[noncentral])))

set.seed(1)
budworm_mpl_ci <- confint(
    budworm_focus,
    method = "mpl",
    level = 0.7,
    nsim = 50,
    grid_size = 2,
    max_level = 0.8,
    interpolation = "cubic"
)
expect_equal(
    budworm_mpl_ci,
    confint(budworm_modified,
            method = "mpl",
            level = 0.7,
            interpolation = "cubic"),
    tolerance = 1e-10
)

set.seed(1)
budworm_rstar_ci <- confint(
    budworm_focus,
    method = "rstar",
    level = 0.7,
    nsim = 50,
    grid_size = 2,
    max_level = 0.8
)
expect_equal(
    budworm_rstar_ci,
    confint(budworm_modified, method = "rstar", level = 0.7),
    tolerance = 1e-10
)

gaussian_fit <- glm(
    mpg ~ wt,
    data = mtcars,
    method = brglm2::brglmFit,
    type = "ML"
)
gaussian_focus <- focus(
    gaussian_fit,
    correction = "no",
    on = function(theta) theta[1],
    on_gradient = function(theta)
        replace(numeric(length(theta)), 1L, 1),
    on_hessian = function(theta)
        matrix(0, length(theta), length(theta))
)

set.seed(2)
gaussian_modified <- modified_profile(
    gaussian_focus,
    nsim = 50,
    grid_size = 2,
    max_level = 0.8
)

expect_true(inherits(gaussian_modified, "modified_profile_focus_list"))
expect_identical(ncol(gaussian_modified$theta),
                 length(coef(gaussian_fit, model = "full")))
gaussian_noncentral <- gaussian_modified$signed_root != 0
expect_true(all(is.finite(gaussian_modified$rstar[gaussian_noncentral])))
