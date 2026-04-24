f_quad <- function(theta, shift = c(0, 0, 0)) {
    z <- theta + shift
    z[1]^2 + 3 * z[1] * z[2] + 4 * z[1] * z[3] +
        5 * z[2]^2 + 6 * z[2] * z[3] + 7 * z[3]^2
}

theta <- c(1, 2, 3)
shift <- c(0.5, -1, 2)
z <- theta + shift

expected_grad <- c(
    2 * z[1] + 3 * z[2] + 4 * z[3],
    3 * z[1] + 10 * z[2] + 6 * z[3],
    4 * z[1] + 6 * z[2] + 14 * z[3]
)
expected_hess <- matrix(
    c(2, 3, 4,
      3, 10, 6,
      4, 6, 14),
    nrow = 3,
    byrow = TRUE
)

gh <- focuson:::grad_hess(f_quad, theta, shift = shift)

expect_equal(gh$grad, expected_grad, tolerance = 1e-5)
expect_equal(gh$hess, expected_hess, tolerance = 1e-4)

gh_default <- focuson:::grad_hess(f_quad, theta, shift = shift)
gh_tuned <- focuson:::grad_hess(
    f_quad,
    theta,
    method.args = list(eps = 1e-4, d = 1e-4, r = 4, v = 2),
    shift = shift
)

expect_equal(gh_default$grad, gh_tuned$grad, tolerance = 1e-4)
expect_equal(gh_default$hess, gh_tuned$hess, tolerance = 1e-3)

expect_identical(names(gh), c("grad", "hess"))
expect_identical(length(gh$grad), length(theta))
expect_identical(dim(gh$hess), c(length(theta), length(theta)))
