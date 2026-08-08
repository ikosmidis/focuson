library("betareg")

data("GasolineYield", package = "betareg")

on_square <- function(theta, j)
    theta[j]^2

on_square_gradient <- function(theta, j)
    replace(numeric(length(theta)), j, 2 * theta[j])

on_square_hessian <- function(theta, j) {
    out <- matrix(0, length(theta), length(theta))
    out[j, j] <- 2
    out
}

on_coordinate <- function(theta, j)
    theta[j]

on_coordinate_gradient <- function(theta, j)
    replace(numeric(length(theta)), j, 1)

on_coordinate_hessian <- function(theta, j)
    matrix(0, length(theta), length(theta))

for (type in c("ML", "BC", "BR")) {
    fit <- betareg(yield ~ temp | temp, data = GasolineYield, type = type)
    theta <- coef(fit, model = "full")
    aux <- enrichwith::get_auxiliary_functions(fit)
    components <- list(
        V = vcov(fit, model = "full"),
        P = aux$Pmat(coefficients = theta),
        Q = aux$Qmat(coefficients = theta)
    )
    estimator <- if (type == "ML") "ML" else "meanBR"

    for (correction in c("no", "mean", "median")) {
        fitted_focus <- focus(
            fit,
            on = on_square,
            correction = correction,
            on_gradient = on_square_gradient,
            on_hessian = on_square_hessian,
            j = 2
        )
        engine_focus <- focus_engine(
            theta,
            components,
            on = on_square,
            correction = correction,
            estimator = estimator,
            on_gradient = on_square_gradient,
            on_hessian = on_square_hessian,
            j = 2
        )

        expect_equal(fitted_focus$estimate, engine_focus$estimate,
                     tolerance = 1e-10)
        expect_equal(fitted_focus$se, engine_focus$se,
                     tolerance = 1e-10)
        expect_identical(fitted_focus$correction, correction)
        expect_identical(fitted_focus$object, fit)
        expect_true(inherits(fitted_focus, "focus_list_betareg"))
    }
}

fit <- betareg(yield ~ temp | temp, data = GasolineYield, type = "ML")
theta <- coef(fit, model = "full")
V <- vcov(fit, model = "full")

coordinate_focus <- focus(
    fit,
    on = on_coordinate,
    correction = "no",
    on_gradient = on_coordinate_gradient,
    on_hessian = on_coordinate_hessian,
    j = 1
)
all_coordinates <- focus_on_all(coordinate_focus)

expect_equal(all_coordinates["estimate", ], theta,
             check.attributes = FALSE)
expect_equal(all_coordinates["se", ], sqrt(diag(V)),
             check.attributes = FALSE)
expect_identical(colnames(all_coordinates), names(theta))
expect_identical(rownames(all_coordinates), c("estimate", "se"))

square_focus <- focus(
    fit,
    on = on_square,
    correction = "median",
    on_gradient = on_square_gradient,
    on_hessian = on_square_hessian,
    j = 2
)
all_corrected <- focus_on_all(square_focus)
se_corrected <- focus_se(square_focus)
aux <- enrichwith::get_auxiliary_functions(fit)

expect_equal(dim(all_corrected), c(2L, length(theta)))
expect_true(all(is.finite(all_corrected)))
expect_true(abs(unname(se_corrected$theta[2]^2) - coef(square_focus)) < 1e-8)
expect_equal(se_corrected$V,
             solve(aux$information(coefficients = se_corrected$theta)),
             tolerance = 1e-10)
expect_equal(length(se_corrected$gradient), length(theta))
expect_identical(se_corrected$replace, 2L)

ci_compatible <- confint(square_focus, se_at = "compatible")
expect_identical(attr(ci_compatible, "se_at"), "compatible")
expect_true(is.list(attr(ci_compatible, "se_info")))

gy_ml <- betareg(yield ~ batch + temp, data = GasolineYield, type = "ML")
gy_focus_no <- focus(gy_ml, on = function(theta) theta[12], correction = "no")
expect_equal(coef(gy_ml, model = "precision"), coef(gy_focus_no), check.attributes = FALSE)

gy_bc <- betareg(yield ~ batch + temp, data = GasolineYield, type = "BC")
gy_focus_mean <- focus(gy_ml, on = function(theta) theta[12], correction = "mean")
expect_equal(coef(gy_bc, model = "precision"), coef(gy_focus_mean), check.attributes = FALSE)

gy_focus_median <- focus(gy_ml, on = function(theta) theta[12], correction = "median")

## Median bias corrected estimates using mbrbetareg with fsmaxit = 1
## starting from ML estimates
## using https://github.com/eulogepagui/mbrbetareg
## 3b047581b9372e9d7d62c66e869277b6ab352572
## with call
## ml <- betareg(yield ~ batch + temp, data = GasolineYield)
## aa <- mbrbetareg(yield ~ batch + temp, data = GasolineYield,
##                 fsmaxit = 1, start = coef(ml))
coefs_mbrb <- c(-6.1489480161980, 1.7252930172155, 1.3202250313036,
                1.5694960561405, 1.0582460707296, 1.1319896689429,
                1.0385832713058, 0.5432673494269, 0.4953643257657,
                0.3849571539474, 0.0109458927842, 279.5435176253078)
ses_mbrb <- c(2.28163611105e-01, 1.26745274805e-01, 1.47570987799e-01,
              1.45343189263e-01, 1.28135865728e-01, 1.29616013392e-01,
              1.32753175470e-01, 1.36542803846e-01, 1.36335753038e-01,
              1.48467057325e-01, 5.16331443936e-04, 6.98427886687e+01)

coefs_focus <- ses_focus <- numeric(12)
for (j in 1:12) {
    fobj <- focus(gy_ml, on = function(theta) theta[j])
    coefs_focus[j] <- coef(fobj)
    ses_focus[j] <- focus_se(fobj)$se
}
expect_equal(coefs_focus, coefs_mbrb)
expect_equal(ses_focus, ses_mbrb)

expect_equal(focus_on_all(gy_focus_median)["estimate", ], coefs_mbrb,
             check.attributes = FALSE)
expect_equal(focus_on_all(gy_focus_no)["estimate", ], coef(gy_ml),
             check.attributes = FALSE)
expect_equal(focus_on_all(gy_focus_mean)["estimate", ], coef(gy_bc),
             check.attributes = FALSE)



