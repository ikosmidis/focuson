library("brglm2")

data("coalition", package = "brglm2")
coalition_fit <- glm(duration ~ fract + numst2,
                     family = Gamma,
                     data = coalition,
                     method = "brglmFit",
                     type = "ML")

focus_square <- focus(coalition_fit,
                      on = function(theta, k) theta[k]^2,
                      correction = "median",
                      k = 2)
se_square <- focus_se(focus_square)
theta_square <- se_square$theta

expect_true(abs(unname(theta_square[2]^2) - unname(coef(focus_square))) < 1e-8)
expect_true(is.numeric(se_square$se))
expect_equal(dim(se_square$V), rep(length(theta_square), 2))
expect_equal(length(se_square$gradient), length(theta_square))
expect_equal(length(se_square$replace), 1L)

all_square <- focus_on_all(focus_square)

expect_equal(rownames(all_square), c("estimate", "se"))
expect_equal(colnames(all_square), names(coef(focus_square$object, model = "full")))
expect_equal(length(theta_square), ncol(all_square))

focus_square_default <- focus(coalition_fit,
                              on = function(theta, k) theta[k]^2,
                              correction = "median",
                              se_at = "object",
                              k = 2)
focus_square_corrected <- focus(coalition_fit,
                                on = function(theta, k) theta[k]^2,
                                correction = "median",
                                se_at = "corrected",
                                k = 2)

expect_equal(focus_square_default$estimate,
             focus_square_corrected$estimate,
             check.attributes = FALSE)
expect_equal(focus_square_default$se_at, "object")
expect_equal(focus_square_corrected$se_at, "corrected")
expect_true(is.null(focus_square_default$se_info))
expect_true(is.list(focus_square_corrected$se_info))
expect_equal(focus_square_corrected$se,
             focus_square_corrected$se_info$se,
             check.attributes = FALSE)
expect_equal(focus_square_corrected$se_info$theta,
             se_square$theta,
             tolerance = 1e-8,
             check.attributes = FALSE)

focus_exp <- focus(coalition_fit,
                   on = function(theta) exp(theta[2]),
                   correction = "median")
se_exp <- focus_se(focus_exp)

expect_true(abs(unname(exp(se_exp$theta[2]) - coef(focus_exp))) < 1e-4)
expect_true(is.numeric(se_exp$se))

warpbreaks_fit <- glm(breaks ~ wool + tension,
                      data = warpbreaks,
                      family = poisson,
                      method = "brglmFit",
                      type = "ML")

focus_warpbreaks <- focus(warpbreaks_fit,
                          on = function(theta, k) theta[k]^2,
                          correction = "median",
                          k = 2)
se_warpbreaks <- focus_se(focus_warpbreaks, tol_opt = 1e-5)

expect_true(is.numeric(se_warpbreaks$se))
expect_equal(dim(se_warpbreaks$V), rep(length(se_warpbreaks$theta), 2))
expect_equal(length(se_warpbreaks$theta), length(coef(focus_warpbreaks$object, model = "mean")))
