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

focus_exp <- focus(coalition_fit,
                   on = function(theta) exp(theta[2]),
                   correction = "median")

expect_error(focus_se(focus_exp),
             pattern = "Could not reconstruct the model parameter vector")

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
