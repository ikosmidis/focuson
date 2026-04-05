# Copyright (C) 2026- Ioannis Kosmidis

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#'
#' focuson: Inference on scalar functions of model parameters.
#'
#' Methods for estimation and inference on user-specified scalar
#' parameters of interest from fitted model objects through [focus()].
#' Estimation can be performed using mean or median bias correction, and
#' inference is based on Wald-type confidence intervals using
#' delta-method standard errors or HulC confidence intervals.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso
#'
#' [focus()]
#'
#' @references
#'
#' Kenne Pagui E C, Salvan A, Sartori N (2017). Median bias
#' reduction of maximum likelihood estimates. *Biometrika*, **104**,
#' 923–938. \doi{10.1093/biomet/asx046}.
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59. \doi{10.1007/s11222-019-09860-6}.
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196. \doi{10.1002/wics.1296}.
#'
#' Kuchibhotla A K, Balakrishnan S, Wasserman L (2024). HulC: high
#' confidence level upper and lower confidence bounds. *Journal of the
#' Royal Statistical Society Series B: Statistical Methodology*, **86**,
#' 586-622. \doi{10.1093/jrsssb/qkad134}.
#'
#' @docType package
#' @aliases focuson-package
#' @name focuson
#' @import brglm2
#' @import enrichwith
#' @importFrom numDeriv grad hessian
#' @importFrom stats coef model.frame qnorm runif update vcov printCoefmat
#' @importFrom utils capture.output
#'
"_PACKAGE"
