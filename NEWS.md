# focuson 0.3

## New functionality

* Added `modified_profile_focus()` for simulation-based construction of the
  modified profile likelihood and the modified signed likelihood-ratio
  statistic $r^*$ for arbitrary scalar functions of model parameters. The
  returned profile includes the nuisance-parameter and information
  adjustments.

* Extended `confint()` for computed focus profiles with
  `method = "pl"`, `"mpl"`, or `"rstar"` for intervals based on the ordinary
  profile likelihood, modified profile likelihood, or modified signed
  likelihood-ratio statistic $r^*$, respectively.

# focuson 0.2

## Bug fixes

* Improved recovery of the original model data for HulC confidence intervals,
  including binomial models with matrix responses.

* Fixed GLM refitting in examples and package use when `brglm2` is imported but
  not attached.

## New functionality

* Added `profile_ci()` for direct computation of profile likelihood
  confidence intervals for arbitrary scalar functions of model parameters
  using the Venzon--Moolgavkar endpoint equations.

* Added `profile_focus()` for constructing complete profile likelihoods for
  arbitrary scalar functions, either on a likelihood-root grid or over a
  user-supplied range of focus values.

* Added profile likelihood support for GLM-based `focus()` results through
  `profile()` and `confint(..., method = "pl")`.

* Added `print()`, `plot()`, and `confint()` methods for computed focus
  profiles, with linear or cubic interpolation.

## Improvements, updates and additions

* Consolidated the common estimation calculations used by `focus()`
  and `focus_engine()`, improving consistency between fitted-model and
  low-level interfaces.

* Standardized the interfaces of `estimate_focus_components()`,
  `estimate_focus_components_fef()`, and
  `estimate_focus_components_iid()`. Additional likelihood and simulation
  arguments are now supplied through `likelihood_args` instead of `...`.
  The separate `data` argument of `estimate_focus_components_fef()` has been
  removed; the observed dataset is now supplied as the named `data` element
  of `likelihood_args`.

* Expanded examples, tests, and documentation for likelihood profiling and
  confidence-interval workflows.

# focuson 0.1

* First public release, along with <https://arxiv.org/abs/2606.28597v1>.
