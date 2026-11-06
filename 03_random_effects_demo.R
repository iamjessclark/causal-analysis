# =============================================================================
# SHOULD VILLAGE AND HERD BE RANDOM EFFECTS, OR COVARIATES IN THE DAG?
# =============================================================================
#
# Short answer: both, and they are not alternatives.
#
#   The DAG decides WHICH variables to condition on. Village and herd cause
#     both the exposure and the outcome, so they are confounders and the DAG
#     says condition on them.
#   Fixed versus random is a separate choice about HOW to condition, made
#     after the DAG has spoken.
#
# The trap is assuming that `(1 | herd)` does the confounding job. It does not,
# because a random intercept assumes the cluster effect is uncorrelated with
# the exposure. When herds that graze on wetland also sit in high-vector areas,
# that assumption fails and a random intercept leaves bias behind.
#
# This script demonstrates the three things you need to know, by simulating
# data where the right answer is known. Run time is about a minute.
#
#   PART 1  A random intercept does not remove cluster-level confounding.
#   PART 2  The fix that keeps both properties (Mundlak / within-between).
#   PART 3  Random-effects odds ratios are not on the same scale as glm ones.
#
# =============================================================================

library(lme4)

inv_logit <- function(x) 1 / (1 + exp(-x))


# ---- PART 1 & 2: confounding ------------------------------------------------

# One herd-level variable U drives BOTH the exposure and the outcome, and is
# never observed. This is the realistic case: unrecorded vector pressure at the
# homestead pushes a herd towards wetland grazing and raises its infection risk
# independently.
#
# TRUE causal log-odds of the exposure = 0.70 throughout.

simulate_herds <- function(n_herd = 60, per_herd = 15, confounding = 1.2) {
  U    <- rnorm(n_herd)
  herd <- rep(seq_len(n_herd), each = per_herd)
  u    <- U[herd]
  X    <- rbinom(length(herd), 1, inv_logit(-0.2 + confounding * u))
  Y    <- rbinom(length(herd), 1, inv_logit(-1.0 + 0.70 * X + confounding * u))
  data.frame(herd = factor(herd), X = X, Y = Y)
}

# Four ways of handling the cluster.
estimates <- function(d) {
  d$X_herdmean <- ave(d$X, d$herd)      # the Mundlak term
  c(
    ignore_cluster = coef(glm(Y ~ X, family = binomial, data = d))["X"],
    fixed_effects  = coef(glm(Y ~ X + herd, family = binomial, data = d))["X"],
    random_int     = fixef(glmer(Y ~ X + (1 | herd),
                                 family = binomial, data = d))["X"],
    mundlak        = fixef(glmer(Y ~ X + X_herdmean + (1 | herd),
                                 family = binomial, data = d))["X"]
  )
}

set.seed(11)
reps <- replicate(60, estimates(simulate_herds()))

cat("\nTRUE causal log-odds = 0.70. Mean estimate over 60 simulated datasets:\n\n")
print(round(data.frame(mean_estimate = rowMeans(reps),
                       bias          = rowMeans(reps) - 0.70), 3))

# WHAT YOU SHOULD SEE
#
#   ignore_cluster  ~1.51   badly confounded, as expected
#   fixed_effects   ~0.78   nearly right
#   random_int      ~1.06   STILL BADLY CONFOUNDED
#   mundlak         ~0.73   nearly right
#
# The random intercept closed less than half the gap. It soaked up variation
# between herds but it did not break the link between U and the exposure,
# because fitting a random intercept assumes that link is absent. Adding
# `(1 | herd)` to a model and calling the result adjusted for herd is the
# mistake this part exists to prevent.
#
# THE MUNDLAK FIX
# Adding the herd mean of the exposure splits its coefficient into a
# within-herd part and a between-herd part. The within-herd part is what the
# fixed-effects model estimates, and the between-herd part absorbs the
# confounding. You keep the random intercept, so you also keep the partial
# pooling and the correct standard errors. Bell, Fairbrother & Jones (2019)
# call this the within-between random effects model.
#
# WHY NOT ALWAYS USE FIXED EFFECTS
# Two reasons, and the second is usually decisive in a serosurvey.
#   1. Logistic fixed effects are biased away from zero when clusters are
#      small, which is the ~0.78 against a truth of 0.70 above. The bias
#      shrinks as clusters grow; rerun with per_herd = 40 to watch it go.
#   2. A fixed effect absorbs any exposure that does not vary within the
#      cluster. If every animal in a herd shares one grazing system, herd
#      fixed effects make the grazing coefficient inestimable. Mundlak has the
#      same problem for a purely herd-level exposure, so for those you fall
#      back on measured herd covariates plus an E-value.


# ---- PART 3: the odds ratio changes meaning ---------------------------------

# Here the exposure is RANDOMISED, so there is no confounding anywhere and
# nothing to adjust for. Any difference between the two columns is not bias.

noncollapsibility <- function(sd_u) {
  n_herd <- 120; per_herd <- 20
  U    <- rnorm(n_herd, 0, sd_u)
  herd <- rep(seq_len(n_herd), each = per_herd)
  X    <- rbinom(length(herd), 1, 0.5)
  Y    <- rbinom(length(herd), 1, inv_logit(-1.0 + 0.70 * X + U[herd]))
  d    <- data.frame(herd = factor(herd), X = X, Y = Y)
  c(herd_SD     = sd_u,
    conditional = unname(fixef(glmer(Y ~ X + (1 | herd),
                                     family = binomial, data = d))["X"]),
    marginal    = unname(coef(glm(Y ~ X, family = binomial, data = d))["X"]))
}

set.seed(7)
cat("\n\nExposure randomised, so there is NO confounding. True log-odds = 0.70.\n\n")
print(round(as.data.frame(t(sapply(c(0.25, 1.0, 2.0, 3.0), noncollapsibility))), 3))

# WHAT YOU SHOULD SEE
# The conditional estimate stays near 0.70 whatever the herd variance. The
# marginal estimate shrinks towards zero as herd variance grows, reaching
# roughly 0.2 at a herd SD of 3.
#
# Neither is wrong. They answer different questions.
#
#   CONDITIONAL (what glmer and brms with (1|herd) report)
#     "Within a given herd, how much does this exposure raise the odds?"
#     This is the number to quote when you are describing a mechanism.
#
#   MARGINAL (what glm reports, and what a control programme needs)
#     "If every herd in the region changed, how much would seroprevalence
#      fall?" This is the counterfactual quantity in section 13 of the
#      tutorial.
#
# The odds ratio is non-collapsible, so these two differ even with no
# confounding at all. Never compare a coefficient from a random-effects model
# with one from a glm and call the difference confounding.
#
# HOW TO GET THE MARGINAL ANSWER FROM A RANDOM-EFFECTS MODEL
# Predict probabilities for every animal, then average them. That is exactly
# what pop_mean() does in section 13, and it marginalises correctly as long as
# you let the random effects vary rather than fixing them at zero:
#
#   posterior_epred(fit, newdata = d)                     # averages over the
#                                                         #   observed herds
#   posterior_epred(fit, newdata = d, re_formula = NA)    # sets every herd to
#                                                         #   zero: this is the
#                                                         #   MEDIAN herd, not
#                                                         #   the average one
#
# For a question about these herds, the first is right. For a question about
# herds in general, give the new data unseen herd labels with
# allow_new_levels = TRUE and sample_new_levels = "gaussian", which integrates
# over the fitted herd distribution.


# ---- WHAT THIS MEANS FOR THE RVF TUTORIAL -----------------------------------
#
# VILLAGE, 5 levels. Keep it a fixed effect, as sections 6 to 14 do.
#   Five clusters is far too few to estimate a variance from. Published
#   guidance asks for 30 to 50 groups before a random intercept is dependable,
#   and with five the variance and its standard error are unreliable. Village
#   is also a confounder you want fully adjusted, not partially pooled.
#
# HERD, around 50 levels. This is where a random effect earns its place.
#   Add (1 | herd_id) for the dependence between animals in a herd, and add
#   the herd mean of any exposure that varies within herds:
#
#     brm(seropositive ~ grazing + grazing_herdmean + village + age_class +
#                        sex + (1 | herd_id), family = bernoulli(), ...)
#
#   For a multi-level factor like grazing, the "herd mean" is the set of herd
#   proportions in each category, so in practice you add one column per
#   non-reference level.
#
# WHAT CHANGES, AND WHAT DOES NOT
#   The adjustment sets do not change. The DAG decides those, and nothing here
#   touches the DAG.
#   The standard errors change, and usually widen, because 900 animals in 50
#   herds carry less information than 900 independent animals.
#   The point estimate changes only if herd-level confounding was present.
#
# HOW THIS RELATES TO THE TUTORIAL DATA
#   rvf_cattle_serology.csv has real herd clustering, with a herd standard
#   deviation of 0.63 on the log-odds scale after village, an intraclass
#   correlation of 0.11. What it does not have is strong herd-level
#   CONFOUNDING, so there the random intercept, the Mundlak model and herd
#   fixed effects all land close to the village-only answer. The parameters in
#   this script are turned up so the difference is visible. Use this script to
#   see the mechanism, and the CSV to practise the analysis.
#
# REFERENCES
#   Bell, Fairbrother & Jones 2019, Qual Quant 53:1051, on within-between
#     random effects and the Mundlak specification.
#   Mundlak 1978, Econometrica 46:69, the original device.
#   Hubbard et al. 2010, Epidemiology 21:467, "To GEE or not to GEE", on
#     marginal versus conditional estimands in epidemiology.
#   Greenland, Robins & Pearl 1999, Statist Sci 14:29, on non-collapsibility.
#
# =============================================================================
