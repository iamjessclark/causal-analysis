# =============================================================================
# CAUSAL ANALYSIS IN R: RVFV SEROLOGY IN CATTLE
#
# The code for the tutorial. The explanation of every section, the function
# reference and how to read each output are in TUTORIAL.md. Read a section
# there, then run the matching section here.
#
# Needs: dagitty, ggdag, tidyverse, brms, fastDummies, CCP, plus a C++
#   toolchain for brms (Xcode command line tools on macOS, Rtools on Windows).
# Data: rvf_cattle_serology.csv, written by 01_simulate_data.R.
# Run time: the first run compiles around a dozen Stan programs and takes 10 to
#   20 minutes. Fits cache in models/, so re-runs are instant. Delete models/
#   if you change the data or a model formula.
#
# Follows the analysis structure of Lim et al. (2026), BMJ Public Health
# 4:e004500, doi:10.1136/bmjph-2025-004500.
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Install missing packages and load the needed ones. Run once, then comment this block out.
needed <- c("dagitty", "ggdag", "tidyverse", "brms", "fastDummies", "CCP")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)

library(dagitty)    # build DAGs, find adjustment sets, test the DAG
library(ggdag)      # draw DAGs
library(tidyverse)  # ggplot2 for every figure, dplyr for mutate() and bind_rows()
library(brms) # what you are using to fit your models

# fastDummies and CCP are installed but never loaded, because dagitty calls
# them inside localTests() in section 5. The tidyverse conflict list printed on
# loading is expected. TUTORIAL.md, section 0.

dir.create("figures", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)   # caches the Bayesian fits

set.seed(1)


# ---- 1. Read the data and set reference levels ------------------------------
# TUTORIAL.md, section 1, for the column list and what the summaries show.

# you can either run
source("01_simulate_data.R") # which creates rvf_cattle_serology.csv
# or if you have already run the 01 script and have the CSV file saved in this
# working directory then you can load the csv directly and # out line 75.

cattle <- read.csv("rvf_cattle_serology.csv", stringsAsFactors = TRUE)

str(cattle)

# R orders factor levels alphabetically, so set the reference to the
# lowest-risk or baseline level yourself.
cattle <- cattle |>
  mutate(
    grazing        = relevel(grazing, ref = "tethered_homestead"),
    water_source   = relevel(water_source, ref = "carried_to_animal"),
    age_class     = factor(age_class, levels = c("calf", "juvenile", "adult")),
    sex            = relevel(sex, ref = "female")
  )

# Describe the data before modelling it. 
# overall seroprevalence
paste0(round(mean(cattle$seropositive) * 100, 1), "%")

# 38 villages is too many to read one by one, so summarise the spread between
# them instead of printing all 38.
village_seroprev <- tapply(cattle$seropositive, cattle$village, mean) * 100
# summaries quantiles
round(summary(village_seroprev), 1)
# lowest and highest
round(sort(village_seroprev)[c(1, 38)], 1)
# grazing habits
round(tapply(cattle$seropositive, cattle$grazing, mean) * 100, 1)
# age structure
round(tapply(cattle$seropositive, cattle$age_class, mean) * 100, 1)
# grazing ~ water access
table(cattle$grazing, cattle$water_source)


# ---- 2. How every model in this tutorial is fitted --------------------------
# TUTORIAL.md, section 2, for the choice of priors and the function reference.

# Weakly informative priors on the log-odds scale, centred on no effect.
pr_bin <- c(prior(normal(0, 1.5), class = "b"),
            prior(normal(0, 2),   class = "Intercept"))

# One wrapper for the sampler settings used everywhere below. file= caches the
# fit, so the compiling only happens once.
fit_bin <- function(formula, name) {
  brm(formula, family = bernoulli(link = "logit"), data = cattle,
      prior = pr_bin,
      chains = 4, iter = 3000, warmup = 1000, cores = 4,
      seed = 1, refresh = 0,
      file = file.path("models", name))
}

# A function to produce log-odds coefficients as odds ratios with a 95% credible interval. 
# keep= is a regular expression matched against the coefficient names.
or_table <- function(fit, keep = NULL) {
  fx  <- fixef(fit)
  out <- data.frame(term = rownames(fx),
                    OR = exp(fx[, "Estimate"]),
                    lo = exp(fx[, "Q2.5"]),
                    hi = exp(fx[, "Q97.5"]))
  if (!is.null(keep)) out <- out[grepl(keep, out$term), ]
  rownames(out) <- NULL
  out[, 2:4] <- round(out[, 2:4], 2)
  out
}

# The proportion of posterior draws above zero, which is the probability that
# the effect raises the odds (another custom-made function).
p_positive <- function(fit, keep) {
  draws <- as.data.frame(as_draws_df(fit))
  round(sapply(draws[, grep(keep, names(draws)), drop = FALSE],
               function(x) mean(x > 0)), 3)
}


# ---- 3. Why a serostatus ~ grazing regression is not sufficient --------------
# TUTORIAL.md, section 3.

# Grazing and nothing else, so every later model can be read as a change from
# this starting point. That is why this is "crude" model. Note using your 
# custom model that you defined above. 
b_crude <- fit_bin(seropositive ~ grazing, "b_crude")
# generate the odds ratio table
or_table(b_crude, keep = "grazing")


# ---- 4. Draw the DAG --------------------------------------------------------
# TUTORIAL.md, section 4, for the syntax, how to decide each node's role and
# how to decide its position.

# "A -> B" means A causes B. Square brackets set the role of a node and its
# position in the drawing. Names must match the column names in your data
# exactly, otherwise the DAG cannot be tested against the data later.

rvf_dag <- dagitty('dag {
  village                [pos="0.0,1.0"]
  age_class              [pos="0.0,2.1"]
  sex                    [pos="0.0,3.0"]
  grazing                [exposure, pos="1.6,2.0"]
  water_source           [pos="3.1,0.3"]
  hours_at_water         [pos="4.4,0.9"]
  wildlife_contact       [pos="3.1,1.7"]
  introduced_recently    [pos="1.6,4.3"]
  seropositive           [outcome, pos="5.8,2.3"]

  village -> grazing
  village -> water_source
  village -> wildlife_contact
  village -> introduced_recently
  village -> seropositive

  age_class -> grazing
  age_class -> seropositive

  sex -> grazing
  sex -> seropositive

  grazing -> water_source
  grazing -> hours_at_water
  grazing -> wildlife_contact
  grazing -> seropositive

  water_source -> hours_at_water
  water_source -> seropositive

  hours_at_water      -> seropositive
  wildlife_contact    -> seropositive
  introduced_recently -> seropositive
}')

# tidy_dagitty() turns the DAG into a data frame of nodes and edges, and
# node_status() flags the exposure and the outcome so the plot can colour them.
p_dag <- rvf_dag |>
  tidy_dagitty() |>
  node_status() |>
  mutate(status = ifelse(is.na(status), "other", as.character(status))) |>
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey60") +
  geom_dag_point(aes(colour = status), size = 5) +
  geom_dag_text(aes(label = name), colour = "grey15", size = 3,
                nudge_y = 0.2) +
  scale_colour_manual(name = NULL,
                      values = c(exposure = "#C0392B", outcome = "#1F618D",
                                 other = "grey40"),
                      breaks = c("exposure", "outcome")) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  theme_dag() +
  theme(legend.position = "bottom") +
  labs(title = "Hypothesised causal structure for RVFV seropositivity in cattle",
       subtitle = "Arrows point from cause to effect")
ggsave("figures/01_dag.png", p_dag, width = 11, height = 7, dpi = 200)
p_dag


# ---- 5. Test the DAG against the data ---------------------------------------
# TUTORIAL.md, section 5, including why a 38-level factor inflates the estimate.

# Every pair of variables with no arrow between them should be independent once
# you condition on the right other variables.
ci <- impliedConditionalIndependencies(rvf_dag)
length(ci)       # how many claims this DAG makes
ci[1:6]          # "A _||_ B | C" reads "A is independent of B given C"

# localTests() checks each claim in the data. "cis.pillai" is the version that
# copes with a mix of categorical and continuous variables.
dag_vars <- names(cattle)[names(cattle) %in% names(rvf_dag)]
tests <- localTests(rvf_dag, cattle[, dag_vars], type = "cis.pillai")
tests <- tests[order(-abs(tests$estimate)), ]
round(head(tests, 10), 3)

sum(abs(tests$estimate) > 0.15)     # how many look like clear violations
sum(tests$p.value < 0.05)           # how many significant, out of nrow(tests)


# ---- 6. Paths, confounders, mediators and colliders -------------------------
# TUTORIAL.md, section 6, for how to read the path list. The three path shapes
# are in section 4.

paths(rvf_dag, from = "grazing", to = "seropositive")$paths |> head(12)
length(paths(rvf_dag, from = "grazing", to = "seropositive")$paths)

# Backdoor paths are the ones that leave the exposure through an arrow pointing
# INTO it, so grepl() picks them out by the arrow direction.
backdoors <- paths(rvf_dag, from = "grazing", to = "seropositive")
length(backdoors$paths[grepl("^grazing <-", backdoors$paths)])
head(backdoors$paths[grepl("^grazing <-", backdoors$paths)], 6)


# ---- 7. Let the DAG choose your covariates ----------------------------------
# TUTORIAL.md, section 7. Record these sets in your methods.

# TOTAL effect: everything grazing does to seropositivity, by any route.
adjustmentSets(rvf_dag, exposure = "grazing", outcome = "seropositive",
               effect = "total")

# DIRECT effect: what is left when you hold the mediators fixed.
adjustmentSets(rvf_dag, exposure = "grazing", outcome = "seropositive",
               effect = "direct")

# The same question for each of the other exposures. One model per exposure,
# not one model per paper.
adjustmentSets(rvf_dag, "water_source", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "hours_at_water", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "wildlife_contact", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "introduced_recently", "seropositive", effect = "total")


# ---- 8. Estimate the total effect -------------------------------------------
# TUTORIAL.md, section 8.

# The exposure plus the adjustment set from section 7, and nothing more.
b_total <- fit_bin(seropositive ~ grazing + village + age_class + sex, "b_total")

# CHECK CONVERGENCE BEFORE READING ANY RESULT. Rhat above 1.01 means the chains
# have not settled on the same answer. Bulk_ESS and Tail_ESS above a few
# hundred are fine.
summary(b_total)
max(rhat(b_total), na.rm = TRUE)

# Posterior predictive check. The outcome is 0/1, so type = "bars": pale bars
# are the observed counts of 0 and 1, dark points are what the model predicts.
p_pp <- pp_check(b_total, type = "bars", ndraws = 500) +
  labs(title = "Posterior predictive check",
       subtitle = "Observed counts of 0 and 1 against what the model predicts")
ggsave("figures/02_ppcheck.png", p_pp, width = 6, height = 4, dpi = 200)
p_pp

# Now read the effect, as odds ratios and as posterior probabilities.
or_table(b_total, keep = "grazing")
p_positive(b_total, "^b_grazing")

# Posterior intervals on the log-odds scale.
# mcmc_plot() fits the x axis to the intervals, which here all lie above zero,
# so the scale_x_continuous() line widens the axis to include the no-effect
# line at 0. suppressMessages() hides ggplot's note that the axis was replaced.
p_eff <- suppressMessages(
  mcmc_plot(b_total, variable = "^b_grazing", regex = TRUE,
            prob = 0.66, prob_outer = 0.95) +
    scale_x_continuous(limits = function(l) range(c(l, 0)))) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "Total effect of grazing system on RVFV seropositivity",
       subtitle = "Posterior log-odds vs tethered cattle; thick 66%, thin 95% credible interval",
       x = "log-odds")
ggsave("figures/03_posterior_grazing.png", p_eff, width = 7, height = 3.5, dpi = 200)
p_eff


# ---- 9. Estimate the direct effect, and read the mediation ------------------
# TUTORIAL.md, section 9. The gap between the two is partial mediation.

b_direct <- fit_bin(seropositive ~ grazing + village + age_class + sex +
                      water_source + hours_at_water + wildlife_contact,
                    "b_direct")

or_table(b_direct, keep = "grazing")

# Put the two side by side.
compare <- merge(or_table(b_total, "grazing"), or_table(b_direct, "grazing"),
                 by = "term", suffixes = c("_total", "_direct"))
compare


# ---- 10. Alternative model specifications -----------------------------------
# TUTORIAL.md, section 10.

# The same outcome fitted with four sets of covariates.

b_all_variables <- fit_bin(seropositive ~ grazing + village + age_class + sex +
                             water_source + hours_at_water + wildlife_contact +
                             introduced_recently,
                           "b_all_variables")

truth <- read.csv("simulation_truth.csv")

wetland_or <- function(fit) {
  round(exp(fixef(fit)["grazingwetland_grazing", "Estimate"]), 2)
}

specifications <- data.frame(
  approach = c("no adjustment", "DAG total effect", "DAG direct effect",
               "everything in the model"),
  wetland_OR = c(wetland_or(b_crude), wetland_or(b_total), wetland_or(b_direct),
                 wetland_or(b_all_variables)),
  answers = c("nothing: confounded by village and age",
              "the total causal effect",
              "the effect not running through the mediators",
              "the direct effect, plus one covariate the DAG leaves out")
)
specifications

# The simulated direct effect that generated these data, for reference:
truth[truth$term == "grazing_wetland_grazing", ]
exp(truth$true_log_odds[truth$term == "grazing_wetland_grazing"])


# ---- 11. A continuous exposure: hours at surface water ----------------------
# TUTORIAL.md, section 11.

# Its adjustment set from section 7 is {grazing, water_source}.
b_hours <- fit_bin(seropositive ~ hours_at_water + grazing + water_source,
                   "b_hours")
or_table(b_hours, keep = "hours_at_water")

# A per-hour OR is hard to interpret, so convert it to the contrasts people
# care about. The effect of k hours is exp(beta * k), not k times the OR.
beta_hours <- fixef(b_hours)["hours_at_water", "Estimate"]
data.frame(extra_hours = c(1, 2, 4),
           OR = round(exp(beta_hours * c(1, 2, 4)), 2))


# ---- 12. Counterfactual interventions: what if we changed practice? ---------
# TUTORIAL.md, section 12, for the recipe and why each step is as it is.

# posterior_epred() returns a draws-by-animals matrix of predicted
# probabilities. Averaging across animals within each draw gives a posterior
# for the population seroprevalence.
pop_mean <- function(fit, newdata, rows = NULL) {
  pp <- posterior_epred(fit, newdata = newdata)
  if (!is.null(rows)) pp <- pp[, rows, drop = FALSE]
  rowMeans(pp)
}

q <- function(x) round(c(mean = mean(x), quantile(x, c(0.025, 0.975))) * 100, 1)

# --- Scenario A: move all wetland grazers to tethered homestead grazing ------
cf_A <- cattle
cf_A$grazing[cf_A$grazing == "wetland_grazing"] <- "tethered_homestead"
wet <- which(cattle$grazing == "wetland_grazing")

obs_all <- pop_mean(b_total, cattle)
cfa_all <- pop_mean(b_total, cf_A)
obs_wet <- pop_mean(b_total, cattle, wet)
cfa_wet <- pop_mean(b_total, cf_A,   wet)

rbind(`whole population, observed`      = q(obs_all),
      `whole population, no wetland`    = q(cfa_all),
      `wetland grazers, observed`       = q(obs_wet),
      `wetland grazers, now tethered`   = q(cfa_wet))

# Report the difference as a distribution, not a subtraction of two means.
q(obs_all - cfa_all)    # percentage points avoided, whole population
q(obs_wet - cfa_wet)    # percentage points avoided, among those changed

# --- Scenario B: remove surface-water contact, do(hours_at_water = 0) -------
cf_B <- cattle
cf_B$hours_at_water <- 0
contact <- which(cattle$hours_at_water > 0)

obsB_all <- pop_mean(b_hours, cattle);  cfB_all <- pop_mean(b_hours, cf_B)
obsB_sub <- pop_mean(b_hours, cattle, contact); cfB_sub <- pop_mean(b_hours, cf_B, contact)

rbind(`whole population, observed`   = q(obsB_all),
      `whole population, no contact` = q(cfB_all),
      `animals with contact, observed`  = q(obsB_sub),
      `animals with contact, set to 0`  = q(cfB_sub))
q(obsB_sub - cfB_sub)

# --- Scenario C: intervene on grazing AND its downstream consequences -------
# Editing the mediators explicitly needs a model that contains them, which is
# the direct-effect model. Scenario A did not, because the total-effect model
# already carries the indirect routes inside its grazing coefficient.
cf_C <- cattle
sw <- cattle$grazing == "wetland_grazing"
cf_C$grazing[sw]          <- "tethered_homestead"
cf_C$water_source[sw]     <- "carried_to_animal"
cf_C$hours_at_water[sw]   <- 0
cf_C$wildlife_contact[sw] <- 0

cfC_wet <- pop_mean(b_direct, cf_C, wet)
obsC_wet <- pop_mean(b_direct, cattle, wet)
rbind(`wetland grazers, observed`               = q(obsC_wet),
      `tethered + mediators changed (direct)`   = q(cfC_wet),
      `tethered only (total-effect model)`      = q(cfa_wet))

# --- Scenario D: does the benefit differ by grazing system? -----------------
# A conditional average treatment effect. One average effect can hide a large
# effect in one subgroup and none in another.
cate <- lapply(levels(cattle$grazing), function(g) {
  rows <- which(cattle$grazing == g & cattle$hours_at_water > 0)
  if (length(rows) < 10) return(NULL)
  d <- pop_mean(b_hours, cattle, rows) - pop_mean(b_hours, cf_B, rows)
  data.frame(grazing = g, n = length(rows),
             median_pp = round(median(d) * 100, 1),
             lo = round(quantile(d, 0.025) * 100, 1),
             hi = round(quantile(d, 0.975) * 100, 1))
})
cate <- bind_rows(cate)
rownames(cate) <- NULL
cate

p_cate <- ggplot(cate, aes(reorder(grazing, median_pp), median_pp)) +
  geom_pointrange(aes(ymin = lo, ymax = hi)) +
  coord_flip() +
  labs(x = NULL, y = "Percentage points of seroprevalence avoided",
       title = "Effect of removing surface-water contact, by grazing system",
       subtitle = "Median and 95% credible interval; animals with contact only") +
  theme_minimal()
ggsave("figures/04_cate.png", p_cate, width = 7, height = 3.5, dpi = 200)
p_cate


# ---- 13. Sensitivity analysis: how wrong could this be? ---------------------
# TUTORIAL.md, section 13, including the three assumptions no test covers.

# (a) UNMEASURED CONFOUNDING: the E-value.
# The minimum association, on the risk-ratio scale, that a hidden confounder
# would need with BOTH exposure and outcome to explain away the effect.
e_value <- function(or, p0) {
  rr <- or / (1 - p0 + p0 * or)         # OR to RR, given baseline risk p0
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}

p0 <- mean(cattle$seropositive[cattle$grazing == "tethered_homestead"])
wet_fx   <- fixef(b_total)["grazingwetland_grazing", ]
or_pt    <- exp(wet_fx["Estimate"])
ci_total <- exp(wet_fx[c("Q2.5", "Q97.5")])

data.frame(baseline_risk = round(p0, 3),
           OR = round(or_pt, 2),
           E_point = round(e_value(or_pt, p0), 2),
           E_lower_CI = round(e_value(ci_total[1], p0), 2))

# The DAG you can draw and the DAG you believe are not the same thing. Mark the
# hidden variable as latent and ask dagitty whether the effect is still
# identifiable from the data you have. An empty result means it is not.
rvf_dag_u <- dagitty('dag {
  U [latent]
  U -> grazing
  U -> water_source
  U -> hours_at_water
  U -> seropositive
  village -> grazing
  village -> seropositive
  age_class -> grazing
  age_class -> seropositive
  sex -> grazing
  sex -> seropositive
  grazing -> water_source
  grazing -> hours_at_water
  grazing -> seropositive
  water_source -> hours_at_water
  water_source -> seropositive
  hours_at_water -> seropositive
  grazing [exposure]
  seropositive [outcome]
}')
adjustmentSets(rvf_dag_u, effect = "total")

# (b) WHAT IF A MEDIATOR IS REALLY A CONFOUNDER?
# You cannot settle this from the data, so show the estimate under every
# combination and let the reader judge.
mediators <- c("water_source", "hours_at_water", "wildlife_contact")
combos <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), 3))
names(combos) <- mediators

grid <- lapply(seq_len(nrow(combos)), function(i) {
  extra <- mediators[unlist(combos[i, ])]
  f <- reformulate(c("grazing", "village", "age_class", "sex", extra),
                   response = "seropositive")
  fit <- fit_bin(f, paste0("b_grid_", i))
  est <- exp(fixef(fit)["grazingwetland_grazing", c("Estimate", "Q2.5", "Q97.5")])
  data.frame(conditioned_on = if (length(extra)) paste(extra, collapse = " + ") else "none (total effect)",
             OR = round(est[1], 2), lo = round(est[2], 2), hi = round(est[3], 2))
})
grid <- bind_rows(grid)
rownames(grid) <- NULL
grid[order(-grid$OR), ]

# (c) MEASUREMENT ERROR IN A SELF-REPORTED EXPOSURE.
# If farmers report only a fraction f of the true hours at water, the fitted
# coefficient under-states the true one by the same fraction.
frac <- c(1.0, 0.8, 0.6, 0.5)
data.frame(reported_fraction = frac,
           OR_per_hour = round(exp(beta_hours / frac), 2),
           OR_over_4h  = round(exp(4 * beta_hours / frac), 2))


# ---- 14. How to write this up -----------------------------------------------
# TUTORIAL.md, section 14, for the methods and results checklists and the
# further reading.

# =============================================================================
