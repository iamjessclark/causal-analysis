# =============================================================================
# CAUSAL ANALYSIS IN R: A BEGINNER'S TUTORIAL
# Worked example: Rift Valley fever virus (RVFV) serology in cattle
# =============================================================================
#
# WHAT THIS TUTORIAL IS FOR
#
# A standard regression tells you which variables move together. It does not
# tell you which variable causes which, and it does not tell you what would
# happen if you changed one of them. Causal analysis adds two things on top of
# regression: a written-down assumption about how the variables are connected
# (a DAG), and a rule that uses that assumption to decide which variables to
# put in the model. Get the DAG right and the regression coefficient answers a
# causal question. Get it wrong, or skip it, and the coefficient answers no
# clear question at all.
#
# The worked example follows the analysis structure of Lim et al. (2026)
# "Local habitual movement as a mechanism for Schistosoma mansoni transmission
# resurgence: a causal analysis of cross-sectional survey data",
# BMJ Public Health 4:e004500, doi:10.1136/bmjph-2025-004500, applied to a
# simulated cattle RVFV serosurvey instead of a human schistosomiasis survey.
#
# THE RESEARCH QUESTION
#
# Does the way cattle are grazed and watered cause RVFV infection, and by how
# much would seroprevalence fall if those practices changed?
#
# HOW TO USE THIS FILE
#
# Run it top to bottom, one block at a time, and read the comment under each
# output before moving on. Sections 1 to 11 use ordinary glm() and run in
# seconds. Sections 12 to 14 use Bayesian models, which take a few minutes the
# first time because the Stan code has to compile, and are instant afterwards
# because the fits are cached in models/.
#
# The whole causal logic is in sections 1 to 11. If brms gives you installation
# trouble, stop at section 11 and you will still have the method.
#
# WHAT YOU NEED FIRST
#
#   source("01_simulate_data.R")   # creates rvf_cattle_serology.csv
#
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Install anything missing. Run once, then comment this block out.
needed <- c("dagitty", "ggdag", "ggplot2", "dplyr", "tidyr",
            "fastDummies", "CCP", "brms", "patchwork")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)

library(dagitty)   # build DAGs, find adjustment sets, test the DAG
library(ggdag)     # draw DAGs
library(ggplot2)
library(dplyr)

dir.create("figures", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)   # caches the Bayesian fits

set.seed(1)


# ---- 1. Read the data and set reference levels ------------------------------

cattle <- read.csv("rvf_cattle_serology.csv", stringsAsFactors = TRUE)

# A single cross-sectional survey: 935 cattle from 52 herds in 5 villages, with
# one blood sample and one management questionnaire per animal.
#
#   seropositive            1 if antibodies to natural RVFV infection detected
#   titre_sp                antibody level (S/P %), recorded only if seropositive
#   grazing                 how the animal is grazed
#   water_source            where its drinking water comes from
#   hours_at_water          hours per day in or at the edge of surface water
#   wildlife_contact        shares grazing or watering points with wildlife
#   night_housing           kraaled or housed overnight
#   vector_control          acaricide / insecticide use
#   introduced_recently     bought into the herd in the last 12 months
#   abortion_reported       abortion or stillbirth recorded for this animal
#   village, age_class, sex, herd_size, herd_id

str(cattle)

# R orders factor levels alphabetically, so the reference category in every
# model would be whichever level happens to come first in the alphabet. Set the
# reference to the lowest-risk or baseline level yourself, otherwise every
# coefficient is a contrast against something arbitrary.
cattle <- cattle |>
  mutate(
    grazing        = relevel(grazing, ref = "tethered_homestead"),
    water_source   = relevel(water_source, ref = "carried_to_animal"),
    vector_control = relevel(vector_control, ref = "none"),
    age_class      = factor(age_class, levels = c("calf", "juvenile", "adult")),
    sex            = relevel(sex, ref = "female")
  )

# Describe the data before modelling it. Anything surprising here is usually a
# data problem, not a finding.
round(mean(cattle$seropositive) * 100, 1)                 # overall seroprevalence
round(tapply(cattle$seropositive, cattle$village, mean) * 100, 1)
round(tapply(cattle$seropositive, cattle$grazing, mean) * 100, 1)
round(tapply(cattle$seropositive, cattle$age_class, mean) * 100, 1)
table(cattle$grazing, cattle$water_source)

# INTERPRETATION
# Seroprevalence is about 36% overall, ranging from 8% in village A to 51% in
# village D, and it rises from about 17% in tethered cattle to about 72% in
# wetland grazers. It also rises with age, from 19% in calves to 45% in adults,
# which is expected for a cumulative marker of past exposure.
#
# The cross-tabulation is far from balanced, because wetland grazers almost
# never have water carried to them. Exposures that travel together like this
# cannot be separated by the data alone, which is the first reason you need a
# causal diagram rather than a kitchen-sink regression.


# ---- 2. Why the crude comparison is not the answer --------------------------

crude <- glm(seropositive ~ grazing, family = binomial, data = cattle)
round(exp(cbind(OR = coef(crude), confint.default(crude))), 2)

# INTERPRETATION
# Wetland grazing comes out at roughly 13 times the odds of seropositivity
# compared with tethered cattle. That number is a description of the survey,
# not an effect of grazing. Two other things differ between those groups.
# Village: wetland grazing happens where there is wetland, and wetland is also
# where the mosquito vectors breed. Age: calves are tethered near the homestead
# and have had less time to be exposed. Both of those push the comparison up on
# their own. The rest of this tutorial is about separating the part caused by
# grazing from the part caused by everything that travels with it.


# ---- 3. Draw the DAG --------------------------------------------------------

# A DAG (directed acyclic graph) is a drawing of what you believe causes what.
# Three rules:
#   1. One node per variable.
#   2. One arrow per direct cause, pointing from cause to effect.
#   3. No loops. If A causes B, nothing downstream of B can cause A.
#
# The arrows come from subject knowledge, published studies and expert opinion,
# never from the data. Missing an arrow that really exists is an assumption that
# can bias your estimate, so it is safer to include a plausible arrow than to
# leave it out.
#
# The arrows below are drawn from the RVF literature: wetland and riverine
# grazing and shared watering points raise exposure to Aedes and Culex vectors;
# cattle not kraaled at night are exposed to night-biting vectors; animals
# bought in to replenish a herd can import infection; acaricide and insecticide
# use lowers vector contact; shared grazing with wildlife raises it.
#
# Syntax: "A -> B" means A causes B. Square brackets set the role of a node and
# its position in the drawing. Names must match the column names in your data
# exactly, otherwise the DAG cannot be tested against the data later.

rvf_dag <- dagitty('dag {
  village                [pos="0.0,1.0"]
  age_class              [pos="0.0,2.1"]
  sex                    [pos="0.0,3.0"]
  grazing                [exposure, pos="1.6,2.0"]
  water_source           [pos="3.1,0.3"]
  hours_at_water         [pos="4.4,0.9"]
  wildlife_contact       [pos="3.1,1.7"]
  night_housing          [pos="3.1,2.6"]
  vector_control         [pos="3.1,3.5"]
  introduced_recently    [pos="1.6,4.3"]
  seropositive           [outcome, pos="5.8,2.3"]
  abortion_reported [pos="7.2,2.3"]

  village -> grazing
  village -> water_source
  village -> wildlife_contact
  village -> night_housing
  village -> vector_control
  village -> introduced_recently
  village -> seropositive

  age_class -> grazing
  age_class -> seropositive

  sex -> grazing
  sex -> seropositive

  grazing -> water_source
  grazing -> hours_at_water
  grazing -> wildlife_contact
  grazing -> night_housing
  grazing -> seropositive
  grazing -> abortion_reported

  water_source -> hours_at_water
  water_source -> seropositive

  hours_at_water      -> seropositive
  wildlife_contact    -> seropositive
  night_housing       -> seropositive
  vector_control      -> seropositive
  introduced_recently -> seropositive

  seropositive -> abortion_reported
}')

# tidy_dagitty() turns the DAG into a data frame of nodes and edges, and
# node_status() flags which node is the exposure and which is the outcome, so
# the plot can colour them. Labels go beside the points rather than inside
# them, because variable names are too long to fit in a circle.
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

# INTERPRETATION
# Read the drawing and check each arrow against what you know about the system.
# Four features of this DAG decide everything that follows.
#   village, age_class and sex sit upstream of both grazing and seropositivity,
#     which makes them confounders of the grazing effect.
#   water_source, hours_at_water, wildlife_contact and night_housing sit
#     between grazing and seropositivity, which makes them mediators. They are
#     part of how grazing causes infection, not a nuisance to be removed.
#   abortion_reported sits downstream of seropositivity, because infection
#     causes abortion rather than the reverse.
#   hours_at_water has no arrow into grazing, because the amount of time an
#     animal spends at water is a consequence of how it is grazed and watered.


# ---- 4. Test the DAG against the data ---------------------------------------

# A DAG makes checkable claims. Every pair of variables with no arrow between
# them should be independent once you condition on the right other variables.
# These claims are the implied conditional independencies.

ci <- impliedConditionalIndependencies(rvf_dag)
length(ci)       # how many claims this DAG makes
ci[1:6]          # "A _||_ B | C" reads "A is independent of B given C"

# localTests() checks each claim in the data. "cis.pillai" is the version that
# copes with a mix of categorical and continuous variables.
dag_vars <- names(cattle)[names(cattle) %in% names(rvf_dag)]
tests <- localTests(rvf_dag, cattle[, dag_vars], type = "cis.pillai")
tests <- tests[order(-abs(tests$estimate)), ]
round(head(tests, 10), 3)

sum(abs(tests$estimate) > 0.15)     # how many clear violations
sum(tests$p.value < 0.05)           # how many significant, out of nrow(tests)

# INTERPRETATION
# "estimate" is the leftover association after conditioning, on a scale where 0
# means the claim holds. Nothing here exceeds 0.11, and 3 of the 47 tests reach
# p < 0.05, which is close to the 2 you would expect from chance at the 5% level.
# The DAG is compatible with the data.
#
# A violation above roughly 0.15 means an arrow you left out probably belongs
# in. Add it, re-test, and say in your methods that you did. What this test
# cannot do is confirm the DAG, because several different DAGs imply exactly the
# same independencies. It only rules out DAGs that the data clearly contradict.


# ---- 5. Paths, confounders, mediators and colliders -------------------------

# Every connection between exposure and outcome runs along a path. Some paths
# carry the causal effect you want. Others carry bias. The job is to block the
# bias paths and leave the causal paths alone.

paths(rvf_dag, from = "grazing", to = "seropositive")$paths |> head(12)
length(paths(rvf_dag, from = "grazing", to = "seropositive")$paths)

# Three path shapes, and what to do with each:
#
#   CHAIN      grazing -> hours_at_water -> seropositive
#              The middle variable is a MEDIATOR. Conditioning on it removes
#              part of the effect you are trying to measure. Leave it out when
#              you want the total effect.
#
#   FORK       grazing <- village -> seropositive
#              The middle variable is a CONFOUNDER. It creates association
#              without causation. Conditioning on it is required.
#
#   COLLIDER   grazing -> abortion_reported <- seropositive
#              Two arrows meet head-on. These variables are already
#              independent through this path, and conditioning on the collider
#              CREATES a false association. Never condition on it.
#
# Confounders are found by looking for open backdoor paths: paths that leave
# the exposure through an arrow pointing INTO it.
backdoors <- paths(rvf_dag, from = "grazing", to = "seropositive")
length(backdoors$paths[grepl("^grazing <-", backdoors$paths)])
head(backdoors$paths[grepl("^grazing <-", backdoors$paths)], 6)

# INTERPRETATION
# 39 paths connect grazing to seropositivity and 9 of them are backdoor paths.
# Every one of those 9 starts at village, age_class or sex, because those are
# the only three variables with an arrow into grazing. All 9 are open, and all 9
# have to be blocked. You do not have to work out how by hand, because the next
# section does it for you. The principle is that paths leaving grazing through a
# backward arrow carry bias, and paths leaving through a forward arrow carry the
# effect.


# ---- 6. Let the DAG choose your covariates ----------------------------------

# This is the step that makes the analysis causal, because adjustmentSets()
# reads the DAG and returns the covariates that block every bias path while
# leaving the causal paths intact. You no longer choose covariates by p-value,
# by stepwise selection or by what was in the last paper, all of which can
# introduce bias of their own.

# TOTAL effect: everything grazing does to seropositivity, by any route.
adjustmentSets(rvf_dag, exposure = "grazing", outcome = "seropositive",
               effect = "total")

# DIRECT effect: what is left when you hold the mediators fixed.
adjustmentSets(rvf_dag, exposure = "grazing", outcome = "seropositive",
               effect = "direct")

# The same question for each of the other exposures.
adjustmentSets(rvf_dag, "water_source", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "hours_at_water", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "night_housing", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "vector_control", "seropositive", effect = "total")
adjustmentSets(rvf_dag, "introduced_recently", "seropositive", effect = "total")

# INTERPRETATION
# Total effect of grazing needs {village, age_class, sex} and nothing else.
# Direct effect needs those three plus all four mediators.
# Each of the other exposures gets its own set, and they are not the same as
# each other: hours at water needs {grazing, water_source}, while vector
# control and recent introductions need {village} alone. One model
# containing every variable at once would answer none of these questions
# correctly, because the mediators belong in some sets and not others. In
# causal analysis you fit one model per exposure, not one model per paper.
#
# Record these sets in your methods, as the paper does in its supplementary
# table S3. They are the formal statement of what your estimates mean.


# ---- 7. Estimate the total effect -------------------------------------------

# Binary outcome, so logistic regression. The exposure plus the adjustment set,
# and nothing more.
m_total <- glm(seropositive ~ grazing + village + age_class + sex,
               family = binomial, data = cattle)

# A small helper to read coefficients as odds ratios.
or_table <- function(model, keep = NULL) {
  out <- data.frame(
    term = names(coef(model)),
    OR   = exp(coef(model)),
    lo   = exp(confint.default(model)[, 1]),
    hi   = exp(confint.default(model)[, 2])
  )
  if (!is.null(keep)) out <- out[grepl(keep, out$term), ]
  rownames(out) <- NULL
  out[, 2:4] <- round(out[, 2:4], 2)
  out
}

or_table(m_total, keep = "grazing")

# INTERPRETATION
# Read each row as: the odds of being seropositive for an animal under this
# grazing system, divided by the odds for a tethered animal of the same age and
# sex in the same village. Wetland grazing comes out near OR 6.5 with a
# confidence interval of 4.0 to 10.6, well clear of 1, so wetland grazing causes
# a large increase in infection risk under this DAG. Communal grazing is around
# OR 2.4. Fenced paddock is 1.2 with an interval spanning 1, meaning no
# detectable effect rather than no effect.
#
# Compare with the crude OR of about 13 in section 2. Half the crude
# association was village and age confounding. An interval that excludes 1 is
# the convention for "supported"; the paper uses the Bayesian equivalent, a 95%
# credible interval excluding zero on the log-odds scale.


# ---- 8. Estimate the direct effect, and read the mediation ------------------

m_direct <- glm(seropositive ~ grazing + village + age_class + sex +
                  water_source + hours_at_water + wildlife_contact + night_housing,
                family = binomial, data = cattle)

or_table(m_direct, keep = "grazing")

# Put the two side by side.
compare <- merge(or_table(m_total, "grazing"), or_table(m_direct, "grazing"),
                 by = "term", suffixes = c("_total", "_direct"))
compare

# INTERPRETATION
# The wetland grazing OR falls from about 6.5 (total) to about 2.5 (direct),
# and its interval drops from 4.0-10.6 to 1.4-4.5. That drop is the answer to a
# scientific question, not a problem with the model. It says most of what wetland grazing
# does to infection risk runs through the mediators: through drinking from
# rivers and dams, through hours spent at the water's edge, and through wildlife
# contact and not being kraaled. This is partial mediation, and the paper
# reports the same pattern for daily travel to Lake Victoria.
#
# Note what this means for control. If you could keep cattle in wetland areas
# but change their watering and housing, you would remove most of the risk. The
# wetland location itself is doing comparatively little.
#
# A warning about reporting: the direct-effect model is not "the better model".
# Its grazing coefficient answers a narrower question, and quoting it as the
# effect of grazing understates the total effect by a factor of nearly three.


# ---- 9. Four ways to get it wrong -------------------------------------------

# Same data, same outcome, four covariate choices. Only one of them answers the
# question "what does grazing do to infection risk".

m_naive     <- glm(seropositive ~ grazing, family = binomial, cattle)
m_kitchen   <- glm(seropositive ~ grazing + village + age_class + sex +
                     water_source + hours_at_water + wildlife_contact +
                     night_housing + vector_control +
                     introduced_recently + abortion_reported,
                   family = binomial, cattle)
m_collider  <- glm(seropositive ~ grazing + village + age_class + sex +
                     abortion_reported, family = binomial, cattle)

truth <- read.csv("simulation_truth.csv")

wrong_ways <- data.frame(
  approach = c("no adjustment", "DAG total effect", "DAG direct effect",
               "everything in the model", "adjusted for a collider"),
  wetland_OR = round(c(
    exp(coef(m_naive)["grazingwetland_grazing"]),
    exp(coef(m_total)["grazingwetland_grazing"]),
    exp(coef(m_direct)["grazingwetland_grazing"]),
    exp(coef(m_kitchen)["grazingwetland_grazing"]),
    exp(coef(m_collider)["grazingwetland_grazing"])
  ), 2),
  answers = c("nothing: confounded by village and age",
              "the total causal effect",
              "the effect not running through the mediators",
              "nothing: mediators and a collider both included",
              "nothing: collider bias added to a correct set")
)
wrong_ways

# The simulated direct effect that generated these data, for reference:
truth[truth$term == "grazing_wetland_grazing", ]
exp(truth$true_log_odds[truth$term == "grazing_wetland_grazing"])

# INTERPRETATION
# The no-adjustment OR of about 12.9 is nearly double the correct total effect
# of 6.5, because confounding by village and age is left in.
#
# The collider model shows why "adjust for more things to be safe" is wrong.
# Starting from a correct adjustment set and adding one variable caused by the
# outcome pulls the estimate from 6.5 down to about 4.6. More covariates is not
# more rigorous. The bias happens to point downwards here, but collider bias
# can go in either direction and you cannot predict which without knowing the
# thing you are trying to estimate.
#
# The kitchen-sink model gives a fourth number again, about 1.7. It conditions
# on the mediators, so it is not the total effect, and it conditions on a
# consequence of infection, so it is not the direct effect either. It is the
# most common specification in the literature and there is no causal quantity
# it corresponds to. Any reader who asks "the effect of grazing on what
# comparison?" cannot be given an answer.
#
# The true simulated direct effect is OR 2.59, and the DAG direct estimate of
# about 2.5 lands close to it. The small gap is the unmeasured confounder that
# section 14 reveals.


# ---- 10. A continuous exposure: hours at surface water ----------------------

# Its adjustment set from section 6 is {grazing, water_source}.
m_hours <- glm(seropositive ~ hours_at_water + grazing + water_source,
               family = binomial, data = cattle)
or_table(m_hours, keep = "hours_at_water")

# A per-hour OR is hard to feel. Convert it to the contrasts people care about.
beta_hours <- coef(m_hours)["hours_at_water"]
data.frame(extra_hours = c(1, 2, 4),
           OR = round(exp(beta_hours * c(1, 2, 4)), 2))

# INTERPRETATION
# Each extra hour per day at surface water multiplies the odds of
# seropositivity by about 1.22. A small per-unit effect accumulates: four hours
# a day, which is normal for wetland-grazed cattle, carries about 2.2 times the
# odds of an animal with no surface-water contact. Report the cumulative
# contrast as well as the coefficient, as the paper does when it converts a 1%
# per-minute effect into 82% higher odds over an hour.
#
# This estimate assumes the effect is linear on the log-odds scale. Check that
# before trusting it, for instance by refitting with splines() or with hours
# cut into bands, and say what you found.


# ---- 11. A second outcome: antibody level among seropositive animals --------

# Infection risk and infection intensity can have different causes. The paper
# models egg burden among egg-positive people; the equivalent here is antibody
# level among seropositive cattle.

positives <- filter(cattle, seropositive == 1)

# Antibody level is positive and right-skewed, so a Gamma family with a log
# link fits better than a Gaussian one.
hist(positives$titre_sp, breaks = 30, main = "Antibody level (S/P %) in seropositives",
     xlab = "S/P %")

m_titre <- glm(titre_sp ~ water_source + grazing + village + age_class + sex,
               family = Gamma(link = "log"), data = positives)
round(exp(cbind(MeanRatio = coef(m_titre),
                confint.default(m_titre)))[2:4, ], 2)

# INTERPRETATION
# Coefficients are mean ratios because of the log link. River or stream water
# comes out near 1.7, meaning animals on that water source have on average
# about 70% higher antibody levels than animals whose water is carried to them.
# Shared dams sit near 1.3 and boreholes near 0.9 with an interval spanning 1.
# Risk of infection and level of antibody respond to the same gradient here,
# which need not be true in general.
#
# READ THIS BEFORE REPORTING ANY "AMONG THE INFECTED" MODEL.
# Restricting to seropositive animals conditions on a variable that sits
# between the exposure and this outcome, which is a form of selection. Any
# unmeasured cause of seropositivity can then show up as a spurious
# exposure-titre association. The estimate is still worth reporting, but it is
# a weaker claim than the seropositivity models and should be labelled as such.
# This caveat applies equally to the paper's egg-burden results.


# ---- 12. The Bayesian version -----------------------------------------------

# Everything so far works with glm(). The paper uses Bayesian regression
# instead, for two concrete reasons. First, it returns a full posterior
# distribution for every effect rather than a point estimate and an interval,
# so you can state how probable an effect of a given size is. Second, those
# posterior draws are what make the counterfactual simulations in section 13
# possible with honest uncertainty attached.
#
# The first run of each model compiles C++ code and takes a minute or two.
# file= caches the fit, so later runs load from disk instantly.

library(brms)

# PRIORS. A prior says what effect sizes are plausible before seeing the data.
# normal(0, 1.5) on the log-odds scale says: centred on no effect, and odds
# ratios beyond about 1/20 to 20 are unlikely. That is weakly informative. It
# rules out absurd values and leaves the data to decide everything else. Run
# the model under a couple of different priors and report that the conclusions
# did not change, as the paper does in its prior sensitivity analysis.
pr_bin <- c(prior(normal(0, 1.5), class = "b"),
            prior(normal(0, 2),   class = "Intercept"))

b_total <- brm(seropositive ~ grazing + village + age_class + sex,
               family = bernoulli(link = "logit"), data = cattle,
               prior = pr_bin,
               chains = 4, iter = 3000, warmup = 1000, cores = 4,
               seed = 1, refresh = 0,
               file = "models/b_total")

b_direct <- brm(seropositive ~ grazing + village + age_class + sex +
                  water_source + hours_at_water + wildlife_contact + night_housing,
                family = bernoulli(link = "logit"), data = cattle,
                prior = pr_bin,
                chains = 4, iter = 3000, warmup = 1000, cores = 4,
                seed = 1, refresh = 0,
                file = "models/b_direct")

b_hours <- brm(seropositive ~ hours_at_water + grazing + water_source,
               family = bernoulli(link = "logit"), data = cattle,
               prior = pr_bin,
               chains = 4, iter = 3000, warmup = 1000, cores = 4,
               seed = 1, refresh = 0,
               file = "models/b_hours")

# CHECK CONVERGENCE BEFORE READING ANY RESULT.
# Rhat compares variation between chains with variation within them. Values
# above 1.01 mean the chains have not settled on the same answer and the
# numbers cannot be trusted. Bulk_ESS and Tail_ESS are effective sample sizes;
# anything above a few hundred is fine.
summary(b_total)
max(rhat(b_total), na.rm = TRUE)

# Posterior predictive check: simulate new datasets from the fitted model and
# compare them with the data you actually have. The outcome is 0/1, so use
# type = "bars": the pale bars are the observed counts of 0 and 1, and the dark
# points with intervals are what the model predicts. The default density
# overlay is designed for continuous outcomes and is unreadable here.
p_pp <- pp_check(b_total, type = "bars", ndraws = 500) +
  labs(title = "Posterior predictive check",
       subtitle = "Observed counts of 0 and 1 against what the model predicts")
ggsave("figures/02_ppcheck.png", p_pp, width = 6, height = 4, dpi = 200)
p_pp

# Posterior intervals on the log-odds scale, the paper's figure 3 layout.
p_eff <- mcmc_plot(b_total, variable = "^b_grazing", regex = TRUE,
                   prob = 0.66, prob_outer = 0.95) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "Total effect of grazing system on RVFV seropositivity",
       subtitle = "Posterior log-odds vs tethered cattle; thick 66%, thin 95% credible interval",
       x = "log-odds")
ggsave("figures/03_posterior_grazing.png", p_eff, width = 7, height = 3.5, dpi = 200)
p_eff

# Probability that each effect is positive, which is a statement a frequentist
# model cannot make.
post <- as.data.frame(as_draws_df(b_total))
sapply(post[, grep("^b_grazing", names(post))], function(x) mean(x > 0))

round(exp(fixef(b_total)[grep("grazing", rownames(fixef(b_total))), c(1, 3, 4)]), 2)

# INTERPRETATION
# Rhat is 1.00 and both effective sample sizes run to several thousand, so the
# chains converged. In the predictive check the observed counts sit inside the
# predicted range, so the model reproduces the observed seroprevalence.
# Odds ratios match the glm results closely, which is what weakly informative
# priors on 935 animals should give. The added value is in the last two
# outputs: the posterior probability that wetland grazing increases risk is
# 1.00, while for fenced paddock it is 0.79, which is a more useful statement
# than "not significant".


# ---- 13. Counterfactual interventions: what if we changed practice? ---------

# A coefficient answers "how much more likely is infection in wetland-grazed
# cattle". A counterfactual answers "what would seroprevalence be if those
# cattle were grazed differently". The second is the question a control
# programme actually asks.
#
# The recipe:
#   1. Copy the observed data.
#   2. Edit the exposure column to the intervened value, plus any variable
#      downstream of it that must change as a consequence.
#   3. Predict from the model that was fitted to the ORIGINAL data.
#   4. Compare predictions on the original and edited data, draw by draw.
#
# Step 3 matters because refitting the model to the edited data would answer a
# different question. The fitted model is the causal mechanism; the edited data
# is the new input.

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

# INTERPRETATION
# Seroprevalence in the whole population falls from about 36% to about 29%, a
# drop of roughly 7 percentage points. Among the wetland-grazed animals
# themselves it falls from about 71% to about 30%, a drop of roughly 41 points.
# That gap is the central finding of the paper restated in cattle: an
# intervention can work very well for the individuals it reaches and still move
# the population number only a little, because those individuals are a minority
# of the herd. Report both numbers. The subgroup figure is what a farmer
# experiences and the population figure is what a control programme achieves.
#
# Taking the difference draw by draw, rather than subtracting two summary
# means, is what gives you an honest interval on the reduction itself. Here it
# runs from 5.6 to 9.0 points and stays clear of zero, so the reduction is
# supported.

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
# When you move an animal from wetland to tethered grazing, its water source
# and its hours at water cannot stay as they were. Scenario A did not need to
# edit them, because the total-effect model deliberately contains no mediators,
# so the indirect routes are already inside its grazing coefficient. To edit
# them explicitly you need a model that contains them, which is the
# direct-effect model.
cf_C <- cattle
sw <- cattle$grazing == "wetland_grazing"
cf_C$grazing[sw]          <- "tethered_homestead"
cf_C$water_source[sw]     <- "carried_to_animal"
cf_C$hours_at_water[sw]   <- 0
cf_C$wildlife_contact[sw] <- 0
cf_C$night_housing[sw]    <- 1

cfC_wet <- pop_mean(b_direct, cf_C, wet)
obsC_wet <- pop_mean(b_direct, cattle, wet)
rbind(`wetland grazers, observed`               = q(obsC_wet),
      `tethered + mediators changed (direct)`   = q(cfC_wet),
      `tethered only (total-effect model)`      = q(cfa_wet))

# INTERPRETATION
# Scenario A left these animals at about 30% and scenario C leaves them at
# about 21%. The two routes disagree because they are not the same
# intervention, which is the lesson of this block.
#
# Scenario A moves the animals to tethered grazing and lets the model supply
# the average consequences of being tethered. Real tethered cattle are a mixed
# group: some still drink at a dam, some are not kraaled, and their mean hours
# at water is above zero. Scenario C instead forces every mediator to
# its most protective value at once, so it describes a far more complete
# intervention than simply changing how an animal is grazed.
#
# Pick the one that matches the policy you want to cost. "Move the herd off the
# wetland" is scenario A. "Move the herd off the wetland AND provide troughs
# AND kraal at night" is scenario C. State in the paper which mediators you
# intervened on, because the number depends entirely on that list.

# --- Scenario D: does the benefit differ by grazing system? -----------------
# This is a conditional average treatment effect, the paper's mechanistic
# counterfactual analysis. One average effect can hide a large effect in one
# subgroup and none in another.
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

# INTERPRETATION
# The reduction is about 15 percentage points for wetland grazers and 11 for
# communal grazers, against about 4 for fenced paddock and 2 for tethered
# cattle, who had little contact to remove in the first place. A programme that
# provided troughs to every herd would spend most of its budget on the last two
# groups for almost no return. Targeting wetland and communal herds captures
# nearly all of the benefit. A single average effect would have hidden that
# completely.


# ---- 14. Sensitivity analysis: how wrong could this be? ---------------------

# Every causal estimate rests on assumptions. Sensitivity analysis asks how far
# off the assumptions would have to be to overturn the conclusion. Three checks
# cover most of what reviewers ask.

# (a) UNMEASURED CONFOUNDING: the E-value.
# The E-value is the minimum strength of association, on the risk-ratio scale,
# that an unmeasured confounder would need with BOTH the exposure and the
# outcome to explain away the observed effect. Large E-value means a strong
# hidden confounder would be needed, so the finding is robust.
e_value <- function(or, p0) {
  rr <- or / (1 - p0 + p0 * or)         # OR to RR, given baseline risk p0
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}

p0 <- mean(cattle$seropositive[cattle$grazing == "tethered_homestead"])
ci_total <- exp(confint.default(m_total))["grazingwetland_grazing", ]
or_pt <- exp(coef(m_total)["grazingwetland_grazing"])

data.frame(baseline_risk = round(p0, 3),
           OR = round(or_pt, 2),
           E_point = round(e_value(or_pt, p0), 2),
           E_lower_CI = round(e_value(ci_total[1], p0), 2))

# INTERPRETATION
# An E-value of about 6.3 for the point estimate means a hidden confounder
# would need roughly a six-fold association with both wetland grazing and
# seropositivity, over and above village, age and sex, to reduce the effect to
# nothing. Few plausible confounders are that strong, so the finding is
# robust. The E-value for the lower confidence limit is about 4.8 and answers
# the more demanding question of what it would take to push the whole interval
# across 1. For comparison, the measured effect of village in this dataset is
# around OR 3 to 5, so the hidden confounder would have to be stronger than
# the strongest thing you did measure.

# WHAT IS ACTUALLY HIDDEN HERE.
# The simulation that made these data contains one variable that is not in the
# CSV: distance from the homestead to standing water. It affects both grazing
# system and infection, so it is a genuine unmeasured confounder, which is why
# the direct-effect estimate in section 9 came in below the simulated truth.
# Real surveys always have at least one of these. Say what you think yours is.

# The DAG you can draw and the DAG you believe are not the same thing. Mark the
# hidden variable as latent and ask dagitty whether the effect is still
# identifiable from the data you have.
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

# INTERPRETATION
# An empty result means no set of measured variables identifies the effect. The
# honest conclusion is that the estimate is conditional on the unmeasured
# confounder being weak, and the E-value quantifies how weak. This is the step
# most papers skip.

# (b) WHAT IF A MEDIATOR IS REALLY A CONFOUNDER?
# The DAG puts water_source, hours_at_water and wildlife_contact between
# grazing and infection. If any of them is in fact a common cause, it belongs
# in the adjustment set. You cannot settle this from the data, so show your
# estimate under every combination and let the reader judge.
mediators <- c("water_source", "hours_at_water", "wildlife_contact")
combos <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), 3))
names(combos) <- mediators

grid <- lapply(seq_len(nrow(combos)), function(i) {
  extra <- mediators[unlist(combos[i, ])]
  f <- reformulate(c("grazing", "village", "age_class", "sex", extra),
                   response = "seropositive")
  fit <- glm(f, family = binomial, data = cattle)
  ci <- exp(confint.default(fit))["grazingwetland_grazing", ]
  data.frame(conditioned_on = if (length(extra)) paste(extra, collapse = " + ") else "none (total effect)",
             OR = round(exp(coef(fit))["grazingwetland_grazing"], 2),
             lo = round(ci[1], 2), hi = round(ci[2], 2))
})
grid <- bind_rows(grid)
rownames(grid) <- NULL
grid[order(-grid$OR), ]

# INTERPRETATION
# The OR ranges from about 2.7 to 6.5 across the eight models. Every interval
# stays above 1, so the direction of the finding survives whichever structural
# assumption you make, and that is the headline of this table. The magnitude
# does not survive: conditioning on water source or on hours at water roughly
# halves the estimate, as partial mediation predicts, while conditioning on
# wildlife contact alone moves it from 6.5 to 5.6.
#
# Report the whole table rather than the single model you prefer. It tells the
# reader that your headline number rests on a structural choice, shows them
# what that choice is, and lets them see how much it is worth. The paper does
# the same across its eight mediator combinations in supplementary table S14.

# (c) MEASUREMENT ERROR IN A SELF-REPORTED EXPOSURE.
# hours_at_water comes from asking the farmer, and farmers under-report time at
# water. If they report only a fraction f of the true hours, the fitted
# coefficient under-states the true one by the same fraction. Correct it across
# a range of f and see whether the conclusion changes.
frac <- c(1.0, 0.8, 0.6, 0.5)
data.frame(reported_fraction = frac,
           OR_per_hour = round(exp(beta_hours / frac), 2),
           OR_over_4h  = round(exp(4 * beta_hours / frac), 2))

# INTERPRETATION
# Under-reporting biases this estimate towards no effect, so the observed
# coefficient is a lower bound. If farmers report half the true time, the
# four-hour contrast rises from OR 2.2 to OR 4.8. The direction of the bias is
# the useful part: it means the real effect is at least as large as reported,
# never smaller.


# ---- 15. How to write this up -----------------------------------------------

# Methods section, in order:
#   1. The DAG, as a figure, with a sentence on where each arrow came from.
#   2. That you tested the implied conditional independencies, and the result.
#   3. A table of the adjustment set used for each exposure-outcome pair.
#   4. The model family and link for each outcome, and the priors if Bayesian.
#   5. Convergence diagnostics, Rhat and effective sample size.
#   6. The sensitivity analyses, including E-values.
#
# Results section:
#   Report total and direct effects separately and label them as such. Give
#   counterfactual results for both the whole population and the targeted
#   subgroup, because they answer different policy questions. Use credible or
#   confidence intervals rather than p-values alone.
#
# Say plainly what the estimates do and do not support. "Under the causal
# structure in figure 1, wetland grazing causes an estimated 6.5-fold increase
# in the odds of RVFV seropositivity" is a claim a reader can check and
# disagree with. "Wetland grazing was significantly associated with
# seropositivity" is not.
#
# Three things to remember:
#   The DAG is an assumption, not a result. Published it so others can argue.
#   Total and direct effects are different quantities. Never report one as the
#     other.
#   More covariates is not more rigorous. Mediators and colliders both make
#     estimates worse.
#
# Further reading:
#   Pearl & Mackenzie, The Book of Why, for the intuition.
#   Hernan & Robins, Causal Inference: What If, free at
#     https://miguelhernan.org/whatifbook, for the formal treatment.
#   Textor et al. 2016, Int J Epidemiol 45:1887, the dagitty paper.
#   VanderWeele & Ding 2017, Ann Intern Med 167:268, on E-values.
#   Lim et al. 2026, BMJ Public Health 4:e004500, the analysis this follows.

# =============================================================================
