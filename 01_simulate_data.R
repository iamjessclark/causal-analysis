# 01_simulate_data.R
# Simulates a cross-sectional RVFV serosurvey in cattle, clustered by herd.
# The structural equations below define the DAG that 02_causal_analysis_tutorial.R
# recovers. Run this once; it writes two CSV files and then you can forget it.
#
# Cattle are sampled from ~52 herds in 5 villages. Management is mostly decided
# at herd level, and two unmeasured herd-level variables drive infection:
# distance from the homestead to standing water (which also drives management,
# so it confounds) and a residual herd risk (which does not, so it only
# induces clustering). Neither is written to the CSV.
#
# Outputs
#   rvf_cattle_serology.csv  the analysis dataset (what a field team would hand you)
#   simulation_truth.csv     the true coefficients, for checking the tutorial's answers

set.seed(2026)

inv_logit <- function(x) 1 / (1 + exp(-x))

softmax_pick <- function(eta, labels) {
  p <- exp(eta)
  p <- p / rowSums(p)
  labels[apply(p, 1, function(pi) sample.int(length(labels), 1L, prob = pi))]
}

# ---- Herd level --------------------------------------------------------------

n_herd <- 52
village_labels <- c("A", "B", "C", "D", "E")

village_h <- factor(sample(village_labels, n_herd, TRUE, prob = c(.22, .18, .20, .22, .18)),
                    levels = village_labels)
wetland_h <- c(A = -0.8, B = 0.9, C = 0.1, D = 0.5, E = -0.4)[as.character(village_h)]
size_h    <- pmax(8, round(rlnorm(n_herd, log(17) + 0.20 * wetland_h, 0.35)))

# Unmeasured herd-level variables.
wprox_h   <- rnorm(n_herd, wetland_h, 0.9)   # homestead distance to standing water
herd_risk_h <- rnorm(n_herd, 0, 0.70)        # how risky this herd is for
                                             # reasons nothing in the survey
                                             # records: local breeding sites,
                                             # neighbouring cattle, and so on

# Herd-level tendencies in management. Kept modest so that village and the
# arrows already in the DAG account for the dependence between these variables.
gshift_h  <- matrix(rnorm(n_herd * 4, 0, 0.5), n_herd, 4)
wshift_h  <- matrix(rnorm(n_herd * 4, 0, 0.5), n_herd, 4)
nhshift_h <- rnorm(n_herd, 0, 0.6)
intro_h   <- rnorm(n_herd, 0, 0.4)

# ---- Expand to animals -------------------------------------------------------

herd <- rep(seq_len(n_herd), times = size_h)
n    <- length(herd)

village   <- village_h[herd]
wetland   <- wetland_h[herd]
herd_size <- size_h[herd]
wprox     <- wprox_h[herd]
herd_id   <- factor(sprintf("%s-%02d", village_h[herd], herd))

age_class <- factor(sample(c("calf", "juvenile", "adult"), n, TRUE, prob = c(.18, .27, .55)),
                    levels = c("calf", "juvenile", "adult"))
age_years <- round(ifelse(age_class == "calf", runif(n, 0.2, 1),
                   ifelse(age_class == "juvenile", runif(n, 1, 3), runif(n, 3, 12))), 1)

sex <- factor(sample(c("female", "male"), n, TRUE, prob = c(.72, .28)),
              levels = c("female", "male"))

is_calf <- as.numeric(age_class == "calf")
is_male <- as.numeric(sex == "male")

# ---- Grazing system ----------------------------------------------------------
# Mostly a herd decision, but calves are tethered near the homestead whatever
# the rest of the herd does, which is what creates within-herd variation.

grazing_labels <- c("tethered_homestead", "fenced_paddock", "communal_grazing", "wetland_grazing")

eta_g <- cbind(
  0.60 + 1.30 * is_calf - 0.70 * wetland - 0.20 * is_male - 0.25 * wprox,
  0.20 + 0.40 * is_calf - 0.40 * wetland - 0.10 * is_male - 0.10 * wprox,
  0.50 - 0.20 * is_calf + 0.50 * wetland + 0.30 * is_male + 0.20 * wprox,
 -0.40 - 0.90 * is_calf + 1.10 * wetland + 0.35 * is_male + 0.45 * wprox
) + gshift_h[herd, ]
grazing <- factor(softmax_pick(eta_g, grazing_labels), levels = grazing_labels)

g_tether   <- as.numeric(grazing == "tethered_homestead")
g_paddock  <- as.numeric(grazing == "fenced_paddock")
g_communal <- as.numeric(grazing == "communal_grazing")
g_wetland  <- as.numeric(grazing == "wetland_grazing")

# ---- Water source ------------------------------------------------------------

water_labels <- c("carried_to_animal", "borehole_trough", "shared_dam", "river_stream")

eta_w <- cbind(
  1.00 + 1.00 * g_tether - 0.30 * wprox - 0.20 * wetland,
  0.60 + 0.50 * g_tether + 0.40 * g_paddock - 0.20 * wprox - 0.10 * wetland,
  0.20 + 0.80 * g_communal + 0.90 * g_wetland + 0.30 * wprox + 0.20 * wetland,
 -0.20 + 0.60 * g_communal + 1.30 * g_wetland + 0.50 * wprox + 0.30 * wetland
) + wshift_h[herd, ]
water_source <- factor(softmax_pick(eta_w, water_labels), levels = water_labels)

w_carried  <- as.numeric(water_source == "carried_to_animal")
w_borehole <- as.numeric(water_source == "borehole_trough")
w_dam      <- as.numeric(water_source == "shared_dam")
w_river    <- as.numeric(water_source == "river_stream")

# ---- Hours per day at surface water ------------------------------------------

mu_h <- exp(-0.80 + 0.30 * g_paddock + 0.90 * g_communal + 1.40 * g_wetland +
              0.50 * w_dam + 0.80 * w_river + 0.15 * wprox)
hours_at_water <- round(rgamma(n, shape = 2.2, rate = 2.2 / mu_h), 2)
hours_at_water[w_carried == 1] <- 0            # water carried to the animal
hours_at_water <- pmin(hours_at_water, 9)

# ---- Remaining management variables ------------------------------------------

wildlife_contact <- rbinom(n, 1, inv_logit(-1.20 + 0.40 * wetland +
                                             0.50 * g_communal + 1.20 * g_wetland))

night_housing <- rbinom(n, 1, inv_logit(1.00 - 0.50 * wetland + nhshift_h[herd] +
                                          0.80 * g_tether + 0.40 * g_paddock -
                                          0.50 * g_wetland))

vc_labels <- c("none", "irregular", "routine")
eta_v <- cbind(
  0.40 + 0.30 * wetland,
  0.30 + 0.10 * wetland,
 -0.10 - 0.20 * wetland
)
vector_control <- factor(softmax_pick(eta_v, vc_labels), levels = vc_labels)

introduced_recently <- rbinom(n, 1, inv_logit(-1.40 + 0.30 * wetland + intro_h[herd] +
                                                0.45 * log(herd_size / 17)))

# ---- Outcome 1: natural-infection seropositivity -----------------------------

truth <- c(
  intercept                  = -2.95,
  village_B                  =  0.95,
  village_C                  =  0.40,
  village_D                  =  0.65,
  village_E                  =  0.20,
  age_juvenile               =  0.70,
  age_adult                  =  1.40,
  sex_male                   =  0.15,
  grazing_fenced_paddock     =  0.20,
  grazing_communal_grazing   =  0.55,
  grazing_wetland_grazing    =  0.95,
  water_borehole_trough      =  0.15,
  water_shared_dam           =  0.60,
  water_river_stream         =  0.85,
  hours_at_water             =  0.22,
  wildlife_contact           =  0.35,
  night_housing              = -0.45,
  vector_control_irregular   = -0.25,
  vector_control_routine     = -0.55,
  introduced_recently        =  0.45,
  unmeasured_water_proximity =  0.35,
  unmeasured_herd_risk_sd    =  0.70
)

eta_i <- truth["intercept"] +
  truth["village_B"] * (village == "B") + truth["village_C"] * (village == "C") +
  truth["village_D"] * (village == "D") + truth["village_E"] * (village == "E") +
  truth["age_juvenile"] * (age_class == "juvenile") +
  truth["age_adult"] * (age_class == "adult") +
  truth["sex_male"] * is_male +
  truth["grazing_fenced_paddock"] * g_paddock +
  truth["grazing_communal_grazing"] * g_communal +
  truth["grazing_wetland_grazing"] * g_wetland +
  truth["water_borehole_trough"] * w_borehole +
  truth["water_shared_dam"] * w_dam +
  truth["water_river_stream"] * w_river +
  truth["hours_at_water"] * hours_at_water +
  truth["wildlife_contact"] * wildlife_contact +
  truth["night_housing"] * night_housing +
  truth["vector_control_irregular"] * (vector_control == "irregular") +
  truth["vector_control_routine"] * (vector_control == "routine") +
  truth["introduced_recently"] * introduced_recently +
  truth["unmeasured_water_proximity"] * wprox +
  herd_risk_h[herd]

seropositive <- rbinom(n, 1, inv_logit(eta_i))

# ---- Outcome 2: antibody level among seropositive animals --------------------

mu_t <- exp(log(45) + 0.20 * w_dam + 0.35 * w_river + 0.030 * hours_at_water +
              0.15 * (age_class == "adult") + 0.10 * g_wetland)
titre_sp <- round(rgamma(n, shape = 3.5, rate = 3.5 / mu_t), 1)
titre_sp[seropositive == 0] <- NA_real_

# ---- A descendant of infection: reported abortion / stillbirth ---------------
# Recorded for this animal, and caused by its own infection, so it is a
# consequence of the outcome and must never be used as a covariate.

abortion_reported <- rbinom(n, 1, inv_logit(-2.20 + 1.80 * seropositive +
                                              1.00 * g_wetland +
                                              0.60 * g_communal))

# ---- Write out ---------------------------------------------------------------

cattle <- data.frame(
  animal_id = sprintf("RVF%04d", seq_len(n)),
  village, herd_id, herd_size,
  age_class, age_years, sex,
  grazing, water_source, hours_at_water,
  wildlife_contact, night_housing, vector_control, introduced_recently,
  seropositive, titre_sp, abortion_reported
)

write.csv(cattle, "rvf_cattle_serology.csv", row.names = FALSE)
write.csv(data.frame(term = names(truth), true_log_odds = unname(truth)),
          "simulation_truth.csv", row.names = FALSE)

cat("Wrote rvf_cattle_serology.csv:", nrow(cattle), "cattle in",
    nlevels(cattle$herd_id), "herds,",
    round(100 * mean(cattle$seropositive)), "% seropositive\n")
