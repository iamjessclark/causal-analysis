# 01_simulate_data.R
# Simulates a cross-sectional RVFV serosurvey in cattle, clustered by herd.
# The structural equations below define the DAG that 02_causal_analysis_tutorial.R
# recovers. Run this once; it writes two CSV files and then you can forget it.
#
# Cattle are sampled from 200 herds in 38 villages. Management is mostly decided
# at herd level, and two unmeasured herd-level variables drive infection:
# distance from the homestead to standing water (which also drives management,
# so it confounds) and a residual herd risk (which does not, so it only
# induces clustering). Neither is written to the CSV as they purposefully 
# represent data we have not collected.
#
# Outputs
#   rvf_cattle_serology.csv  the analysis dataset (what a field team would hand you)
#   simulation_truth.csv     the true coefficients, for checking the tutorial's answers

# nothing on a computer is truly random, it follows an algorithm so this pins 
# where the algorithm starts so that you can always reproduce the same output
set.seed(19) 

# Turns a number on the log-odds scale into a probability between 0 and 1. Every
# binary variable below is built by writing its log-odds as a sum of causes,
# passing that sum through here, and drawing from rbinom() with the result.
inv_logit <- function(x) 1 / (1 + exp(-x))

# The same idea for a variable with more than two categories, such as grazing
# system. `scores` is a matrix with one row per animal and one column per
# category, holding that animal's score for each category which is then 
# modulated later by the simulated effect sizes and `labels` names
# the categories in the same column order.
softmax_pick <- function(scores, labels) {
  p <- exp(scores)          # make every score positive
  p <- p / rowSums(p)       # divide by the row total so each row sums to 1
  # each row is now a set of probabilities, so draw one category per row
  labels[apply(p, 1, function(pi) sample.int(length(labels), 1L, prob = pi))]
}

# ---- Herd levåel --------------------------------------------------------------

# taken from the real study design but there were more than 200 households/herd
n_village <- 38
n_herd    <- 200
# this is just creating V01, V02 to name the villages
village_labels <- sprintf("V%02d", seq_len(n_village))

# Each village contributes 5 or 6 herds. Villages differ in how much standing
# water they sit near, which is the root of everything downstream.
wetland_v <- rnorm(n_village, 0, 0.65)

village_h <- factor(sample(rep(village_labels, length.out = n_herd)),
                    levels = village_labels)
v_idx     <- as.integer(village_h)
wetland_h <- wetland_v[v_idx]
size_h    <- pmax(6, round(rlnorm(n_herd, log(10) + 0.20 * wetland_h, 0.35)))

# A village's effect on infection is mostly its wetland position, plus a
# village-specific residual. With 38 villages this is a variance you could
# estimate, which is why village is a candidate random effect in real data.
vill_slope <- 0.90
vill_sd    <- 0.30
vill_eff_v <- vill_slope * wetland_v + rnorm(n_village, 0, vill_sd)

## Unmeasured herd-level variables.

# Boma distance to standing water
wprox_h   <- rnorm(n_herd, wetland_h, 0.9)  
# how risky this herd is for reasons nothing in the survey records: local 
# breeding sites, neighbouring cattle, and so on
herd_risk_h <- rnorm(n_herd, 0, 0.70)        
                                        
## Herd-level tendencies in management. Kept modest so that village and the
# arrows already in the DAG account for the dependence between these variables.

# Each herd's preference for each of the four grazing systems: a 200 x 4 matrix,
# one row per herd and one column per grazing system. A positive entry means
# this herd favours that system more than average, a negative entry less.
gshift_h  <- matrix(rnorm(n_herd * 4, 0, 0.5), n_herd, 4)

# Each herd's preference for each of the four water sources, same layout.
wshift_h  <- matrix(rnorm(n_herd * 4, 0, 0.5), n_herd, 4)

# Each herd's tendency to buy animals in, again one value per herd.
intro_h   <- rnorm(n_herd, 0, 0.4)

# ---- Expand to animals -------------------------------------------------------

herd <- rep(seq_len(n_herd), times = size_h)
n    <- length(herd)

village   <- village_h[herd]
wetland   <- wetland_h[herd]
vill_eff  <- vill_eff_v[v_idx[herd]]
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
# Grazing is mostly a herd-level decision, so without something to break it up
# every animal in a herd would share one system and grazing could not be told
# apart from herd. Calves are tethered near the homestead whatever the rest of
# the herd does, which supplies that within-herd variation. 198 of the 200
# herds end up with at least two grazing systems.
grazing_labels <- c("tethered_homestead", "fenced_paddock", "communal_grazing", "wetland_grazing")

# One column per grazing system, in the order of grazing_labels, and one row per
# animal. Each number says how much that animal is pulled towards that system:
# calves towards tethering (+1.30) and away from wetland (-0.90), animals in
# wetland villages the other way. Adding gshift_h[herd, ] tilts all four columns
# by that animal's herd preference.
grazing_score <- cbind(
  0.60 + 1.30 * is_calf - 0.70 * wetland - 0.20 * is_male - 0.25 * wprox,
  0.20 + 0.40 * is_calf - 0.40 * wetland - 0.10 * is_male - 0.10 * wprox,
  0.50 - 0.20 * is_calf + 0.50 * wetland + 0.30 * is_male + 0.20 * wprox,
 -0.40 - 0.90 * is_calf + 1.10 * wetland + 0.35 * is_male + 0.45 * wprox
) + gshift_h[herd, ]

# softmax_pick() converts each animal's four scores into four probabilities that
# sum to 1, then draws one grazing system using them. A score 1 point above
# another gives about 2.7 times the probability, so the highest-scoring system
# is the most likely without being certain. This is what puts a few wetland
# grazers in dry villages and a few tethered animals in wet ones, which is what
# the tutorial needs in order to separate grazing from village.
grazing <- factor(softmax_pick(grazing_score, grazing_labels), levels = grazing_labels)

g_tether   <- as.numeric(grazing == "tethered_homestead")
g_paddock  <- as.numeric(grazing == "fenced_paddock")
g_communal <- as.numeric(grazing == "communal_grazing")
g_wetland  <- as.numeric(grazing == "wetland_grazing")

# ---- Water source ------------------------------------------------------------

water_labels <- c("carried_to_animal", "borehole_trough", "shared_dam", "river_stream")

water_score <- cbind(
  1.00 + 1.00 * g_tether - 0.30 * wprox - 0.20 * wetland,
  0.60 + 0.50 * g_tether + 0.40 * g_paddock - 0.20 * wprox - 0.10 * wetland,
  0.20 + 0.80 * g_communal + 0.90 * g_wetland + 0.30 * wprox + 0.20 * wetland,
 -0.20 + 0.60 * g_communal + 1.30 * g_wetland + 0.50 * wprox + 0.30 * wetland
) + wshift_h[herd, ]
water_source <- factor(softmax_pick(water_score, water_labels), levels = water_labels)

w_carried  <- as.numeric(water_source == "carried_to_animal")
w_borehole <- as.numeric(water_source == "borehole_trough")
w_dam      <- as.numeric(water_source == "shared_dam")
w_river    <- as.numeric(water_source == "river_stream")

# ---- Hours per day at surface water ------------------------------------------

# Hours depend only on grazing system and water source, so the DAG's claim that
# village is independent of hours once both are held fixed is true.
mu_h <- exp(-0.80 + 0.30 * g_paddock + 0.90 * g_communal + 1.40 * g_wetland +
              0.50 * w_dam + 0.80 * w_river)
hours_at_water <- round(rgamma(n, shape = 2.2, rate = 2.2 / mu_h), 2)
hours_at_water[w_carried == 1] <- 0            # water carried to the animal
hours_at_water <- pmin(hours_at_water, 9)

# ---- Remaining management variables ------------------------------------------

wildlife_contact <- rbinom(n, 1, inv_logit(-1.20 + 0.40 * wetland +
                                             0.50 * g_communal + 1.20 * g_wetland))

introduced_recently <- rbinom(n, 1, inv_logit(-1.40 + 0.30 * wetland + intro_h[herd] +
                                                0.45 * log(herd_size / 10)))

# ---- Outcome: natural-infection seropositivity -------------------------------

truth <- c(
  intercept                  = -3.20,
  village_wetland_slope      = vill_slope,
  village_residual_sd        = vill_sd,
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
  introduced_recently        =  0.45,
  unmeasured_water_proximity =  0.35,
  unmeasured_herd_risk_sd    =  0.70
)

infection_logodds <- truth["intercept"] +
  vill_eff +
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
  truth["introduced_recently"] * introduced_recently +
  truth["unmeasured_water_proximity"] * wprox +
  herd_risk_h[herd]

seropositive <- rbinom(n, 1, inv_logit(infection_logodds))

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
  wildlife_contact, introduced_recently,
  seropositive, abortion_reported
)

write.csv(cattle, "rvf_cattle_serology.csv", row.names = FALSE)
write.csv(data.frame(term = names(truth), true_log_odds = unname(truth)),
          "simulation_truth.csv", row.names = FALSE)

cat("Wrote rvf_cattle_serology.csv:", nrow(cattle), "cattle in",
    nlevels(cattle$herd_id), "herds in", nlevels(cattle$village), "villages,",
    round(100 * mean(cattle$seropositive)), "% seropositive\n")
