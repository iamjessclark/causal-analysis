# Causal analysis in R: a beginner's tutorial

A worked, runnable introduction to structural causal modelling in R, using a
simulated Rift Valley fever virus (RVFV) serosurvey in cattle.

The analysis follows the structure of Lim et al. (2026), *Local habitual
movement as a mechanism for* Schistosoma mansoni *transmission resurgence: a
causal analysis of cross-sectional survey data*, BMJ Public Health 4:e004500,
[doi:10.1136/bmjph-2025-004500](https://doi.org/10.1136/bmjph-2025-004500).
Where that paper models human travel to Lake Victoria and schistosome
infection, this tutorial models cattle grazing and watering practice and RVFV
seropositivity. The method is identical: DAG, adjustment sets, total and direct
effects, Bayesian estimation, counterfactual simulation, sensitivity analysis.

## Files

| File | What it is |
|---|---|
| `01_simulate_data.R` | Generates the dataset. Run once. Not part of the tutorial. |
| `02_causal_analysis_tutorial.R` | The tutorial. Read and run it section by section. |
| `03_random_effects_demo.R` | Optional appendix: should village and herd be random effects? Needs `lme4`. |
| `rvf_cattle_serology.csv` | The analysis dataset, written by `01`. |
| `simulation_truth.csv` | The true coefficients, so you can check the answers. |
| `figures/` | Plots written by the tutorial. |
| `models/` | Cached Bayesian fits, so re-runs are fast. |

## Running it

```r
install.packages(c("dagitty", "ggdag", "ggplot2", "dplyr", "tidyr",
                   "fastDummies", "CCP", "brms", "patchwork"))
```

```r
source("01_simulate_data.R")
```

Then open `02_causal_analysis_tutorial.R` and work through it. Sections 1 to 11
use `glm()` and run in seconds. Sections 12 to 14 use `brms`, which compiles
Stan code on first use and takes a few minutes. After that the `file=` cache in
`models/` makes re-runs instant. Delete `models/` if you change the data or a
model formula.

`brms` needs a working C++ toolchain. On macOS install Xcode command line tools
with `xcode-select --install`; on Windows install Rtools. If the install gives
you trouble, sections 1 to 11 stand on their own and cover the whole causal
logic.

## What the tutorial covers

1. Reading the data and setting reference levels
2. Why the crude comparison is not a causal estimate
3. Writing a DAG in `dagitty`
4. Testing the DAG against the data with implied conditional independencies
5. Paths, confounders, mediators and colliders
6. Letting the DAG choose the covariates, with `adjustmentSets()`
7. Estimating the total effect
8. Estimating the direct effect, and reading mediation from the gap
9. Four ways to get the answer wrong, scored against the simulated truth
10. A continuous exposure, and converting a per-unit effect into a useful contrast
11. A second outcome, with a note on selection bias in "among the infected" models
12. The Bayesian version: priors, convergence, posterior probabilities
13. Counterfactual interventions, population versus subgroup, and effects by subgroup
14. Sensitivity analysis: E-values, mediator-versus-confounder structure, measurement error
15. How to write it up

## The simulated data

935 cattle from 52 herds across 5 villages, one blood sample and one management
questionnaire each. Overall seroprevalence is 36%, matching the source paper.
Management is decided mostly at herd level, so animals in the same herd
resemble each other.

Six management exposures, each chosen because the RVF literature reports it as
a risk factor:

| Variable | Levels | Literature basis |
|---|---|---|
| `grazing` | tethered at homestead / fenced paddock / communal / wetland | Grazing in park and farm areas raised seropositivity in Ugandan cattle (AOR 4.9) |
| `water_source` | carried to the animal / borehole trough / shared dam / river or stream | Rivers and streams at grazing fields carried an OR near 10.8 in Nigerian pastoral herds; evidence on water-point type is mixed, with a Madagascar study finding no effect |
| `night_housing` | kraaled overnight or not | South African surveys found higher seroprevalence in cattle not kraaled at night |
| `vector_control` | none / irregular / routine | Acaricide and insecticide use reduces vector contact |
| `wildlife_contact` | shares grazing or watering with wildlife | Shared grazing with wildlife is reported from Ethiopia and Cameroon |
| `introduced_recently` | bought into the herd in the last 12 months | Owners who buy animals to renew the herd had higher seroconversion in Madagascar |

Plus `hours_at_water`, a continuous duration of surface-water contact that acts
as the mediator, and `abortion_reported`, which is caused by infection and is in
the data specifically so the tutorial can show what happens when you adjust for
a collider.

The simulation also contains two variables that are **not** written to the CSV.
Distance from the homestead to standing water causes both grazing system and
infection, so it is an unmeasured confounder, and section 14 uses it to show
why the recovered direct effect sits below the simulated truth and how an
E-value quantifies the gap. A residual herd risk raises or lowers every
animal's risk in a herd without touching management, so it produces clustering
without confounding.

## What the simulation simplifies

**Clustering.** Animals cluster in herds, with a herd standard deviation of
0.63 on the log-odds scale after village is accounted for, which is an
intraclass correlation of 0.11. The tutorial treats village as a fixed effect
and leaves herd alone, which is what the source paper does.

Because adding `(1 | herd_id)` is the obvious next thought, and because it does
not do what most people expect, `03_random_effects_demo.R` works the question
through on data where the right answer is known. The summary: a random
intercept does **not** control cluster-level confounding, and in that script it
leaves a log-odds estimate at 1.06 against a truth of 0.70, where fixed effects
give 0.78 and the Mundlak specification gives 0.73. Village stays a fixed effect
here because five clusters is too few to estimate a variance from. Clustering
changes the standard errors and, where cluster-level confounding exists, the
point estimate; it never changes the adjustment sets, because those come from
the DAG.

## Further reading

- Pearl & Mackenzie, *The Book of Why*, for the intuition
- Hernán & Robins, *Causal Inference: What If*, free at <https://miguelhernan.org/whatifbook>
- Textor et al. 2016, *Int J Epidemiol* 45:1887, the `dagitty` paper
- VanderWeele & Ding 2017, *Ann Intern Med* 167:268, on E-values
- Lim et al. 2026, *BMJ Public Health* 4:e004500, the analysis this follows
