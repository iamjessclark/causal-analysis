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
| `rvf_cattle_serology.csv` | The analysis dataset, written by `01` - can just load this and not run 1 |
| `02_causal_analysis_tutorial.R` | The tutorial code. Read and run it section by section. |
| `TUTORIAL.md` | The written tutorial: the explanation for each section, the function reference and how to read each output. |
| `03_random_effects_demo.R` | Optional appendix: should village and herd be random effects? Needs `lme4`. |
| `simulation_truth.csv` | The true coefficients, so you can check the answers. |
| `figures/` | Plots written by the tutorial. |
| `models/` | Cached Bayesian fits, so re-runs are fast. |

## Running it

```r
install.packages(c("dagitty", "ggdag", "tidyverse", "brms",
                   "fastDummies", "CCP"))
```

```r
source("01_simulate_data.R")
```

Then open `02_causal_analysis_tutorial.R` and work through it. Every model is
fitted in `brms`, as in the source paper, so the first run compiles around a
dozen Stan programs and takes 10 to 20 minutes. After that the `file=` cache in
`models/` makes re-runs instant. Delete `models/` if you change the data or a
model formula.

`brms` needs a working C++ toolchain. On macOS install Xcode command line tools
with `xcode-select --install`; on Windows install Rtools.

## What the tutorial covers

1. Reading the data and setting reference levels
2. How every model is fitted: `brms`, priors, and the helpers used throughout
3. Why the crude comparison is not a causal estimate
4. Writing a DAG in `dagitty`
5. Testing the DAG against the data with implied conditional independencies
6. Paths, confounders, mediators and colliders
7. Letting the DAG choose the covariates, with `adjustmentSets()`
8. Estimating the total effect, with convergence checks and posterior probabilities
9. Estimating the direct effect, and reading mediation from the gap
10. Alternative model specifications, compared against the simulated truth
11. A continuous exposure, and converting a per-unit effect into a useful contrast
12. Counterfactual interventions, population versus subgroup, and effects by subgroup
13. Sensitivity analysis: E-values, mediator-versus-confounder structure, measurement error
14. How to write it up

## The simulated data

2220 cattle from 200 herds across 38 villages, with one blood sample per animal
and the management of each animal recorded from its owner's questionnaire,
matching the design of the CEPI study. Overall seroprevalence is 36%, as in the
source paper, with a village median of 26% and a range from 2% to 79%. Every
management variable varies within a herd, because the youngest animals are
tethered at the homestead while the rest go out, and 197 of the 200 herds hold
two or more grazing systems.

Four management exposures, each chosen because the RVF literature reports it as
a risk factor. Herd-level practices such as night kraaling and acaricide use
are left out, because they take one value for a whole herd and testing a DAG
that contains them needs multilevel methods beyond a beginner tutorial.

| Variable | Levels | Literature basis |
|---|---|---|
| `grazing` | tethered at homestead / fenced paddock / communal / wetland | Grazing in park and farm areas raised seropositivity in Ugandan cattle (AOR 4.9) |
| `water_source` | carried to the animal / borehole trough / shared dam / river or stream | Rivers and streams at grazing fields carried an OR near 10.8 in Nigerian pastoral herds; evidence on water-point type is mixed, with a Madagascar study finding no effect |
| `wildlife_contact` | shares grazing or watering with wildlife | Shared grazing with wildlife is reported from Ethiopia and Cameroon |
| `introduced_recently` | bought into the herd in the last 12 months | Owners who buy animals to renew the herd had higher seroconversion in Madagascar |

Plus `hours_at_water`, a continuous duration of surface-water contact that acts
as the mediator.

The simulation also contains two variables that are **not** written to the CSV.
Distance from the homestead to standing water causes both grazing system and
infection, so it is an unmeasured confounder, and section 13 uses it to show
why the recovered direct effect misses the simulated truth and how an
E-value quantifies the gap. A residual herd risk raises or lowers every
animal's risk in a herd without touching management, so it produces clustering
without confounding.

## What the simulation simplifies

**Clustering.** Animals cluster at two levels. After village is accounted for,
the herd standard deviation is 0.76 on the log-odds scale, an intraclass
correlation of 0.15; fit village and herd together as nested random intercepts
and the village standard deviation is 0.97 against 0.91 for herd within
village, so the two levels contribute about equally. The tutorial treats village as a
fixed effect and leaves herd alone, which is what the source paper does.

Because adding `(1 | herd_id)` is the obvious next thought, and because it does
not do what most people expect, `03_random_effects_demo.R` works the question
through on data where the right answer is known. The summary: a random
intercept does **not** control cluster-level confounding, and in that script it
leaves a log-odds estimate at 1.06 against a truth of 0.70, where fixed effects
give 0.78 and the Mundlak specification gives 0.73.

With 38 villages there are enough clusters to estimate a village variance, so
village could be a random effect here. It stays a fixed effect for a different
reason: 37 dummy variables remove every village-level confounder without
assuming anything about how the village effects relate to grazing, and a
random intercept does not. The cost is 37 coefficients you have to scroll past
in `summary()` and no estimate of the village variance itself. Clustering
changes the standard errors and, where cluster-level confounding exists, the
point estimate. It never changes the adjustment sets, because those come from
the DAG.

## Further reading

- Pearl & Mackenzie, *The Book of Why*, for the intuition
- Hernán & Robins, *Causal Inference: What If*, free at <https://miguelhernan.org/whatifbook>
- Textor et al. 2016, *Int J Epidemiol* 45:1887, the `dagitty` paper
- VanderWeele & Ding 2017, *Ann Intern Med* 167:268, on E-values
- Lim et al. 2026, *BMJ Public Health* 4:e004500, the analysis this follows
