# COVID-19 Vaccination - Epidemiological Simulations

## Introduction

This repository includes the scripts used to generate the results in "Testing fractional doses of COVID-19 vaccines".

We take the model in "Neutralizing antibody levels are highly predictive of immune protection from symptomatic SARS-CoV-2 infection" (Khoury et al., 2021) and extend the analysis to immunogenicity levels observed for doses lower than the ones currently approved for distribution.

We implement a SEIR (Susceptible, Exposed, Infectious, Recover) model and estimate the health outcomes under different vaccination scenarios.

In the baseline scenarios, individuals receive the standard dose size of vaccines that are 70% and 95% against death and infection.

The most important alternative strategy analyzed considers a scenario where a lower dose regimen allows for an increase in vaccination speed. We project the potential reduction in deaths and infections and the number of extra doses that would be available each month under this strategy.

We also compare various starting points for the vaccination campaign in relation to the peak of the epidemic, and different efficacies against transmission while keeping efficacy against death constant.

All strategies are compared according to the total number of infections and deaths.

## Packages

* **R**: 4.0.5
* **tidyverse**: 1.3.1
* **odin**: 1.0.8
* **socialmixr**: 0.1.8
* **MASS**: 7.3-53.1

## Setup

All results can be replicated through the [run-all-cases.R](run-all-cases.R) script.

The main scripts called are: 

* **[project_setup.R](project_setup.R)**: Loads epidemiological models (ordinary differential equations implemented with odin), auxiliary functions, and initialize parameters.
* * **[immune_response.R](cases/immune_response.R)**: Creates the plot that extends the model in Khoury et al. (2021). The final figure (Figure 1) presented in the paper was generated outside of R, the script here uses the same data for an exploration exercise.
* **[prep-results.R](cases/prep-results.R)**: Estimates the model using the baseline strategy for different vaccination rates and vaccine efficacies
* **[general-example.R](cases/general-example.R)**: Uses estimates from prep-results.R to analyze outcomes under different vaccination rates (Figure S2 and S5).
* **[lower_efficacy_baseline_grid.R](cases/lower_efficacy_baseline_grid.R)**: Uses estimates from prep-results.R to compare the baseline strategies (70% and 95% efficacy) with an alternative that uses a lower dose-regimen at a faster rate but at a potentially lower efficacy (Figure S3).
* * * * **[extra_doses.R](cases/extra_doses.R)**: Auxiliary script used to estimate the additional vaccine supply generated through the adoption of lower doses. This contains only auxiliary calculations and does not generate the final estimates shown in Table S3, which are compiled in a separate spreadsheet.
* * **[dev-pdeath.R](cases/dev-pdeath.R)**: Compares scenarios where vaccination leads to varying efficacy against infection/transmission but fixed efficacy against death, allowing for the comparison of indirect and direct benefits (Figure S4).
* **[delay-impact.R](cases/delay-impact.R)**: Compares scenarios with different starting points for vaccination relative to the peak of transmissions (Figure S6)
* **[generate-figures.R](cases/generate-figures.R)**: Plots the results of the previous scripts.

Three variables from preparation script [project-setup.R](project-setup.R) deserve special attention.

First, the vector `pop` stores the demographic distribution according to age group used in the simulations. The default case uses values compatible with high-income countries. The same applies to the vector `default_pdeath` which contains the probability of death given infection for each age group.

Lastly, `default_supply_ceiling` defines the maximum percentage of the population that can get vaccines due to supply constraints, the default value is 1.
