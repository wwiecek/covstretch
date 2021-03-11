# COVID-19 Vaccination - Epidemiological Simulations

## Introduction

This repository includes the scripts used to generate the results in "Stretching Vaccines: Optimal Vaccine Use in a Pandemic".

We implement a SEIR (Susceptible, Exposed, Infectious, Recover) model and estimate the health outcomes under different vaccination scenarios.

The most important scenarios studied include a baseline strategy, a first-doses-first strategy, and a scenario where a low-efficacy vaccine/dose regimen allows for an increase in vaccination speed.

In the baseline strategy, individuals that received the first dose take priority over unvaccinated ones. Alternatively, for the first-doses-first strategy, the priority is to vaccinate everyone with the first dose before moving to the second.

All strategies are compared according to the final number of infections and deaths.

## Packages

* **R**: 4.0.2
* **tidyverse**: 1.3.0
* **odin**: 1.0.8
* **socialmixr**: 0.1.8

## Setup

All results can be replicated through the [run-all.R](run-all.R) script.

The main scripts called are, in order: 

* **[setup.R](setup.R)**: Loads epidemiological models (ordinary differential equations implemented with odin), auxiliary functions, and initialize parameters.
* **[fdf-prep-delta.R](cases/fdf-prep-delta.R)**: Calibrates the vaccination speed in the first-dose-first scenario to match the total number of doses from baseline.
* **[fdf-results.R](cases/fdf-results.R)**: Estimates the model for the first-doses-first and a hybrid strategy and compares the results with baseline (article section 7).
* **[prep-results.R](cases/prep-results.R)**: Estimates the model using the baseline strategy for different vaccination rates and vaccine efficacies
* **[general-example.R](cases/general-example.R)**: Uses estimates from prep-results.R to analyze the effects of different vaccination rates and vaccine efficacies (article section 5).
* **[lower-efficacy.R](cases/lower-efficacy.R)**: Uses estimates from prep-results.R to compare the baseline strategy with an alternative that uses a lower-efficacy vaccine at a faster rate.
* **[lower-efficacy-delay.R](cases/lower-efficacy.R)**: Estimates the model for the case where a lower-efficacy vaccine is immediately available and compares it to the baseline scenario with a delay on vaccine availability.
* **[kappa-impact.R](cases/kappa-impact.R)**: Estimates the model allowing for loss of immunity for those vaccinated at a given rate and compares the results with the case without loss of immunity.
* **[generate-figures.R](cases/generate-figures.R)**: Plots the results of the previous scripts.

Three variables from [run-all.R](run-all.R) deserve special attention.

First, the vector `pop` stores the demographic distribution according to age group used in the simulations. The default case uses values compatible with high-income countries.

Second, `kappa_default` defines the rate at which loss of immunity happens, the default value is 0.

Lastly, `default_supply_ceiling` defines the maximum percentage of the population that can get vaccines due to supply constraints, the default value is 1.