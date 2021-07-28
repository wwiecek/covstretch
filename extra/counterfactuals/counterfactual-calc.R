library(readr)
CovidData_UK_final <- read_csv("counterfactuals/CovidData_UK_final.csv")
CovidData_US_final <- read_csv("counterfactuals/CovidData_US_final.csv")

head(CovidData_US_final)
summary(CovidData_US_final)


# UNITED STATES
# Values extracted from CDC source on 24 February 2021
# https://covid.cdc.gov/covid-data-tracker/#vaccination-demographic
# this is distribution across different age groups
c(18, 30, 40, 50, 65, 75)
c(0.1, 7, 10, 10.5, 18.3, 25.9, 28.3)/100
c(01., 6.6, 9.5, 10.1, 18.9, 28, 26.7)/100
sum(ifr*c(0, 0, 6.6, 9.5, 10.1, 18.9/1.5, 18.9/3 + 28/2, 28/2+26.7/2, 26.7/2)/100)


# 55% of vaccines went to over 65

# ENGLAND:
# Excel spreadsheet from NHS, extracted from 
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/
# this is % of age group vaccinates
c(70, 75, 80)
c(5.7, 2.56, 1.93, 2.65)/c(37.9, 2.78, 1.94, 2.83) #the denominator for 1st group is without under 16s
# 55% of vaccines went to over 70's:
c(5.7, 2.56+1.93+2.65)
c(5.7, 2.56, 1.93, 2.65)/sum(c(5.7, 2.56, 1.93, 2.65))



library(tidyverse)
CovidData_UK_final %>% 
  filter(!is.na(days_apprvl)) %>%
  select(country, date_td, first_dose, second_dose) %>% View

share_people_vaccinated_covid <- read_csv("counterfactuals/share-people-vaccinated-covid.csv")
people_fully_vaccinated_covid <- read_csv("counterfactuals/share-people-fully-vaccinated-covid.csv")
filter(share_people_vaccinated_covid, Entity %in% c("United Kingdom", "United States"))
filter(people_fully_vaccinated_covid, Entity %in% c("United Kingdom", "United States"))

# If there is no transmission, then
# Mortality risk = (daily infection risk) x IFR x (length of time unprotected)

load("data/default_inputs.Rdata")
us_pop <- as.numeric(pbc_spread["840",])
uk_pop <- as.numeric(pbc_spread["826",])
# https://github.com/1DaySooner/RiskModel/blob/master/ifr-model/bayesian_ifr_model.pdf
ifr <- exp(1.12*(1:9 - 3) - 8.80)
us_ifr <- sum(ifr*us_pop)/sum(us_pop)
uk_ifr <- sum(ifr*uk_pop)/sum(uk_pop)
ifr_all <- .01

efficacy1 <- .8 #efficacy following 1st dose
efficacy2 <- .95 #efficacy following 2nd dose

# IFR reweighted due to actual vaccination policy, likely an underestimate:
ifr_rew <- .55*sum(ifr[8:9]*uk_pop[8:9]) / sum(uk_pop[8:9]) + 
  .45*sum(ifr[3:7]*uk_pop[3:7]) / sum(uk_pop[3:7])

# Estimated number of people who had Covid at different points in time
owid_covid_data <- read_csv("counterfactuals/owid-covid-data.csv")
df <- owid_covid_data %>% 
  filter(location %in% c("United Kingdom", "United States")) %>% 
  filter(date <= "2021-02-23") %>%
  select(location, date, new_deaths_smoothed_per_million, new_cases_smoothed_per_million) %>%
  left_join(share_people_vaccinated_covid %>% rename(location = Entity, date = Date), by = c("location", "date")) %>%
  left_join(people_fully_vaccinated_covid %>% rename(location = Entity, date = Date), by = c("location", "date")) %>%
  select(-Code.x, -Code.y) %>% 
  group_by(location) %>%
  # Count number of infections to date
  mutate(new_infections_per_milion = lead(new_deaths_smoothed_per_million, 30)/uk_ifr) %>% 
  mutate(case_to_inf = new_infections_per_milion[date == "2021-01-24"]/new_cases_smoothed_per_million[date == "2021-01-24"]) %>%
  mutate(new_infections_per_milion = ifelse(is.na(new_infections_per_milion), 
                                            case_to_inf*new_cases_smoothed_per_million, 
                                            new_infections_per_milion)) %>%
  # mutate(new_infections_per_milion = zoo::na.locf(new_infections_per_milion, na.rm=F)) %>% 
  mutate(infected_to_date_per_milion = cumsum(new_infections_per_milion)) %>% 
  # linear approximation for fraction vaccinated:
  mutate(people_vaccinated_per_hundred = zoo::na.approx(people_vaccinated_per_hundred, na.rm=FALSE)) %>%
  mutate(people_vaccinated_per_hundred = 
           ifelse(is.na(people_vaccinated_per_hundred), 0, people_vaccinated_per_hundred)) %>%
  mutate(people_fully_vaccinated_per_hundred = zoo::na.approx(people_fully_vaccinated_per_hundred, na.rm=FALSE)) %>%
  mutate(people_fully_vaccinated_per_hundred = 
           ifelse(is.na(people_fully_vaccinated_per_hundred), 0, people_fully_vaccinated_per_hundred)) %>%
  # Susceptible proportion
  mutate(susceptible_proportion = 1 - 
           efficacy1*(people_vaccinated_per_hundred-people_fully_vaccinated_per_hundred)/100 -
           efficacy2*(people_fully_vaccinated_per_hundred)/100 -
           infected_to_date_per_milion/1e06) %>%
  mutate(new_1st_doses = people_vaccinated_per_hundred/100 - lag(people_vaccinated_per_hundred, 1, default = 0)/100) %>%
  mutate(new_2nd_doses = people_fully_vaccinated_per_hundred/100 - lag(people_fully_vaccinated_per_hundred, 1, default = 0)/100) %>%
  mutate(infection_risk_from_t = (max(infected_to_date_per_milion, na.rm=T) - infected_to_date_per_milion) / 1e06) %>%
  mutate(deaths_averted_by_vaccinations_per_milion = 
           1e06*efficacy1*new_1st_doses*ifr_rew*infection_risk_from_t,
         deaths_averted_by_2nd_doses_per_milion = 
           1e06*(efficacy2 - efficacy1)*new_2nd_doses*ifr_rew*infection_risk_from_t,
         deaths_averted_by_2nd_doses_ctfl_per_milion = 
           1e06*efficacy1*new_2nd_doses*ifr_rew*infection_risk_from_t)

df %>% ggplot(aes(x = date, y = new_infections_per_milion, color = location)) + geom_line()
df %>% ggplot(aes(x = date, y = new_deaths_smoothed_per_million, color = location)) + geom_line()

df %>% 
  ggplot(aes(x = date, y = infection_risk_from_t, color = location)) + geom_line()

df %>% 
  ggplot(aes(x = date, y = deaths_averted_by_vaccinations_per_milion, color = location)) + geom_line()

df %>%
  summarise(
    d1 = sum(new_1st_doses),
    d2 = sum(new_2nd_doses),
    v1 = sum(deaths_averted_by_vaccinations_per_milion, na.rm=T),
    b1 = sum(deaths_averted_by_2nd_doses_per_milion, na.rm=T),
    b2 = sum(deaths_averted_by_2nd_doses_ctfl_per_milion, na.rm=T)
    ) %>%
  mutate(pop = c(55, 331)) %>%
  mutate(v1 = pop*v1, b1 = pop*b1, b2 = pop*b2) %>%
  mutate(total_b = v1+b1, total_b_ast = v1+b2, rr = (v1+b2)/(v1+b1))
  
