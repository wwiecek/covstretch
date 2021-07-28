myfile <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"
csv <- read_csv(myfile, col_types = cols(
  total_vaccinations = col_number(),
  people_vaccinated = col_number(),
  total_vaccinations_per_hundred         = col_number(),
  people_vaccinated_per_hundred          = col_number(),
  new_vaccinations_smoothed              = col_number(),
  new_vaccinations_smoothed_per_million  = col_number(),
  new_vaccinations                       = col_number(),
  people_fully_vaccinated                = col_number(),
  people_fully_vaccinated_per_hundred    = col_number()
)) %>%
  select(location, date, new_deaths_smoothed_per_million, new_cases_smoothed_per_million,
         people_vaccinated_per_hundred, people_fully_vaccinated_per_hundred)

saveRDS(csv, file = "counterfactuals/owin-data-formatted.rds")