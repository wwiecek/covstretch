library(tidyverse)

df <- read_csv('data/071221-unicef_vaccines.csv')
df <- df[!is.na(df$vaccine),]
df[ df == "(Blank)" ] <- NA

unique(df$vaccine)

df[df$vaccine %in% c("SII - Covishield","AstraZeneca - Vaxzevria","Vaxzevria"), "vaccine"] <- "AstraZeneca"

df %>% select(vaccine,total) %>% group_by(vaccine) %>% mutate(total=sum(total)) %>% unique()

df %>% select(country,total) %>% group_by(country) %>% mutate(total=sum(total)) %>% unique()

vax_rate <- read_csv('data/071221-world_vaccinations.csv')
head(vax_rate)

vax_rate <- vax_rate %>% mutate(date= as.Date(date,format = "%m/%d/%y"))
head(vax_rate)

vax_rate["year_month"] <- format(vax_rate$date,"%y-%m")

total_vax <- sum(vax_rate$vaccinated)

vax_rate <- vax_rate %>% select(year_month,vaccinated) %>% group_by(year_month) %>% mutate(vax_mean=sum(vaccinated)/total_vax) %>% 
  select(vax_mean,year_month) %>% unique()

vaccines <- c("Pfizer BioNTech - Comirnaty","Moderna - mRNA-1273","Janssen - Ad26.COV 2.S","AstraZeneca")
df_vax <- df %>% filter(vaccine %in% vaccines) %>% select(vaccine,total) %>% group_by(vaccine) %>% mutate(total=sum(total)) %>% unique()

df_vax_month <- data.frame()
for (ym in vax_rate$year_month){
  prop <- vax_rate[vax_rate$year_month==ym,]$vax_mean
  df.tmp <- (df_vax %>% mutate(total=total*prop))
  df.tmp["year_month"] <- ym
  df.tmp<-spread(df.tmp,vaccine,total)
  df_vax_month <- rbind(df_vax_month,df.tmp)
}

(df_vax_month[df_vax_month$year_month=="21-06",] %>% select(-year_month))*6
df_vax_month[df_vax_month$year_month=="21-06",] %>% select(-year_month)
