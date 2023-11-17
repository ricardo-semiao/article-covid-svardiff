# Packages and Functions --------------------------------------------------
library(stargazer)
library(vars)
library(varutils)
library(tidyverse)

theme_set(theme_bw())

source("functions.R")

# Data Import -------------------------------------------------------------
#path <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"
#data_raw <- read_csv(path) %>% filter(location == "Brazil") %>% filter(row_number() >= which(!is.na(new_deaths))[1])

data_daily <- data_raw %>%
  transmute(
    date = date,
    cases_s = new_cases_smoothed,
    deaths_s = new_deaths_smoothed,
    vaccines_s = new_vaccinations_smoothed,
    cases_p = new_cases_smoothed_per_million,
    deaths_p = new_deaths_smoothed_per_million,
    vaccines_p = new_vaccinations_smoothed_per_million
  ) %>%
  filter(
    !if_all(-date, is.na) & date < as.Date("2023-03-23") & cases_s + deaths_s > 0
  ) %>%
  na_approx() %>%
  mutate(across(ends_with("_s"), log, .names = "{str_remove(.col, '_s')}_l"))

data_weekly <- data_daily %>%
  group_by(week = format(date, "%Y-%U")) %>%
  summarise(
    date = median(date),
    across(-date, sum)
  ) %>%
  select(-week)

start_vaccine <- which(!is.na(data_weekly$vaccines_s))[1]

#start_vaccine <- list(
#  day = which(!is.na(data_daily$vaccines_s))[1],
#  week = which(!is.na(data_weekly$vaccines_s))[1]
#)
#
#date_deaths <- list(
#  day = which(data_daily$deaths_s > 0)[1]],
#  week = which(data_weekly$deaths_s > 0)[1]
#)

n <- list(day = nrow(data_daily), week = nrow(data_weekly))
vars_names <- c("Casos", "Mortes", "Vacinas")

# Exploratory Analisys ----------------------------------------------------
ggvar_values(
  select(data_weekly, ends_with("l")), index = data_weekly$date,
  args_facet = list(scales = "free_y")
)

ggvar_acf(data_weekly %>% select(ends_with("l")) %>% na.omit())
ggvar_acf(data_weekly %>% select(ends_with("l")) %>% na.omit(), type = "partial")
ggvar_ccf(data_weekly %>% select(ends_with("l")) %>% na.omit())

table_adf(data_weekly)

gap_size <- 12

data_ba <- list(
  Before = data_weekly[1:(start_vaccine - 1), ],
  After = data_weekly[(start_vaccine):nrow(data_weekly), ],
  Afterp = data_weekly[(start_vaccine + gap_size):nrow(data_weekly), ]
)

#data_ba$Before <- select(data_ba$Before, !starts_with("vaccines"))

ccf_ba(data_ba, c("Antes", "Depois", "Depois+"))

map(data_ba, ~ select(.x, ends_with("_s")) %>% VARselect(lag.max = 8))
map(data_ba, ~ select(.x, ends_with("_s"))) %>% granger_ba(p = 2, test = "Instant", type = "text")

data_weekly$cases_l %>% ts(frequency = 4*2) %>% decompose("additive") %>% plot()
data_weekly$deaths_p %>% ts(frequency = 4) %>% decompose() %>% plot()
data_weekly$vaccines_p %>% ts(frequency = 4) %>% decompose() %>% plot()



# Modeling ----------------------------------------------------------------
exogen <- data_weekly %>%
  slice(start_vaccine:n$week) %>%
  transmute(
    season_high = format(date, "%m") %in% c("11", "12", "01", "02", "03", "04")
  )

VARselect(data_weekly[start_vaccine:n$week,] %>% select(ends_with("s")),
          exogen = exogen, type = "both")

mod_main <- VAR(
  data_weekly[start_vaccine:n$week,] %>% select(ends_with("s")),
  p = 2, exogen = exogen
)

ggvar_fit(mod_main, args_facet = list(scales = "free_y", ncol = 1))
ggvar_values(mod_main, args_facet = list(scales = "free_y", ncol = 1))
logLik(mod_main)

ggvar_acf(mod_main)
ggvar_ccf(mod_main)

ggvar_dispersion(mod_main, args_facet = list(scales = "free"))
ggvar_distribution(mod_main, args_facet = list(scales = "free"))
ggvar_stability(mod_main)

serial.test(mod_main)
serial.test(mod_main, type = "ES")
residuals(mod_main) %>% as_tibble() %>% map(tseries::adf.test)

arch.test(mod_main)
normality.test(mod_main)

roots(mod_main) %>%
  map_chr(~ paste0(round(., 4), " (", ifelse(.<1, "<1", ">=1"), ")")) %>%
  cat(sep = ",\t")

summary(mod_main)$varresult %>%
  map(~ list(`R squared` = .$r.squared, `F test` =.$fstatistic)) %>%
  transpose()

#------------------  SVAR           ------------------
amat <- length(vars_names) %>% matrix(NA, ., .)
amat[upper.tri(amat)] <- 0

mod_svar <- SVAR(mod_main, "scoring", Amat = amat, lrtest = FALSE)

# -- Coefficients and statistics
summary(mod_main)$corres %>% stargazer(summary = FALSE, rownames = FALSE)
mod_svar$A %>% stargazer(summary = FALSE, rownames = FALSE)
#summary(mod_svar)
#logLik(mod_svar)

