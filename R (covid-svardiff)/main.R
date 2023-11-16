# Packages and Functions --------------------------------------------------
library(stargazer)
library(vars)
library(varutils)
library(tidyverse)

theme_set(theme_bw())

pvalue.ast = function(x) {
  c("", "*", "**", "***")[pmap_dbl(tibble(x), ~ sum(. < c(0.01, 0.05, 0.10)) + 1)]
}

na_approx <- function(data) {
  numeric_cols <- sapply(data, is.numeric)
  data[numeric_cols] <- map_dfc(data[numeric_cols], function(col) {
    index <- (which(!is.na(col))[1]):length(col)
    col[index] <- zoo::na.approx(col[index], na.rm = FALSE)
    col
  })
  data
}

table_adf <- function(data) {
  data %>%
    filter(date >= date_vaccine$week) %>%
    select(ends_with("l")) %>%
    map_dfr(function(col) {
      col %>%
        tseries::adf.test() %>%
        `[`(c("statistic", "parameter", "p.value"))
    }) %>%
    set_names(c("Estatística", "Lag", "P-valor")) %>%
    mutate(across(-Lag, ~ formatC(.x, digits = 2, format = "f"))) %>%
    mutate(Variável = vars_names, .before = 1) %>%
    stargazer(summary = FALSE, title = "Teste Dickey-Fuller Aumentado", label = "tb:testadf")
}

ccf_ba <- function(data_list, labels) {
  data_ccf <- imap_dfr(data_list, function(df, name) {
    confInterval <- qnorm((1 - 0.95)/2)/sqrt(nrow(df))
    ccf_values <- ccf(df$cases_s, df$deaths_s, lag.max = (start_vaccine), plot = FALSE)
    tibble(
      name = as.factor(name), Lag = ccf_values$lag[,,1],
      value = ccf_values$acf[,,1], ci = confInterval
    )
    })
  
  ggplot(data_ccf, aes(x = Lag, y = value)) +
    geom_area() + 
    geom_ribbon(aes(ymin = -ci, ymax = ci), linetype = 2, fill = NA) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    facet_wrap(vars(name), ncol = 3, labeller = labeller(
      name = set_names(labels, c("Before", "After", "Afterp"))
    )) +
    labs(y = "Valor", x = "Lag (casos)")
}

granger_ba = function(data, p, test, ...){
  data <- map(data, ~ select(.x, !starts_with("vaccines")))
  mods <- map(data, ~VAR(.x, p))
  mods_vcov <- map(mods, ~vcovHC(.x, "HC1"))
  
  data_granger <- map2(mods, mods_vcov, function(mod, vcov) {
    causality(mod, cause = "cases_s", vcov. = vcov)[[test]] %>% 
      `[`(c("statistic", "parameter", "p.value")) %>%
      reduce(c)
  }) %>%
    reduce(rbind) %>%
    as_tibble()
  
  data_granger %>%
    setNames(c("Estatística", "GL 1", "GL 2", "P-valor")) %>%
    mutate(
      Estatística = round(Estatística, 2),
      across(contains("GL"), round),
      `P-valor` = paste(round(`P-valor`, 4), pvalue.ast(`P-valor`))
    ) %>%
    mutate(Variável = c("Antes", "Depois", "Depois+"), .before = 1) %>%
    stargazer(
      summary = FALSE, title = "Causalidade de Granger",
      label = "tb:grangerba", notes = "GL: graus de liberdade",
      rownames = FALSE, ...
    )
}

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




