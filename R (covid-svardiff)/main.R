# Packages and Functions --------------------------------------------------
library(stargazer)
library(forecast)
library(vars)
library(svars)
library(varutils)
library(tidyverse)
library(rlang)

theme_set(theme_bw())

source("functions.R")
source("functions (package).R")



# Parameters --------------------------------------------------------------
week <- TRUE
trans <- "level"
include <- "season" #season|intramonth

chol_order <- paste0(c("cases", "vaccines", "deaths"), "_", "s")
gap_size <- 4*2 + 2

n.ahead <- 20
runs <- 100
method <- "direct"



# Data Import -------------------------------------------------------------
#path <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"
#data_raw <- read_csv(path) %>% filter(location == "Brazil") %>% filter(row_number() >= which(!is.na(new_deaths))[1])
#saveRDS(data_raw, "data_raw.RDS")

data_raw <- readRDS("data_raw.RDS")

data_basic <- get_basic_data(data_raw, trans)

if (week) {
  data <- data_basic %>%
    group_by(week = format(date, "%Y-%U")) %>%
    summarise(
      date = median(date),
      across(-date, sum)
    ) %>%
    select(-week)
} else {
  data <- data_basic
}

data <- data %>% select(c("date", chol_order))

start_vac <- select(data, starts_with("vaccines")) %>% {which(!is.na(.))[1]}
n <- nrow(data)
vars_names <- c("Casos", "Mortes", "Vacinas")



# Exploratory Analisys ----------------------------------------------------
# ---- Historic values and stacionarity ----
# Historic values:
ggvar_values(
  select(data, -date), index = data$date,
  args_facet = list(scales = "free_y", ncol = 1)
)

# Correlations:
ggvar_acf(select(data, -date) %>% na.omit())
ggvar_acf(select(data, -date) %>% na.omit(), type = "partial")
ggvar_ccf(select(data, -date) %>% na.omit())

vars_sds <- data %>%
  select(-date) %>%
  map_dbl(~ sd(., na.rm = TRUE) * ifelse(trans == "log", 100, 1))

# ADF tests:
table_adf(data, start_vac, vars_names, type = "text")

# Covid waves and seasons:
breaks_waves <- as.Date(c(
  start1 = "2021-02-21", end1 = "2021-06-22",
  start2 = "2022-01-03", end2 = "2022-03-12"
))

ggwaves <- imap(breaks_waves, function(x, name) {
  geom_vline(
    xintercept = as.Date(x), linetype = 2, color = c("green", "red")[1 + grepl("2", name)]
  )
})

ggvar_values(
  select(data, -date), index = data$date,
  args_facet = list(scales = "free_y", ncol = 1)
) + ggwaves
  


# ---- Seasonality ----
ts(format(data$date, "%U"))
as_ts <- function(x, freq = 13, start = 8) ts(x, frequency = freq, start = start)

# Season plots:
walk(list(decompose, ggseasonplot, ggsubseriesplot), function(fun){
  walk(select(data, -date), ~ as_ts(.x) %>% fun() %>% plot())
})

# Cycle averages:
walk(
  exprs(
    cycle(as_ts(date, 8, 4)),
    format(date, "%b"),
    as.numeric(format(data$date, "%U")) %% 2
    ), function(expr) {
      data %>%
        group_by(!!expr) %>%
        summarise(across(where(is.numeric),
          list(mean = mean, sd = sd), na.rm = TRUE, .names = "{.fn}_{.col}")
        ) %>%
        print()
})

# Season tests
walk(list(ocsb.test, uroot::hegy.test), function(fun){ #uroot::ch.test
  walk(select(data, -date), ~ as_ts(.x) %>% fun() %>% print)
})

# Auto arima
walk(select(data, -date), function(x) {
  auto.arima(as_ts(x)) %>%
    `[[`("arma") %>%
    set_names(c("AR", "MA", "SAR", "SMA", "P", "D", "SD")) %>%
    print()
})


# ---- Before and after analisys ----
data_ba <- list(
  Before = data[1:(start_vac - 1), ],
  After  = data[start_vac:n, ],
  Afterp = data[(start_vac + gap_size):n, ]
)

data_ba$Before <- select(data_ba$Before, !starts_with("vaccines"))

# CCF B&A:
ccf_ba(data_ba, c("Antes", "Depois", "Depois+"), lag.max = 12)

# Granger tests B&A:
map(data_ba, ~ select(.x, -date) %>% VARselect(lag.max = 10)) %>% transpose() %>% pluck("selection")

granger_ba(map(data_ba, ~ select(.x, -date)), p = 2, test = "Granger", type = "text")
granger_ba(map(data_ba, ~ select(.x, -date)), p = 2, test = "Instant", type = "text")



# Modeling (VAR) ----------------------------------------------------------
exogen <- create_exogen(data %>% slice(start_vac:n), include)

VARselect(select(data[start_vac:n,], -date), exogen = exogen)
p <- 2

main_var <- create_var(data, p, start_vac:n, "date", include)



# Model Diagnostics -------------------------------------------------------
summary(main_var)

# Fitted and residual values:
ggvar_fit(main_var, index = data$date[-(1:(start_vac - 1))],
  args_facet = list(scales = "free_y", ncol = 1),
) + ggwaves

ggvar_values(main_var, index = data$date[-(1:(start_vac - 1 + p))],
  args_facet = list(scales = "free_y", ncol = 1)
) + ggwaves

logLik(main_var)

# Residual correlations
ggvar_acf(main_var)
ggvar_acf(main_var, type = "partial")
ggvar_ccf(main_var)

serial.test(main_var)
serial.test(main_var, type = "ES")
residuals(main_var) %>% as_tibble() %>% map(tseries::adf.test)

# Residual distribution and stability
ggvar_dispersion(main_var, args_facet = list(scales = "free"))
arch.test(main_var)

ggvar_distribution(main_var, args_facet = list(scales = "free"))
normality.test(main_var)

ggvar_stability(main_var)

# General VAR statistics
roots(main_var) %>%
  map_chr(~ paste0(round(., 4), " (", ifelse(.<1, "<1", ">=1"), ")")) %>%
  cat(sep = ",\t")

summary(main_var)$varresult %>%
  map(~ list(`R squared` = .$r.squared, `F test` =.$fstatistic)) %>%
  transpose()



# Modeling (SVAR) ---------------------------------------------------------
# Cholesky:
amat <- length(vars_names) %>% matrix(NA, ., .)
amat[upper.tri(amat)] <- 0

main_svar <- SVAR(main_var, Amat = amat, estmethod = method, lrtest = FALSE)
#main_svar <- id.chol(main_var)

# Residual correlations:
summary(main_var)$corres %>% stargazer(summary = FALSE, rownames = FALSE)
mod_svar$A %>% stargazer(summary = FALSE, rownames = FALSE)

logLik(mod_svar)



# Results: single model ---------------------------------------------------
# ---- Impulse response functions ----
main_svar_irf <- irf(main_svar, n.ahead = n.ahead)
ggvar_irf(mod_svar_irf, facet = "ggh4x", args_facet = list(scales = "free_y", independent = "y"))
#plot(main_svar_irf)

main_svar_irf[1:3] <-  mod_svar_irf[1:3] %>%
  map(function(line){imap(line, ~ .x / vars_sds[.y])})

ggvar_irf(main_svar_irf, facet = "ggh4x", args_facet = list(scales = "free_y", independent = "y"))
#my.ggsave("pictures/varIRF VCM.png", 16, 10)

# ---- Forecast error variance decomposition ----
main_svar_fevd <- fevd(main_svar, n.ahead = n.ahead)
ggvar_fevd(main_svar_fevd)
plot(main_svar_fevd)
#my.ggsave("varFEVD.png")


# Results: before & after -------------------------------------------------
exogen2 <- create_exogen(data %>% slice(1:(start_vac - 1)), include) #to solve bug in vars

ba_var <- list(
  before = data %>%
    slice(1:(start_vac - 1)) %>%
    select(-matches("date|vaccines")) %>%
    VAR(p = p, exogen = exogen2), #to solve bug in vars
  after_none = create_var(data, p, start_vac:n, "date|vaccines", include),
  after_vac = create_var(data, p, start_vac:n, "date", include)
)

amat2 <- matrix(c(NA, NA, 0, NA), nrow = 2)

max.iter <- 60
ba_svar <- list(
  SVAR(ba_var$before, Amat = amat2, max.iter = max.iter, estmethod = method, lrtest = FALSE),
  SVAR(ba_var$after_none, Amat = amat2, estmethod = method, lrtest = FALSE),
  SVAR(ba_var$after_vac, Amat = amat, estmethod = method, lrtest = FALSE)
) %>%
  set_names(names(ba_var))

#ba_svar <- list(
#  id.chol(ba_var$before),
#  id.chol(ba_var$after_none),
#  id.chol(ba_var$after_vac)
#) %>%
#  set_names(names(ba_var))

# Correlations:
map(ba_var, ~ summary(.x)$corres %>% round(3))
walk(ba_svar, ~ .$A %>% stargazer(summary = FALSE, rownames = FALSE))

# IRF's:
max.iter <- 8 #to solve bug in vars
ba_svar_irf <- map(ba_svar, function(mod) {
  ba_irf <- irf(mod, n.ahead = n.ahead, cumulative = TRUE, runs = runs)
  ba_irf[1:3] <-  ba_irf[1:3] %>% map(function(line){imap(line, ~ .x / vars_sds[.y])})
})

ba_svar_irf_data <- ba_svar_irf %>% 
  imap(function(irf, name) {
    irf_values <- irf %>%
      transpose() %>%
      pluck(grep("cases", names(.), value = TRUE)) %>%
      map_dfc(~ .x[,grepl("deaths", colnames(.x))])
    tibble(irf_values, model = name, horizon = 1:nrow(irf_values))
  }) %>%
  bind_rows() %>%
  mutate(model = factor(model, c("before", "after_none", "after_vac")))

ggvar_irf2(ba_svar_irf_data)

ba_svar_irf_data[1:3] <- ba_svar_irf_data[1:3]/vars_sds[grepl("cases", names(vars_sds))]
ggvar_irf2(ba_svar_irf_data)

#map_dfc(ba_svar, ~ irf(.x)$irf["epsilon[ cases_s ] %->% deaths_s"]) %>%
#  mutate(Horizon = 1:20) %>%
#  pivot_longer(-Horizon) %>%
#  ggplot(aes(Horizon, value)) +
#  geom_line() +
#  facet_wrap(vars(name), nrow = 1)

# Results: counterfactual -------------------------------------------------
