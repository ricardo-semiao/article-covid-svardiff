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
include <- "season|intramonth" #season|intramonth|waves|intraweek

chol_order <- paste0(c("cases", "vaccines", "deaths"), "_", "s")
gap_size <- if (week) 4*2.5 else 30*2.5

n.ahead <- 52*2
lag.max <- if (week) 4*2.5 else 30*2.5



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

data <- data %>% select(all_of(c("date", chol_order)))

start_vac <- select(data, starts_with("vaccines")) %>% {which(!is.na(.))[1]}
n <- nrow(data)

vars_names <- set_names(
  c("Casos", "Mortes", "Vacinas"),
  data %>% select(-date) %>% colnames()
)

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



# Exploratory Analisys ----------------------------------------------------
# ---- Historic values and stacionarity ----
# Historic values:
ggvar_values(
  select(data, -date), index = data$date,
  args_facet = list(scales = "free_y", ncol = 1, labeller = labeller(serie = vars_names))
) + labs(title = "Valores históricos das séries", x = "Lags", y = "Correlação") +
  ggwaves

# Correlations:
ggvar_acf(select(data, -date) %>% na.omit(),
  args_facet = list(labeller = labeller(serie = vars_names))
) + labs(title = "Auto-correlação das séries", x = "Lags", y = "Correlação")

ggvar_acf(select(data, -date) %>% na.omit(), type = "partial",
  args_facet = list(labeller = labeller(serie = vars_names))
) + labs(title = "Auto-correlação parcial das séries", x = "Lags", y = "Correlação")

ggvar_ccf(select(data, -date) %>% na.omit(),
  args_facet = list(labeller = labeller(var_row = vars_names, var_col = vars_names))
) + labs(title = "Correlação cruzada das séries", x = "Lags", y = "Correlação")

vars_sds <- data %>%
  select(-date) %>%
  map_dbl(~ sd(., na.rm = TRUE) * ifelse(trans == "log", 100, 1))

# ADF tests:
table_adf(data, start_vac, vars_names, type = "text")


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
          list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
          .names = "{.fn}_{.col}"
        )) %>%
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
names_ba <- c("Antes", "Depois", "Depois+")
ccf_ba(data_ba, names_ba, lag.max = lag.max)

# Granger tests B&A:
map(data_ba, ~ VARselect(select(.x, -date), lag.max = lag.max, exogen = create_exogen(.x, include, week))) %>%
  transpose() %>%
  pluck("selection") %>%
  bind_rows() %>%
  mutate(Modelos = names_ba, .before = 1)

granger_ba(map(data_ba, ~ select(.x, -date)), p = 2, test = "Granger", type = "text")
granger_ba(map(data_ba, ~ select(.x, -date)), p = 2, test = "Instant", type = "text")



# Modeling (VAR) ----------------------------------------------------------
VARselect(select(data[start_vac:n,], -date),
  lag.max = lag.max,
  exogen = create_exogen(data %>% slice(start_vac:n), include, week)
)

p <- if (week) 2 else 16

modav_var <- create_var(data, p, start_vac:n, "date", include)

summary(modav_var)



# Diagnostics (VAR) -------------------------------------------------------
# Fitted and residual values:
ggvar_fit(modav_var, index = data$date[-(1:(start_vac - 1))],
  args_facet = list(scales = "free_y", ncol = 1, labeller = labeller(serie = vars_names)),
) + labs(title = "VAR fit", linetype = "Variável") +
  ggwaves

ggvar_values(modav_var, index = data$date[-(1:(start_vac - 1 + p))],
  args_facet = list(scales = "free_y", ncol = 1, labeller = labeller(serie = vars_names))
) + labs(title = "VAR resíduos", y = "Resíduos") +
  ggwaves

logLik(modav_var)

# Residual correlations
ggvar_acf(modav_var, args_facet = list(labeller = labeller(serie = vars_names))) +
  labs(title = "Auto-correlação dos resíduos",y = "Correlação")
ggvar_acf(modav_var, type = "partial", args_facet = list(labeller = labeller(serie = vars_names))) +
  labs(title = "Auto-correlação parcial dos resíduos",y = "Correlação")
ggvar_ccf(modav_var, args_facet = list(labeller = labeller(serie = vars_names))) +
  labs(title = "Correlação cruzada dos resíduos",y = "Correlação")

serial.test(modav_var)
serial.test(modav_var, type = "ES")
residuals(modav_var) %>% as_tibble() %>% map(tseries::adf.test)

# Residual distribution and stability
ggvar_dispersion(modav_var, args_facet = list(scales = "free", labeller = labeller(serie = vars_names))) +
  labs(title = "VAR resíduos (dispersão)", x = "Fit", y = "Resíduos")
arch.test(modav_var)

ggvar_distribution(modav_var, args_facet = list(scales = "free", labeller = labeller(serie = vars_names))) +
  labs(title = "VAR resíduos (distribuição)", x = "Resíduos", y = "Densidade")
normality.test(modav_var)

ggvar_stability(modav_var) #????

# General VAR statistics
roots(modav_var) %>%
  map_chr(~ paste0(round(., 4), " (", ifelse(.<1, "<1", ">=1"), ")")) %>%
  cat(sep = ",\t")

summary(modav_var)$varresult %>%
  map(~ list(`R squared` = .$r.squared, `F test` =.$fstatistic)) %>%
  transpose()



# Modeling (SVAR) ---------------------------------------------------------
modav_svar <- id.chol(modav_var)

# Residual correlations:
summary(modav_var)$corres %>% stargazer(summary = FALSE, rownames = FALSE, type = "text")
modav_svar$B %>% `colnames<-`(rownames(.)) %>% stargazer(summary = FALSE, rownames = FALSE, type = "text")



# Results: single model ---------------------------------------------------
# ---- Impulse response functions ----
modbv_svar_irf <- custom_svars_irf(modav_svar, n.ahead = n.ahead)
plot(modbv_svar_irf)

modbv_svar_irf$irf[,-1] <- modbv_svar_irf$irf[,-1] %>% map2_dfc(rep(vars_sds, 3), ~ .x/.y)
plot(modbv_svar_irf)
#my.ggsave("pictures/varIRF VCM.png", 16, 10)

# ---- Forecast error variance decomposition ----
modbv_svar_fevd <- custom_svars_fevd(modav_svar, n.ahead = n.ahead)
plot(modbv_svar_fevd)
#my.ggsave("varFEVD.png")



include = "season"
mod_var <- list(
  modbn  = create_var(data, 1, 1:(start_vac - 1), "date|vaccines", include, type = "both"),
  modan  = create_var(data, p, start_vac:n, "date|vaccines", include, type = "both"),
  modav  = create_var(data, p, start_vac:n, "date", include, type = "both"),
  modanp = create_var(data, p, (start_vac + gap_size):n, "date|vaccines", include, type = "both"),
  modavp = create_var(data, p, (start_vac + gap_size):n, "date", include, type = "both")
)
# Results: before & after -------------------------------------------------
VARselect(select(data[1:(start_vac - 1),], -matches("date|vaccines")),
  lag.max = lag.max,
  exogen = create_exogen(data %>% slice(1:(start_vac - 1)), include, week)
)



mod_svar <- map(mod_var, id.chol)

# Correlations:
map(mod_var, ~ summary(.x)$corres %>% round(3))
walk(mod_svar, ~ .$B %>% `colnames<-`(rownames(.)) %>% stargazer(summary = FALSE, rownames = FALSE, type = "text"))

# IRF's:
mod_svar_irf <- map(mod_svar, ~ custom_svars_irf(.x, n.ahead = n.ahead))

mod_sds <- list(vars_sds[-2], vars_sds[-2], vars_sds,vars_sds[-2], vars_sds)
mod_svar_irf <- map2(mod_svar_irf, mod_sds, function(x, sds) {
  x$irf[,-1] <- x$irf[,-1] %>% map2_dfc(rep(sds, sqrt(ncol(.))), ~ .x/.y)
  x
})

#plot(mod_svar_irf$modbn); plot(mod_svar_irf$modan); plot(mod_svar_irf$modav)

data_irf <- map_dfc(mod_svar_irf, ~ .x$irf["epsilon[ cases_s ] %->% deaths_s"]) %>%
  set_names(c("Antes", "Depois s/", "Depois c/", "Depois+ s/", "Depois+ c/"))

data_irf %>%
  mutate(Horizon = 1:n.ahead) %>%
  pivot_longer(-Horizon) %>%
  ggplot(aes(Horizon, value)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(name), nrow = 1)

apply(data_irf[1:40,], 2, sum)
apply(data_irf, 2, `[`, 1)

# Results: counterfactual -------------------------------------------------
mod_var$modbn

modpred <- mod_var$modavp
modpred$varresult$Vacinas$coefficients <- length(modpred$varresult$Vacinas$coefficients) %>% numeric()
modpred$datamat <- mutate(modpred$datamat, across(matches("Vacinas"), ~ 0))


# Var select --------------------------------------------------------------
lag.max = 8

mod_select <- list(
  modbn  = select_var(data, lag.max, 1:(start_vac - 1), "date|vaccines", include, type = "both", season = 4),
  modan  = select_var(data, lag.max, start_vac:n, "date|vaccines", include, type = "both", season = 4),
  modav  = select_var(data, lag.max, start_vac:n, "date", include, type = "both", season = 4),
  modanp = select_var(data, lag.max, (start_vac + gap_size):n, "date|vaccines", include, type = "both", season = 4),
  modavp = select_var(data, lag.max, (start_vac + gap_size):n, "date", include, type = "both", season = 4)
)

transpose(mod_select)$selection

mod_select %>%
  imap_dfr(function(sel, mod) {
    crit <- as_tibble(apply(t(sel$criteria), 2, \(x) x/x[1]))
    tibble(lag = 1:lag.max, mod = mod, select(crit, -`FPE(n)`))
  }) %>%
  pivot_longer(-c(lag, mod)) %>%
  ggplot(aes(lag, value, color = name)) +
  geom_line() +
  facet_wrap(vars(mod), scales = "free_y")
