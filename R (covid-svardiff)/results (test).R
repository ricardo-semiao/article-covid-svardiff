# Packages and Functions --------------------------------------------------
library(stargazer)
library(forecast)
library(vars)
library(varutils)
library(tidyverse)
library(rlang)

theme_set(theme_bw())

source("functions.R")



# Parameters --------------------------------------------------------------
week <- TRUE
trans <- "level"
use_exog <- "season|intramonth"

n.ahead <- 100
runs <- 10
method <- "scoring"



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

start_vaccine <- select(data, starts_with("vaccines")) %>% {which(!is.na(.))[1]}
n <- nrow(data)
vars_names <- c("Casos", "Mortes", "Vacinas")

vars_sds <- data %>%
  select(-date) %>%
  map_dbl(~ sd(., na.rm = TRUE) * ifelse(trans == "log", 100, 1))

gap_size <- 4*2 + 2



# Modeling (VAR) ----------------------------------------------------------
exogen <- create_exogen(data %>% slice(start_vaccine:n), use_exog)

VARselect(select(data[start_vaccine:n,], -date), exogen = exogen)
p <- 3

mod_main <- create_var(data, p, start_vaccine:n, "date", use_exog)



# Modeling (SVAR) ---------------------------------------------------------
# Cholesky:
amat <- length(vars_names) %>% matrix(NA, ., .)
amat[upper.tri(amat)] <- 0

mod_svar <- SVAR(mod_main, Amat = amat, estmethod = method, lrtest = FALSE)

# Residual correlations:
summary(mod_main)$corres %>% stargazer(summary = FALSE, rownames = FALSE, type = "text")
mod_svar$A %>% stargazer(summary = FALSE, rownames = FALSE, type = "text")



# Results: single model ---------------------------------------------------
# ---- Impulse response functions ----
#mod_svar_irf <- irf(mod_svar, ortho = TRUE, n.ahead = n.ahead, cumulative = TRUE, runs = runs)
#ggvar_irf(mod_svar_irf, facet = "ggh4x", args_facet = list(scales = "free_y", independent = "y"))
#
#mod_svar_irf[1:3] <-  mod_svar_irf[1:3] %>%
#  map(function(line){imap(line, ~ .x / vars_sds[.y])})
#
#ggvar_irf(mod_svar_irf, facet = "ggh4x", args_facet = list(scales = "free_y", independent = "y"))
#my.ggsave("pictures/varIRF VCM.png", 16, 10)

# ---- Forecast error variance decomposition ----
#mod_svar_fevd <- fevd(mod_svar, n.ahead = n.ahead)
#ggvar_fevd(mod_svar_fevd)
#my.ggsave("varFEVD.png")



# Results: before & after -------------------------------------------------
exogen2 <- create_exogen(data %>% slice(1:(start_vaccine - 1)), use_exog) #to solve bug in vars

ba_var <- list(
  before = data %>%
    slice(1:(start_vaccine - 1)) %>%
    select(-matches("date|vaccines")) %>%
    VAR(p = p, exogen = exogen2), #to solve bug in vars
  after_none = create_var(data, p, start_vaccine:n, "date|vaccines", use_exog),
  after_vac = create_var(data, p, start_vaccine:n, "date", use_exog)
)

amat2 <- matrix(c(NA, NA, 0, NA), nrow = 2)

max.iter <- 78
ba_svar <- list(
  SVAR(ba_var$before, Amat = amat2, max.iter = max.iter, estmethod = method, lrtest = FALSE),
  SVAR(ba_var$after_none, Amat = amat2, estmethod = method, lrtest = FALSE),
  SVAR(ba_var$after_vac, Amat = amat, estmethod = method, lrtest = FALSE)
) %>%
  set_names(names(ba_var))

# Correlations:
map(ba_var, ~ summary(.x)$corres %>% round(3))
walk(ba_svar, ~ .$A %>% stargazer(summary = FALSE, rownames = FALSE, type = "text"))

# IRF's:
max.iter <- 22 #to solve bug in vars
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

#ggvar_irf2(ba_svar_irf_data)

ba_svar_irf_data[1:3] <- ba_svar_irf_data[1:3]/vars_sds[grepl("cases", names(vars_sds))]
ggvar_irf2(ba_svar_irf_data, scales = "free_y")


# Results: counterfactual -------------------------------------------------
