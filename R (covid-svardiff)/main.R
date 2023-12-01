# Packages and Functions --------------------------------------------------
library(stargazer)
library(varutils)
library(patchwork)
library(forecast)
library(vars)
library(svars)
library(rlang)
library(tidyverse)

theme_set(theme_bw())

source("functions.R")
source("functions (package).R")

sci <- function(x) format(x, digits = 1, scientific = TRUE)

custom_ggsave = function(filename, w, h, ...){
  ggsave(paste0("../Latex/Figuras/", filename),
         device = "png", width = w, height = h, units = "cm", ...
  )
}


# Parameters --------------------------------------------------------------
week <- TRUE #TRUE, FALSE
trans <- "level" #level, log, prop
include <- "waves" #season|waves|intramonth|intraweek
type <- "const" #none, const, trend, both
p <- if (week) list(1, 2, 2, 2, 2) else c(8, 16, 16, 16, 16)

vars_order <- paste0(
  c("cases", "vaccines", "deaths"), "_",
  c(level = "s", log = "l", prop = "p")[trans]
)
vars_breaks <- list(
  as.Date(c(sqr = "2021-02-23", none = "2021-08-09")),
  as.Date(c(sqr = "2022-01-03", none = "2022-03-30")),
  as.Date(c(sqr = "2022-05-25", none = "2022-09-03")),
  as.Date(c(sqr = "2022-11-15", none = "2023-01-30"))
)

n_start <- NULL #default to ignore start
n_end <- as.Date("2023-03-20")

gap_size <- if (week) 4*2.5 else 30*2.5
lag.max <- if (week) 2*2.5 else 15*2.5
n.ahead <- round(if (week) 53*0.75 else 365*0.75)

pal <- c("#1B9E77", "#D95F02") #RColorBrewer::brewer.pal(3, "Dark2")
w <- 15.2
h <- 7


# Data Import -------------------------------------------------------------
#path <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"
#data_raw <- read_csv(path) %>% filter(location == "Brazil") %>% filter(row_number() >= which(!is.na(new_deaths))[1])
#saveRDS(data_raw, "data_raw.RDS")

if (week) {
  data <- get_basic_data(readRDS("data_raw.RDS"), trans) %>%
    group_by(week = format(date, "%Y-%U")) %>%
    summarise(
      date = median(date),
      across(-date, sum)
    ) %>%
    select(-week)
} else {
  data <- get_basic_data(readRDS("data_raw.RDS"), trans)
}

data <- data %>%
  filter(date <= n_end) %>%
  select(all_of(c("date", vars_order)))

n_vac <- select(data, starts_with("vaccines")) %>% {which(!is.na(.))[1]}
n_start <- n_start %||% which.min(abs(data$date - as.Date("2020-06-01")))
n <- nrow(data)

namings <- list(
  vars = set_names(
    c("Casos", "Vacinas", "Mortes"),
    vars_order
  ),
  ba_expl = set_names(c("Antes", "Depois", "Depois+")),
  ba_mods = set_names(
    c("Antes (1)", "Depois s/ (2)", "Depois c/ (3)", "Depois+ s/ (4)", "Depois+ c/ (5)"),
    c("bf_no", "af_no", "af_vac", "afp_no", "afp_vac")
  )
)

# Covid waves:
ggwaves <- function(ind = FALSE, date = data$date, breaks = vars_breaks, gap = n_vac + gap_size) {
  map(breaks, function(x) {
    xs <- if (ind) map_vec(x, ~ which.min(abs(date[gap:n] - .x))) else x
    annotate("rect",
      fill = pal[1], alpha = 0.2,
      ymin = -Inf, ymax = Inf, xmin = xs[1], xmax = xs[2]
    )
  })
}

# Variables S.D.'s:
vars_sds <- map_dbl(select(data, -date), ~ sd(., na.rm = TRUE))



# Exploratory Analysis ----------------------------------------------------
# ---- Historic values and stacionarity ----
# Historic values:
ggvar_values(
  select(data, -date), index = data$date,
  args_facet = list(scales = "free_y", ncol = 1, labeller = labeller(serie = namings$vars))
) +
  labs(title = NULL, x = "Lags", y = "Correlação") + #Valores históricos das séries
  scale_x_date(date_breaks = "6 months", date_labels = "%Y/%m") +
  scale_y_continuous(n.breaks = 3, labels = sci) +
  ggwaves()
custom_ggsave("stat_historic.png", w = w, h = h, scale = 2)

# Correlations:
list(
  ggvar_ccf(select(data, -date) %>% na.omit(),
    args_facet = list(labeller = labeller(var_row = namings$vars, var_col = namings$vars)),
    args_ribbon = list(linetype = 2, color = pal[1], fill = NA)
  ) +
    labs(title = NULL, x = NULL, y = "Correlação"),
  ggvar_acf(select(data, -date) %>% na.omit(), type = "partial",
    args_facet = list(labeller = labeller(serie = namings$vars)),
    args_ribbon = list(linetype = 2, color = pal[1], fill = NA)
  ) +
    labs(title = NULL, x = NULL, y = "Correlação")
) %>%
  wrap_plots(nrow = 2, heights = c(3, 1)) +
  plot_annotation(title = NULL) #Correlação cruzada e parcial das séries
custom_ggsave("stat_acf.png", w = w, h = h*1.25, scale = 2)

# ADF tests:
table_adf(data, n_vac, namings$vars, type = "text")


# ---- Seasonality ----
as_ts <- function(x,
                  freq = if (week) 4 else 30, #experiment with different values
                  start = format(data$date, if (week) "%U" else "%d")[1]) {
  ts(x, frequency = freq, start = start)
} 

# Season plots:
iwalk(select(data, -date), function(var, name) {
  graphs <- map(
    list(\(x) autoplot(decompose(x)), ggseasonplot, ggsubseriesplot),
    ~ as_ts(var) %>% .x() + theme(legend.position = "none")
  )
  g <- wrap_plots(graphs) + plot_annotation(title = paste("Seasonality analisys for", name))
  plot(g)
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

data %>%
  mutate(
    season_high = format(date, "%m") %in% c("12", "01", "02"),
    season_mid = format(date, "%m") %in% c("10", "11", "03", "04")
  ) %>% 
  lm(cases_s ~ season_high + season_mid, data = .) %>%
  summary()

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


# ---- Before and after analysis ----
ba_data <- list(
  Before = data[1:(n_vac - 1), ],
  After  = data[n_vac:n, ],
  Afterp = data[(n_vac + gap_size):n, ]
)

ba_data$Before <- select(ba_data$Before, !starts_with("vaccines"))

# CCF B&A:
ccf_ba(ba_data, namings$ba_expl, lag.max = lag.max*3) +
  labs(title = NULL) #Correlação mortes-casos pré e pós vacinas
custom_ggsave("stat_ba.png", w = w, h = h*0.75, scale = 2)

# Granger tests B&A:
map(ba_data, ~ VARselect(select(.x, -date), lag.max = lag.max, exogen = create_exogen(.x, include, week))) %>%
  transpose() %>%
  pluck("selection") %>%
  bind_rows() %>%
  mutate(Modelos = namings$ba_expl, .before = 1)

granger_ba(map(ba_data, ~ select(.x, -date)), p = 2, test = "Granger", type = "text")
granger_ba(map(ba_data, ~ select(.x, -date)), p = 2, test = "Instant", type = "text")



# Modeling ----------------------------------------------------------------
# ---- Lag selection ----
mods_args <- list( #bf = before, af = after, afp = after + gap
  bf_no = list(slice = n_start:(n_vac - 1), remove_match = "date|vaccines"),
  af_no = list(slice = n_vac:n, remove_match = "date|vaccines"),
  af_vac = list(slice = n_vac:n, remove_match = "date"),
  afp_no = list(slice = (n_vac + gap_size):n, remove_match = "date|vaccines"),
  afp_vac = list(slice = (n_vac + gap_size):n, remove_match = "date")
)

mods_select <- map(mods_args, function(args) {
  inject(select_var(data, lag.max, include = include, !!!args, type = type)) #season = 4
})

transpose(mods_select)$selection
ggvar_select(mods_select, lag.max)

mods_args <- map2(mods_args, p, ~ c(.x, p = .y))


# ---- VAR and SVAR ----
mods_var <- map(mods_args, function(args) {
  inject(create_var(data, include = include, !!!args, type = type)) #season = 4
})
walk(mods_var, ~ print(summary(.x)$corres %>% .[1, ncol(.)]))

mods_svar <- map(mods_var, id.chol)
#walk(mods_svar, ~ print(.x$B))



# Diagnostics -------------------------------------------------------------
mod <- mods_var$afp_vac #run this section model by model

# Fitted and residual values:
list(
  ggvar_fit(mod, #index = data$date[-(1:(n_vac - 1))]
            args_facet = list(scales = "free_y", ncol = 1, labeller = labeller(serie = namings$vars)),
  ) +
    labs(title = NULL, x = NULL, y = "Valores", linetype = "Variável") +
    ggwaves(TRUE) +
    scale_y_continuous(n.breaks = 3, labels = sci),
  ggvar_values(mod, #index = data$date[-(1:(n_vac - 1 + p))]
               args_facet = list(scales = "free_y", ncol = 1, labeller = labeller(serie = namings$vars))
  ) +
    labs(title = NULL, y = "Resíduos") +
    ggwaves(TRUE) +
    scale_y_continuous(n.breaks = 3, labels = sci)
) %>%
  wrap_plots(nrow = 2, heights = c(2, 1)) +
  plot_annotation(title = NULL) #Fit e resíduos do VAR
custom_ggsave("diag_fit-res.png", w = w, h = h*1.25, scale = 2)


logLik(mod)

# Residual correlations
list(
  ggvar_ccf(mod,
    args_facet = list(labeller = labeller(var_row = namings$vars, var_col = namings$vars)),
    args_ribbon = list(linetype = 2, color = pal[1], fill = NA)
  ) +
    labs(title = NULL, y = "Correlação"),
  ggvar_acf(mod, type = "partial",
    args_facet = list(labeller = labeller(serie = namings$vars)),
    args_ribbon = list(linetype = 2, color = pal[1], fill = NA)
  ) +
    labs(title = NULL, y = "Correlação")
) %>%
  wrap_plots(nrow = 2, heights = c(3, 1)) +
  plot_annotation(title = NULL) #Correlação cruzada e parcial dos resíduos
custom_ggsave("diag_acf.png", w = w, h = h*1.25, scale = 2)


serial.test(mod, type = "PT.adjusted")
serial.test(mod, type = "ES")
residuals(mod) %>% as_tibble() %>% map(~ tseries::adf.test(.x, k = mod$p)) #urca::ur.df

# Residual distribution and stability
list(
  ggvar_dispersion(mod,
    args_facet = list(scales = "free", labeller = labeller(serie = namings$vars))
  ) +
    labs(title = NULL, x = "Fit", y = "Resíduos"),
  ggvar_distribution(mod,
    args_facet = list(scales = "free", labeller = labeller(serie = namings$vars))
  ) +
    labs(title = NULL, x = "Resíduos", y = "Densidade")
) %>%
  wrap_plots(nrow = 2) +
  plot_annotation(title = NULL) #Dispersão e distribuição dos resíduos
custom_ggsave("diag_disp-dist.png", w = w, h = h, scale = 2)

arch.test(mod)
normality.test(mod)

ggvar_stability(stability(mod, type = "OLS-CUSUM"),
  args_hline = list(linetype = 2, color = pal[1]),
  args_facet = list(labeller = labeller(equation = namings$vars))
) +
  labs(title = NULL, x = "Índice", y = "EFP") #Estabilidade dos resíduos
custom_ggsave("diag_stab.png", w = w, h = h*0.75, scale = 2)

# General VAR statistics
roots(mod) %>%
  map_chr(~ paste0(round(., 4), " (", ifelse(.<1, "<1", ">=1"), ")")) %>%
  cat(sep = ",\t")

summary(mod)$varresult %>%
  map(~ list(`R squared` = .$r.squared, `F test` =.$fstatistic)) %>%
  transpose()



# Results: single model ---------------------------------------------------
# IRF's:
single_mod <- mods_svar$afp_vac
single_irfs <- custom_svars_fun(irf, single_mod, n.ahead = n.ahead)

single_irfs$irf[,-1] <- single_irfs$irf[,-1] %>%
  map2_dfc(rep(vars_sds[names(single_mod$VAR$varresult)], single_mod$K), ~ .x*1000/.y)
plot(single_irfs)

#boot <- custom_svars_fun(wild.boot, single_mod, nboot = 1000)# %>% ba.boot()
#plot(boot)


# FEVD's:
single_fevds <- custom_svars_fun(fevd, single_mod, n.ahead = n.ahead)
plot(single_fevds)



# Results: before & after -------------------------------------------------
# Correlations:
walk2(mods_var, mods_svar, function(var, svar) {
  B_var <- summary(var)$corres %>%
    round(3) %>%
    `colnames<-`(namings$vars[colnames(.)])
  B_svar <- paste(
    cbind(as.vector(svar$B), diag(svar$B)) %>% 
      {`dim<-`(.[,1] / .[,2], dim(B_var))} %>%
      round(3),
    paste0("(", svar$B %>% round(), ")")
  ) %>%
    `dim<-`(dim(B_var)) %>%
    `colnames<-`(colnames(B_var))
  
  cbind(`VAR:` = "", B_var, "|", `SVAR:` = "", B_svar) %>%
    stargazer(summary = FALSE, rownames = FALSE, type = "text")
})

# IRF's:
ba_irfs <- map(mods_svar, ~ custom_svars_fun(irf, .x, n.ahead = n.ahead))

ba_irfs <- map2(
  ba_irfs,
  list(vars_sds[-2], vars_sds[-2], vars_sds, vars_sds[-2], vars_sds),
  function(x, sds) {
    x$irf[,-1] <- x$irf[,-1] %>% map2_dfc(rep(sds, sqrt(ncol(.))), ~ .x*1000/.y)
    x
})

ba_irfs_deaths <- map_dfc(ba_irfs, ~ .x$irf %>% .[grep(".+cases.+deaths.+", names(.))]) %>%
  set_names(namings$ba_mods)

ba_irfs_deaths %>%
  mutate(Horizon = 1:n.ahead) %>%
  pivot_longer(-Horizon) %>%
  mutate(name = factor(name, unique(name))) %>%
  ggplot(aes(Horizon, value)) +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(name), nrow = 1) +
  labs(title = NULL) #IRF casos em mortes - modelos (1) a (5)
custom_ggsave("res1_irfs.png", w = w, h = h*0.75, scale = 2)

ba_sums <- apply(ba_irfs_deaths[1:n.ahead,], 2, sum) %>% round(4)
ba_peak <- apply(ba_irfs_deaths[1:n.ahead,], 2, `[`, 1) %>% round(4)
print(ba_sums)

#
#mods_svar$bf_no %>% wild.boot(design="recursive", distr = "gaussian") %>% plot()



# Results: counterfactual -------------------------------------------------
cf_mods <- list(cf1a3 = mods_var$bf_no, cf2a4 = mods_var$af_vac)

for (x in c("deaths", "cases")) {
  coefs <- cf_mods$cf2a4$varresult[[name(x)]]$coefficients
  coefs[name("vaccines", coefs)] <- 0
  cf_mods$cf2a4$varresult[[name(x)]]$coefficients <- coefs
  
  cf_mods$cf1a3$varresult[[name(x)]]$coefficients <- c(
    cf_mods$cf1a3$varresult[[name(x)]]$coefficients,
    cf_mods$cf2a4$varresult[[name(x)]]$coefficients %>%
      .[!grepl("deaths|cases|vaccines|const|trend", names(.))] #%>% {`if`(x == "deaths", ., .)}
  )
}

cf_mods$cf1a3$datamat <- cbind(
  cf_mods$cf1a3$datamat,
  create_exogen(data, include, week, slice = (n_start + cf_mods$cf1a3$p):(n_vac - 1))
)

cf_preds <- list(
  cf1 = create_cf_1(cf_mods$cf1a3, ci = 0.8),
  cf2 = create_cf_1(cf_mods$cf2a4, ci = 0.8, nnn = 1),
  cf3 = create_cf_2(cf_mods$cf1a3, ci = 0.8),
  cf4 = create_cf_2(cf_mods$cf2a4, ci = 0.8)
)

cf_plots <- map(cf_preds, function(pred) {
  ggvar_predict(pred, data[n_vac:n,],
    args_facet = list(scales = "free_y", labeller = labeller(serie = namings$vars)),
    args_ribbon = list(linetype = 2, color = pal[1], fill = NA)
  ) +
    labs(title = NULL, linetype = "Legenda") +
    scale_linetype_manual(labels = c("Baseline", "Contrafactual"), values = c(1, 2))
})

wrap_plots(cf_plots[1:2], nrow = 2, byrow = TRUE, guides = "collect") +
  plot_annotation(title = NULL) #Contrafactuais (i) e (ii)
custom_ggsave("res2_cf1-2.png", w = w, h = h, scale = 2)

wrap_plots(cf_plots[3:4], nrow = 2, byrow = TRUE, guides = "collect") +
  plot_annotation(title = NULL) #Contrafactuais (i*) e (ii*)
custom_ggsave("res2_cf4-5.png", w = w, h = h, scale = 2)

map(cf_preds, function(pred) {
  pred$fcst[[name("deaths")]][,-4] %>%
    as.data.frame() %>%
  sapply(\(x)  {sum(x - data[[name("deaths")]][n_vac:n])})
})

apply(cf_preds$cf2$fcst$cases_s[,-4], 2, \(x) sum(x - data$cases_s[n_vac:n]))

cf_irfs <- bind_rows(
  tibble(
    index = 1:length(n_vac:n),
    base_true_actual = data$deaths_s[n_vac:n],
    base_true_cf = cf_preds$cf2$fcst[[name("deaths")]][,"fcst"],
    base_irf_actual = data$cases_s[n_vac:n] * ba_sums[5]/1000,
    base_irf_cf = cf_preds$cf2$fcst[[name("cases")]][,"fcst"] * ba_sums[5]/1000
  ) %>%
    pivot_longer(-index, names_sep = "_", names_to = c("line", "base", "type")),
  tibble(
    index = 1:length(n_vac:n),
    cf_actual = data$cases_s[n_vac:n] * ba_sums[1]/1000,
    cf_cf = cf_preds$cf2$fcst[[name("cases")]][,"fcst"] * ba_sums[1]/1000
  ) %>%
    pivot_longer(-index, names_sep = "_", names_to = c("line", "type")) %>%
    map_dfc(~ rep(.x, 2)) %>%
    mutate(base = rep(c("true", "irf"), each = nrow(.)/2))
)

cf_irfs_diff <- cf_irfs %>%
  group_by(base, type) %>%
  reframe(
    diff = value[line == "cf"] - value[line == "base"],
    index = unique(index),
    linetype = 3
  )
 
ggplot(cf_irfs, aes(index, value, linetype = line)) +
  geom_line() +
  geom_segment(
    aes(yend = diff, y = 0, xend = index), color = pal[1],
    linetype = 1, data = cf_irfs_diff
  ) +
  facet_grid(base ~ type,
    labeller = labeller(
      type = c(actual = "Realizado", cf = "Contrafactual"),
      base = c(irf = "IRF", true = "Real")
    )
  ) +
  scale_linetype_manual(labels = c("Baseline", "Contrafactual"), values = c(1, 2)) +
  varutils:::create_sec_axis(xseclab = "Casos", yseclab = "Baseline") +
  labs(title = NULL, x = "Índice", y = "Valores", linetype = "Legenda") #Contrfactuais (iii) e (iii*)
custom_ggsave("res2_cf3-6.png", w = w, h = h, scale = 2)

cf_irfs_diff %>%
  group_by(base, type) %>%
  summarise(sum(diff))

#custom_svars_fun(cf, mods_svar$bf_no, 2) %>% plot()
#custom_svars_fun(cf, mods_svar$afp_vac, 3) %>% plot()

