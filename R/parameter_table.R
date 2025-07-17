library(gt)
library(dplyr)
library(knitr)
library(gtExtras) # for two-column table
library(tidyverse)

# label
human_parms <- data.frame(
  label = c("nu_H", "nu_TH", "beta_HA", "beta_HB", "epsilon_H", "alpha_H",
            "rho_RH", "rho_VH", "b_H", "mu_H", "gamma_H", "theta_H"),

  description = c("human vaccination rate", "human vaccine effectiveness",
                  "transmission rate from humans to Culex",
                  "transmission rate from humans to Aedes",
                  "human incubation period", "human infectious period",
                  "waning immunity from natural infection",
                  "waning immunity from vaccination",
                  "vaccine booster rate", "death rate", "birth rate",
                  "delay to full vaccine effectiveness"),

  point_estimate = c("--", "90%", "0.81*", "0.89*", "4 days*", "4 days*",
                     "900 days*", ">5 years", "5 years", "70 years^",
                     "1.5 children/person^", "14 days"),

  range = c("[0 - 100]%", "[85-95]%", "0.78-0.84*", "0.77-1.00*",
            "[2-6] days^", "[2-7] days^", "", "", "",
            "40-90 years^", "", ""))


human_parms <- rbind(human_parms, NA)


livestock_parms <- data.frame(
  label = c("nu_M", "nu_TM", "beta_MA", "beta_MB", "epsilon_M", "alpha_M",
            "rho_RM", "rho_VM", "b_M", "mu_M", "gamma_1", "gamma2", "theta_M"),

  description = c("livestock vaccination rate",
                  "livestock vaccine effectiveness",
                  "transmission rate from livestock to Culex",
                  "transmission rate from livestock to Aedes",
                  "livestock incubation period", "livestock infectious period",
                  "waning immunity from natural infection",
                  "waning immunity from vaccination",
                  "vaccine booster rate", "death rate",
                  "birth rate (noninfected livestock)",
                  "birth rate (infected livestock)",
                  "delay to full vaccine effectiveness"),

  point_estimate = c("--", "90%", "0.81*", "0.89*", "42 hours*", "4 days*",
                     "900 days",
                     "", "", "5 years^", "1.5/lifetime*", "0.1*1.5/lifetime*",
                     "14 days"),

  range = c("[0-100]%", "[85-100]%", "[78-84]%*", "[77-100]%*", "12-72 hours*",
            "1-7 days^", "", "", "", "1-10 years^", "0-2.3/lifetime^", "", ""))

culex_parms <- data.frame(
  label = c("beta_AH", "beta_AM", "tau_A", "zeta_A", "omega_A", "omegaA11",
            "kappaA", "mu_PA", "mu_SA", "mu_IA", "mu_QA"),

  description = c("transmission from Culex to humans",
                  "transmission from Culex to livestock", "Culex hatch rate",
                  "proportion of eggs infected by vertical transmission",
                  "Culex egg lay rate", "??", "Culex carrying capacity",
                  "mortality rate of noninfected eggs",
                  "mortality rate of noninfected adults",
                  "mortality rate of infectious adults",
                  "mortality rate of infected eggs"),

  point_estimate = c("0.07*", "0.07*", "0.2*", "0.25*", "25/day*", "", "",
                     "0.002*", "10 days*", "10 days*", "0.002*"),

  range = c("0-0.10^", "0.001-0.14^", "0.01-0.33^", "0.01-0.5^", "", "",
            "1,750-17,500*", "", "3-60 days^", "3-60 days^", ""))

aedes_parms <- data.frame(
  label = c("beta_BH", "beta_BM", "tau_B", "zeta_B", "omega_B", "omegaB11",
            "kappa_B", "mu_PB", "mu_SB", "mu_IB", "mu_QB"),

  description = c("transmission from Aedes to humans",
                  "transmission from Aedes to livestock", "Aedes hatch rate",
                  "proportion of eggs infected by vertical transmission",
                  "Aedes egg lay rate", "??", "Aedes carrying capacity",
                  "mortality rate of noninfected eggs",
                  "mortality rate of noninfected adults",
                  "mortality rate of infectious adults",
                  "mortality rate of infected eggs"),

  point_estimate = c("0.01*", "0.01*", "0.2*", "0.25*", "15*", "", "175,000*",
                     "1 year*", "4 days*", "4 days*", "0.0025*"),

  range = c("0.000046-0.02^", "0.001-0.54^", "0.01-0.33^", "0.0-0.5^",
            "0.05-200^", "", "15,000-1e9^", "0.00001-0.005*", "3-60 days^",
            "3-60 days^", ""))


tab1 <- rbind(human_parms, culex_parms)
gt_table1 <-
  tab1 |>
  select(label, description, point_estimate, range) |>

  gt(rowname_col = "label") |>
  cols_label(label = "Parameter",
             description = "Description",
             point_estimate = "Point Estimate",
             range = "Range") |>

  tab_row_group(label = "Culex",
                rows = 14:23) |>
  tab_row_group(label = "Humans",
                rows = 1:13) |>

  tab_stubhead(label = "Parameter")

tab2 <- rbind(livestock_parms, aedes_parms)
gt_table2 <-
  tab2 |>
  select(label, description, point_estimate, range) |>

  gt(rowname_col = "label") |>
  cols_label(label = "Parameter",
             description = "Description",
             point_estimate = "Point Estimate",
             range = "Range") |>

  tab_row_group(label = "Aedes",
                rows = 14:23) |>
  tab_row_group(label = "Livestock",
                rows = 1:13) |>

  tab_stubhead(label = "Parameter")

listed_tables <- list(gt_table1, gt_table2)

gt_two_column <- gt_two_column_layout(listed_tables, output = "save",
                                      filename = "gt_two_column.png",
                                      vwidth = 1600, vheight = 600)


