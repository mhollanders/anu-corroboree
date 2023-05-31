# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, janitor, lubridate, nimble, ggdist)

# plot theme
theme_set(theme_classic(base_size = 10))
theme_update(axis.ticks = element_line(color = "#333333", linewidth = 1/4),
             axis.line = element_line(color = "#333333", linewidth = 1/4),
             axis.title = element_text(color = "#333333"),
             axis.text = element_text(color = "#333333"),
             legend.title = element_text(color = "#333333"),
             legend.text = element_text(color = "#333333"),
             legend.position = "none")

# load Micalong and Warogong-Cotterill (2018--2020) data
dat_mwc <- read_csv("data/NCF_demographics-Oct-22.csv") |>
  select(-Site) |>
  clean_names() |>
  filter(site != "Museum",
         !is.na(bd_load)) |> 
  mutate(site = factor(site, levels = c("Micalong", "Warogong-Cotterill")),
         infected = ifelse(bd_load > 0, 1, 0),
         year = case_when(
           date == "2018" ~ 2018,
           date == "7/03/2018" ~ 2018,
           date == "3/03/2019" ~ 2019,
           date == "3/03/2020" ~ 2020
         )) |>
  select(site, infected, year) |> 
  glimpse()

# load Micalong 2012 data
dat_m12 <- readxl::read_xlsx("data/P_pengilleyi_Micalong_2012.xlsx") |>
  clean_names() |>
  mutate(site = factor(site),
         infected = ifelse(!is.na(zse), 1, 0),
         year = year(ymd(date))) |>
  select(site, infected, year) |>
  glimpse()

# join tibbles
dat <- bind_rows(dat_mwc, dat_m12) |>
  mutate(year = factor(year, levels = c(2012, 2018:2020))) |>
  glimpse()

# model code
code <- nimbleCode({
  
  # priors
  alpha[1] ~ dlogis(0, 1)
  alpha[2] ~ dlogis(0, 1)
  sigma ~ dexp(1)
  
  # random year effects (noncentered)
  for (j in 1:n_year) {
    z[j] ~ dnorm(0, 1)
    epsilon[j] <- z[j] * sigma
  } # j
  
  # observation process
  for (i in 1:n) {
    psi[i] <- expit(alpha[site[i]] + epsilon[year[i]])
    y[i] ~ dbern(psi[i])
  } # i
  
  # derived quantities
  
  # yearly prevalence
  for (j in 1:n_year) {
    psi_year[j] <- expit(alpha[1] + epsilon[j])  # Micalong
  } # j
  psi_year[5] <- expit(alpha[2] + epsilon[2])    # Warogong-Cotterill
  
  # marginal mean Micalong
  psi_mu[1] <- mean(psi_year[1:n_year])
  psi_mu[2] <- psi_year[5]
  
})

# data, constants, parameters to monitor
data <- list(y = dat$infected) |> glimpse()
constants <- list(n = nrow(dat),
                  year = as.numeric(dat$year),
                  site = as.numeric(dat$site),
                  n_year = nlevels(dat$year)) |> glimpse()

# model and MCMC
Cmodel <- nimbleModel(code, constants, data) |> 
  compileNimble()
Cmcmc <- configureMCMC(Cmodel, monitors = c("alpha", "sigma", "psi_year", "psi_mu")) |>
  buildMCMC() |>
  compileNimble(project = Cmodel, resetFunctions = T)
fit <- runMCMC(Cmcmc, niter = 6e4, nburnin = 1e4, nchains = 4)
draws <- do.call(rbind, fit) |> as_tibble()
write.csv(draws, "mcmc/prevalence.csv")

# summary
draws |>
  # odds ratio
  mutate(or = (`psi_mu[1]` / (1 - `psi_mu[1]`)) / (`psi_mu[2]` / (1 - `psi_mu[2]`))) |>
  select(contains("psi_mu"), or) |>
  median_hdci()

# plot
draws |>
  select(contains("psi_year")) |>
  pivot_longer(everything(), names_to = "psi", values_to = "iteration") |>
  mutate(site = c(rep("Micalong", 4), "Warogong-Cotterill") |>
           rep(n() / 5) |>
           factor(levels(dat$site)),
         year = c(levels(dat$year), 2018) |>
           rep(n() / 5) |>
           factor()) |>
  ggplot(aes(year, iteration, fill = site)) +
  geom_jitter(data = dat,
              aes(y = if_else(infected == 1, 0.9, 0.1),
                  color = site),
              shape = 16,
              height = 1/15,
              width = 1/3,
              alpha = 1/2,
              show.legend = F) +
  stat_halfeye(point_interval = "median_hdci",
               .width = 0.95,
               interval_size_range = c(1/3, 2/3),
               slab_alpha = 0.382,
               position = position_dodge(width = 1/3)) +
  scale_color_manual(values = c("#333333", "darkred")) +
  scale_fill_manual(values = c("#333333", "darkred")) +
  scale_y_continuous(breaks = seq(0.2, 1, 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Infection prevalence",
       fill = "Site") +
  theme(legend.position = c(0.95, 0.75),
        legend.justification = c("right", "top"))
ggsave("figs/fig-psi.png", width = 5, height = 4, dpi = 300)
