# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, nimble, MCMCvis, ggdist)

# plot theme
theme_set(theme_classic(base_size = 10))
theme_update(axis.ticks = element_line(color = "#333333", linewidth = 1/4),
             axis.line = element_line(color = "#333333", linewidth = 1/4),
             axis.title = element_text(color = "#333333"),
             axis.text = element_text(color = "#333333"),
             legend.title = element_text(color = "#333333"),
             legend.text = element_text(color = "#333333"),
             legend.position = "none")

# load data
dat <- read_csv("data/NCF_demographics-Oct-22.csv") |>
  mutate(site = factor(site, levels = c("Micalong", "Warogong-Cotterill", "Museum")))

# model code
code <- nimbleCode({
  
  # priors
  for (j in 1:n_site) {
    p[j] ~ dbeta(1, 1)
    r[j] ~ dexp(0.5)
  }
  
  # likelihood
  for (i in 1:n_ind) {
    age[i] ~ dnegbin(prob = p[site[i]], size = r[site[i]])
  }
  
  # derived
  mean[1:n_site] <- r[1:n_site] * (1 - p[1:n_site]) / p[1:n_site]
  
})

# data, constants, and parameters to monitor
data <- list(age = dat$age) |> glimpse()
consts <- list(n_ind = nrow(dat),
               n_site = length(levels(dat$site)),
               site = as.numeric(dat$site)) |> glimpse()
mons <- c("p", "r", "mean")

# model and mcmc
Cmodel <- nimbleModel(code, consts, data) |> compileNimble()
conf <- configureMCMC(Cmodel, monitors = mons)
Cmcmc <- buildMCMC(conf) |> compileNimble(project = Cmodel, resetFunctions = T)
n_chains <- 4 ; n_iter <- 6e4 ; n_burnin <- 1e4 ; n_thin <- 10
samples <- runMCMC(Cmcmc, nchains = n_chains, niter = n_iter, nburnin = n_burnin, thin = n_thin)

# inspect chain histories and calculate summaries
MCMCtrace(samples, Rhat = T, n.eff = T, ind = T, pdf = F)
MCMCsummary(samples, round = 3)

# save samples
draws <- do.call(rbind, samples)  |> as_tibble()
write.csv(draws, "mcmc/age.csv")

# extract mean age samples and prepare for plotting
age_mean <- draws |>
  select(starts_with("mean")) |>
  set_names(levels(dat$site)) |>
  pivot_longer(everything(), names_to = "site", values_to = "age") |>
  mutate(site = factor(site, levels = levels(dat$site)))

# plot
colors <- c("#559e83", "#ae5a41")
dat |>
  ggplot(aes(site)) +
  # observed ages (slightly spread for aesthetics)
  geom_dots(aes(y = age + rnorm(nrow(dat), 0, 0.1),
                color = site == "Museum",
                fill = site == "Museum"),
            normalize = "groups") +
  stat_halfeye(data = age_mean,
               aes(y = age),
               point_interval = "median_hdci",
               .width = 0.95, 
               interval_size_range = c(1/3, 2/3),
               slab_alpha = 0.382,
               color = "#333333",
               fill = "#333333") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = seq(1, 8, 1),
                     limits = c(0, 8.5),
                     expand = c(0, 0)) +
  labs(x = "Site",
       y = "Age")
ggsave("figs/fig-age.png", width = 5, height = 4, units = "in", dpi = 600)
