# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, nimble, nimbleNoBounds, MCMCvis, ggdist)

# load museum data
dat1 <- read_csv("data/NCF_demographics-Oct-22.csv") |>
  filter(site == "Museum",  # museum data only
         age > 3            # remove low detection ages
  )

# load Kosch et al. (2019) data and add censoring info
n_days <- 103    # experiment duration
n_control <- 76  # number of control individuals (not included in data)
dat2 <- read_csv("data/PhenotypeData.csv") |>
  mutate(surv = if_else(DaysSurvived == n_days, NA, DaysSurvived),
         censored = if_else(is.na(surv), 1, 0))

# model code
code <- nimbleCode({
  
  # priors
  alpha ~ dgamma(2, 1)       # annual mortality rate
  beta ~ dnorm(0, sd = 1.5)  # effect of Bd on mortality
  
  # likelihood 1: age distribution of museum specimens
  phi_y <- 1 - exp(-alpha)
  for (i in 1:n_ind[1]) {
    y[i] ~ dnegbin(prob = phi_y, size = 1)
  } # i
  
  # likelihood 2: Bd effect
  for (j in 1:n_ind[2]) {
    phi_x[j] <- 1 - exp(-exp(log(alpha / 365) 
                             + beta * inf[j]
                             )
                        )
    x[j] ~ dnegbin(prob = phi_x[j], size = 1)
    
    # right-censored survivors
    censored[j] ~ dinterval(x[j], c[j])
  } # j
  
})

# data, constants, initial values, and parameters to monitor
data <- list(y = dat1$age - 4,
             x = c(dat2$surv - 1, rep(NA, n_control)),
             censored = c(dat2$censored, rep(1, n_control))) |> glimpse()
consts <- list(n_ind = c(nrow(dat1), nrow(dat2) + n_control),
               c = if_else(data$censored == 1, n_days - 1, Inf),
               inf = c(rep(1, nrow(dat2)), rep(0, n_control))) |> glimpse()
inits <- list(x = if_else(is.na(data$x), 1000, NA),
              alpha = 0.5,
              beta = 1) |> glimpse()
mons <- c("alpha", "beta")

# model and mcmc
Cmodel <- nimbleModel(code, consts, data, inits) |> compileNimble()
conf <- configureMCMC(Cmodel, monitors = mons)
Cmcmc <- buildMCMC(conf) |> compileNimble(project = Cmodel, resetFunctions = T)
n_chains <- 4 ; n_iter <- 6e4 ; n_burnin <- 1e4 ; n_thin <- 10
samples <- runMCMC(Cmcmc, nchains = n_chains, niter = n_iter, nburnin = n_burnin, thin = n_thin)

# inspect chain histories and calculate summaries
MCMCtrace(samples, Rhat = T, n.eff = T, ind = T, pdf = F)
MCMCsummary(samples, round = 3)

# save samples
draws <- do.call(rbind, samples) |> as_tibble()
write.csv(draws, "mcmc/survival.csv")

# hazard ratio of infected:uninfected
hist(exp(draws$beta),
     yaxt = "n",
     main = NULL,
     xlab = "Hazard ratio",
     ylab = NULL,
     breaks = "fd")