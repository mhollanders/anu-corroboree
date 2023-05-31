# load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(popbio, tidyverse, reshape2, ggdist, patchwork)

# plot theme
theme_set(theme_classic(base_size = 10))
theme_update(axis.ticks = element_line(color = "#333333", linewidth = 1/4),
             axis.line = element_line(color = "#333333", linewidth = 1/4),
             axis.title = element_text(color = "#333333"),
             axis.text = element_text(color = "#333333"),
             legend.title = element_text(color = "#333333"),
             legend.text = element_text(color = "#333333"),
             legend.position = "none",
             strip.background = element_blank())

# load mcmc samples from previous analysis
draws <- read_csv("mcmc/survival.csv")

# number of age classes
n_age <- 7

# survival probabilities
n_psi <- 4  # prevalence classes
n_F <- 4    # fecundity classes (2, 3, 4, and 5 year sexual maturity)

# survival of uninfected and infected
psi <- c(0, 0.2, 0.4, 0.6)
P <- rep(NA, n_psi)
for (j in 1:n_psi) {
  P[j] <- 
    # uninfected individuals
    median(exp(-draws$alpha)) * (1 - psi[j]) +
    # infected individuals
    median(exp(-exp(log(draws$alpha / 365) + draws$beta) * 365)) * psi[j]
}

# fecundity that yields lambda = 1
f <- matrix(NA, n_psi, n_F)

# prevalence 0 (uninfected), age at maturation 2
f[1, 1] <- 0.83
# uninfected, age at maturation 3
f[1, 2] <- 1.55
# uninfected, age at maturation 4
f[1, 3] <- 2.92
# uninfected, age at maturation 5
f[1, 4] <- 5.75

# infected, first prevalence, age at maturation 2
f[2, 1] <- 1.27
# age at maturation 3
f[2, 2] <- 2.91
# age at maturation 4
f[2, 3] <- 6.71
# age at maturation 5
f[2, 4] <- 16.0

# infected, second prevalence, age at maturation 2
f[3, 1] <- 2.01
# age at maturation 3
f[3, 2] <- 6.08
# age at maturation 4
f[3, 3] <- 18.5
# age at maturation 5
f[3, 4] <- 57.1

# infected, third prevalence, age at maturation 2
f[4, 1] <- 3.51
# age at maturation 3
f[4, 2] <- 15.8
# age at maturation 4
f[4, 3] <- 71.4
# age at maturation 5
f[4, 4] <- 325

# elasticity analysis
m <- e <- array(NA, c(n_psi, n_F, n_age, n_age))
e_f <- e_P <- lambda <- matrix(NA, n_psi, n_F)
for (i in 1:n_psi) {
  for (j in 1:n_F) {
    
    # Leslie matrices
    m[i, j, , ] <- matrix(c(rep(0, j), rep(f[i, j], n_age - j),
                            P[i], rep(0, n_age - 1),
                            rep(0, 1), P[i], rep(0, n_age - 2),
                            rep(0, 2), P[i], rep(0, n_age - 3),
                            rep(0, 3), P[i], rep(0, n_age - 4),
                            rep(0, 4), P[i], rep(0, n_age - 5),
                            rep(0, 5), P[i], 0),
                          n_age, byrow = T)
    
    # population growth rates (to confirm equilibrium)
    lambda[i, j] <- lambda(m[i, j, , ])
    
    # elasticities
    e[i, j, , ] <- elasticity(m[i, j, , ])
    
    # averaged for fecundity and survival for each "treatment"
    e_f[i, j] <- mean(e[i, j, 1, j:n_age])
    e_P[i, j] <- mean(e[i, j, 2:n_age, ][e[i, j, 2:n_age, ] > 0])
    
    
  }
}

# convert f and P to long format
e_f_long <- melt(e_f, value.name = "e_f", varnames = c("prevalence", "age to maturation"))
e_P_long <- melt(e_P, value.name = "e_P", varnames = c("prevalence", "age to maturation"))

# plot
dat <- left_join(e_f_long, e_P_long) |>
  mutate(f = c(f)) |>
  pivot_longer(c(e_f, e_P), values_to = "e", names_to = "rate") |>
  mutate(prevalence = factor(prevalence, labels = psi),
         age = factor(`age to maturation`, labels = c(1:n_F) + 1),
         rate = factor(rate, levels = c("e_P", "e_f"), labels = c("Survival", "Fecundity")))

# elasticities
fig_e <- dat |>
  ggplot(aes(age, e)) +
  geom_point(aes(shape = rate, alpha = prevalence),
             size = 3,
             color = "#6b200c") +
  scale_shape_manual(values = c(16, 17)) +
  scale_alpha_manual(values = seq(0.4, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0.03, 0.12, 0.03),
                     limits = c(0.01, 0.14)) +
  theme(legend.position = "right") +
  labs(x = "Age to maturation",
       y = "Elasticity",
       shape = "Vital rate",
       alpha = "Prevalence")

# fecundity
fig_f <- dat |>
  distinct(prevalence, age, f) |>
  ggplot(aes(age, f)) +
  geom_hline(yintercept = 15,
             linetype = "dashed",
             color = "#333333") +
  geom_point(aes(alpha = prevalence),
             color = "#133e7e",
             shape = 16,
             size = 3) +
  scale_alpha_manual(values = seq(0.4, 1, 0.2)) +
  theme(legend.position = "right") +
  scale_y_log10() +
  labs(x = "Age to maturation",
       y = expression(paste("Fecundity at ", lambda, " = 1")),
       alpha = "Prevalence") +
  annotate("text", x = 1, y = 20, label = expression(paste(italic(f)[max], " = 15")))

# save
(fig_e | fig_f) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(face = "bold"))
ggsave("figs/fig-e-f.png", width = 9, height = 4)
