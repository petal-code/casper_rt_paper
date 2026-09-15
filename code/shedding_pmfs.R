# Derive shedding pmfs

source(here::here("code/00_config.R"))


# SARS-CoV-2 ----
## Wannigamma et al ----
# Medians extracted using plot digitizer
wannigamma <- read_csv(here("documents/Wannigama_2024_figA_digitized.csv"))

wannigamma$variant <- factor(wannigamma$variant, 
                                     levels = c(
                                       "XBB.1.5",
                                       "EG.5.1",
                                       "HV.1",
                                       "JD.1.1",
                                       "BA.2.86",
                                       "JN.1"
                                     ))

# Plot
ggplot(wannigamma,aes(day,cp_mg,fill=variant,col=variant)) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "loess") +
  scale_y_log10() +
  scale_fill_viridis_d(option = "C", begin = 0.1, end = 0.9) +
  scale_color_viridis_d(option = "C", begin = 0.1, end = 0.9) +
  labs(title = "Wannigamma et al. (2024), digitized medians from fig A",
       x = "Days after symptom onset",
       y = "Fecal concentration of SARS-CoV-2 viral RNA\n(cp/mg)")
# magnitudes differ across variants but overall shape is similar

# Compute pmf 
# use linear interpolation to fill missing days
wannigamma_pmf <- wannigamma %>% 
  group_by(day) %>% 
  summarise(mean_conc = mean(cp_mg)) %>% 
  complete(day = -3:21) %>%
  mutate(mean_conc = na.approx(mean_conc, x = day)) %>% 
  mutate(pmf = mean_conc/sum(mean_conc))
sum(wannigamma_pmf$pmf) # verify sums to 1

plot(NonParametric(wannigamma_pmf$pmf))

## Ando et al ----
# Parameters (posterior medians)
A     <- 1550
alpha <- 0.20
beta  <- 25.4

# Evaluate at days 0:25
days <- 0:25
C_t  <- exp(-alpha * days) * (1 - exp(-(beta - alpha) * days))

# Normalize to PMF
ando_pmf <- C_t / sum(C_t)

ando_pmf
sum(ando_pmf)
plot(NonParametric(ando_pmf))

# Zheng et al ----
gamma_func <- function(t, A, B, C) {
  A * (t^B) * exp(-C * t)
}

logistic_func <- function(t, A, B, C) {
  A / (1 + exp(B * (t - C)))
}

# log10(Concentration + 1) -> Concentration, clipped at 0
to_conc <- function(y) pmax(10^y - 1, 0)

## Noro ----
noro_y <- gamma_func(t = 0:30,
                        A=5.10127,
                        B=0.47654,
                        C=0.09760)
noro_pred <- to_conc(noro_y)
noro_pmf <- noro_pred/sum(noro_pred)
sum(noro_pmf)
plot(NonParametric(noro_pmf)) # days after infection

## Rota ----
rota_y <- gamma_func(t = 0:30,
                     A=5.06567, 
                     B=0.16051, 
                     C=0.05349)
rota_pred <- to_conc(rota_y)
rota_pmf <- rota_pred/sum(rota_pred)
sum(rota_pmf)
plot(NonParametric(rota_pmf)) # days after infection
