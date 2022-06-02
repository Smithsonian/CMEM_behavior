# Morris parabola bayesian

library(rjags)
library(ggmcmc)
library(tidyverse)
load.module("glm")
library(rCMEM)

# prep data space
dataFile <- read_csv("data/Morris_2013/derivative/Morris_et_al_2013_plant_plot.csv")

dataFileJags <- dataFile %>% 
  mutate(index = 1:nrow(dataFile),
         X = (elevation - MSL) / (MLHW-MSL),
         Y = AGB_total,
         obs.prec = 1/(AGB_total_se^2)) %>%
  rename(n = AGB_total_n) %>% 
  select(index, X, Y, obs.prec, n) 

dataFileJags[dataFileJags == Inf] <- NA

parabolic_bmass <- "model{
  
  emin ~ dnorm(0, 0.001)
  emax ~ dnorm(0, 0.001) T(emin, )
  epeak ~ dnorm(0, 0.001) T(emin, emax)
  bmax ~ dnorm(0, 0.001) T(0, )

  prec.y ~ dgamma(0.001, 0.001)

  # parabola function
  dummy_emin_high <- epeak-(emax-epeak)
  dummy_emax_low <- epeak+(epeak-emin)

  a_up <- -((-dummy_emin_high * bmax - emax * bmax) / 
    ((dummy_emin_high - epeak) * (-emax + epeak)))
  b_up <- -(bmax / 
    ((dummy_emin_high - epeak) * (-emax + epeak)))
  c_up <- (dummy_emin_high * emax * bmax) / 
    ((dummy_emin_high - epeak) * (emax - epeak))

  # Solve for the parametrs of the lower curve.
  a_low <- -((-emin * bmax - dummy_emax_low * bmax) / 
    ((emin - epeak) * (-dummy_emax_low + epeak)))
  b_low <- -(bmax / 
    ((emin - epeak) * (-dummy_emax_low + epeak)))
  c_low <- (emin * dummy_emax_low * bmax) / 
    ((emin - epeak) * (dummy_emax_low - epeak))
  
  # Iterate through points
  for (i in 1:14) {

    Z[i] <- ifelse(X[i] <= emin, 0,
        ifelse(X[i]>=emax, 0,
          ifelse(X[i]>epeak, 
            a_up*X[i] + b_up*X[i]^2 + c_up, 
            a_low*X[i] + b_low*X[i]^2 + c_low)))
    
    tau.y[i] <- prec.y*n[i]
    u1[i] <- n[i]/2                             
    u2[i] <- n[i]/(2*prec.y)
    obs.prec[i] ~ dgamma(u1[i], u2[i])
    
    Y[i] ~ dnorm(Z[i], tau.y[i]) T(0,)
    
  }
}"

data <- list(Y = dataFileJags$Y,
             X = dataFileJags$X,
             obs.prec = dataFileJags$obs.prec,
             n = dataFileJags$n)

j.model <- jags.model(file = textConnection(parabolic_bmass),
                      data=data,
                      n.chains = 3)

jags.out <- coda.samples(model=j.model, variable.names=c("emin",
                                                         "emax",
                                                         "epeak",
                                                         "bmax",
                                                         "prec.y"),
                         n.iter=10000)
tidyJags<- ggs(jags.out)

ggs_traceplot(tidyJags)

tidyJagsDQs <- tidyJags %>%
  group_by(Parameter) %>%
  summarise(median = median(value),
            UpperCI=quantile(value, 0.975),
            LowerCI=quantile(value, 0.025),
            sd = sd(value))

write_csv(tidyJagsDQs, "data/Morris_2013/tidyJagsSAparabloa.csv")

morris_vis <- dataFile %>% 
  mutate(zStar = (elevation-MSL)/(MLHW-MSL))

target_x <- seq(-0.5, 2.15, by = 0.05)
target_y <- predictBiomass(z = target_x,
                           bMax = 866,
                           zVegMax = 2.08, 
                           zVegMin = -0.470 , 
                           zVegPeak = 0.831)

ggplot(data = morris_vis, aes(x=zStar, y=AGB_total)) +
  geom_point() +
  geom_segment(aes(xend = zStar, y = AGB_total+AGB_total_se, yend = AGB_total-AGB_total_se)) +
  geom_line(data = data.frame(target_x, target_y), aes(x = target_x, y = target_y),
             lty = 2, size = 2, color = "grey") +
  theme_minimal()

