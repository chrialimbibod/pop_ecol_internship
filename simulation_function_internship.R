#############################################################################################
# SCRIPT WITH FUNCTION FOR A SIMULATION OF A HYPOTHETICAL POPULATION OF WILD BOAR
#
# Idea: Take a population vector (6 states = 3 weight classes x 2 sexes) and simulate
#       population
#
# Christoph Imboden
#
# Created: 27.05.2024 CAI
#
#############################################################################################

# clear R's brain
rm(list=ls())

# Space to load libraries
library(jagsUI)
library(IPMbook)
library(tidyverse) # not strictly required

# Set working directory
setwd("~/.../...")

##### 1) HUNTED POPULATION #####

# INPUT:                      1)  N0 = population vector of length 6 -> number of individuals in the stages
#                                 ATTENTION: too low numbers will end in collapse of population, and difficulties for JAGS
#                                 causing it to return an error (most presumably because I tries to divide by 0 [invalid parent node])
#                             2)  t_span = number of timesteps to simulate
#                             3)  hunting_regime = 1, 2, 3, 4, 5, or 6 (must be one of those!)
#                             4)  params_to_monitor = used for JAGS, by default pretty much all of interest, but can be specified -> must be a vector with string/char content
#                                 Quite loaded by default - can be reduced when running simulations!

# OUTPUT:                     1)  JAGS_output = list output from JAGS (contains the results for all the "params_to_monitor")
#                             2)  mean_HR_df = mean hunting bag, number of hunted animals per sex- and weight-class
#                             3)  rand.realization = random realization that was drawn for first visualizations
#                             4)  N.tot.t.random.realization = total N in random realization (vector)
#                             5)  N.hunted.t.random.realization = total N of hunted animals (vector)
#                             6)  HR_df.random.realization = random realization hunting bag, number of hunted animals per sex- and weight-class


sim_population_hunted_JAGS <- function(N0, t_span, hunting_regime, params_to_monitor = c("N","Ntot", "extinct", "r", "lambda", "stage.distr", "stable.stage.distr", "hunted", "C.hunted", "C.Ntot", "new_recruits_f_s", "new_transition_f_m", "new_recruits_f_m",
                                                                                         "new_transition_ML_f_l", "new_transition_SL_f_l", "new_recruits_m_s", "new_transition_m_m", "new_recruits_m_m", "new_transition_ML_m_l", "new_transition_SL_m_l",
                                                                                         "surv_f_s", "surv_f_m", "surv_f_l", "surv_m_s", "surv_m_m", "surv_m_l", "hunt_f_s", "hunt_f_m", "hunt_f_l",
                                                                                         "hunt_m_s", "hunt_m_m", "hunt_m_l", "mean.BPs", "mean.BPm", "mean.BPl", "mean.Sns", "mean.Snm", "mean.Snl",
                                                                                         "mean.Ssm", "mean.Smm", "mean.Slm", "mean.LSs", "mean.LSm", "mean.LSl", "opfs", "opfm", "opfl", "OPFS", "OPFM", "OPFL", "OPF",
                                                                                         "nat_mort_f_s", "nat_mort_f_m", "nat_mort_f_l", "nat_mort_m_s", "nat_mort_m_m", "nat_mort_m_l", "p.hunted.recov", "p.Ntot.detection", "r.annual", "lambda.annual")) {
  
  # we need JAGS for simulation
  require(jagsUI)
  # IPMbood is not strictly needed in original output, but is used for getting alpha and beta parameters of a beta distribution
  require(IPMbook)
  
  # get the current working directory (path may be needed for the running of JAGS)
  curr_working_directory <- getwd()
  
  # quickly checks whether the initial population vecotr is of length 6 and otherwise spits out an error
  if (length(N0) != 6) {
    stop("Initial population size vector must be of length 6!")
  }
  
  # Split the initial population vector into components and let the numbers vary by using Poisson distribution
  N1_1 <- rpois(1, N0[1])
  N1_2 <- rpois(1, N0[2])
  N1_3 <- rpois(1, N0[3])
  N1_4 <- rpois(1, N0[4])
  N1_5 <- rpois(1, N0[5])
  N1_6 <- rpois(1, N0[6])
  
  # assign number of time steps to T (just because using T is nicer during development)
  T <- t_span
  
  # Hunting regime 1: more or less random against all size classes and both sexes
  # Hunting regime 2: preferred against small and medium males and females
  # Hunting regime 3: preferred against large females
  # Hunting regime 4: preferred and directed against small and medium females ("young" & most important reproductive units)
  # Hunting regime 5: preferred against small females
  # Hunting regime 6: preferred against medium females
  
  # We take hXX as proportion (= success in Binomial trials) of population in this stage and sex being killed by hunting
  # Furthermore, as we can assume the hunting pressure stays about the same, an increase in one hXX.mean does lead to a +/- decrease
  # in another hXX.mean
  # This is very much a simulation and based on own opinion rather than supported by data! Quite static and unrealistic!
  # For nomenclature and explanation see report.
  
  if (hunting_regime == 1) {
    # female hunting
    hsf.mean <- 0.5
    hsf.sd.t <- 0.35          # this is a guesstimate
    hmf.mean <- 0.4
    hmf.sd.t <- 0.29          # this is a guesstimate
    hlf.mean <- 0.5
    hlf.sd.t <- 0.45          # this is a guesstimate
    
    # male hunting
    hsm.mean <- 0.45
    hsm.sd.t <- 0.35          # this is a guesstimate
    hmm.mean <- 0.45
    hmm.sd.t <- 0.35          # this is a guesstimate
    hlm.mean <- 0.6
    hlm.sd.t <- 0.51          # this is a guesstimate
    
  } else if (hunting_regime == 2) { # following is own idea/guesses/opinions
    # female hunting
    hsf.mean <- 0.65
    hsf.sd.t <- 0.325
    hmf.mean <- 0.65
    hmf.sd.t <- 0.325
    hlf.mean <- 0.17
    hlf.sd.t <- 0.1
    
    # male hunting
    hsm.mean <- 0.65
    hsm.sd.t <- 0.325
    hmm.mean <- 0.65
    hmm.sd.t <- 0.325
    hlm.mean <- 0.17
    hlm.sd.t <- 0.1
    
  } else if (hunting_regime == 3) {
    # female hunting
    hsf.mean <- 0.2
    hsf.sd.t <- 0.15
    hmf.mean <- 0.2
    hmf.sd.t <- 0.15
    hlf.mean <- 0.8
    hlf.sd.t <- 0.4
    
    # male hunting
    hsm.mean <- 0.2
    hsm.sd.t <- 0.15
    hmm.mean <- 0.2
    hmm.sd.t <- 0.15
    hlm.mean <- 0.25
    hlm.sd.t <- 0.15
    
  } else if (hunting_regime == 4) {
    # female hunting
    hsf.mean <- 0.75
    hsf.sd.t <- 0.412
    hmf.mean <- 0.75
    hmf.sd.t <- 0.412
    hlf.mean <- 0.15
    hlf.sd.t <- 0.145
    
    # male hunting
    hsm.mean <- 0.4
    hsm.sd.t <- 0.35
    hmm.mean <- 0.4
    hmm.sd.t <- 0.35
    hlm.mean <- 0.3
    hlm.sd.t <- 0.145
    
  } else if (hunting_regime == 5) {
    # female hunting
    hsf.mean <- 0.868
    hsf.sd.t <- 0.5
    hmf.mean <- 0.36
    hmf.sd.t <- 0.32
    hlf.mean <- 0.3
    hlf.sd.t <- 0.3
    
    # male hunting
    hsm.mean <- 0.32
    hsm.sd.t <- 0.31
    hmm.mean <- 0.34
    hmm.sd.t <- 0.31
    hlm.mean <- 0.34
    hlm.sd.t <- 0.145
    
  } else if (hunting_regime == 6) {
    # female hunting
    hsf.mean <- 0.36
    hsf.sd.t <- 0.32
    hmf.mean <- 0.82
    hmf.sd.t <- 0.5
    hlf.mean <- 0.24
    hlf.sd.t <- 0.145
    
    # male hunting
    hsm.mean <- 0.32
    hsm.sd.t <- 0.31
    hmm.mean <- 0.34
    hmm.sd.t <- 0.31
    hlm.mean <- 0.34
    hlm.sd.t <- 0.145
    
  } else {
    stop("Please use 1, 2, 3, 4, or 5 as a hunting regime!")
  }
  
  # Definition of mean and temporal variability of the demographic 
  # parameters that are not assumed to be constant
  # For nomenclature and explanation see report.
  
  mean.BPs <- 0.139           # on probability scale (from Touzot et al 2020)
  sd.BPs.t <- 0.017
  
  mean.BPm <- 0.639           # on probability scale (from Touzot et al 2020)
  sd.BPm.t <- 0.022
  
  mean.BPl <- 0.801           # on probability scale (from Touzot et al 2020)
  sd.BPl.t <- 0.025
  
  mean.Sns <- 0.978           # on probability scale (from Touzot et al 2020)
  sd.Sns.t <- 0.002
  
  mean.Snm <- 0.920           # on probability scale (from Touzot et al 2020)
  sd.Snm.t <- 0.007
  
  mean.Snl <- 0.917           # on probability scale (from Touzot et al 2020)
  sd.Snl.t <- 0.006
  
  mean.Ssm <- 0.962           # on probability scale (from Gamelon et al 2012)
  sd.Ssm.t <- 0.002           # assumed to be the same as in females
  
  mean.Smm <- 0.777           # on probability scale (from Gamelon et al 2012)
  sd.Smm.t <- 0.007           # assumed to be the same as in females
  
  mean.Slm <- 0.904           # on probability scale (from Gamelon et al 2012)
  sd.Slm.t <- 0.006           # assumed to be the same as in females
  
  mean.LSs <- 3.918           # from Touzot et al 2020
  sd.LSs.t <- 0.250
  
  mean.LSm <- 4.707           # from Touzot et al 2020
  sd.LSm.t <- 0.104
  
  mean.LSl <- 6.096           # from Touzot et al 2020
  sd.LSl.t <- 0.121
  
  # Other data from the transition matrix that are deemed +/- constant over time (could be made time dependent)
  # For nomenclature and explanation see report.
  
  mpimm <- 0.322              # from Gamelon et al 2012
  mpSS <- 0.253               # from Gamelon et al 2012
  Spn <- 0.750                # from expert opinion Gamelon et al. 2012
  mpiOs <- 0.6                # from expert opinion Gamelon et al. 2012
  piOs <- 0.6                 # from expert opinion Gamelon et al. 2012
  pSS <- 0.137                # from Touzot et al 2020 (variation omitted due to missing data in males)
  pSM <- 0.234                # from Touzot et al 2020 (variation omitted due to missing data in males)
  pSL <- 0.556                # from Touzot et al 2020 (variation omitted due to missing data in males)
  pML <- 0.490                # from Touzot et al 2020 (variation omitted due to missing data in males)
  
  # Detection/Recovery probability means and sd
  # Not in report because it was not reviewed!
  mean.p.hunted.recov <- 0.98
  sd.p.hunted.recov <- 0.002
  mean.p.N.detection <- 0.20
  sd.p.N.detection <- 0.1
  
  # Bundle data
  jags.data <- list(mean.BPs = mean.BPs, mean.BPm = mean.BPm, mean.BPl = mean.BPl,
                    mean.Sns = mean.Sns, mean.Snm = mean.Snm, mean.Snl = mean.Snl,
                    mean.Ssm = mean.Ssm, mean.Smm = mean.Smm, mean.Slm = mean.Slm,
                    mean.LSs = mean.LSs, mean.LSm = mean.LSm, mean.LSl = mean.LSl,
                    sd.BPs.t = sd.BPs.t, sd.BPm.t = sd.BPm.t, sd.BPl.t = sd.BPl.t, 
                    sd.Sns.t = sd.Sns.t, sd.Snm.t = sd.Snm.t, sd.Snl.t = sd.Snl.t, 
                    sd.Ssm.t = sd.Ssm.t, sd.Smm.t = sd.Smm.t, sd.Slm.t = sd.Slm.t,
                    sd.LSs.t = sd.LSs.t, sd.LSm.t = sd.LSm.t, sd.LSl.t = sd.LSl.t,
                    mpimm = mpimm, mpSS = mpSS, Spn = Spn, mpiOs = mpiOs, piOs = piOs, 
                    pSS = pSS, pSM = pSM, pSL = pSL, pML = pML,
                    hsf.mean = hsf.mean, hsf.sd.t = hsf.sd.t, hmf.mean = hmf.mean, hmf.sd.t = hmf.sd.t,
                    hlf.mean = hlf.mean, hlf.sd.t = hlf.sd.t, hsm.mean = hsm.mean, hsm.sd.t = hsm.sd.t,
                    hmm.mean = hmm.mean, hmm.sd.t = hmm.sd.t, hlm.mean = hlm.mean, hlm.sd.t = hlm.sd.t,
                    T = T, N1_1 = N1_1, N1_2 = N1_2, N1_3 = N1_3, N1_4 = N1_4, N1_5 = N1_5, N1_6 = N1_6,
                    alpha.p.hunted = getBeta2Par(mean.p.hunted.recov, sd.p.hunted.recov)[1], beta.p.hunted = getBeta2Par(mean.p.hunted.recov, sd.p.hunted.recov)[2],
                    alpha.p.Ntot = getBeta2Par(mean.p.N.detection, sd.p.N.detection)[1], beta.p.Ntot = getBeta2Par(mean.p.N.detection, sd.p.N.detection)[2]) 
  
  # Write JAGS model file 
  cat(file = "sim_model_wildboar_hunted.txt", "
      model {
      # Calculate precision of the temporal variability in parameters/demographic rates
      tau.logit.BPs.t <- pow(sd.BPs.t, -2)
      tau.logit.BPm.t <- pow(sd.BPm.t, -2)
      tau.logit.BPl.t <- pow(sd.BPl.t, -2)
      tau.logit.Sns.t <- pow(sd.Sns.t, -2)
      tau.logit.Snm.t <- pow(sd.Snm.t, -2)
      tau.logit.Snl.t <- pow(sd.Snl.t, -2)
      tau.logit.Ssm.t <- pow(sd.Ssm.t, -2)
      tau.logit.Smm.t <- pow(sd.Smm.t, -2)
      tau.logit.Slm.t <- pow(sd.Slm.t, -2)
      tau.logit.hsf.t <- pow(hsf.sd.t, -2)
      tau.logit.hmf.t <- pow(hmf.sd.t, -2)
      tau.logit.hlf.t <- pow(hlf.sd.t, -2)
      tau.logit.hsm.t <- pow(hsm.sd.t, -2)
      tau.logit.hmm.t <- pow(hmm.sd.t, -2)
      tau.logit.hlm.t <- pow(hlm.sd.t, -2)
      tau.LSs.t <- pow(sd.LSs.t, -2)
      tau.LSm.t <- pow(sd.LSm.t, -2)
      tau.LSl.t <- pow(sd.LSl.t, -2)
      
      # Use of RNG to accomodate temporal variability of parameters/demographic rates +
      # Use of RNG to accomodate temporal variability of hunting
      for (t in 1:T) {
        BPs[t] <- ilogit(logit.BPs[t])
        logit.BPs[t] ~ dnorm(logit(mean.BPs), tau.logit.BPs.t)
        BPm[t] <- ilogit(logit.BPm[t])
        logit.BPm[t] ~ dnorm(logit(mean.BPm), tau.logit.BPm.t)
        BPl[t] <- ilogit(logit.BPl[t])
        logit.BPl[t] ~ dnorm(logit(mean.BPl), tau.logit.BPl.t)
        Sns[t] <- ilogit(logit.Sns[t])
        logit.Sns[t] ~ dnorm(logit(mean.Sns), tau.logit.Sns.t)
        Snm[t] <- ilogit(logit.Snm[t])
        logit.Snm[t] ~ dnorm(logit(mean.Snm), tau.logit.Snm.t)
        Snl[t] <- ilogit(logit.Snl[t])
        logit.Snl[t] ~ dnorm(logit(mean.Snl), tau.logit.Snl.t)
        Ssm[t] <- ilogit(logit.Ssm[t])
        logit.Ssm[t] ~ dnorm(logit(mean.Ssm), tau.logit.Ssm.t)
        Smm[t] <- ilogit(logit.Smm[t])
        logit.Smm[t] ~ dnorm(logit(mean.Smm), tau.logit.Smm.t)
        Slm[t] <- ilogit(logit.Slm[t])
        logit.Slm[t] ~ dnorm(logit(mean.Slm), tau.logit.Slm.t)
        hsf[t] <- ilogit(logit.hsf[t])
        logit.hsf[t] ~ dnorm(logit(hsf.mean), tau.logit.hsf.t)
        hmf[t] <- ilogit(logit.hmf[t])
        logit.hmf[t] ~ dnorm(logit(hmf.mean), tau.logit.hmf.t)
        hlf[t] <- ilogit(logit.hlf[t])
        logit.hlf[t] ~ dnorm(logit(hlf.mean), tau.logit.hlf.t)
        hsm[t] <- ilogit(logit.hsm[t])
        logit.hsm[t] ~ dnorm(logit(hsm.mean), tau.logit.hsm.t)
        hmm[t] <- ilogit(logit.hmm[t])
        logit.hmm[t] ~ dnorm(logit(hmm.mean), tau.logit.hmm.t)
        hlm[t] <- ilogit(logit.hlm[t])
        logit.hlm[t] ~ dnorm(logit(hlm.mean), tau.logit.hlm.t)
        LSs[t] ~ dnorm(mean.LSs, tau.LSs.t)T(0,)                  # Truncation to positive values
        LSm[t] ~ dnorm(mean.LSm, tau.LSm.t)T(0,)                  # Truncation to positive values
        LSl[t] ~ dnorm(mean.LSl, tau.LSl.t)T(0,)                  # Truncation to positive values
      }
      
      # Not in report because it was not reviewed!
      # Also use RNG to accomodate tempral variability in detection/recovery probabilities
      for (t in 1:T) {
        p.hunted.recov[t] ~ dbeta(alpha.p.hunted, beta.p.hunted)
        p.Ntot.detection[t] ~ dbeta(alpha.p.Ntot, beta.p.Ntot)
      }
      
      # Model (not really model...) for initial state
      N[1,1] <- N1_1
      N[2,1] <- N1_2
      N[3,1] <- N1_3
      N[4,1] <- N1_4
      N[5,1] <- N1_5
      N[6,1] <- N1_6
      
      Ntot[1] <- N1_1 + N1_2 + N1_3 + N1_4 + N1_5 + N1_6
      
      # Loop over time
      for (t in 1:T) {
        
        # Population model using state specific population equations and distributions to account for demographic stochasticity
        # Cannot use matrix multiplication in JAGS (Nimble would support this)
        
        # Changes in number of small females:
        new_rec_contribution_f_s_f_s[t] ~ dpois((N[1,t]*BPs[t]*LSs[t]*0.5*Spn*piOs))
        new_rec_contribution_f_s_f_m[t] ~ dpois((N[2,t]*BPm[t]*LSm[t]*0.5*Spn*piOs))
        new_rec_contribution_f_s_f_l[t] ~ dpois((N[3,t]*BPl[t]*LSl[t]*0.5*Spn*piOs))
        new_recruits_f_s[t] ~ sum(new_rec_contribution_f_s_f_s[t] + new_rec_contribution_f_s_f_m[t] + new_rec_contribution_f_s_f_l[t])
        surv_f_s[t] ~ dbin((pSS*Sns[t]), N[1,t])
        hunt_f_s[t] ~ dbin(hsf[t], (surv_f_s[t] + new_recruits_f_s[t]))
        nat_mort_f_s[t] <- N[1,t] - surv_f_s[t]                                 # only looks at the one that have already been in this stage
        N[1,t+1] <- (surv_f_s[t] + new_recruits_f_s[t]) - hunt_f_s[t]
        
        # Changes in number of medium females:
        new_transition_f_m[t] ~ dbin(pSM*Snm[t], N[1,t])
        new_rec_contribution_f_m_f_s[t] ~ dpois(N[1,t]*BPs[t]*LSs[t]*0.5*Spn*(1-piOs))
        new_rec_contribution_f_m_f_m[t] ~ dpois(N[2,t]*BPm[t]*LSm[t]*0.5*Spn*(1-piOs))
        new_rec_contribution_f_m_f_l[t] ~ dpois(N[3,t]*BPl[t]*LSl[t]*0.5*Spn*(1-piOs))
        new_recruits_f_m[t] ~ sum(new_rec_contribution_f_m_f_s[t] + new_rec_contribution_f_m_f_m[t] + new_rec_contribution_f_m_f_l[t])
        surv_f_m[t] ~ dbin((1-pML)*Snm[t], N[2,t])
        hunt_f_m[t] ~ dbin(hmf[t], (new_transition_f_m[t] + new_recruits_f_m[t] + surv_f_m[t]))
        nat_mort_f_m[t] <- N[2,t] - surv_f_m[t]                                 # only looks at the one that have already been in this stage
        N[2,t+1] <- (new_transition_f_m[t] + new_recruits_f_m[t] + surv_f_m[t]) - hunt_f_m[t]
        
        # Changes in number of large females:
        new_transition_ML_f_l[t] ~ dbin(pML*Snl[t], N[2,t])
        new_transition_SL_f_l[t] ~ dbin(pSL*Snl[t], N[1,t])
        surv_f_l[t] ~ dbin(Snl[t], N[3,t])
        hunt_f_l[t] ~ dbin(hlf[t], (new_transition_ML_f_l[t] + new_transition_SL_f_l[t] + surv_f_l[t]))
        nat_mort_f_l[t] <- N[3,t] - surv_f_l[t]                                 # only looks at the one that have already been in this stage
        N[3,t+1] <- (new_transition_ML_f_l[t] + new_transition_SL_f_l[t] + surv_f_l[t]) - hunt_f_l[t]
        
        # Changes in number of small males:
        new_rec_contribution_m_s_f_s[t] ~ dpois((N[1,t]*BPs[t]*LSs[t]*0.5*Spn*mpiOs))
        new_rec_contribution_m_s_f_m[t] ~ dpois((N[2,t]*BPm[t]*LSm[t]*0.5*Spn*mpiOs))
        new_rec_contribution_m_s_f_l[t] ~ dpois((N[3,t]*BPl[t]*LSl[t]*0.5*Spn*mpiOs))
        new_recruits_m_s[t] ~ sum(new_rec_contribution_m_s_f_s[t] + new_rec_contribution_m_s_f_m[t] + new_rec_contribution_m_s_f_l[t])
        surv_m_s[t] ~ dbin(mpSS*Ssm[t], N[4,t])
        hunt_m_s[t] ~ dbin(hsm[t], (new_recruits_m_s[t] + surv_m_s[t]))
        nat_mort_m_s[t] <- N[4,t] - surv_m_s[t]                                 # only looks at the one that have already been in this stage
        N[4,t+1] <- (new_recruits_m_s[t] + surv_m_s[t]) - hunt_m_s[t]
        
        # Changes in number of medium males:
        new_transition_m_m[t] ~ dbin((1-mpSS)*Smm[t], N[4,t])
        new_rec_contribution_m_m_f_s[t] ~ dpois((N[1,t]*BPs[t]*LSs[t]*0.5*Spn*(1-mpiOs)))
        new_rec_contribution_m_m_f_m[t] ~ dpois((N[2,t]*BPm[t]*LSm[t]*0.5*Spn*(1-mpiOs)))
        new_rec_contribution_m_m_f_l[t] ~ dpois((N[3,t]*BPl[t]*LSl[t]*0.5*Spn*(1-mpiOs)))
        new_recruits_m_m[t] ~ sum(new_rec_contribution_m_m_f_s[t] + new_rec_contribution_m_m_f_m[t] + new_rec_contribution_m_m_f_l[t])
        surv_m_m[t] ~ dbin(mpimm*Smm[t], N[5,t])
        hunt_m_m[t] ~ dbin(hmm[t], (new_transition_m_m[t] + new_recruits_m_m[t] + surv_m_m[t]))
        nat_mort_m_m[t] <- N[5,t] - surv_m_m[t]                                 # only looks at the one that have already been in this stage
        N[5,t+1] <- (new_transition_m_m[t] + new_recruits_m_m[t] + surv_m_m[t]) - hunt_m_m[t]
        
        # Changes to number of large males:
        new_transition_ML_m_l[t] ~ dbin((1-mpimm)*Slm[t], N[5,t])
        new_transition_SL_m_l[t] ~ dbin(pSL*Slm[t], N[4,t])
        surv_m_l[t] ~ dbin(Slm[t], N[6,t])
        hunt_m_l[t] ~ dbin(hlm[t], (new_transition_ML_m_l[t] + new_transition_SL_m_l[t] + surv_m_l[t]))
        nat_mort_m_l[t] <- N[6,t] - surv_m_l[t]                                 # only looks at the one that have already been in this stage
        N[6,t+1] <- (new_transition_ML_m_l[t] + new_transition_SL_m_l[t] + surv_m_l[t]) - hunt_m_l[t]
        
        # Extinction parameter
        extinct[t] <- equals(N[1,t+1] + N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1], 0) # to determine whether population went extinct or not!
        # 0 = NOT extinct and 1 = went extinct
        
        # Total population per time step
        Ntot[t+1] <- N[1,t+1] + N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1]
        
        # Total number of hunted individuals per time step
        hunted[t] <- hunt_f_s[t] + hunt_f_m[t] + hunt_f_l[t] + hunt_m_s[t] + hunt_m_m[t] + hunt_m_l[t]
        # hunted[t] <- sum(hunt_f_s[t], hunt_f_m[t], hunt_f_l[t], hunt_m_s[t], hunt_m_m[t], hunt_m_l[t])
        
        # Not in report because it was not reviewed!
        # Offspring per female
        opfs[t] <- ifelse(N[1,t] > 0, (new_rec_contribution_f_s_f_s[t] + new_rec_contribution_f_m_f_s[t] + new_rec_contribution_m_s_f_s[t] + new_rec_contribution_m_m_f_s[t])/N[1,t], 0)
        # opfs[t] <- (new_rec_contribution_f_s_f_s[t] + new_rec_contribution_f_m_f_s[t] + new_rec_contribution_m_s_f_s[t] + new_rec_contribution_m_m_f_s[t])/N[1,t]
        opfm[t] <- ifelse(N[2,t] > 0, (new_rec_contribution_f_s_f_m[t] + new_rec_contribution_f_m_f_m[t] + new_rec_contribution_m_s_f_m[t] + new_rec_contribution_m_m_f_m[t])/N[2,t], 0)
        # opfm[t] <- (new_rec_contribution_f_s_f_m[t] + new_rec_contribution_f_m_f_m[t] + new_rec_contribution_m_s_f_m[t] + new_rec_contribution_m_m_f_m[t])/N[2,t]
        opfl[t] <- ifelse(N[3,t] > 0, (new_rec_contribution_f_s_f_l[t] + new_rec_contribution_f_m_f_l[t] + new_rec_contribution_m_s_f_l[t] + new_rec_contribution_m_m_f_l[t])/N[3,t], 0)
        # opfl[t] <- (new_rec_contribution_f_s_f_l[t] + new_rec_contribution_f_m_f_l[t] + new_rec_contribution_m_s_f_l[t] + new_rec_contribution_m_m_f_l[t])/N[3,t]
        
        # Annual growth rate on log scale
        r.annual[t] <- log(N[1,t+1] + N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1]) - log(N[1,t] + N[2,t] + N[3,t] + N[4,t] + N[5,t] + N[6,t])
        lambda.annual[t] <- (N[1,t+1] + N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])/(N[1,t] + N[2,t] + N[3,t] + N[4,t] + N[5,t] + N[6,t])
        
        # Not in report because it was not reviewed!
        # Scaled stage distributions
        stage.distr[1,t] <- N[1,t+1] / (N[1,t+1] +N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])
        stage.distr[2,t] <- N[2,t+1] / (N[1,t+1] +N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])
        stage.distr[3,t] <- N[3,t+1] / (N[1,t+1] +N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])
        stage.distr[4,t] <- N[4,t+1] / (N[1,t+1] +N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])
        stage.distr[5,t] <- N[5,t+1] / (N[1,t+1] +N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])
        stage.distr[6,t] <- N[6,t+1] / (N[1,t+1] +N[2,t+1] + N[3,t+1] + N[4,t+1] + N[5,t+1] + N[6,t+1])
        
      }
      
      # Derived quantities 
      
      # Not in report because it was not reviewed!
      # Mean offspring per female (split per size class)
      OPFS <- mean(opfs[1:T])
      OPFM <- mean(opfm[1:T])
      OPFL <- mean(opfl[1:T])
      # Overall mean offspring per female (not split per size class)
      OPF <- (OPFS + OPFM + OPFL)/3
      
      # r and lambda
      r <- mean(r.annual[5:T]) # assumes stable stage/age distribution = equilibrium is reached after 5 timesteps...
      lambda <- exp(r)
      
      # Not in report because it was not reviewed!
      # Stable stage distribution
      stable.stage.distr <- stage.distr[,T]
      
      # Not in report because it was not reviewed!
      # Observation process
      for (t in 1:T) {
        C.hunted[t] ~ dbin(p.hunted.recov[t], hunted[t])
        C.Ntot[t] ~ dbin(p.Ntot.detection[t], Ntot[t])
      }
      
      }") # end of JAGS model file
  
  # Parameters monitored (for sure "N" needed!)
  parameters <- params_to_monitor
  # MCMC settings (we are not estimating anything and therefor do not need multiple chains, and also no burn in period)
  ni <- 50000; nt <- 1; nb <- 0; nc <- 1; na <- 0
  
  # Call JAGS from R (ART < 1 min)
  output_sim <- jags(data = jags.data, inits = NULL, parameters.to.save = parameters, model.file = paste(curr_working_directory,"sim_model_wildboar_hunted.txt", sep = "/"),
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC = FALSE)
  
  # Get average/mean of hunted individuals per time step
  mean.hunted.small.female <- output_sim$mean$hunt_f_s
  mean.hunted.medium.female <- output_sim$mean$hunt_f_m
  mean.hunted.large.female <- output_sim$mean$hunt_f_l
  mean.hunted.small.male <- output_sim$mean$hunt_m_s
  mean.hunted.medium.male <- output_sim$mean$hunt_m_m
  mean.hunted.large.male <- output_sim$mean$hunt_m_l
  
  mean_HR_df <- rbind(mean.hunted.small.female, mean.hunted.medium.female, mean.hunted.large.female,
                      mean.hunted.small.male, mean.hunted.medium.male, mean.hunted.large.male)
  
  # PICK RANDOM REALIZATION
  # set a +/- random seed -> that is most of the time different for people
  # used to pick random number/realization
  set.seed(Sys.time())
  rand.int <- sample(1:ni, 1)
  
  # extract information on Ntot of random realization/process
  N_total <- output_sim$sims.list$Ntot[rand.int[1], ]
  
  # extract information on number of hunted of random realization/process
  N_hunted <- output_sim$sims.list$hunted[rand.int[1], ]
  
  # extract the HR dataframe for random realization
  rand.hunted.f.s <- output_sim$sims.list$hunt_f_s[rand.int,]
  rand.hunted.f.m <- output_sim$sims.list$hunt_f_m[rand.int,]
  rand.hunted.f.l <- output_sim$sims.list$hunt_f_l[rand.int,]
  rand.hunted.m.s <- output_sim$sims.list$hunt_m_s[rand.int,]
  rand.hunted.m.m <- output_sim$sims.list$hunt_m_m[rand.int,]
  rand.hunted.m.l <- output_sim$sims.list$hunt_m_l[rand.int,]
  
  rand_HR_df <- rbind(rand.hunted.f.s, rand.hunted.f.m, rand.hunted.f.l, rand.hunted.m.s, rand.hunted.m.m, rand.hunted.m.l)
  
  # Prepare list for return
  return(list(JAGS_output = output_sim, mean_HR_df = mean_HR_df, rand.realization = rand.int, N.tot.t.random.realization = N_total,
              N.hunted.t.random.realization = N_hunted, HR_df.random.realization = rand_HR_df))
  
}
