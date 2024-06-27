#############################################################################################
# "Utils" for visualizing results from wild boar simulation
#
# Christoph Imboden
#
# Created: 26.06.2024 CAI
#
#############################################################################################

# clear R's brain
rm(list=ls())

# Function to calculate the number of available for hunting individuals prior to hunting!

# Input:        1) sim.output = output from function "sim_population_hunted_JAGS" 
#                               MUST be the whole simulation function output from my own function!

# Output:       1) avail =      vector with total available individuals per time step
#                               prior to the hunting spike!

calc.avail <- function(sim.output) {
  
  # determine the random realization that was used in the output
  r <- sim.output$rand.realization
  
  # somewhat awkward calculations...
  # just tallies up the number of transitions, recruits, survivors...
  avail <- sim.output$JAGS_output$sims.list$new_recruits_f_s[r,] +
    sim.output$JAGS_output$sims.list$new_transition_f_m[r,] +
    sim.output$JAGS_output$sims.list$new_recruits_f_m[r,] +
    sim.output$JAGS_output$sims.list$new_transition_ML_f_l[r,] +
    sim.output$JAGS_output$sims.list$new_transition_SL_f_l[r,] +
    sim.output$JAGS_output$sims.list$new_recruits_m_s[r,] +
    sim.output$JAGS_output$sims.list$new_transition_m_m[r,] +
    sim.output$JAGS_output$sims.list$new_recruits_m_m[r,] +
    sim.output$JAGS_output$sims.list$new_transition_ML_m_l[r,] +
    sim.output$JAGS_output$sims.list$new_transition_SL_m_l[r,] +
    sim.output$JAGS_output$sims.list$surv_f_s[r,] +
    sim.output$JAGS_output$sims.list$surv_f_m[r,] +
    sim.output$JAGS_output$sims.list$surv_f_l[r,] +
    sim.output$JAGS_output$sims.list$surv_m_s[r,] +
    sim.output$JAGS_output$sims.list$surv_m_m[r,] +
    sim.output$JAGS_output$sims.list$surv_m_l[r,]
  
  return(avail)
}