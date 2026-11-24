## Trial_function
library(TrialSimulator)
library(MyPackageV2)

Threshold_prob = function(lamda, gamma, n_total, n_current, symmetric = FALSE){
  ## Function responsible for the threshold probabilities for early stopping
  # gamma and lamda are the same variables in the BOP2 paper
  # n_total is the total sample size (of the arm)
  # n_current is the current sample size (of the arm)
  
  if (lamda > 1){
    if (n_total == n_current){
      return(c(1,1))
    }
    return(c(0,1))
  }
  
  if (symmetric){
    
    if (n_total == n_current){
      return(c(1,1))
    }
    return(c(0.1,0.9))
  }
  
  Futility = lamda * (n_current/n_total)^gamma
  Super = 2 * pnorm ( (qnorm ( (1 + lamda)/2))/sqrt(n_current/n_total)) - 1
  return(c(Futility, Super))
}

Trial_function = function(RAR = 0, Test = 0, n = 80, IAn = c(20,40,60,80), p0 = 0.2, p1 = 0.4, fixed_ratio = c(1,1), lambda = 0.9, gamma = 0.9, time_trend = FALSE, time_trend_amount = 0.2, n_sim = 1000){
  # Uses TrailSimulator to run presets of different trials, for easy comparison
  # RAR determines which RAR is used
    # RAR = 0, ER
    # RAR = 1, BRAR
    # RAR = 2, Urn
    # RAR = 3, RSHIR, SMLE
    # RAR = 4, RSHIR, ERADE
  # Test determines which test is used (type one error = 10% if necessary to be specified)
    # Test = 0, Logistic
    # Test = 1, BOP2 (calibrated for ddefalt paremters)
    # Test = 2, Randomisation 
    # Test = 3, Wald
  # n is max number of patients
  # IAn are the IA placements
  # p0 and p1 and the control and experimental success rates respectively
  # fixed ratio is the starting (and continuing if ER) allocation
  # lambda and gamma are parameters for BOP2
  # time trend is a boolean representing whether there are time trends
  # time trend amount indiciates how much the arms success rates will change by
  # n_sim is the number of simulations
  
  # n is the maximum number of patients
  # IAn is the timing of the IA, note that this example function doesn't support more or less than 3 IA, then the final IA
  # p0 and p1 are the control and experimental success rates respectively
  
  # Control
  eff_pbo = endpoint(name = "orr", type = "non-tte", 
                     readout = c(orr = 0), generator  = rbinom, size = 1, prob = p0)
  pbo = arm(name = "placebo")
  pbo$add_endpoints(eff_pbo)
  
  # Experimental 
  eff_trt = endpoint(name = "orr", type = "non-tte",  
                     readout = c(orr = 0), generator = rbinom, size = 1, prob = p1)
  trt = arm(name = "treatment")
  trt$add_endpoints(eff_trt)
  
  trial = trial(
    name = "demo", n_patients = n,
    enroller = StaggeredRecruiter,
    accrual_rate = data.frame(end_time = Inf,
                              piecewise_rate = 10) # This example function doesn't support modified accrual rate
  )
  
  trial$add_arms(sample_ratio = fixed_ratio, pbo, trt) # Starts at 1:1 for the burn-in
  trial$save(RAR, name = 'RAR')
  trial$save(Test, name = 'Test')
  
  #IA 1
  source("IA_1.R")
  interim_1 = milestone(name = "IA_1",
                        when = eventNumber("orr", n = IAn[1]),
                        action = IA_1,
                        RAR = RAR, Test = Test, fixed_ratio = fixed_ratio, lambda = lambda, gamma = gamma)
  
  # IA 2
  source("IA_2.R")
  interim_2 = milestone(name = "IA_2",
                        when = eventNumber("orr", n = IAn[2]),
                        action = IA_2,
                        RAR = RAR, Test = Test, fixed_ratio = fixed_ratio, lambda = lambda, gamma = gamma, time_trend = time_trend, time_trend_amount = time_trend_amount, p0 = p0, p1 = p1)
  # IA 3
  source("IA_3.R")
  interim_3 = milestone(name = "IA_3",
                        when = eventNumber("orr", n = IAn[3]),
                        action = IA_3,
                        RAR = RAR, Test = Test, fixed_ratio = fixed_ratio, lambda = lambda, gamma = gamma)
  
  # IA 4/ Final IA
  source("final_action.R")
  final = milestone(name = "final",
                    when = eventNumber("orr", n = n),
                    action = final_action,
                    RAR = RAR, Test = Test, lambda = lambda, gamma = gamma)
  
  listener = listener(silent = TRUE)
  listener$add_milestones(interim_1,
                          interim_2,
                          interim_3,
                          final)
  controller = controller(trial, listener)
  controller$run(n = n_sim, silent = TRUE)
  out = controller$get_output()
  
  return(out)
}