## IA 3

IA_3 = function(trial, RAR = 0, Test = 0, fixed_ratio = c(1,1), lambda = 0.9, gamma = 0.9){
  locked_data = trial$get_locked_data("IA_3")
  #browser()
  
  # BOP2 test
  num_pbo = sum(locked_data$arm == "placebo")
  num_trt = sum(locked_data$arm == "treatment")
  
  succ_pbo = sum(locked_data$orr[locked_data$arm == "placebo"])
  succ_trt = sum(locked_data$orr[locked_data$arm == "treatment"])
  
  fail_pbo = num_pbo - succ_pbo
  fail_trt = num_trt - succ_trt
  
  if (Test == 0){
    # Logistic
    
    fit = fitLogistic(orr ~ arm, placebo = "placebo", data = locked_data, alternative = "greater", scale = "coefficient")
    
    p_val = fit$p
    
    threshold = c(0.05,0.5) # magic numbers, adjust as needed
    
    # For consistency with the BOP2 test
    post_prob = 1 - p_val
    threshold = 1 - threshold 
    
  }
  else if (Test == 1){
    #BOP2
    
    # Lambda, Gamma = 0.9
    post_prob = 1 - rcpp_exact_beta(1 + succ_trt, 1 + fail_trt,1 + succ_pbo, 1 +  fail_pbo)
    
    threshold = Threshold_prob(lamda = lambda, gamma = gamma, 80,num_pbo + num_trt)
    #print(threshold)
    
  }
  else if (Test == 2){
    # Radomisation test - swap the labels and output and see if its significant
    
    cur_diff = succ_trt - succ_pbo # Current difference
    
    num_succ = succ_trt + succ_pbo
    num_fail = fail_trt + fail_pbo
    
    to_reorder = c(rep(1,num_succ),rep(0,num_fail))
    diff_vec = numeric(10000)
    
    for (i in 1:10000){
      temp_reorder = sample(to_reorder, num_succ + num_fail)
      temp_succ_trt = sum(temp_reorder[1:num_trt])
      #temp_succ_pbo = num_succ - temp_succ_trt
      
      diff_vec[i] = 2*temp_succ_trt - num_succ
      
    }
    
    p_val = mean(diff_vec >= cur_diff)
    
    threshold = c(0.05,0.5) # magic numbers, adjust as needed
    
    # For consistency with the BOP2 test
    post_prob = 1 - p_val
    threshold = 1 - threshold 
    
  }
  else if (Test == 3){
    # Wald 
    pbo_mle = succ_pbo/(num_pbo)
    trt_mle = succ_trt/(num_trt)
    
    var_pbo = pbo_mle*(1-pbo_mle)*num_pbo
    var_trt = trt_mle*(1-trt_mle)*num_trt
    
    Z_val = (trt_mle - pbo_mle) /(var_pbo + var_trt)
    
    p_val = 1 - pnorm(Z_val)
    
    threshold = c(0.05,0.5) # magic numbers, adjust as needed
    
    # For consistency with the BOP2 test
    post_prob = 1 - p_val
    threshold = 1 - threshold 
    
  }
  
  if (post_prob > threshold[2]){
    # End early for supe
    
    now = trial$get_current_time()
    trial$update_milestone("final", when  = calendarTime(time = now)) # end the trial now
    trial$save(1, name = 'reject_IA_3')
    
  }
  else if (post_prob < threshold[1]){
    # end early for fute
    
    now = trial$get_current_time()
    trial$update_milestone("final", when  = calendarTime(time = now)) # end the trial now
    trial$save(0, name = 'reject_IA_3')
    
  }
  else{
    # Continue
    # Update the RAR
    if (RAR == 0){
      # ER
      
      r = fixed_ratio
    }
    else if (RAR == 1) {
      # BRAR
      
      c = (num_pbo + num_trt)/(2*n) # tuning parameter
      
      post_prob = 1 - rcpp_exact_beta(1 + succ_trt, 1 + fail_trt,1 + succ_pbo, 1 +  fail_pbo)
      
      allo_prob = post_prob^c/(post_prob^c + (1 - post_prob)^c)
      allo_prob = floor((allo_prob*20))/20 # So rounding always favours the null
      r = c(allo_prob, 1 - allo_prob)
    }
    else if (RAR == 2) {
      # Urn
      
      # Start with a ball in each container
      control_ball = 1 + succ_pbo + fail_trt
      exp_ball = 1 + succ_trt + fail_pbo
      
      r = c(exp_ball, control_ball)/(control_ball + exp_ball)
      
    }
    else if (RAR == 3) {
      # RSHIRH, SMLE
      
      pbo_mle = succ_pbo/(num_pbo)
      trt_mle = succ_trt/(num_trt)
      
      target_allo = sqrt(trt_mle)/(sqrt(pbo_mle) + sqrt(trt_mle))
      
      r = c(target_allo, 1 - target_allo) # Don't need to do anything in addition for SMLE
      
    }
    else if (RAR == 4) {
      
      # RSHIRH, ERADE
      
      pbo_mle = succ_pbo/(num_pbo)
      trt_mle = succ_trt/(num_trt)
      
      target_allo = sqrt(trt_mle)/(sqrt(pbo_mle) + sqrt(trt_mle))
      
      a = 0.5
      
      current_prop = num_trt/(num_trt + num_pbo)
      
      if (current_prop > target_allo){
        # Prop is too high, need to go down
        
        r = c(a*target_allo, 1 - a*target_allo)
        
      }
      else if (current_prop < target_allo){
        # Prop is too low, need to go up
        
        r = c(1 - a*(1 - target_allo) , a*(1 - target_allo))
        
      }
      else{
        # On target, continue
        
        r = c(target_allo, 1 - target_allo) # Don't need to do anything in addition for SMLE
      }
      
      
    }
    
    
    
    
    #r = c(1,1)
    if (r[1] == 0 | r[2] == 0){
      
      r = (r + 1/3) # clipping, means the minmium allocation is 0.2 if it's too extreme
    }
    
    trial$update_sample_ratio(
      arm_names = c("treatment", "placebo"),
      sample_ratios = r)
    
  }
}