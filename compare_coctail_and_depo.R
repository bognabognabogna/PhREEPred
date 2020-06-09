library(deSolve)
library(ggplot2)
library(dplyr)


######### Helper functions ################################
# # our model of growth (using differential euqutions)
two_phages_and_bacteria = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G); 
    if (bf==0) {
      burst_factor = 1
    } else {
      burst_factor = 10*Jg/Vh #bacteria will not burst and release new phages when they don't grow
    }
    # depolymerase activity
    Jdepo = E*V_depo*B1/(K_depo+B1); # bacteria with capsule that will become capsuleless because of the dpeolymerase
    #equations during the grow
    Edot    = -depo_decay_rate*E + P1*c;
    Gdot    = -Jg*(B1 + B2 + B3);
    B1dot   = a*Jg*B1 - phi_depo*B1*P1 - epsilonB1toB2*B1 + epsilonB2toB1*B2 - epsilonB1toB3*B1 + epsilonB3toB1*B3 - Jdepo; #those can be infected by p1 only
    B2dot   = a*Jg*B2 - phi_non_depo*B2*P2 - epsilonB2toB1*B2 + epsilonB1toB2*B1 -epsilonB2toB3*B2 + epsilonB3toB2*B3 + Jdepo; #those can be infected by p2 only
    B3dot   = a*Jg*B3 + epsilonB2toB3*B2 - epsilonB3toB2*B3 + epsilonB1toB3*B1 - epsilonB3toB1*B3;  #completely resistant strain
    P1dot  = beta_depo*phi_depo*B1*P1*burst_factor - decay*P1;
    P2dot =  beta_non_depo*phi_non_depo*B2*P2*burst_factor - decay*P2;
    batch  = list(c(Edot, Gdot,P1dot,P2dot,B1dot,B2dot,B3dot))
  })
}



Simulate_Phage_Coctail = function(Vh, Kh,a, beta_non_depo, beta_depo, phi_non_depo, phi_depo, decay, V_depo, K_depo, c, depo_decay_rate, epsilonB2toB1, epsilonB1toB2,epsilonB2toB3, epsilonB3toB2,epsilonB1toB3, epsilonB3toB1, B0, 
                                  e0=1, MOI =1, G0=29, Tmax = 24, title_plot = "", colors = NULL, bf = 1, propPP =0.5) {
  

  if (is.null(colors)) {
    colors = c("All bacteria: no phage" = "black",
               "All bacteria: external depo" = "gray",
               "All bacteria: with P_P only" = "green",
               "All bacteria: with P_N only" = "blue",
               "All bacteria: with P_N plus depo" = "magenta",
               "All bacteria: P_N & P_P" = "red")
  }
  
  time = seq(0,Tmax,1)
  e_0 = 0
  pars =  c(a = a, 
            Vh = Vh, 
            Kh = Kh,
            beta_non_depo = beta_non_depo, 
            beta_depo = beta_depo, 
            phi_non_depo = phi_non_depo, 
            phi_depo = phi_depo,
            decay = decay, 
            V_depo = V_depo, 
            K_depo = K_depo, 
            depo_decay_rate =depo_decay_rate, 
            c=c,
            bf = bf,
            epsilonB2toB1 = epsilonB2toB1, 
            epsilonB1toB2 = epsilonB1toB2,
            epsilonB2toB3 = epsilonB2toB3,
            epsilonB3toB2 = epsilonB3toB2,
            epsilonB1toB3 = epsilonB1toB3,
            epsilonB3toB1 = epsilonB3toB1)


  P0=B0*MOI;
  
  
  simulated_data = data.frame(time = numeric(0), 
                              bacteria_capsule = numeric(0), 
                              bacteria_no_capcule = numeric(0),
                              description = character(0), stringsAsFactors = FALSE)
  
  # No external depolimerase: different combinations of phages
  decriptions = c("All bacteria: with P_N only", 
                  "All bacteria: P_N & P_P",
                  "All bacteria: with P_P only")
  i=0
  for (freq in c(0, propPP, 1)) {
  i=i+1
  inits = c(E= 0,
            G= G0, 
            P1=freq*P0,
            P2=(1-freq)*P0,
            B1= B0,
            B2 = 0,
            B3 = 0)
  

  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(description = decriptions[i])
  
  simulated_data = rbind(simulated_data, simulation)
  }
  
  # only pahges without depo and some additoonal depo
  inits['P1']=0 
  inits['P2']=P0
  inits['E'] =e0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(description = paste0("All bacteria: with P_N plus depo"))
  simulated_data = rbind(simulated_data, simulation) 

  
  # only bacteria and depo
  inits['P1'] =  0 
  inits['P2'] =  0
  inits['E'] =e0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(description = "All bacteria: external depo")
  simulated_data = rbind(simulated_data, simulation) 
  
  
  # only bacteria 
  inits['P1'] =  0 
  inits['P2'] =  0
  inits['E'] = 0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(description = "All bacteria: no phage")
  simulated_data = rbind(simulated_data, simulation) %>%
    mutate(bacteria = bacteria_capsule + bacteria_no_capsule + bacteria_resistant)
  
  
  g1=ggplot(simulated_data,
         aes(x = time, y = bacteria, col = description)) +
    scale_color_manual(values=colors) +
    geom_line() +
    theme_bw() +
    ggtitle(title_plot)
  
  print(g1)
  
  g2=ggplot(simulated_data,
            aes(x = time, col = description)) +
    geom_line(aes(y = bacteria_capsule), linetype = "solid") +
    geom_line(aes(y = bacteria_no_capsule), linetype = "dashed") +
    geom_line(aes(y = bacteria_resistant), linetype = "dotted") +
    theme_bw() +
    ggtitle(title_plot)
  
  #print(g2)
  
  g3=ggplot(simulated_data,
            aes(x = time, col = description)) +
    geom_line(aes(y = phage_depo), linetype = "solid") +
    geom_line(aes(y = phage_no_depo), linetype = "dashed") +
    theme_bw() +
    ggtitle(title_plot)
  #print(g3)
  return(invisible(NULL))
}
