library(deSolve)
library(ggplot2)
library(dplyr)


######### Helper functions ################################
# # our model of growth (using differential euqutions)
two_phages_and_bacteria_with_delay = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G); 
    growth_factor = 5*Jg/Vh #ifelse(G > 5, 1, 0)
    # depolymerase activity
    Jdepo = (e_0 + P1*c)*V_depo*B1/(K_depo+B1); # bacteria with capsule that will become capsuleless because of the dpeolymerase
    #equations during the grow
    Gdot   = -Jg*(B1 + B2);
    B1dot   = a*Jg*B1 - phi_depo*B1*P1free - epsilonB1toB2*B1 + epsilonB2toB1*B2 - epsilonB1toB3*B1 + epsilonB3toB1*B3 - Jdepo; #those can be infected by p1 only
    B2dot   = a*Jg*B2 - phi_non_depo*B2*P2free - epsilonB2toB1*B2 + epsilonB1toB2*B1 -epsilonB2toB3*B2 + epsilonB3toB2*B3 + Jdepo; #those can be infected by p2 only
    B3dot   = a*Jg*B3 + epsilonB2toB3*B2 - epsilonB3toB2*B3 + epsilonB1toB3*B1 - epsilonB3toB1*B3;  #completely resistant strain
    P1freedot   = beta_depo*P1*growth_factor - decay*P1free - phi_depo*B1*P1free; # free phages: burst if bacteria divide: that's why we multiply by Jg
    P1dot       = phi_depo*B1*P1free - P1*growth_factor;   # phages locked inside bacteria
    P2freedot   = beta_non_depo*P2*growth_factor - decay*P2free - phi_non_depo*B2*P2free; # free phages: burst if bacteria divide: that's why we multiply by Jg
    P2dot       = phi_non_depo*B2*P2free - P2*growth_factor;   # phages locked inside bacteria
    batch       = list(c(Gdot,B1dot,B2dot,B3dot,P1dot,P2dot,P1freedot, P2freedot))
  })
}



Simulate_Phage_Coctail_with_delay = function(Vh, Kh,a, beta_non_depo, beta_depo, phi_non_depo, phi_depo, decay, V_depo, K_depo, c, epsilonB2toB1, epsilonB1toB2,epsilonB2toB3, epsilonB3toB2,epsilonB1toB3, epsilonB3toB1, B0, MOI =1, G0=29, Tmax = 24, title_plot = "") {
  time = seq(0,Tmax,1)
  e_0 = 0
  pars =  c(a = a, 
            Vh = Vh, 
            Kh = Kh,
            c=c,
            beta_non_depo = beta_non_depo, 
            beta_depo = beta_depo, 
            phi_non_depo = phi_non_depo, 
            phi_depo = phi_depo,
            decay = decay, 
            V_depo = V_depo, 
            K_depo = K_depo, 
            e_0 =e_0, 
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
                  "All bacteria: 50% P_N & 50% P_P",
                  "All bacteria: with P_P only")
  i=0
  for (freq in c(0, 0.5, 1)) {
  i=i+1
  pars['e_0'] = 0
  inits = c(G= G0, 
            B1= B0,
            B2 = 0,
            B3 = 0,
            P1 = 0,
            P2 = 0,
            P1free = freq*P0,
            P2free = (1-freq)*P0)
  
  
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria_with_delay, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1free, phage_no_depo = P2free) %>%
    mutate(description = decriptions[i])
  
  simulated_data = rbind(simulated_data, simulation)
  }
  
  # only pahges without depo and some additoonal depo
  inits['P1free']=0 
  inits['P2free']= P0
  pars['e_0'] = 1
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1free, phage_no_depo = P2free) %>%
    mutate(description = paste0("All bacteria: P_N plus depo"))
  simulated_data = rbind(simulated_data, simulation)
  
  # only bacteria and dpeo
  inits['P1free'] =  0 
  inits['P2free'] =  0
  pars['e_0'] = 0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1free, phage_no_depo = P2free) %>%
    mutate(description = "All bacteria: no phage")
  simulated_data = rbind(simulated_data, simulation) %>%
    mutate(bacteria = bacteria_capsule + bacteria_no_capsule + bacteria_resistant)
  
  
  g1=ggplot(simulated_data,
         aes(x = time, y = bacteria, col = description)) +
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


