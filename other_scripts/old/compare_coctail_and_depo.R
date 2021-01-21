library(deSolve)
library(ggplot2)
library(dplyr)


SetDefaultParameters = function(mg_count = 10^12) {
  # PARAMETERS
  #mutation rate
  epsilon = 10^(-3) #10^(-6)
  #epsilon2= 0.00005
  epsilon2 = 10^(-15) #10^(-20)
  epsilon3 = 0
  
  default_params = 
    list(
      Tmax = 24,
  
  
      # initial values
      # started with 10^6 CFU/mL i.e. 10^9 CFU/L i.e. ~10^(-3) g/L

      B0  = 10^9/mg_count, #[g/L]
      MOI = 10, # so  both phages and bacteria are now counted in units of (10^12) 
      G0 = 13.9,#55.6
      
      # Weitz: K= 4*10^(-6)g/mL = 4*10^(-3)g/L (4miligram per liter ~ 0.01 mmol/L)
      # according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC98929/
      # K for E-Coli is in the range of 40-100 000 miligrams/L i.e. 0.15-500 mmol/L
      
      
      
      # fitted to KP77, first 24 h
      Vh=323.92,
      Kh=500,
      a=0.17,
      
      # fitted to KP77, ffirst 4 h
      #Vh=597.39
      #Kh=500
      #a=0.17#0.04
      
      
      # fitted to Kp486, first 4h
      #Kh=500
      #Vh=809.5
      #a=0.039
      
      # fagi z depolimerazą są mniejsze więc krocej się namnażają
      # burst size
      # Weitz 2005: burst size = 71
      # Da Paeppe: burst size = 50-3600
      # Cairns 2009:: burst size ~ 2
      beta_non_depo = 15,
      beta_depo = 50,
      
      
      
      # adsororption rate to the specific phage
      # Cairns 2009: adsorption ~ 8*10^(-8) ml/(CFU*h). After converting CFU to mg we get: 10^9 * 8*10^(-8) ~ 80
      # Weitz 2005: adsorrption ~ 6.24*10^(-8) ml/[hr*CFU].  i.e. 6.24*10^(-11) L/[hr*CFU] i.e. ~ 62.4 L/[hr * 10^12 phage cells]
      #phi_non_depo = 0.25; # phage 1 (with depolymerase) to bacteria with polisaccharide
      #phi_depo = 0.25; # phage 2 (without depolimerase) to bacteria without polissacharide
      # Da Paeppe 10^(-11) - 10^(-9) [1/min] = 10^(-9)-10^(-7) [1/min] i.e 1-100
      phi_non_depo = 20,
      phi_depo = 20,
      
      # phage decay rate
      #Cairns 2009 decay = 0.0106 [1/h], This constant doesn't depend on the unit of phages. 
      # Levin 2004 decay: 0.1 [1/h]
      # Da Paeppe 0.07-0.5 [1/24h] = 0.003-0.2 [1/h]
      #So it is ok if we count phages in units of 10^12
      decay = 0.0106, 
      
      # depolymerase activity
      depo_decay_rate = 0.2,
      e0 = 0.05,#0.9; # when added depolymerase
      c=0,
      V_depo = 50, #15; # depolymerase max activity
      K_depo = 1, # 5 depolymerase Michaelis Menten constant
      # as in Michaelis Menten kinetics V = e_0 * k_2 where e_0 is the enzyme concentration
      

      epsilonB2toB1 = epsilon,
      epsilonB1toB2 = epsilon,
      epsilonB2toB3 = epsilon2,
      epsilonB3toB2 = epsilon2,
      epsilonB1toB3 = epsilon3,
      epsilonB3toB1 = epsilon3
    )
  
  return(default_params)
}


######### Helper functions ################################

two_simultaneous_phages_and_bacteria_with_two_receptors = function(Time, State, Pars) {
  # B1 - same as B_P bacteria with both receptors sentitive to Pdepo (attaches to receptor A) and P2 (attaches to receptor B)
  # B2 - bacteria with receptor A, and lost receptor B, sensitive to Pdepo (attaches to receptor A)
  # B3 - bacteria with receptor B, and lost receptor A, sensitive to Pnondepo (attaches to the lost receptor A) and P2 (attaches to receptor B)
  # B4 - bacteria without either:lost receptor A and B, sensitive to Pnondepo (attaches to the lost receptor A)
  
  # B2 and B3 sum up to B_N in our main model
  # Pdepo (attaches to receptor A and encodes depolymerase) therefore it makes some of the bacteria loose that receptor and still keep alive if we assume c > 0!
  # P2 (attaches to receptor B)
  # Pnondepo (attaches to the lost receptor A and doesn;t encode depolymerase) 
  
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G); 
    if (bf==0) {
      burst_factor = 1
    } else {
      #burst_factor = 10*Jg/Vh #bacteria will not burst and release new phages when they don't grow
      burst_factor = G/(1+G)
    }
    Gdot    = -Jg*(B1 + B2 + B3 + B4);
    B1dot   = a*Jg*B1 - phi_depo*B1*Pdepo - phi2*B1*P2 - epsilonB1toB2*B1 - epsilonB1toB3*B1 + epsilonB2toB1*B2 + epsilonB3toB1*B3; 
    B2dot   = a*Jg*B2 - phi_depo*B2*Pdepo - epsilonB2toB1*B2 -epsilonB2toB4*B2  + epsilonB1toB2*B1+ epsilonB4toB2*B4; 
    B3dot   = a*Jg*B3 - phi_non_depo*B3*Pnondepo - phi2*B3*P2 - epsilonB3toB4*B3 - epsilonB3toB1*B3 + epsilonB4toB3*B4 + epsilonB1toB3*B1; 
    B4dot   = a*Jg*B4 - phi_non_depo*B4*Pnondepo - epsilonB4toB3*B4 - epsilonB4toB2*B4  + epsilonB2toB4*B2 + epsilonB3toB4*B3; 
    
    Pdepodot  = beta_depo*phi_depo*(B1 +B2)*Pdepo*burst_factor - decay*Pdepo;
    P2dot =  beta2*phi2*(B1 + B3)*P2*burst_factor - decay*P2;
    Pnondepodot =  beta_non_depo*phi_non_depo*(B3 + B4)*Pnondepo*burst_factor - decay*Pnondepo;
    
    batch  = list(c(Gdot,Pdepodot,P2dot,Pnondepodot, B1dot,B2dot,B3dot,B4dot))
  })
}





# # our model of growth (using differential euqutions)
two_phages_and_bacteria = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G); 
    if (bf==0) {
      burst_factor = 1
    } else {
      burst_factor = G/(1+G)#10*Jg/Vh #bacteria will not burst and release new phages when they don't grow
    }
    # depolymerase activity
    Jdepo = E*V_depo*B1/(K_depo+B1); # bacteria with capsule that will become capsuleless because of the dpeolymerase
    #equations during the grow
    Edot    = -depo_decay_rate*E # if we allow some free depolymerase then:+ P1*c;
    Gdot    = -Jg*(B1 + B2 + B3);
    B1dot   = a*Jg*B1 - phi_depo*B1*P1 - epsilonB1toB2*B1 + epsilonB2toB1*B2 - epsilonB1toB3*B1 + epsilonB3toB1*B3 - Jdepo; #those can be infected by p1 only
    B2dot   = a*Jg*B2 - phi_non_depo*B2*P2 - epsilonB2toB1*B2 + epsilonB1toB2*B1 -epsilonB2toB3*B2 + epsilonB3toB2*B3 + Jdepo; #those can be infected by p2 only
    B3dot   = a*Jg*B3 + epsilonB2toB3*B2 - epsilonB3toB2*B3 + epsilonB1toB3*B1 - epsilonB3toB1*B3;  #completely resistant strain
    P1dot  = beta_depo*phi_depo*B1*P1*burst_factor - decay*P1;
    P2dot =  beta_non_depo*phi_non_depo*B2*P2*burst_factor - decay*P2;
    batch  = list(c(Edot, Gdot,P1dot,P2dot,B1dot,B2dot,B3dot))
  })
}




Simulate_Phage_Coctail = function(Vh, Kh,a, beta_non_depo, beta_depo, phi_non_depo, phi_depo, decay, V_depo, K_depo, depo_decay_rate, epsilonB2toB1, epsilonB1toB2,epsilonB2toB3, epsilonB3toB2,epsilonB1toB3, epsilonB3toB1, B0, 
                                  e0=1, MOI =1, G0=29, Tmax = 24,  bf = 1, propPP =0.5, propB2init = 0, ode_method = "ode45", mg_count = (10^12), tspan = 0.2) {
  
  


  
  time = seq(0,Tmax,tspan)
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
                              Treatment = character(0), stringsAsFactors = FALSE)
  
  # No external depolimerase: different combinations of phages
  Treatments = c(           "Phage with no depolymerase",
                             "Phage cocktail" ,
                             "Phage with depolymerase")
  i=0
  for (freq in c(0, propPP, 1)) {
  i=i+1
  inits = c(E= 0,
            G= G0, 
            P1=freq*P0,
            P2=(1-freq)*P0,
            B1= (1-propB2init)*B0,
            B2 = propB2init*B0,
            B3 = 0)
  

  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(Treatment = Treatments[i])
  
  simulated_data = rbind(simulated_data, simulation)
  }
  
  # only pahges without depo and some additoonal depo
  inits['P1']=0 
  inits['P2']=P0
  inits['E'] =e0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(Treatment = paste0("Phage with no depolymerase + external depolymerase"))
  simulated_data = rbind(simulated_data, simulation) 

  
  # only bacteria and depo
  inits['P1'] =  0 
  inits['P2'] =  0
  inits['E'] =e0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(Treatment = "External depolymerase")
  simulated_data = rbind(simulated_data, simulation) 
  
  
  # only bacteria 
  inits['P1'] =  0 
  inits['P2'] =  0
  inits['E'] = 0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2) %>%
    mutate(Treatment = "No phage")
  simulated_data = rbind(simulated_data, simulation) %>%
    mutate(bacteria_capsule = bacteria_capsule*mg_count,
           bacteria_no_capsule = bacteria_no_capsule*mg_count,
           bacteria_resistant = bacteria_resistant*mg_count) %>%
    mutate(all_bacteria_CFU_per_L = bacteria_capsule + bacteria_no_capsule + bacteria_resistant)
    
  


  return(simulated_data)
}

PlotSimulatedPhageAndBacteria = function( simulated_data,
                                          title_plot = "", 
                                          colors = NULL,
                                          ymin=-1*10^11,
                                          ymax = 2.5*10^12,
                                          tmax = 24,
                                          minCFU = 6.3*10^9,
                                          text_size = 6,
                                          linetypes = NULL) {
  
  my_theme = theme(
    panel.grid.major = element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text_size),
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black")
  )
  
  if (is.null(colors)) {
    colors = c("No phage" = "black",
               "External depolymerase" = "gray",
               "Phage with depolymerase" = "darkgreen",
               "Phage with no depolymerase" = "darkblue",
               "Phage with no depolymerase + external depolymerase" = "darkmagenta",
               "Phage cocktail" = "darkred")
  }
  if (is.null(linetypes)) {
    linetypes = c("No phage" = "solid",
               "External depolymerase" = "solid",
               "Phage with depolymerase" = "dotted",
               "Phage with no depolymerase" = "dashed",
               "Phage with no depolymerase + external depolymerase" = "solid",
               "Phage cocktail" = "longdash")
  }
  descriptions_to_show = names(colors)
  g1=ggplot(simulated_data %>% 
              filter(Treatment %in% descriptions_to_show) %>%
              mutate(all_bacteria_CFU_per_L = ifelse(all_bacteria_CFU_per_L > minCFU, all_bacteria_CFU_per_L, minCFU)),
            aes(x = time, y = all_bacteria_CFU_per_L, col = Treatment, linetype = Treatment)) +
    scale_color_manual(values=colors) +
    scale_linetype_manual(values = linetypes) +
    geom_line() +
    #theme_bw() +
    ggtitle(title_plot) +
    xlab("time [h]") +
    ylab("bacteria [CFU/L]") + 
    scale_y_log10(limits = c(ymin, ymax)) +
    my_theme +
    xlim(c(0,tmax))
  return(g1)
  }
  
