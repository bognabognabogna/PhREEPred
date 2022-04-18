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

      B0  = 10^9/mg_count, #[g/L] = [mg/mL]
      MOI = 10, # so  both phages and bacteria are now counted in units of (10^12 / L) 
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
      # Da Paeppe 10^(-11) - 10^(-9) [1/min] = 10^(-9)-10^(-7) [1/min] i.e 1-100
      phi_non_depo = 2*20,
      phi_depo = 2*20,
      
      # phage decay rate
      #Cairns 2009 decay = 0.0106 [1/h], This constant doesn't depend on the unit of phages. 
      # Levin 2004 decay: 0.1 [1/h]
      # Da Paeppe 0.07-0.5 [1/24h] = 0.003-0.2 [1/h]
      #So it is ok if we count phages in units of 10^12
      decay = 0.0106, 
      
      # depo activity
      depo_decay_rate = 0.2,
      e0 = 0.05,#0.9; # when added depo
      c=0,
      V_depo = 50, #15; # depo max activity
      K_depo = 1, # 5 depo Michaelis Menten constant
      # as in Michaelis Menten kinetics V = e_0 * k_2 where e_0 is the enzyme concentration
      

      epsilonB2toB1 = epsilon,
      epsilonB1toB2 = epsilon,
      epsilonB2toB3 = epsilon2,
      epsilonB3toB2 = epsilon2,
      epsilonB1toB3 = epsilon3,
      epsilonB3toB1 = epsilon3,
      latent_period_depo = 15/60,
      latent_period_non_depo = 25/60,
      propB2init = 0
    )
  
  return(default_params)
}

SetMaxParameters = function(mg_count = 10^12) {
  Epsilon = 10^15
  max_params = 
    list(
      Tmax = 240,
      B0  = 10^15/mg_count, 
      MOI = 10^6, # so  both phages and bacteria are now counted in units of (10^12 / L) 
      G0 = 1000,
      Vh = Epsilon,
      Kh = Epsilon,
      a = 1,
      beta_non_depo = 10^4,
      beta_depo = 10^4,
      phi_non_depo = 10^4,
      phi_depo = 10^4,
      decay = 10^4, 
      depo_decay_rate = 10^4,
      e0 = 10^4,
      c=10^4,
      V_depo = Epsilon,
      K_depo = Epsilon,
      epsilonB2toB1 = 1,
      epsilonB1toB2 = 1,
      epsilonB2toB3 = 1,
      epsilonB3toB2 = 1,
      epsilonB1toB3 = 1,
      epsilonB3toB1 = 1
    )
  
  return(max_params)
}


SetMinParameters = function(mg_count = 10^12) {

  # some min value so that we don't get a numerical error
  epsilon = 10^(-16);
  min_params = 
    list(
      Tmax = 0,
      B0  = 0, #[g/L] = [mg/mL]
      MOI = 0, # so  both phages and bacteria are now counted in units of (10^12 / L) 
      G0 = 0,
      Vh = 0,
      Kh = epsilon,
      a = 0,
      beta_non_depo = 0,
      beta_depo = 0,
      phi_non_depo = 0,
      phi_depo = 0,
      decay = 0, 
      depo_decay_rate = 0,
      e0 = 0,
      c=0,
      V_depo = 0,
      K_depo = epsilon,
      epsilonB2toB1 = 0,
      epsilonB1toB2 = 0,
      epsilonB2toB3 = 0,
      epsilonB3toB2 = 0,
      epsilonB1toB3 = 0,
      epsilonB3toB1 = 0
    )
  
  return(min_params)
}


######### Helper functions ################################

two_simultaneous_phages_and_bacteria_with_two_receptors = function(Time, State, Pars) {
  # B1 - same as B_P bacteria with both receptors sentitive to Pdepo (attaches to receptor A) and P2 (attaches to receptor B)
  # B2 - bacteria with receptor A, and lost receptor B, sensitive to Pdepo (attaches to receptor A)
  # B3 - bacteria with receptor B, and lost receptor A, sensitive to Pnondepo (attaches to the lost receptor A) and P2 (attaches to receptor B)
  # B4 - bacteria without either:lost receptor A and B, sensitive to Pnondepo (attaches to the lost receptor A)
  
  # B2 and B3 sum up to B_N in our main model
  # Pdepo (attaches to receptor A and encodes depo) therefore it makes some of the bacteria loose that receptor and still keep alive if we assume c > 0!
  # P2 (attaches to receptor B)
  # Pnondepo (attaches to the lost receptor A and doesn;t encode depo) 
  
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



# # our model of growth (using differential euqutions)
two_phages_and_bacteria = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G); 
    if (bf==0) {
      burst_factor = 1
    } else {
      burst_factor = G/(1+G)#10*Jg/Vh #bacteria will not burst and release new phages when they don't grow
    }
    # depo activity
    Jdepo = E*V_depo*B1/(K_depo+B1); # bacteria with capsule that will become capsuleless because of the dpeolymerase
    #equations during the grow
    Edot    = -depo_decay_rate*E # if we allow some free depo then:+ P1*c;
    Gdot    = -Jg*(B1 + B2 + B3);
    B1dot   = a*Jg*B1 - phi_depo*B1*P1 - epsilonB1toB2*B1 + epsilonB2toB1*B2 - epsilonB1toB3*B1 + epsilonB3toB1*B3 - Jdepo; #those can be infected by p1 only
    B2dot   = a*Jg*B2 - phi_non_depo*B2*P2 - epsilonB2toB1*B2 + epsilonB1toB2*B1 -epsilonB2toB3*B2 + epsilonB3toB2*B3 + Jdepo; #those can be infected by p2 only
    B3dot   = a*Jg*B3 + epsilonB2toB3*B2 - epsilonB3toB2*B3 + epsilonB1toB3*B1 - epsilonB3toB1*B3;  #completely resistant strain
    P1dot  = beta_depo*phi_depo*B1*P1*burst_factor - decay*P1;
    P2dot =  beta_non_depo*phi_non_depo*B2*P2*burst_factor - decay*P2;
    Gdot_used_by_decapsulated_only = -Jg*(B2+B3);
    Gdot_used_by_capsulated_only = -Jg*B1;
    killed_bacteria_dot =  phi_depo*B1*P1 + phi_non_depo*B2*P2;
    batch  = list(c(Edot, Gdot,P1dot,P2dot,B1dot,B2dot,B3dot,Gdot_used_by_decapsulated_only, Gdot_used_by_capsulated_only, killed_bacteria_dot))
  })
}

# # our model of growth (using differential euqutions)
two_phages_and_bacteria_delayed = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    if (Time < latent_period_depo) {
      B1L =(1-propB2init)*B0
      P1L =freq*P0
    }
    else {
      B1L = lagvalue(Time - latent_period_depo,5)
      P1L = lagvalue(Time - latent_period_depo,3)
    }
    
    if (Time < latent_period_non_depo) {
      P2L =(1-freq)*P0
      B2L =propB2init*B0
    }
    else {
      B2L = lagvalue(Time - latent_period_non_depo,6)
      P2L = lagvalue(Time - latent_period_non_depo,4)
    }
    
    Jg =Vh*G/(Kh+G); 
    if (bf==0) {
      burst_factor = 1
    } else {
      burst_factor = G/(1+G)#10*Jg/Vh #bacteria will not burst and release new phages when they don't grow
    }
    # depo activity
    Jdepo = E*V_depo*B1/(K_depo+B1); # bacteria with capsule that will become capsuleless because of the dpeolymerase
    #equations during the grow
    Edot    = -depo_decay_rate*E # if we allow some free depo then:+ P1*c;
    Gdot    = -Jg*(B1 + B2 + B3);
    P1dot  = -phi_depo*B1*P1 + beta_depo*phi_depo*B1L*P1L*burst_factor - decay*P1;
    P2dot =  -phi_non_depo*B2*P2 + beta_non_depo*phi_non_depo*B2L*P2L*burst_factor - decay*P2;
    B1dot   = a*Jg*B1 - phi_depo*B1*P1 - epsilonB1toB2*B1 + epsilonB2toB1*B2 - epsilonB1toB3*B1 + epsilonB3toB1*B3 - Jdepo; #those can be infected by p1 only
    B2dot   = a*Jg*B2 - phi_non_depo*B2*P2 - epsilonB2toB1*B2 + epsilonB1toB2*B1 -epsilonB2toB3*B2 + epsilonB3toB2*B3 + Jdepo; #those can be infected by p2 only
    B3dot   = a*Jg*B3 + epsilonB2toB3*B2 - epsilonB3toB2*B3 + epsilonB1toB3*B1 - epsilonB3toB1*B3;  #completely resistant strain
    Gdot_used_by_decapsulated_only = -Jg*(B2+B3);
    Gdot_used_by_capsulated_only = -Jg*(B1);
    killed_bacteria_dot =  phi_depo*B1*P1 + phi_non_depo*B2*P2;
    batch  = list(c(Edot, Gdot,P1dot,P2dot,B1dot,B2dot,B3dot, Gdot_used_by_decapsulated_only, Gdot_used_by_capsulated_only, killed_bacteria_dot))
  })
}


Simulate_Phage_Coctail = function(Vh, Kh,a, beta_non_depo, beta_depo, phi_non_depo, phi_depo, decay, V_depo, K_depo, depo_decay_rate, epsilonB2toB1, epsilonB1toB2,epsilonB2toB3, epsilonB3toB2,epsilonB1toB3, epsilonB3toB1, B0, 
                                  e0=1, MOI =1, G0=29, Tmax = 24,  bf = 1, propPP =0.5, propB2init = 0, ode_method = "ode45", mg_count = (10^12), tspan = 0.2, 
                                  model = "delayed", latent_period_depo = 15/60, latent_period_non_depo = 25/60,
                                  Treatments.included = c(1,1,1,1,1,1)) {
  

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
            epsilonB3toB1 = epsilonB3toB1,
            latent_period_depo = latent_period_depo,
            latent_period_non_depo = latent_period_non_depo,
            propB2init = propB2init)


  P0=B0*MOI;
  
  
  simulated_data = data.frame(time = numeric(0), 
                              bacteria_capsule = numeric(0), 
                              bacteria_no_capcule = numeric(0),
                              Treatment = character(0), 
                              Depolymerase = numeric(0), 
                              Glucose = numeric(0),
                              stringsAsFactors = FALSE)
  
  # No external depolimerase: different combinations of phages
  Treatments = c(           "capsule-independent phage (KP15/KP27)",
                             "phage cocktail (KP15/KP27 + KP34)" ,
                             "depo-equipped phage (KP34)")
  Treatments.data = data.frame(treatment = Treatments) 
  Treatments.data$freq = c(0, propPP, 1)

  for (i in 1:nrow(Treatments.data)) {
  #for (freq in c(0, propPP, 1)) {
  if (Treatments.included[i] == 1){
  freq = Treatments.data$freq[i]
  inits = c(E= 0,
            G= G0, 
            P1=freq*P0,
            P2=(1-freq)*P0,
            B1= (1-propB2init)*B0,
            B2 = propB2init*B0,
            B3 = 0,
            G_for_decapsulated = G0,
            G_for_capsulated = G0,
            killed_bacteria = 0)
  

  if (model == "standard") {
    yout = ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)   
  } else if (model == "delayed") {
    parsdede = pars
    parsdede$freq = freq
    parsdede$propB2init = propB2init
    parsdede$P0 = P0
    parsdede$B0 = B0
    yout = dede(inits, time, two_phages_and_bacteria_delayed, parsdede)
    print("evaluating delayed model")
  }
  simulation <- as.data.frame(yout) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2, Depolymerase = E, Glucose = G, G_for_decapsulated = G_for_decapsulated, G_for_capsulated = G_for_capsulated, killed_bacteria = killed_bacteria) %>%
    mutate(Treatment = Treatments[i])
  
  simulated_data = rbind(simulated_data, simulation)
  }
  }
  
  # only pahges without depo and some additoonal depo
  if (Treatments.included[4] == 1){
  inits['P1']=0 
  inits['P2']=P0
  inits['E'] =e0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2, Depolymerase = E, Glucose = G, G_for_decapsulated = G_for_decapsulated, G_for_capsulated = G_for_capsulated, killed_bacteria = killed_bacteria) %>%
    mutate(Treatment = paste0("capsule-independent phage (KP15/KP27) + depo"))
  simulated_data = rbind(simulated_data, simulation) 
}
  
  # only bacteria and depo
  if (Treatments.included[5] == 1){
  inits['P1'] =  0 
  inits['P2'] =  0
  inits['E'] =e0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2, Depolymerase = E, Glucose = G, G_for_decapsulated = G_for_decapsulated, G_for_capsulated = G_for_capsulated, killed_bacteria = killed_bacteria) %>%
    mutate(Treatment = "External depo")
  simulated_data = rbind(simulated_data, simulation) 
  }
  
  # only bacteria 
  if (Treatments.included[6] == 1){
  inits['P1'] =  0 
  inits['P2'] =  0
  inits['E'] = 0
  simulation <- as.data.frame(ode(inits, time, two_phages_and_bacteria, pars, method = ode_method)) %>% 
    select(time, bacteria_capsule = B1, bacteria_no_capsule = B2, bacteria_resistant = B3, phage_depo = P1, phage_no_depo = P2, Depolymerase = E, Glucose = G, G_for_decapsulated = G_for_decapsulated, G_for_capsulated = G_for_capsulated, killed_bacteria = killed_bacteria) %>%
    mutate(Treatment = "no phage")
  simulated_data = rbind(simulated_data, simulation) 
  }
  
  simulated_data = simulated_data %>%
    mutate(bacteria_capsule = bacteria_capsule*mg_count/1000,
           bacteria_no_capsule = bacteria_no_capsule*mg_count/1000,
           bacteria_resistant = bacteria_resistant*mg_count/1000,
           killed_bacteria = killed_bacteria*mg_count/1000) %>%
    mutate(all_bacteria_CFU_per_mL = bacteria_capsule + bacteria_no_capsule + bacteria_resistant)
    
  


  return(simulated_data)
}


GetVisualisationConstants2 = function(text_size=12, 
                                      no_phage_name = "no phage",
                                      depo_phage_name = "depo-equipped phage (KP34)" ,
                                      no_depo_phage_name = "capsule-independent phage (KP15/KP27)",
                                      no_depo_phage_plus_depo_name = "capsule-independent phage (KP15/KP27) + depo" ,
                                      phage_cocktail_name = "phage cocktail (KP15/KP27 + KP34)") {
  my_theme = theme(
    panel.grid.major = element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text_size),
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black")
  )
  
  colors = list()
  
  
  
  linetypes = c(no_phage_name = "solid",
                #"External depo" = "solid",
                depo_phage_name = "solid",
                no_depo_phage_name= "dashed",
                no_depo_phage_plus_depo_name = "solid",
                phage_cocktail_name = "solid")
  
  
  linesizes = c(no_phage_name = 0.5,
                #"External depo" = 1,
                depo_phage_name = 0.5,
                no_depo_phage_name = 0.5,
                no_depo_phage_plus_depo_name = 0.5,
                phage_cocktail_name = 0.5)
  
  return(list(my_theme = my_theme, 
              colors = colors, 
              linetypes = linetypes, 
              linesizes = linesizes))
}



PlotSimulatedPhageAndBacteria = function( simulated_data,
                                          title_plot = "", 
                                          colors = NULL,
                                          ymin=-1*10^8,
                                          ymax = 2.5*10^9,
                                          tmax = 24,
                                          minCFU = 6.3*10^6,
                                          text_size = 12,
                                          linetypes = NULL,
                                          linesizes = NULL) {
  
  VisualisationConstants = GetVisualisationConstants(text_size)
  my_theme = VisualisationConstants$my_theme
  if (is.null(colors)) {colors = VisualisationConstants$colors}
  if (is.null(linetypes)) {linetypes = VisualisationConstants$linetypes}
  if (is.null(linesizes)) {linesizes = VisualisationConstants$linesizes}
  descriptions_to_show = names(colors)
  
  data_polished = simulated_data %>% 
    filter(Treatment %in% descriptions_to_show) %>%
    mutate(all_bacteria_CFU_per_mL = ifelse(all_bacteria_CFU_per_mL > minCFU, all_bacteria_CFU_per_mL, minCFU))
  
  g1=ggplot(data_polished,
            aes(x = time, 
                y = all_bacteria_CFU_per_mL, 
                col = Treatment, 
                linetype = Treatment, 
                size = Treatment)) +
    geom_line() +
    ggtitle(title_plot) +
    xlab("time [h]") +
    ylab("bacteria [CFU/mL]") + 
    scale_color_manual(values=colors) +
    scale_linetype_manual(values = linetypes) +
    scale_size_manual(values = linesizes) +
    scale_y_log10(limits = c(ymin, ymax)) +
    my_theme +
    xlim(c(0,tmax))
  
  return(g1)
}


PlotSimulatedDataByBacteriaType = function( simulated_data,
                                          title_plot = "", 
                                          colors = NULL,
                                          ymin=-1*10^8,
                                          ymax = 2.5*10^9,
                                          tmax = 24,
                                          minCFU = 6.3*10^6,
                                          text_size = 12,
                                          linetypes = NULL,
                                          linesizes = NULL) {
  
  VisualisationConstants = GetVisualisationConstants(text_size)
  my_theme = VisualisationConstants$my_theme
  if (is.null(colors)) {colors = VisualisationConstants$colors}
  if (is.null(linetypes)) {linetypes = VisualisationConstants$linetypes}
  if (is.null(linesizes)) {linesizes = VisualisationConstants$linesizes}
  descriptions_to_show = names(colors)
  
  data_polished = simulated_data %>% 
    filter(Treatment %in% descriptions_to_show) %>%
    select(-all_bacteria_CFU_per_mL, -bacteria_resistant) %>%
    tidyr::gather(key = "Bacteria Type", value = CFU_per_ml, 
                  bacteria_capsule, bacteria_no_capsule) %>%
    mutate(CFU_per_ml = ifelse(CFU_per_ml > minCFU, CFU_per_ml, minCFU))
    
  
  g1=ggplot(data_polished,
            aes(x = time, 
                y = CFU_per_ml, 
                col = Treatment, 
                linetype = Treatment, 
                size = Treatment)) +
    geom_line() +
    ggtitle(title_plot) +
    xlab("time [h]") +
    ylab("bacteria [CFU/mL]") + 
    scale_color_manual(values=colors) +
    scale_linetype_manual(values = linetypes) +
    scale_size_manual(values = linesizes) +
    scale_y_log10(limits = c(ymin, ymax)) +
    my_theme +
    xlim(c(0,tmax)) +
    facet_grid(`Bacteria Type` ~.)
  
  return(g1)
}



PlotResultsByBacteriumType =function( simulated_data,
                                     title_plot = "", 
                                     ymin=-1*10^8,
                                     ymax = 2.5*10^9,
                                     tmax = 24,
                                     minCFU = 6.3*10^6,
                                     text_size = 12) {
  VisualisationConstants = GetVisualisationConstantsApp(text_size)
  my_theme = VisualisationConstants$my_theme
  data_polished = simulated_data %>% 
  mutate(all_bacteria_CFU_per_mL = ifelse(all_bacteria_CFU_per_mL > minCFU, all_bacteria_CFU_per_mL, minCFU))

  
  data_long = data_polished %>%
    select(-all_bacteria_CFU_per_mL, -bacteria_resistant) %>%
    tidyr::gather(key = "Bacteria Type", value = CFU_per_ml, 
                  bacteria_capsule, bacteria_no_capsule) %>%
    mutate(CFU_per_ml = ifelse(CFU_per_ml > minCFU, CFU_per_ml, minCFU))

    

g1=ggplot(data_long,
          aes(x = time, 
              y = CFU_per_ml, 
              col = `Bacteria Type`)) +
  geom_line() +
  ggtitle(title_plot) +
  xlab("time [h]") +
  ylab("bacteria [CFU/mL]") + 
  #scale_color_manual(values=colors) +
  my_theme +
  xlim(c(0,tmax)) +
  scale_y_log10(limits = c(ymin, ymax)) +
  facet_wrap(facets = "Treatment", ncol =2)
return(g1)
}



GetVisualisationConstantsApp = function(text_size=12) {
  my_theme = theme(
    panel.grid.major = element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text_size),
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black")
  )
  
  
  colors = c("no phage" = "black",
             #"External depo" = "gray",
             "depo-equipped phage" = "darkgreen",
             "capsule-independent phage" = "darkblue",
             "capsule-independent phage + depo" = "darkmagenta",
             "phage cocktail" = "darkred")
  
  
  linetypes = c("no phage" = "solid",
                #"External depo" = "solid",
                "depo-equipped phage" = "solid",
                "capsule-independent phage" = "dashed",
                "capsule-independent phage + depo" = "solid",
                "phage cocktail" = "solid")
  
  
  linesizes = c("no phage" = 0.5,
                #"External depo" = 1,
                "depo-equipped phage" = 0.5,
                "capsule-independent phage" = 0.5,
                "capsule-independent phage + depo" = 0.5,
                "phage cocktail" = 0.5)
  
  return(list(my_theme = my_theme, 
              colors = colors, 
              linetypes = linetypes, 
              linesizes = linesizes))
}


GetVisualisationConstants = function(text_size=12) {
  my_theme = theme(
    panel.grid.major = element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text_size),
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black")
  )
  

    colors = c("no phage" = "black",
               #"External depo" = "gray",
               "depo-equipped phage (KP34)" = "darkgreen",
               "capsule-independent phage (KP15/KP27)" = "darkblue",
               "capsule-independent phage (KP15/KP27) + depo" = "darkmagenta",
               "phage cocktail (KP15/KP27 + KP34)" = "darkred")
  

    linetypes = c("no phage" = "solid",
                  #"External depo" = "solid",
                  "depo-equipped phage (KP34)" = "dotted",
                  "capsule-independent phage (KP15/KP27)" = "dashed",
                  "capsule-independent phage (KP15/KP27) + depo" = "dashed",
                  "phage cocktail (KP15/KP27 + KP34)" = "solid")
  

    linesizes = c("no phage" = 1.2,
                  #"External depo" = 1,
                  "depo-equipped phage (KP34)" = 1.2,
                  "capsule-independent phage (KP15/KP27)" = 1.2,
                  "capsule-independent phage (KP15/KP27) + depo" = 1.2,
                  "phage cocktail (KP15/KP27 + KP34)" = 1.2)
  
  return(list(my_theme = my_theme, 
              colors = colors, 
              linetypes = linetypes, 
              linesizes = linesizes))
}




PlotSimulatedPhageAndBacteria = function( simulated_data,
                                          title_plot = "", 
                                          colors = NULL,
                                          ymin=-1*10^8,
                                          ymax = 2.5*10^9,
                                          tmax = 24,
                                          minCFU = 6.3*10^6,
                                          text_size = 12,
                                          linetypes = NULL,
                                          linesizes = NULL) {
  
  VisualisationConstants = GetVisualisationConstants(text_size)
  my_theme = VisualisationConstants$my_theme
  if (is.null(colors)) {colors = VisualisationConstants$colors}
  if (is.null(linetypes)) {linetypes = VisualisationConstants$linetypes}
  if (is.null(linesizes)) {linesizes = VisualisationConstants$linesizes}
  descriptions_to_show = names(colors)
  
  data_polished = simulated_data %>% 
    dplyr::filter(Treatment %in% descriptions_to_show) %>%
    mutate(all_bacteria_CFU_per_mL = ifelse(all_bacteria_CFU_per_mL > minCFU, all_bacteria_CFU_per_mL, minCFU))
  
  g1=ggplot(data_polished,
            aes(x = time, 
                y = all_bacteria_CFU_per_mL, 
                col = Treatment, 
                linetype = Treatment, 
                size = Treatment)) +
    geom_line(alpha = 0.5) +
    ggtitle(title_plot) +
    xlab("time [h]") +
    ylab("bacteria [CFU/mL]") + 
    scale_color_manual(values=colors) +
    scale_linetype_manual(values = linetypes) +
    scale_size_manual(values = linesizes) +
    scale_y_log10(limits = c(ymin, ymax)) +
    my_theme +
    xlim(c(0,tmax))
  
  return(g1)
  }
  



PlotSimulatedPhageAndBacteriaMultipleScenarios = function( simulated_data,
                                          title_plot = "", 
                                          colors = NULL,
                                          ymin=-1*10^8,
                                          ymax = 2.5*10^9,
                                          tmax = 24,
                                          minCFU = 6.3*10^6,
                                          text_size = 12,
                                          linetypes = NULL,
                                          linesizes = NULL,
                                          ncol = 1,
                                          legend.position="right",
                                          strip.background.color = "gray") {
   n = n_distinct(simulated_data$scenario)
   nrow = ceiling(n/ncol)
   g1 = PlotSimulatedPhageAndBacteria(simulated_data = simulated_data,
                                       title_plot = title_plot,
                                      colors = colors,
                                      ymin = ymin, 
                                      ymax = ymax, 
                                      tmax = tmax,
                                      minCFU = minCFU,
                                      text_size = text_size,
                                      linetypes = linetypes,
                                      linesizes = linesizes)
   g2 = g1 + 
     facet_wrap('scenario', nrow = nrow, ncol = ncol) +
     theme(legend.position=legend.position,
           strip.background = element_rect(fill=strip.background.color))
   return(g2)
  }

# for sensitivity parameter plots
Max.Bacreria.In.Time.Window = function(data, max.time = 24) {
  data.till.max.time = data %>% filter(time <= max.time)
  max.biomass.till.max.time = max(data.till.max.time$all_bacteria_CFU_per_mL)
  return(max.biomass.till.max.time)
}

First.Time.When.Bacteria.Increase = function(data, increase.by = 10) {
  initial.bacteria =  data %>% filter(time == 0) %>% pull(all_bacteria_CFU_per_mL)
  bactera.increased = data %>% filter(all_bacteria_CFU_per_mL >= increase.by*initial.bacteria)
  first.time.bacteria.increased = min(bactera.increased$time)
  return(first.time.bacteria.increased)
}
linspace = function(a,b,n) {
  seq(a, b, by = (b-a)/(n-1))
}
Make.break.labels = function(breaks) {
  labels = paste0(breaks, c(" (low)", " (default)", " (high)"))
}

Make.Sensitivity.Data.To.Adsorption = function(PHI_NON_DEPO_RANGE, PHI_DEPO_RANGE, model.params, minCFU){
  simulated.data.list = list()
  phage.cocktail.success.data = expand.grid(phi.n = PHI_NON_DEPO_RANGE, 
                                            phi.d = PHI_DEPO_RANGE) %>%
    mutate(max.bacteria = NA, 
           min.time = NA)
  
  i=0
  for (PHI_N in PHI_NON_DEPO_RANGE) {
    for (PHI in PHI_DEPO_RANGE) {
      i = i+1
      simulated_data = Simulate_Phage_Coctail(
        model = model.params$model, 
        Vh=model.params$Vh, 
        Kh=model.params$Kh,
        a=model.params$a, 
        beta_non_depo=model.params$beta_non_depo, 
        beta_depo=model.params$beta_depo, 
        phi_non_depo=PHI_N, 
        phi_depo=PHI, 
        decay=model.params$decay, 
        V_depo=model.params$V_depo, 
        K_depo=model.params$K_depo, 
        depo_decay_rate=model.params$depo_decay_rate,
        epsilonB2toB1=model.params$epsilonB2toB1, 
        epsilonB1toB2=model.params$epsilonB1toB2,
        epsilonB2toB3=model.params$epsilonB2toB3, 
        epsilonB3toB2=model.params$epsilonB3toB2,
        epsilonB1toB3=model.params$epsilonB1toB3, 
        epsilonB3toB1=model.params$epsilonB3toB1,  
        B0=model.params$B0, 
        e0=model.params$e0, 
        MOI =model.params$MOI, 
        G0=model.params$G0, 
        Tmax = model.params$Tmax,
        Treatments.included = c(0,1,0,0,0,0)) %>%
        select(time, all_bacteria_CFU_per_mL, Treatment) %>%
        mutate(scenario = "varying_phi",
               phi_depo = PHI,
               phi_non_depo = PHI_N) %>% 
        filter(Treatment == "phage cocktail (KP15/KP27 + KP34)")
      simulated.data.list[[i]] = simulated_data
      
      
      index = which(phage.cocktail.success.data$phi.n == PHI_N & phage.cocktail.success.data$phi.d == PHI)
      max.bact = simulated_data %>% Max.Bacreria.In.Time.Window()
      phage.cocktail.success.data$max.bacteria[index] = max.bact
      
      min.t = simulated_data %>% First.Time.When.Bacteria.Increase()
      phage.cocktail.success.data$min.time[index] = min.t
    }
  }
  
  #simulated.data.all = do.call(rbind, simulated.data.list)
  data = phage.cocktail.success.data %>%
    mutate(min.time = if_else(min.time == Inf, 24, min.time),
           max.bacteria = if_else(max.bacteria >= minCFU,max.bacteria, minCFU),
           phi.n = factor(phi.n, levels = PHI_NON_DEPO_RANGE),
           phi.d = factor(phi.d, levels = PHI_DEPO_RANGE)
    ) 
  return(data)
}

Make.Sensitivity.Plot.To.Adsorption.Max.Bacteria = function(data, BREAKS_NON_DEPO,BREAKS_DEPO,minCFU) {
  ggplot(data) +
    geom_tile(aes(x = phi.n, y = phi.d, fill = max.bacteria)) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                         midpoint = log10(minCFU),
                         name  = "Max bacterial load\nwithin 24 hours",
                         trans = "log10") + 
    scale_x_discrete(breaks = factor(BREAKS_NON_DEPO),
                     labels = Make.break.labels(BREAKS_NON_DEPO))+
    scale_y_discrete(breaks = factor(BREAKS_DEPO),
                     labels = Make.break.labels(BREAKS_DEPO))+
    xlab(latex2exp::TeX("Adsorption rate  $\\phi_N \\,\\, \\lbrack L/(gram \\, bacteria * h) \\rbrack $")) +
    ylab(latex2exp::TeX("Adsorption rate $\\phi_P  \\,\\, \\lbrack L/(gram \\, bacteria * h) \\rbrack $"))
}

Make.Sensitivity.Plot.To.Adsorption.Min.Time = function(data, BREAKS_NON_DEPO,BREAKS_DEPO) {
  ggplot(data ) +
    geom_tile(aes(x = phi.n, y = phi.d, fill = min.time)) +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "darkred", 
                         midpoint = 24,
                         name  = "Time to 10 fold \nbacterial growth") + 
    scale_x_discrete(breaks = factor(BREAKS_NON_DEPO),
                     labels = Make.break.labels(BREAKS_NON_DEPO))+
    scale_y_discrete(breaks = factor(BREAKS_DEPO),
                     labels = Make.break.labels(BREAKS_DEPO))+
    xlab(latex2exp::TeX("Adsorption rate  $\\phi_N \\,\\, \\lbrack L/(gram \\, bacteria * h) \\rbrack $")) +
    ylab(latex2exp::TeX("Adsorption rate $\\phi_P  \\,\\, \\lbrack L/(gram \\, bacteria * h) \\rbrack $"))
}


Make.Sensitivity.Data.To.Decay = function(DECAY_RATE_RANGE, EPSILON_NR_RANGE,model.params,minCFU) {
  simulated.data.list = list()
  phage.cocktail.success.data2 = expand.grid(decay.rate = DECAY_RATE_RANGE, 
                                             epsilon.to.robust = EPSILON_NR_RANGE) %>%
    mutate(max.bacteria = NA, 
           min.time = NA)
  
  i=0
  for (decay.rate in DECAY_RATE_RANGE) {
    for (epsilon.to.robust in EPSILON_NR_RANGE) {
      i = i+1
      simulated_data = Simulate_Phage_Coctail(
        model = model.params$model, 
        Vh=model.params$Vh, 
        Kh=model.params$Kh,
        a=model.params$a, 
        beta_non_depo=model.params$beta_non_depo, 
        beta_depo=model.params$beta_depo, 
        phi_non_depo=model.params$phi_non_depo, 
        phi_depo=model.params$phi_depo, 
        decay=decay.rate, 
        V_depo=model.params$V_depo, 
        K_depo=model.params$K_depo, 
        depo_decay_rate=model.params$depo_decay_rate,
        epsilonB2toB1=model.params$epsilonB2toB1, 
        epsilonB1toB2=model.params$epsilonB1toB2,
        epsilonB2toB3=epsilon.to.robust, 
        epsilonB3toB2=epsilon.to.robust,
        epsilonB1toB3=model.params$epsilonB1toB3, 
        epsilonB3toB1=model.params$epsilonB3toB1,  
        B0=model.params$B0, 
        e0=model.params$e0, 
        MOI =model.params$MOI, 
        G0=model.params$G0, 
        Tmax = model.params$Tmax,
        Treatments.included = c(0,1,0,0,0,0)) %>%
        select(time, all_bacteria_CFU_per_mL, Treatment) %>%
        mutate(scenario = "varying_phi",
               epsilon.to.robust = epsilon.to.robust,
               decay.rate = decay.rate) %>% 
        filter(Treatment == "phage cocktail (KP15/KP27 + KP34)")
      simulated.data.list[[i]] = simulated_data
      
      
      index = which(phage.cocktail.success.data2$epsilon.to.robust == epsilon.to.robust & phage.cocktail.success.data2$decay.rate == decay.rate)
      max.bact = simulated_data %>% Max.Bacreria.In.Time.Window()
      phage.cocktail.success.data2$max.bacteria[index] = max.bact
      
      min.t = simulated_data %>% First.Time.When.Bacteria.Increase()
      phage.cocktail.success.data2$min.time[index] = min.t
    }
  }
  
  #simulated.data.all = do.call(rbind, simulated.data.list)
  data = phage.cocktail.success.data2 %>%
    mutate(min.time = if_else(min.time == Inf, 24, min.time),
           max.bacteria = if_else(max.bacteria >= minCFU,max.bacteria, minCFU),
           decay.rate = factor(decay.rate, levels = DECAY_RATE_RANGE),
           epsilon.to.robust = factor(epsilon.to.robust, levels = EPSILON_NR_RANGE)
    ) 
  return(data)
}

Make.Sensitivity.Data.To.Epsilon = function(EPSILON_PN_RANGE, EPSILON_NR_RANGE,model.params,minCFU) {
  simulated.data.list = list()
  phage.cocktail.success.data = expand.grid(epsilon.to.mutant = EPSILON_PN_RANGE, 
                                             epsilon.to.robust = EPSILON_NR_RANGE) %>%
    mutate(max.bacteria = NA, 
           min.time = NA)

  i=0
  for (epsilon.to.mutant in EPSILON_PN_RANGE) {
    for (epsilon.to.robust in EPSILON_NR_RANGE) {
      i = i+1
      simulated_data = Simulate_Phage_Coctail(
        model = model.params$model, 
        Vh=model.params$Vh, 
        Kh=model.params$Kh,
        a=model.params$a, 
        beta_non_depo=model.params$beta_non_depo, 
        beta_depo=model.params$beta_depo, 
        phi_non_depo=model.params$phi_non_depo, 
        phi_depo=model.params$phi_depo, 
        decay=model.params$decay, 
        V_depo=model.params$V_depo, 
        K_depo=model.params$K_depo, 
        depo_decay_rate=model.params$depo_decay_rate,
        epsilonB2toB1=epsilon.to.mutant, 
        epsilonB1toB2=epsilon.to.mutant,
        epsilonB2toB3=epsilon.to.robust, 
        epsilonB3toB2=epsilon.to.robust,
        epsilonB1toB3=model.params$epsilonB1toB3, 
        epsilonB3toB1=model.params$epsilonB3toB1,  
        B0=model.params$B0, 
        e0=model.params$e0, 
        MOI =model.params$MOI, 
        G0=model.params$G0, 
        Tmax = model.params$Tmax,
        Treatments.included = c(0,1,0,0,0,0)) %>%
        select(time, all_bacteria_CFU_per_mL, Treatment) %>%
        mutate(scenario = "varying_phi",
               epsilon.to.robust = epsilon.to.robust,
               epsilon.to.mutant = epsilon.to.mutant) %>% 
        filter(Treatment == "phage cocktail (KP15/KP27 + KP34)")
      simulated.data.list[[i]] = simulated_data
      
      
      index = which(phage.cocktail.success.data$epsilon.to.robust == epsilon.to.robust & phage.cocktail.success.data$epsilon.to.mutant == epsilon.to.mutant)
      max.bact = simulated_data %>% Max.Bacreria.In.Time.Window()
      phage.cocktail.success.data$max.bacteria[index] = max.bact
      
      min.t = simulated_data %>% First.Time.When.Bacteria.Increase()
      phage.cocktail.success.data$min.time[index] = min.t
    }
  }
  
  #simulated.data.all = do.call(rbind, simulated.data.list)
  data = phage.cocktail.success.data %>%
    mutate(min.time = if_else(min.time == Inf, 24, min.time),
           max.bacteria = if_else(max.bacteria >= minCFU,max.bacteria, minCFU),
           epsilon.to.mutant = factor(epsilon.to.mutant, levels = EPSILON_PN_RANGE),
           epsilon.to.robust = factor(epsilon.to.robust, levels = EPSILON_NR_RANGE)
    ) 
  return(data)
}

Make.Sensitivity.Plot = function(data,
                                 BREAKS_X,
                                   BREAKS_Y,
                                   legend.name,
                                   xlab.name,
                                   ylab.name,
                                   fill.variable,
                                   x.variable,
                                   y.variable,
                                   midpoint,
                                   trans = 'identity') {
  data$x = data %>% pull(x.variable)
  data$y = data %>% pull(y.variable)
  data$fill = data %>% pull(fill.variable)
  
  ggplot(data) +
    geom_tile(aes(x = x, y = y, fill = fill)) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                         midpoint = midpoint,
                         name  = legend.name,
                         trans = trans) + 
    scale_x_discrete(breaks = factor(BREAKS_X), 
                     labels = Make.break.labels(BREAKS_X))+
    scale_y_discrete(breaks = factor(BREAKS_Y),
                     labels = Make.break.labels(BREAKS_Y)) +
    xlab(xlab.name) +
    ylab(ylab.name)
}




Make.Sensitivity.Plot.To.Decay.Min.Time = function(data, BREAKS_DECAY,BREAKS_EPSILON) {
  Make.Sensitivity.Plot(data,
                        BREAKS_X = BREAKS_DECAY,
                        BREAKS_Y = BREAKS_EPSILON,
                        legend.name = "Time to 10 fold \nbacterial growth",
                        xlab.name = latex2exp::TeX("Decay rate $\\gamma \\, \\, \\lbrack 1/h \\rbrack $"),
                        ylab.name = latex2exp::TeX("rate of change $\\epsilon_N^R \\, \\, \\lbrack 1/h \\rbrack $"),
                        fill.variable = "min.time",
                        x.variable = "decay.rate",
                        y.variable = "epsilon.to.robust",
                        midpoint = 24)
}



Make.Sensitivity.Plot.To.Decay.Max.Bacteria = function(data, BREAKS_DECAY,BREAKS_EPSILON, minCFU) {
  Make.Sensitivity.Plot(data,
                        BREAKS_X = BREAKS_DECAY,
                          BREAKS_Y = BREAKS_EPSILON,
                          legend.name = "Max bacterial load\nwithin 24 hours",
                          xlab.name = latex2exp::TeX("Decay rate $\\gamma \\, \\, \\lbrack 1/h \\rbrack $"),
                          ylab.name = latex2exp::TeX("rate of change $\\epsilon_N^R \\, \\, \\lbrack 1/h \\rbrack $"),
                          fill.variable = "max.bacteria",
                          x.variable = "decay.rate",
                          y.variable = "epsilon.to.robust",
                          midpoint = log10(minCFU),
                          trans = "log10" )
}

Make.Sensitivity.Plot.To.Epsilon.Max.Bacteria = function(data, BREAKS_EPSILON_PN,BREAKS_EPSILON, minCFU) {
  Make.Sensitivity.Plot(data,
                        BREAKS_X = BREAKS_EPSILON_PN,
                          BREAKS_Y = BREAKS_EPSILON,
                          legend.name = "Max bacterial load\nwithin 24 hours",
                          xlab.name = latex2exp::TeX("rate of change $\\epsilon_P^N \\, \\, \\lbrack 1/h \\rbrack $"),
                          ylab.name = latex2exp::TeX("rate of change $\\epsilon_N^R \\, \\, \\lbrack 1/h \\rbrack $"),
                          fill.variable = "max.bacteria",
                          x.variable = "epsilon.to.mutant",
                          y.variable = "epsilon.to.robust",
                          midpoint = log10(minCFU),
                         trans = "log10")
}

Make.Sensitivity.Plot.To.Epsilon.Min.Time = function(data, BREAKS_EPSILON_PN,BREAKS_EPSILON) {
  Make.Sensitivity.Plot(data,
                        BREAKS_X = BREAKS_EPSILON_PN,
                          BREAKS_Y = BREAKS_EPSILON,
                          legend.name = "Time to 10 fold \nbacterial growth",
                          xlab.name = latex2exp::TeX("rate of change $\\epsilon_P^N \\, \\, \\lbrack 1/h \\rbrack $"),
                          ylab.name = latex2exp::TeX("rate of change $\\epsilon_N^R \\, \\, \\lbrack 1/h \\rbrack $"),
                          fill.variable = "min.time",
                          x.variable = "epsilon.to.mutant",
                          y.variable = "epsilon.to.robust",
                          midpoint = 24)
}
