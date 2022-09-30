suppressWarnings(suppressPackageStartupMessages(library(deSolve)))
suppressWarnings(suppressPackageStartupMessages(library(DEoptim)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))


###################### helper functions ##############################


bacteria_growth_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G) 
    Gdot   = -Jg*N 
    Ndot   = a*Jg*N
    return(list(c(Gdot, Ndot)))
  })
}




simulate_growth_single_strain = function(a,Vh,Kh,H0,N0, time) {
  pars <- c(a = a, Vh = Vh, Kh = Kh)
  inits = c(G=H0, N=N0)
  simulation <- as.data.frame(ode(inits, time, bacteria_growth_model, pars))
  return(simulation$N)
}


sumLeastSquaresFitGrowthToDeoptim = function(param, data,a, H0,N0) {
  Vh = param[1]
  Kh = param[2]
  simulatedN = simulate_growth_single_strain(a,Vh,Kh,H0,N0, data$time)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}






Fit.Bacterial.Growth.Parameters.To.Data = function(N0,H0, max_time, data.path) {
  deopticontrol = DEoptim.control(itermax = 200, reltol = 10^(-8), trace = 100)
  data = Get.Data(data.path, N0)

  # find initial biomass
  #N0 = data$biomass[1]
  
  # find final biomass
  Nend = data$biomass[nrow(data)]
  N1=data$biomass[1]
  # find the parameter a i.e. how many grams of proteins can be created per mmol of glucose
  a = mean((Nend-N1)/(H0)) # 0.03 * 10^9 cells per mmol of glucose
  data =    data %>%
    filter(time <= max_time)

  
  # fit params to fresh cells using Deoptim
  #Parameters data, a,N0 and H0 are taken from the global environment
  deoptim_out = DEoptim(fn = sumLeastSquaresFitGrowthToDeoptim, 
                        lower = c(0,0), 
                        upper=c(5000,500),
                        control = deopticontrol,
                        data,
                        a,
                        H0,
                        N1)
  Vh=deoptim_out$optim$bestmem[1] %>% as.numeric()
  Kh=deoptim_out$optim$bestmem[2] %>% as.numeric()
  
  
  
  # Use MLE
  
  #free = c(Vh=50, Kh=50, sigma=1)
  #fixed = c(a=a, mu = 0,  N0=N0, H0=H0)
  
  # This is a function used for MLE fitting
  #nll = function(a, Vh, Kh, mu, sigma, N0, H0) {
  #  simulatedN = simulate_growth_single_strain(a, Vh, Kh,H0, N0, data$time)
  #  difference = data$biomass - simulatedN
  #  logf = suppressWarnings(dnorm(difference[-1],mu,sigma,log = TRUE))
  #  ll = sum(logf)
  #  return(-ll)
  #}
  
  #fit <- mle2(nll, start = as.list(free), fixed = as.list(fixed), method = "Nelder-Mead")
  return(list(a=a, Vh = Vh, Kh = Kh))
  }
  
Get.Data = function(data.path, N0 = NULL) {
  # we will work on the scale of mg of bacteria 
  # https://bionumbers.hms.harvard.edu/bionumber.aspx?&id=103905&ver=7
  # assuming that 1 bacterial cells weighs ~ 0.95*10^(-12) g;
  becterial_mass = 10^(-12);
  # in other words 1mg is equal to that many cells
  mg_count = (10^9) 
  
  data = read.csv(data.path,
                  header = FALSE,
                  dec=".", 
                  stringsAsFactors=FALSE)
  
  names(data) = c("time", "CFUperML")
  if (is.na(as.numeric(data[1,1]))) {
    data = data[-1,]
  }
  data = data %>%
    mutate(time = as.numeric(time),
           CFUperML = as.numeric(CFUperML)) %>%
    filter(CFUperML > 0) %>%
    # biomass is now in g/L
    mutate(biomass = CFUperML/mg_count) %>%
    select(time, biomass)
  
  # if the user knows what should be the first read correct for that in the data
  if (!is.null(N0)) {
    data$biomass[1] = N0/mg_count
  }
  return(data)
} 

Plot.Fitted.Bacterial.Growth.Parameters.To.Data = function(N0, H0,params, data.path) {
  data = Get.Data(data.path, N0)
  N0 = data$biomass[1]
  # Here data is taken from the global environment
  data$predicted = simulate_growth_single_strain(params$a,params$Vh,params$Kh,H0,N0,data$time)
  data_mean = data %>%
    rename(real= biomass) %>%
    tidyr::gather(key="biomass.data.type", value="biomass", real, predicted)
  p1=ggplot(data_mean) + 
    geom_point(aes(time,biomass, color = biomass.data.type)) +
    geom_line(aes(time,biomass, color = biomass.data.type)) +
    theme_bw() +
    ylab("Biomass mg/mL")
  return(p1)
}


##################### example usage ###################
#N0 = 10^6
# molar mass of glucose: 180g / mol = 0.180 g per mmol
# 0.25% glucose i.e. 2.5g per 1L = 2.5[g/L] * 1/0.18 [mmol/g] = 2.5/0.18 [mmol/L]
#H0 = 13.9 #mmol / L
#max_time = 24 # By max_time we had around 40% of bacterial biomass! so there must be some glcuose left after that time. 

# first column: time,
# second column: mean CFU per ml
#data.path = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Z_Drulis_Kawa/phage_cocktail_efficiency/shiny_app/PhREEPred/R/example_data.csv"


#params =Fit.Bacterial.Growth.Parameters.To.Data(N0,H0, max_time, data.path) 
#Plot.Fitted.Bacterial.Growth.Parameters.To.Data(N0,H0,params, data.path)

