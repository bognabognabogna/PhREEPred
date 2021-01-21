suppressWarnings(suppressPackageStartupMessages(library(deSolve)))
suppressWarnings(suppressPackageStartupMessages(library(DEoptim)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

deopticontrol = DEoptim.control(itermax = 1000, reltol = 10^(-8), trace = 100)

###################### helper functions ##############################

# per ML
CFUFromOD = function(OD) {
  #CFU = as.numeric(mod$coefficients[1]) + as.numeric(mod$coefficients[2])*exp(OD)
  # this will give us CFU/mL
  CFU = -1079035517 + 1014525658*exp(OD)
  # this will give us CFU/L
  #CFU = 1000*CFU
  return(CFU)
}

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

sumLeastSquaresFitGrowthToDeoptim = function(param) {
  Vh = param[1]
  Kh = param[2]
  #a = param[3]
  simulatedN = simulate_growth_single_strain(a,Vh,Kh,H0,N0, data$time)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}

##################### load data and find best params ###################
# we will work on the scale of mg of bacteria 
# https://bionumbers.hms.harvard.edu/bionumber.aspx?&id=103905&ver=7
# assuming that 1 bacterial cells weighs ~ 0.95*10^(-12) g;
becterial_mass = 10^(-12);
# in other words 1mg is equal to that many cells
mg_count = (10^9) 
# cultures were started with 10^6 cells/mL i.e. N0 mg/mL
N0=10^6/mg_count
# molar mass of glucose: 180g / mol = 0.180 g per mmol
# 0.25% glucose i.e. 2.5g per 1L = 2.5[g/L] * 1/0.18 [mmol/g] = 2.5/0.18 [mmol/L]
H0 = 13.9 #mmol / L
max_time = 24 # By max_time we had around 40% of bacterial biomass! so there must be some glcuose left after that time. 



data486depo=readxl::read_excel(path = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Z_Drulis_Kawa/data/2020_06/Dane dla Bogny2.xlsx", 
                               sheet = "KP486 do papieru",
                               range = "AL121:AU170")
data77depo=readxl::read_excel(path = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Z_Drulis_Kawa/data/2020_06/Dane dla Bogny2.xlsx", 
                              sheet = "KP77 do papieru ",
                              range = "BG126:BP175")


names(data486depo)[1] = "time"
names(data77depo)[1] = "time"
data486 = data486depo %>% 
  select(time,Kp486) %>%
  mutate(CFU = CFUFromOD(Kp486))  %>%
  select(time,CFU) %>% mutate(strain = "Kp486")
data77 = data77depo %>% 
  mutate(CFU = CFUFromOD(Kp77)) %>%
  select(time,CFU)  %>% mutate(strain = "Kp77")
dataAll = data486 %>% rbind(data77) %>% tidyr::spread(strain, CFU)




freshCells = dataAll %>%
  select(time, meanBiomass = Kp77) %>% 
  filter(meanBiomass > 0) %>%
  mutate(biomass = meanBiomass/mg_count) %>%
  select(time, biomass)
# find initial biomass
#N0 = freshCells$biomass[1]
# find final biomass
Nend = freshCells$biomass[nrow(freshCells)]
# find the parameter a i.e. how many grams of proteins can be created per mmol of glucose
a = mean((Nend-N0)/(H0)) # 0.03 * 10^9 cells per mmol of glucose
data =    freshCells %>%
  filter(time <= max_time)




# fit params to fresh cells using Deoptim
#Parameters data, a,N0 and H0 are taken from the global environment
deoptim_out = DEoptim(fn = sumLeastSquaresFitGrowthToDeoptim, lower = c(0,0), upper=c(5000,500),
                      control = deopticontrol)
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



# Here data is taken from the global environment
freshCells$predicted = simulate_growth_single_strain(a,Vh,Kh,H0,freshCells$biomass[1], freshCells$time)
data_mean = freshCells %>%
  tidyr::gather(key="selected_type", value="biomass", biomass, predicted)
p1=ggplot(data_mean) + 
  geom_point(aes(time,biomass, color = selected_type)) +
  theme_bw() +
  ylab("Biomass mg/mL")
print(p1)
print(a)
print(Vh)
print(Kh)


# K [mmol / L] not dependent on units f bacteria
# according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC98929/
# K for E-Coli is in the range of 40-100 000 miligrams/L i.e. 0.15-500 mmol/L

# V [mmol / mg biomass * h]
# a [mg biomass / mmol glucose]
