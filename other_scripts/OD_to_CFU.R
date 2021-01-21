library(readxl)
library(dplyr)
library(ggplot2)
OD_calibration_data = readxl::read_excel(path = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Z_Drulis_Kawa/data/2020_02/Dane dla KP77.xlsx",
                   sheet = "OD vs CFU_ml",
                   range = "B5:C18",
                   col_names = TRUE)


data = OD_calibration_data %>%
  mutate(expOD = exp(`OD średnia`),
         OD2 = `OD średnia`^2,
         logCFU = log(`CFU/ml`))


ggplot(data = data, aes(x = expOD, y = `CFU/ml`)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x)


mod = lm(`CFU/ml` ~ expOD,  data = data)
summary(mod)

CFUFromOD = function(OD) {
  #CFU = as.numeric(mod$coefficients[1]) + as.numeric(mod$coefficients[2])*exp(OD)
  CFU = -1079035517 + 1014525658*exp(OD)
  return(CFU)
}

data_predicted = data %>%
  mutate(CFU.predicted = CFUFromOD(`OD średnia`)) %>%
  rename(`CFU.real` = `CFU/ml`) %>%
  tidyr::gather(key = "data.type", value = "CFU", `CFU.real`, CFU.predicted)

# recalculate to CFU/L
fig = ggplot(data = data_predicted %>% mutate(CFU = CFU), 
       aes(x =`OD średnia` , y = CFU, col = data.type)) + 
  geom_line() +
  scale_y_log10() +
  geom_point() + 
  theme_bw() +
  ylab("CFU/mL") +
  xlab("Mean OD")

experimental_figures_path =  "/Users/bognasmug/MGG Dropbox/Bogna Smug/Z_Drulis_Kawa/manuskrypt/experimental_figures/"
ggsave(filename = paste0(experimental_figures_path, "Fig13_OD_calibration",  ".jpg"), plot = fig, width = 20, height = 10, units = "cm")
print(fig)

