#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)
library(deSolve)
library(latex2exp)

source("R/helpers.R")


source("R/fit_bacterial_growth_parameters.R")
mg_count = (10^12) 
min.cfu = 6.3*10^6
default_params = SetDefaultParameters(mg_count)
min_params = SetMinParameters(mg_count)
max_params = SetMaxParameters(mg_count)


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    withMathJax(),
    # Application title
    titlePanel("PhREEPred: Phage Resistance EmergencE Prediction"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Calculating...",id="loadmessage")),
            h2("Experimental setup:"),
            numericInput("Tmax",
                         label = "T: Experiment time length [h]",
                         min = 0,
                         max = 240,
                         default_params$Tmax),
            numericInput("MOI",
                         label = "MOI: multiplicity of infection",
                         min = 0,
                         max = 100,
                         default_params$MOI),
            #numericInput("B0",
            #             label = "\\(B_0\\): Initial density of bacteria [CFU/mL]",
            #             min = 0,
            #             step = 1000,
            #             default_params$B0*(mg_count/1000)),
            selectInput("B0",
                        label = "\\(B_0\\): Initial density of bacteria [CFU/mL]",
                        c(10, 10^2,10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10),
                        selected = default_params$B0*(mg_count/1000)),
            numericInput("G0",
                         label = "\\(G_0\\): Initial concentration of glucose [mmol/L] (13.9 mmol/L for standard TSB medium)",
                         min = 0,
                         max = 140,
                         default_params$G0),
            br(),
            br(),
            h2("Model Parameters:"),
            h4("Phage virulence:"),
            numericInput("beta_depo",
                         label = "\\( \\beta_P \\): burst size of the phage \\(P_P\\) [phage particles per bacterial cell]",
                         min = 0,
                         max = 10^4,
                         default_params$beta_depo),
            numericInput("beta_non_depo",
                         label = "\\( \\beta_N \\): burst size of the phage \\(P_N\\) [phage particles per bacterial cell]",
                         min = 0,
                         max = 10^4,
                         default_params$beta_non_depo),
            numericInput("latent_period_depo",
                         label = "\\( T^L_P \\): Latent period of \\(P_P\\)  [min]",
                         min = 0,
                         max = 10^4,
                         60*default_params$latent_period_depo),
            numericInput("latent_period_non_depo",
                         label = "\\( T^L_N \\):  Latent period of  \\(P_N\\) [min]",
                         min = 0,
                         max = 10^4,
                         60*default_params$latent_period_non_depo),
            br(),
            br(),
            h4("Bacterial growth:"),
            numericInput("Vh",
                         label = "\\( V^G \\): maximal glucose uptake rate [mmol glucose / g. bacteria x h])",
                         min = 0,
                         default_params$Vh),
            numericInput("Kh",
                         label ="\\( K^G \\): Michaelis-Menten constant for glucose uptake rate [mmol glucose]",
                         min = 0,
                         default_params$Kh),
            numericInput("a",
                         label = " \\( \\eta \\): Efficiency of glucose uptake [g. bacteria/mmol glucose])",
                         min = 0,
                         max = 1,
                         step = 0.01,
                         default_params$a),
            br(),
            br(),
          h2("Optional parameters:"),
          h4("Phage characteristics:"),
          numericInput("phi_depo",
                       label = "\\( \\Phi_P \\): Adsorption rate of \\(P_P\\) to \\(B_P\\) [L/(g. bacteria*h]. This is therate at which a bacterium and phage encounter each other and result in an infection.",
                       min = 0,
                       max = 10^4,
                       default_params$phi_depo),
          numericInput("phi_non_depo",
                       label = "\\( \\Phi_N \\): Adsorption rate of \\(P_N\\) to \\(B_N\\) [L/(g. bacteria*h]. This is the rate at which a bacterium and phage encounter each other and result in an infection.",
                       min = 0,
                       max = 10^4,
                       default_params$phi_non_depo),
          numericInput("decay",
                       label = "\\( \\gamma \\): phage decay rate [1/h]",
                       min = 0,
                       step = 0.1,
                       max = 0.2,
                       default_params$decay),
            #h4("Externally added depolymerase:"),
            #numericInput("e0",
            #             label = "\\(E_0\\): Initial depolymerase concentration [g/L]",
            #             min = 0,
            #             step = 0.1,
            #             default_params$e0),
            #numericInput("V_depo",
            #             label = "\\(V^D\\): Maximal depolymerase activity rate (scaled) [g/L*h]",
            #             min = 0,
            #             default_params$V_depo),
            #numericInput("K_depo",
            #             label = "\\(K^D\\): Michaelis-Menten constant for depolymerase activity [g bacteria / L]",
            #             min = 0,
            #             default_params$K_depo),
            #numericInput("depo_decay_rate",
            #             label = "\\(\\gamma^E\\): Depolymerase decay rate[1/h]",
            #             min = 0,
            #             step = 0.1,
            #             default_params$depo_decay_rate),
            #br(),
            #br(),
            #h4("Phenotypical changes between bacterial strains:"),
            #numericInput("epsilonB1toB2",
            #             label = "\\(\\epsilon_P^N\\): switch rate from \\(B_P\\) to \\(B_N\\) [1/h]",
            #             min = 0,
            #             step = default_params$epsilonB1toB2,
            #             default_params$epsilonB1toB2),
            #numericInput("epsilonB2toB1",
            #             label = "\\(\\epsilon_N^P\\): switch rate from \\(B_N\\) to \\(B_P\\) [1/h]",
            #             min = 0,
            #             step = default_params$epsilonB2toB1,
            #             default_params$epsilonB2toB1),
            #numericInput("epsilonB1toB3",
            #             label = "\\(\\epsilon_P^R\\): switch rate from \\(B_P\\) to \\(B_R\\) [1/h]",
            #             min = 0,
            #             step = default_params$epsilonB1toB3,
            #             default_params$epsilonB1toB3),
            #numericInput("epsilonB3toB1",
            #             label = "\\(\\epsilon_R^P\\): switch rate from \\(B_R\\) to \\(B_P\\) [1/h]",
            #             min = 0,
            #             step = default_params$epsilonB3toB1,
            #             default_params$epsilonB3toB1),
            #numericInput("epsilonB3toB2",
            #             label = "\\(\\epsilon_R^N\\): switch rate from \\(B_R\\) to \\(B_N\\) [1/h]",
            #             min = 0,
            #             step = default_params$epsilonB3toB2,
            #             default_params$epsilonB3toB2),
            #numericInput("epsilonB2toB3",
            #             label = "\\(\\epsilon_N^R\\): switch rate from \\(B_N\\) to \\(B_R\\) [1/h]",
            #             min = 0,
            #             step = default_params$epsilonB2toB3,
            #             default_params$epsilonB2toB3)
        width = 4),
        
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type = "tabs",
                      tabPanel("Bacteria growth dynamics", 
                               h4("Phage therapy effects:"),
                               textOutput("errortext"),
                               plotOutput("timeplot"),
                               br(),
                               br(),
                               h4("Bacterial growth by bacteria type:"),
                               plotOutput("bacteriatypeplot"),
                               #plotOutput("bacteriatypeplot2")
                               ),
                      tabPanel("Phage cocktail: sensitivity to parameters",
                               h4("Efficiency of the phage cocktail treatments depending on optional parameters which may be difficult to measure.
                                  The efficiency is measured by: 
                                  (left) max bacterial load within 24 hours,
                                  (right) min. time at which bacterial load increases 10 times from the initial concentration."),
                               fluidRow(column(12, br(),h4("1) depending on phage adsorption rates."))),
                               fluidRow(
                                 column(6,
                                        plotOutput("sensitivity.to.adsorption.max.bacteria.plot")),
                                 column(6,
                                        plotOutput("sensitivity.to.adsorption.min.time.plot")),
                               ),
                               fluidRow(column(12, br(),h4("2) depending on phage decay rate and switch rate to the completely resistant mutant."))),
                               fluidRow(
                                 column(6,
                                        plotOutput("sensitivity.to.decay.max.bacteria.plot")),
                                 column(6,
                                        plotOutput("sensitivity.to.decay.min.time.plot"))
                                 ),
                               fluidRow(12, br(),h4("3) depending on switch rates to the decpasulated mutant and to the completely resistant mutant.")),
                               fluidRow(
                                 column(6,
                                        plotOutput("sensitivity.to.epsilon.max.bacteria.plot")),
                                 column(6,
                                        plotOutput("sensitivity.to.epsilon.min.time.plot"),
                                 )),
     
                      ),
                      tabPanel("Phage + depo: sensitivity to parameters",
                               h4("Efficiency of the cpasule independent phage and depolymerase combination treatments depending on optional parameters which may be difficult to measure.
                                  The efficiency is measured by: 
                                  (left) max bacterial load within 24 hours,
                                  (right) min. time at which bacterial load increases 10 times from the initial concentration."),
                               fluidRow(column(12, br(),h4("1) depending on strength and decay rate of the externally added depolymerase."))),
                               fluidRow(
                                 column(6,
                                        plotOutput("sensitivity.to.depo.max.bacteria.plot")
                                        ),
                                 column(6,
                                        plotOutput("sensitivity.to.depo.min.time.plot")
                                        ),
                               ),
                      ),
                      tabPanel("Get Bacterial Growth Parameters", 
                               h4("The acterial growth parameters will be calculated upon input of the growth curve data (bacteria in absence of phages).
                                  Please format your data as a .csv file where the first column is time [hours], and the second column is bacterial load [CFU/mL]"),
                               br(),
                               fileInput("growth.curve.file", "Choose CSV File (first column = time, second column = CFU/mL)",
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv"),
                                         placeholder = "No file selected, using example data:",
                                         width = '80%'),
                               numericInput("H0bacteria",
                                            label = "Initial glucose concentration [mmol / L]",
                                            min = 0,
                                            step = 0.1,
                                            max = 0.2,
                                            default_params$G0),
                               numericInput("N0bacteria",
                                            label = "Initial bacterial load [CFU / mL]",
                                            min = 0,
                                            step = 0.1,
                                            max = 0.2,
                                            10^6),
                               numericInput("max_time",
                                            label = "Time when to stop fitting (bacteria in stationary phase) [h]",
                                            min = 0,
                                            step = 0.1,
                                            max = 200,
                                            24),
                               tableOutput("fitted.params"),
                               plotOutput("fitted.params.plot")
                      )
          )

          
        )
    )
))


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
    isValid_input <- reactive({
        #!is.null(input$Vh) && input$Vh>=min_params$Vh && input$Vh<=max_params$Vh &&
         #   !is.null(input$Kh) && input$Kh>=min_params$Kh && input$Kh<=max_params$Kh &&
         #   !is.null(input$a) && input$a>=min_params$a && input$a<=max_params$a &&
         #   !is.null(input$beta_non_depo) && input$beta_non_depo>=min_params$beta_non_depo && input$beta_non_depo<=max_params$beta_non_depo &&
         #   !is.null(input$beta_depo) && input$beta_depo>=min_params$beta_depo && input$beta_depo<=max_params$beta_depo &&
         #   !is.null(input$phi_non_depo) && input$phi_non_depo>=min_params$phi_non_depo && input$phi_non_depo<=max_params$phi_non_depo &&
         #   !is.null(input$phi_depo) && input$phi_depo>=min_params$phi_depo && input$phi_depo<=max_params$phi_depo &&
         #   !is.null(input$decay) && input$decay>=min_params$decay && input$decay<=max_params$decay &&
         #   !is.null(input$V_depo) && input$V_depo>=min_params$V_depo &&  input$V_depo<=max_params$V_depo &&
         #   !is.null(input$K_depo) && input$K_depo>=min_params$K_depo && input$K_depo<=max_params$K_depo &&
         #   !is.null(input$depo_decay_rate) && input$depo_decay_rate>=min_params$depo_decay_rate && input$depo_decay_rate<=max_params$depo_decay_rate &&
            #!is.null(input$epsilonB2toB1) && input$epsilonB2toB1>=0 &&
            #!is.null(input$epsilonB1toB2) && input$epsilonB1toB2>=0 &&
            #!is.null(input$epsilonB2toB3) && input$epsilonB2toB3>=0 &&
            #!is.null(input$epsilonB3toB2) && input$epsilonB3toB2>=0 &&
            #!is.null(input$epsilonB1toB3) && input$epsilonB1toB3>=0 &&
            #!is.null(input$epsilonB3toB1) && input$epsilonB3toB1>=0 &&
            !is.null(input$B0) && as.numeric(input$B0)/(mg_count/1000)>=min_params$B0 &&  as.numeric(input$B0)/(mg_count/1000)<=max_params$B0 &&
            #!is.null(input$e0) && input$e0>=min_params$e0 && input$e0<=max_params$e0 &&
            !is.null(input$MOI) && input$MOI>=min_params$MOI && input$MOI<=max_params$MOI &&
            !is.null(input$G0) && input$G0>=min_params$G0 && input$G0<=max_params$G0 &&
            !is.null(input$Tmax) && input$Tmax>=min_params$Tmax && input$Tmax<=max_params$Tmax })
    
    
    model.params = reactive({
      list(
      Vh=input$Vh, 
      Kh=input$Kh,
      a=input$a, 
      beta_non_depo=input$beta_non_depo, 
      beta_depo=input$beta_depo, 
      phi_non_depo=input$phi_non_depo, 
      phi_depo=input$phi_depo, 
      decay=input$decay, 
      V_depo=default_params$V_depo, 
      K_depo=default_params$K_depo, 
      depo_decay_rate=default_params$depo_decay_rate,
      epsilonB2toB1=default_params$epsilonB2toB1,#input$epsilonB2toB1, 
      epsilonB1toB2=default_params$epsilonB1toB2,#input$epsilonB1toB2,
      epsilonB2toB3=default_params$epsilonB2toB3,#input$epsilonB2toB3, 
      epsilonB3toB2=default_params$epsilonB3toB2,#input$epsilonB3toB2,
      epsilonB1toB3=default_params$epsilonB1toB3,#input$epsilonB1toB3, 
      epsilonB3toB1=default_params$epsilonB3toB1,#input$epsilonB3toB1,  
      B0=as.numeric(input$B0)/(mg_count/1000), # so that it is in g/L
      e0=default_params$e0, 
      MOI =input$MOI, 
      G0=input$G0, 
      Tmax = input$Tmax,
      latent_period_depo = input$latent_period_depo/60,
      latent_period_non_depo = input$latent_period_non_depo/60,
      model = default_params$model,
      bf = default_params$bf,
      propB2init = default_params$propB2init)
    })
    
    output$errortext <- renderText({
        if(isValid_input()){ invisible(NULL) 
        }else{ #
            "Parameter values are outside a reasonable range.
             Note: They must be numeric and non-negative"
        }
    })
    
    fitted.params = reactive({
      if (is.null(input$growth.curve.file)) {
          data.path = "R/simulated_data.csv"
        } else {
          data.path =  input$growth.curve.file$datapath
        }
        fitted.params = Fit.Bacterial.Growth.Parameters.To.Data(input$N0bacteria,input$H0bacteria, input$max_time, data.path) 
    })
    
    output$fitted.params <- renderTable({
      #req(input$growth.curve.file)
      params = fitted.params()
      #growth.curve.data = Get.Data(, N0 = NULL)
      params.to.display = data.frame(
        #parameter = c(HTML("$$\\mbox{Maximal glucose uptake rate: }V^G$$"),
        #              HTML("$$\\mbox{Michaelis-Menten constant for glucose uptake rate: }K^G$$"),
        #              HTML("$$\\mbox{Efficiency of glucose uptake: }\\eta$$")),
        parameter = c(HTML("Maximal glucose uptake rate: V^G"),
                      HTML("Michaelis-Menten constant for glucose uptake rate: K^G"),
                      HTML("Efficiency of glucose uptake: eta")),
        value = c(params$Vh,params$Kh,params$a))
      return(params.to.display)
    }, sanitize.text.function = function(x) x)
    
    
    
    output$fitted.params.plot = renderPlot({
      #req(input$growth.curve.file)
      params = fitted.params()
      if (is.null(input$growth.curve.file)) {
        data.path = "R/simulated_data.csv"
      } else {
        data.path =  input$growth.curve.file$datapath
      }
        fig = Plot.Fitted.Bacterial.Growth.Parameters.To.Data(input$N0bacteria,input$H0bacteria,params, data.path)
      return(fig)
    })
    

    simulated_data =  reactive({Simulate_Phage_Coctail(model.parameters = model.params()) %>% 
        mutate(Treatment = plyr::revalue(Treatment,
                                         c( "depo-equipped phage (KP34)" =  "depo-equipped phage", 
                                            "capsule-independent phage (KP15/KP27)"  = "capsule-independent phage" ,
                                            "capsule-independent phage (KP15/KP27) + depo" = "capsule-independent phage + depo",
                                            "phage cocktail (KP15/KP27 + KP34)" = "phage cocktail" ))) %>%
            filter(Treatment != "External depo" & Treatment != "capsule-independent phage + depo")

    })
    
    num.in.range = 21
    PHI_NON_DEPO_RANGE = reactive({model.params()$phi_non_depo*c(10^linspace(-1,1, n = num.in.range))})
    PHI_DEPO_RANGE = reactive({model.params()$phi_depo*c(10^linspace(-1,1, n =num.in.range))})
    DECAY_RATE_RANGE = reactive({model.params()$decay*c(10^linspace(-1,1, n = num.in.range))})
    EPSILON_NR_RANGE = reactive({model.params()$epsilonB2toB3*c(10^linspace(-10,10, n = num.in.range))})
    EPSILON_PN_RANGE = reactive({model.params()$epsilonB1toB2*c(10^linspace(-3,3, n = num.in.range))})
    DEPO_DECAY_RANGE = reactive({model.params()$depo_decay_rate*c(10^linspace(-3,3, n = num.in.range))})
    DEPO_STRENGTH_RANGE = reactive({model.params()$V_depo*c(10^linspace(-3,3, n = num.in.range))})
    
    BREAKS_NON_DEPO =  reactive({c(PHI_NON_DEPO_RANGE()[1], model.params()$phi_non_depo, PHI_NON_DEPO_RANGE()[length(PHI_NON_DEPO_RANGE())])})
    BREAKS_DEPO =  reactive({c(PHI_DEPO_RANGE()[1], model.params()$phi_depo, PHI_DEPO_RANGE()[length(PHI_DEPO_RANGE())])})
    BREAKS_DECAY =  reactive({c(DECAY_RATE_RANGE()[1], model.params()$decay, DECAY_RATE_RANGE()[length(DECAY_RATE_RANGE())])})
    BREAKS_EPSILON = reactive({c(EPSILON_NR_RANGE()[1], model.params()$epsilonB2toB3, EPSILON_NR_RANGE()[length(EPSILON_NR_RANGE())])})
    BREAKS_EPSILON_PN = reactive({c(EPSILON_PN_RANGE()[1], model.params()$epsilonB1toB2, EPSILON_PN_RANGE()[length(EPSILON_PN_RANGE())])})
    BREAKS_DEPO_DECAY = reactive({c(DEPO_DECAY_RANGE()[1], model.params()$depo_decay_rate, DEPO_DECAY_RANGE()[length(DEPO_DECAY_RANGE())])})
    BREAKS_DEPO_STRENGTH = reactive({c(DEPO_STRENGTH_RANGE()[1], model.params()$V_depo, DEPO_STRENGTH_RANGE()[length(DEPO_STRENGTH_RANGE())])})
    
    sensitivity.data.to.adsorption =  reactive({
      sensitivity.Data.To.Adsorption = Make.Sensitivity.Data.To.Adsorption(PHI_NON_DEPO_RANGE=PHI_NON_DEPO_RANGE(), 
                                                                           PHI_DEPO_RANGE=PHI_DEPO_RANGE(),
                                                                           model.params=model.params(),
                                                                           minCFU = min.cfu)
      return(sensitivity.Data.To.Adsorption)
    })
    

    sensitivity.data.to.decay = reactive({
      sensitivity.Data.To.Decay = Make.Sensitivity.Data.To.Decay(DECAY_RATE_RANGE=DECAY_RATE_RANGE(), 
                                                                 EPSILON_NR_RANGE=EPSILON_NR_RANGE(), 
                                                                 model.params = model.params(),
                                                                 minCFU = min.cfu)
      return(sensitivity.Data.To.Decay)
    })
    
    sensitivity.data.to.epsilon = reactive({
      sensitivity.Data.To.epsilon = Make.Sensitivity.Data.To.Epsilon(EPSILON_PN_RANGE=EPSILON_PN_RANGE(), 
                                                                 EPSILON_NR_RANGE=EPSILON_NR_RANGE(), 
                                                                 model.params = model.params(),
                                                                 minCFU = min.cfu)
      return(sensitivity.Data.To.epsilon)
    })
    
    sensitivity.data.to.depo = reactive({
      sensitivity.Data.To.depo = Make.Sensitivity.Data.To.Depo(DEPO_DECAY_RANGE = DEPO_DECAY_RANGE(), 
                                                               DEPO_STRENGTH_RANGE = DEPO_STRENGTH_RANGE(), 
                                                               model.params = model.params(),
                                                               minCFU = min.cfu)
      return(sensitivity.Data.To.depo)
    })
    
    output$sensitivity.to.adsorption.max.bacteria.plot <- renderPlot({
      if(isValid_input()){
        fig = Make.Sensitivity.Plot.To.Adsorption.Max.Bacteria(sensitivity.data.to.adsorption(), BREAKS_NON_DEPO(),BREAKS_DEPO(), minCFU=min.cfu)
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    output$sensitivity.to.adsorption.min.time.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Adsorption.Min.Time(sensitivity.data.to.adsorption(),BREAKS_NON_DEPO(),BREAKS_DEPO())
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })

    
    output$sensitivity.to.decay.max.bacteria.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Decay.Max.Bacteria(sensitivity.data.to.decay(), BREAKS_DECAY(), BREAKS_EPSILON(), minCFU=min.cfu)
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    output$sensitivity.to.decay.min.time.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Decay.Min.Time(sensitivity.data.to.decay(), BREAKS_DECAY(), BREAKS_EPSILON())
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    output$sensitivity.to.epsilon.max.bacteria.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Epsilon.Max.Bacteria(data=sensitivity.data.to.epsilon(), BREAKS_EPSILON_PN(), BREAKS_EPSILON(), minCFU=min.cfu)
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    output$sensitivity.to.epsilon.min.time.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Epsilon.Min.Time(data=sensitivity.data.to.epsilon(), BREAKS_EPSILON_PN(), BREAKS_EPSILON())
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    
    
    output$sensitivity.to.depo.max.bacteria.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Depo.Max.Bacteria(data=sensitivity.data.to.depo(),BREAKS_DEPO_DECAY(), BREAKS_DEPO_STRENGTH(), minCFU=min.cfu)
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    output$sensitivity.to.depo.min.time.plot <- renderPlot({
      if(isValid_input()){ 
        fig = Make.Sensitivity.Plot.To.Depo.Min.Time(data=sensitivity.data.to.depo(), BREAKS_DEPO_DECAY(), BREAKS_DEPO_STRENGTH())
        return(fig)
      }else{ 
        return(invisible(NULL))
      }
    })
    
    remove_x <- function(vec, x) {
      vec[setdiff(names(vec), x)]
    }
    
    output$timeplot <- renderPlot({
        if(isValid_input()){ 
            
            VisualisationConstants = GetVisualisationConstantsApp()
            my_theme = VisualisationConstants$my_theme
            colors = VisualisationConstants$colors %>% remove_x("capsule-independent phage + depo")
            linetypes = VisualisationConstants$linetypes %>% remove_x("capsule-independent phage + depo")
            linesizes = VisualisationConstants$linesizes %>% remove_x("capsule-independent phage + depo")
            
        
        fig =PlotSimulatedPhageAndBacteria(simulated_data(),
                                           title_plot = "", 
                                           colors = colors,
                                           ymin = NULL,
                                           ymax = NULL,
                                           tmax = input$Tmax,
                                           minCFU =min.cfu,
                                           text_size = 12,
                                           linetypes = linetypes,
                                           linesizes = linesizes) 
        
        return(fig)
        }else{ #
           return(invisible(NULL))
        }
        
        
    })
  
    output$bacteriatypeplot <- renderPlot({
        if(isValid_input()){ 
            fig =PlotResultsByBacteriumType(simulated_data(),
                                               title_plot = "", 
                                               ymin = NULL,
                                               ymax = NULL,
                                               tmax = input$Tmax,
                                               minCFU = min.cfu,
                                               text_size = 12) 
            
            return(fig)
        }else{ #
            return(invisible(NULL))
        }
        
    })
    
    output$bacteriatypeplot2 <- renderPlot({
        if(isValid_input()){ 
            VisualisationConstants = GetVisualisationConstantsApp()
            my_theme = VisualisationConstants$my_theme
            colors = VisualisationConstants$colors
            linetypes = VisualisationConstants$linetypes
            linesizes = VisualisationConstants$linesizes
            fig =PlotSimulatedDataByBacteriaType(simulated_data(),
                                                 title_plot = "", 
                                                 ymin = NULL,
                                                 ymax = NULL,
                                                 tmax = input$Tmax,
                                                 minCFU = min.cfu,
                                                 text_size = 12,
                                                 linetypes = linetypes,
                                                 linesizes = linesizes,
                                                 colors = colors) 
            
            return(fig)
        }else{ #
            return(invisible(NULL))
        }
    })
})

# Run the application 
shinyApp(ui = ui, server = server)
