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

source("R/helpers.R")
mg_count = (10^12) 
default_params = SetDefaultParameters(mg_count)
min_params = SetMinParameters(mg_count)
max_params = SetMaxParameters(mg_count)


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    withMathJax(),
    # Application title
    titlePanel("PhARTNER: PhAge ResisTaNce EmeRgence"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h2("Model Parameters:"),
            br(),
            br(),
            h4("Experimental setup:"),
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
            h4("Phage virulence:"),
            numericInput("decay",
                         label = "\\( \\gamma \\): phage decay rate [1/h]",
                         min = 0,
                         step = 0.1,
                         max = 0.2,
                         default_params$decay),
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
            h4("Externally added depolymerase:"),
            numericInput("e0",
                         label = "\\(E_0\\): Initial depolymerase concentration [g/L]",
                         min = 0,
                         step = 0.1,
                         default_params$e0),
            numericInput("V_depo",
                         label = "\\(V^D\\): Maximal depolymerase activity rate (scaled) [g/L*h]",
                         min = 0,
                         default_params$V_depo),
            numericInput("K_depo",
                         label = "\\(K^D\\): Michaelis-Menten constant for depolymerase activity [g bacteria / L]",
                         min = 0,
                         default_params$K_depo),
            numericInput("depo_decay_rate",
                         label = "\\(\\gamma^E\\): Depolymerase decay rate[1/h]",
                         min = 0,
                         step = 0.1,
                         default_params$depo_decay_rate),
            br(),
            br(),
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
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            h4("Phage therapy effects:"),
            textOutput("errortext"),
            plotOutput("timeplot"),
            br(),
            br(),
            h4("Bacterial growth by bacteria type:"),
            plotOutput("bacteriatypeplot"),
            #plotOutput("bacteriatypeplot2")
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
            !is.null(input$e0) && input$e0>=min_params$e0 && input$e0<=max_params$e0 &&
            !is.null(input$MOI) && input$MOI>=min_params$MOI && input$MOI<=max_params$MOI &&
            !is.null(input$G0) && input$G0>=min_params$G0 && input$G0<=max_params$G0 &&
            !is.null(input$Tmax) && input$Tmax>=min_params$Tmax && input$Tmax<=max_params$Tmax })
    
    
    output$errortext <- renderText({
        if(isValid_input()){ invisible(NULL) 
        }else{ #
            "Parameter values are outside a reasonable range.
             Note: They must be numeric and non-negative"
        }
    })
    
    simulated_data =  reactive({Simulate_Phage_Coctail(Vh=input$Vh, 
                                             Kh=input$Kh,
                                             a=input$a, 
                                             beta_non_depo=input$beta_non_depo, 
                                             beta_depo=input$beta_depo, 
                                             phi_non_depo=input$phi_non_depo, 
                                             phi_depo=input$phi_depo, 
                                             decay=input$decay, 
                                             V_depo=input$V_depo, 
                                             K_depo=input$K_depo, 
                                             depo_decay_rate=input$depo_decay_rate,
                                             epsilonB2toB1=default_params$epsilonB2toB1,#input$epsilonB2toB1, 
                                             epsilonB1toB2=default_params$epsilonB1toB2,#input$epsilonB1toB2,
                                             epsilonB2toB3=default_params$epsilonB2toB3,#input$epsilonB2toB3, 
                                             epsilonB3toB2=default_params$epsilonB3toB2,#input$epsilonB3toB2,
                                             epsilonB1toB3=default_params$epsilonB1toB3,#input$epsilonB1toB3, 
                                             epsilonB3toB1=default_params$epsilonB3toB1,#input$epsilonB3toB1,  
                                             B0=as.numeric(input$B0)/(mg_count/1000), # so that it is in g/L
                                             e0=input$e0, 
                                             MOI =input$MOI, 
                                             G0=input$G0, 
                                             Tmax = input$Tmax) %>% 
        mutate(Treatment = plyr::revalue(Treatment,
                                         c("no-depo-phage (KP15/KP27)" ="no-depo-phage", 
                                           "no-depo-phage+depo (KP15/KP27 + KP34p57)" ="no-depo-phage+depo",
                                           "phage cocktail (KP15/KP27 + KP34)" = "phage cocktail",
                                           "depo-phage (KP34)" = "depo-phage" ))) %>%
            filter(Treatment != "External depo")
            
        
    })
    
    output$timeplot <- renderPlot({
        if(isValid_input()){ 
            
            VisualisationConstants = GetVisualisationConstantsApp()
            my_theme = VisualisationConstants$my_theme
            colors = VisualisationConstants$colors
            linetypes = VisualisationConstants$linetypes
            linesizes = VisualisationConstants$linesizes
            
        
        fig =PlotSimulatedPhageAndBacteria(simulated_data(),
                                           title_plot = "", 
                                           colors = colors,
                                           ymin = NULL,
                                           ymax = NULL,
                                           tmax = input$Tmax,
                                           minCFU = 6.3*10^6,
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
                                               minCFU = 6.3*10^6,
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
                                                 colors = colors,
                                                 ymin = NULL,
                                                 ymax = NULL,
                                                 tmax = input$Tmax,
                                                 minCFU = 6.3*10^6,
                                                 text_size = 12,
                                                 linetypes = linetypes,
                                                 linesizes = linesizes) 
            
            return(fig)
        }else{ #
            return(invisible(NULL))
        }
    })
})

# Run the application 
shinyApp(ui = ui, server = server)
