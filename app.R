#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load the required packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    set.seed(333) # reproducible
    library(shiny)
    require(shinydashboard)
    library(ggplot2)
    library(dplyr)
    library(directlabels)
    library(shiny) 
    library(shinyalert)
    library(Hmisc)
    library(rms)
    library(ggplot2)
    library(tidyverse)
    library(plotly)
    library(survminer)
    library(rms)
    library(scales) # For the trans_format function
    options(max.print=1000000)    

   # https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
    formatz <- function(x){
        
        if (!is.na(x)  ) {
                  
                formatC(signif(x,digits=5), digits=5,format="fg", flag="#",big.mark=",")
             
         }
        
    }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logit <- function(p) log(1/(1/p-1))
    expit <- function(x) 1/(1/exp(x) + 1)
    inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
    is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
    
    options(width=200)
    options(scipen=999)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# function to create minor lines to match log tick values https://r-graphics.org/recipe-axes-axis-log-ticks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    coxdata <- function(n, allocation, hr) { 
      
      trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
      
      cens <- 15*runif(n)
      
      h <- exp(log(hr)*(trt==1))
      
      dt <- -log(runif(n))/h
      
      label(dt) <- 'Follow-up Time'
      
      e <- ifelse(dt <= cens,1,0)
      
      dt <- pmin(dt, cens)
      
      units(dt) <- "Year"
      
      d <- data.frame(cbind(dt, e, trt=trt))
      
      d <- plyr::arrange(d, dt)
      
      dd <- datadist(d)
      options(datadist='dd')
      
      S <- Surv(dt,e)
      f <- cph(S ~  trt, x=TRUE, y=TRUE,d)
      
      f1 <- survfit(S ~ trt, data = d)
      
      return(list(f=f,d=d, f1=f1))
      
    }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Dashboard header carrying the title of the dashboard
header <- dashboardHeader(title = "Transformations")  

#Sidebar content of the dashboard
sidebar <- dashboardSidebar(width=300,
                            br(),
                            tags$head(
                              tags$style(HTML('#resample{background-color:palegreen}'))
                            ),
                            actionButton("resample"," Hit to sample another data set", icon = icon("th"),  width =250  ),
                            
                            sidebarMenu(
                                id = "tabs",
                             #  br(),
                             
                                
                          
                                
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             menuItem("Define parameters ", icon = icon("bar-chart-o"),
                                splitLayout(
                                    
                                    tags$div(
                                        textInput(inputId="n", label='N', width = '90%' , value="100"),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='allocation', label='allocation', width = '90%' , ".5"),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='hr', label='Slrope', width = '90%' , "2"),
                                    )
                                    
                                ) #,
                                
                                
                             ),
                                 
                             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             
                             
                             
                             
                             menuItem("Code & link to explanation", icon = icon("bar-chart-o"),
                                      menuSubItem("Shiny",  
                                                  icon = icon("send",lib='glyphicon'), 
                                                  href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/dashboard1/app.R"),
                                      
                                      
                                      menuSubItem("R",  
                                                  icon = icon("send",lib='glyphicon'), 
                                                  href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/Rcode.R") ,
                                      
                                      
                                      
                                      menuSubItem("Click for bells and whistles main app.",  
                                                  icon = icon("send",lib='glyphicon'), 
                                                  href = "https://eamonn3.shinyapps.io/LoQs/")
                                      
                                      
                                      
                                      
                                      
                             )
                             
                             
                             
                          
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            )
                            
                                                 
)

    frow1 <- fluidRow(
        valueBoxOutput("value1")
        ,valueBoxOutput("value2")
        ,valueBoxOutput("value3")
     )
    
    frow2 <- fluidRow(
        
        box(
            title = "Fitted Analysis Model"
            ,status = "primary"
            ,solidHeader = TRUE 
            ,collapsible = TRUE 
            ,plotlyOutput("plot1", height = "750px")
        )
        
        ,box(
            title = "Fitted Analysis Model with logarithmic transformation"
            ,status = "primary"
            ,solidHeader = TRUE 
            ,collapsible = TRUE 
            ,plotlyOutput("plot2", height = "750px")
        ) 
        
    )


    # combine the two fluid rows to make the body
    body <- dashboardBody(frow1, frow2 
                          
                          )
    #completing the ui part with dashboardPage
    
    ui <- dashboardPage(title = 'Examples of transforming data for analysis', header, sidebar, body, skin='blue')


# create the server functions for the dashboard  
server <- function(input, output) { 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # https://stackoverflow.com/questions/55043092/r-shinydashboard-display-sum-of-selected-input-in-a-valuebox
    output$value1 <- renderValueBox({
      
        valueBox( 
         formatz(setUpByName())
         ,subtitle = tags$p("Sum of squares of residuals (smaller the better)", style = "font-size: 150%;")
         ,icon = icon("server")
         ,color = "purple")
        
    })
    
    output$value2 <- renderValueBox({
        
        valueBox(
            formatz(setUpByName2())
            ,subtitle = tags$p('Model sigma', style = "font-size: 150%;")
            ,icon = icon("stats",lib='glyphicon')
            ,color = "green")
        
    })
    
    output$value3 <- renderValueBox({
        
        valueBox(
            value =  tags$p(paste0(formatz(setUpByName3())," ( ",formatz(setUpByName3()),"; ",formatz(setUpByName3())," )")
            ,style = "font-size: 100%;")
            ,subtitle = tags$p(paste0("Prediction at ",formatz(setUpByName3())," with 95% confidence"), style = "font-size: 150%;")
            ,icon = icon("education",lib='glyphicon')
            ,color = "yellow")
        
    })

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is where a new sample is instigated and inputs converted to numeric
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    random.sample <- reactive({
        
        foo <- input$resample
        
        n <- as.numeric(input$n)
        
        allocation <-as.numeric(input$allocation)

        hr <- as.numeric(input$hr)
      
        return(list(  
            n=n,
            allocation =allocation,
            hr=hr
        ))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GENERATE THE DATA Execute analysis
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dat <- reactive({
        
        sample <- random.sample()
        
        n=         sample$n
        allocation=sample$allocation
        hr=        sample$hr
        
        res <- coxdata(n, allocation, hr)
        
        
        return(list(  d=res$d, f=res$f, f1=res$f1))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    setUpByName <- reactive ({
        f <-dat()$f  # Get the  data
        y <- as.numeric(as.character(f$coefficients))
        return(y)
    })
    
    setUpByName2 <- reactive ({
      f <- dat()$f  # Get the  data
      y <- as.numeric(as.character(f$coefficients))
      return(y)
    })

    setUpByName3 <- reactive ({
      f <- dat()$f  # Get the  data
      y <- as.numeric(as.character(f$coefficients))
      return(y)
    })
    
   
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAIN PLOT! updated with log transformation  option
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot1 <- renderPlotly({
      
      f <- dat()$f1  # Get the  data
      
     
      p1 <- ggsurvplot(f, main = "Kaplan-Meier Curve for the NCCTG Lung Cancer Data")
      ggplotly(p1[[1]])
      
     # plot_ly(iris, x = ~get(input$choice), y = ~Sepal.Length, type = 'scatter', mode = 'markers')
      
      
    })
    
    
    output$plot2<-renderPlotly({     
        
      f <- dat()$f1  # Get the  data
      
      p1 <- ggsurvplot(f, main = "Kaplan-Meier Curve for the NCCTG Lung Cancer Data")
      plotly::ggplotly(p1[[1]])
        
    })

    
}


shinyApp(ui, server)