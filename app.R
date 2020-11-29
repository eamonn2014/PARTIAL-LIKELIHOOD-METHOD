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
    library(DT)
    library(survival)
    options(max.print=1000000)    

   # https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
    formatz <- function(x){
        
        if (!is.na(x)  ) {
                  
                formatC(signif(x,digits=5), digits=5,format="fg", flag="#",big.mark=",")
             
         }
        
    }
    
    formatz0 <- function(x){
      sprintf(x, fmt = '%s')  
    }
    formatz1 <- function(x){
      sprintf(x, fmt = '%#.1f')  
    }
    formatz2 <- function(x){
      sprintf(x, fmt = '%#.2f')  
    }
    formatz00 <- function(x){
      round(x,0) 
    }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logit <- function(p) log(1/(1/p-1))
    expit <- function(x) 1/(1/exp(x) + 1)
    inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
    is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
    
    options(width=200)
    options(scipen=999)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# function to create data and analyse
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    coxdata <- function(n, allocation, hr, baseline) { 
      
      #n=1000; allocation =.5; hr=2
      
      trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
      
      cens <- 15*runif(n)
   
      h <- baseline*exp(log(hr)*(trt==1))
      
      dt <- -log(runif(n))/h
      label(dt) <- 'Follow-up Time'
      
      e <- ifelse(dt <= cens,1,0)
      
      dt <- pmin(dt, cens)
      
      units(dt) <- "Year"
      
      d <<- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?
      
      dd <<- datadist(d)
      options(datadist='dd')
      
     # S <- Surv(dt,e)
      f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
      
      LL1 <- f$loglik[2]
      LL0 <- f$loglik[1]
      
      sf <- summary(f)
      
      f1 <- survfit(Surv(dt,e) ~ trt, data = d)
      
      np <- npsurv(Surv(dt,e) ~ trt,d)

      return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0))
      
    }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ui <- dashboardPage( 
#Dashboard header carrying the title of the dashboard,
      dashboardHeader(title = "COX PH")  ,
      
      #Sidebar content of the dashboard
      sidebar <- dashboardSidebar(width=300,
                            br(),
                            tags$head(
                              tags$style(HTML('#resample{background-color:palegreen}'))
                            ),
                            actionButton("resample"," Hit to sample another data set", icon = icon("th"),  width =250  ),
                            
                            sidebarMenu(
                                id = "tabs",
                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             menuItem("Define parameters ", icon = icon("bar-chart-o"),
                                splitLayout(
                                    
                                    tags$div(
                                        textInput(inputId="n", label='Total sample size', width = '90%' , value="800"),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='allocation', label='random allocation', width = '90%' , ".5"),
                                    ) 
                                ) ,
                                
                                splitLayout(
                                  
                                  tags$div(
                                    textInput(inputId='baseline', label='baseline hazard', width = '90%' , ".4"),
                                  ),
                                  
                                  tags$div(
                                    textInput(inputId='hr', label='Hazard ratio', width = '90%' , "2"),
                                  ) 
                                  
                                ) 
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
                             ),

                             menuItem("Survival analysis", tabName = "OVERVIEW",  icon = icon("bar-chart-o"), selected = FALSE),
                             
                             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             menuItem("3 Partial likelihood",  startExpanded = FALSE,  icon = icon("bar-chart-o"),
                                      
                                      menuSubItem("i Listing & Likelihood",                 tabName = "RESULTS1"),
                                      menuSubItem("ii xxxxxxxxxx",                 tabName = "RESULTS2")
                                      
                             ))),
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        
           dashboardBody(
            
              fluidRow(
                   valueBoxOutput("value1")
                  ,valueBoxOutput("value2")
                  ,valueBoxOutput("value3")
               ),
            
            tabItems(
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            tabItem("OVERVIEW",
                fluidRow(
                    box(
                        title = "Kaplan-Meier Curve"
                        ,status = "primary"
                        ,solidHeader = TRUE 
                        ,collapsible = TRUE 
                        ,plotlyOutput("plot1", height = "750px")
                    )
                    
                    ,box(
                        title='Difference in two Kaplan-Meier estimates with approximate confidence bands for differences'
                        ,status = "primary"
                        ,solidHeader = TRUE 
                        ,collapsible = TRUE 
                        ,plotOutput("plot2", height = "750px")
                    ))),               
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     tabItem("RESULTS1",
        fluidRow(        
           box(
              title = "Data listing"
                 ,status = "primary"
                   ,solidHeader = TRUE 
                    ,collapsible = TRUE 
              , DT::dataTableOutput("mytable")
           )
                                    
            ,box(
               title = "Partial Likelihood by hand! Here we show how to calculate the log likelihood given the model HR. In reality a starting HR is supplied and an iterative process used to search for the value that maximises the log partial likelihood function."
                     ,status = "primary"
                          ,solidHeader = TRUE 
                             ,collapsible = TRUE 
                                      #,plotOutput("PF", height = "500px")
              , DT::dataTableOutput("mytable2")
                                ) 
                                  )
                          ),
                          #~~~~~~~~~~~~~
                          
                          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                          tabItem("RESULTS2",
                                  box(" ",
                                     # DT::dataTableOutput("mytable2")
                                  )  )  )
   ))
   


# create the server functions for the dashboard  
server <- function(input, output) { 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # https://stackoverflow.com/questions/55043092/r-shinydashboard-display-sum-of-selected-input-in-a-valuebox
    output$value1 <- renderValueBox({
      
      valueBox(
        value =  tags$p(paste0(formatz0(setUpByName())," ; ",formatz0(setUpByNamea()) ," ; ",formatz00(setUpByNameb()) ," ; ",formatz1(setUpByNamec()) ," ; ",formatz2(setUpByNamea()/setUpByNameb()  )    )
                        ,style = "font-size: 100%;")
        ,subtitle = tags$p('Treat. 0 : N; Events; Exposure; Median survival; Hazard', style = "font-size: 150%;")
        ,icon = icon("stats",lib='glyphicon')
        ,color = "yellow" )
        
    })
    
    output$value2 <- renderValueBox({
        
        valueBox(
          value =  tags$p(paste0(formatz0(setUpByName2())," ; ",formatz0(setUpByName2a()) ," ; ",formatz00(setUpByName2b()) ," ; ",formatz1(setUpByName2c()) ," ; ",formatz2(setUpByName2a()/setUpByName2b()  )    )
          ,style = "font-size: 100%;")
          ,subtitle = tags$p('Treat. 1 : N; Events; Exposure; Median survival; Hazard', style = "font-size: 150%;")
          ,icon = icon("stats",lib='glyphicon')
            ,color = "purple")
        
    })
    
    output$value3 <- renderValueBox({
        
        valueBox(
            value =  tags$p(paste0(formatz2(setUpByName4())," ( ",formatz2(setUpByName5()),", ",formatz2(setUpByName6())," ) " ," ; ",formatz1(setUpByNameLL()))
            ,style = "font-size: 100%;")
            ,subtitle = tags$p(paste0("Hazard ratio trt 1 v 0 with 95% conf. ; log likelihood"), style = "font-size: 150%;")
            ,icon = icon("education",lib='glyphicon')
            ,color = "green")
        
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
      
        baseline <- as.numeric(input$baseline)
        
        return(list(  
            n=n,
            allocation =allocation,
            hr=hr,
            baseline=baseline
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
        baseline=  sample$baseline
        
        res <- coxdata(n, allocation, hr, baseline)
  
        return(list(  d=res$d, f=res$f, f1=res$f1, sf=res$sf, np=res$np , LL1=res$LL1, LL0=res$LL0))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    setUpByName <- reactive ({
        f <-dat()$np  # Get the  data
        f <-  f$n[1]
        y <- as.numeric(as.character(f))
        return(y)
    })

    setUpByNamea <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  f$numevents[1]
      y <- as.numeric(as.character(f))
      return(y)
    })
    
    setUpByNameb <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  f$exposure[1]
      y <- as.numeric(as.character(f))
      return(y)
    })
    
    setUpByNamec <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  summary(f)$table[,'median'][1]
      y <- as.numeric(as.character(f))
      return(y)
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    setUpByName2 <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  f$n[2]
      y <- as.numeric(as.character(f))
      return(y)
    })
    
    setUpByName2a <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  f$numevents[2]
      y <- as.numeric(as.character(f))
      return(y)
    })
    

    setUpByName2b <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  f$exposure[2]
      y <- as.numeric(as.character(f))
      return(y)
    })
    
    setUpByName2c <- reactive ({
      f <-dat()$np  # Get the  data
      f <-  summary(f)$table[,'median'][2]
      y <- as.numeric(as.character(f))
      return(y)
    })
    
    setUpByNameLL <- reactive ({
      f <-dat()$LL1  
      y <- as.numeric(as.character(f))
      return(y)
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    setUpByName3 <- reactive ({
      f <- dat()$f  # Get the  data
      y <- as.numeric(as.character(f$coefficients))
      return(y)
    })
    
    
    setUpByName4 <- reactive ({
      f <-dat()$sf  # Get the  data
      y <- as.numeric(as.character(f[2,c("Effect")]))
      return(y)
    })
    
  
    setUpByName5 <- reactive ({
      f <-dat()$sf  # Get the  data
      y <- as.numeric(as.character(f[2,c("Lower 0.95")]))
      return(y)
    })
    
    setUpByName6 <- reactive ({
      f <-dat()$sf  # Get the  data
      y <- as.numeric(as.character(f[2,c("Upper 0.95")]))
      return(y)
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAIN PLOT! updated with log transformation  option
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot1 <- renderPlotly({
      
      f <- dat()$f1  # Get the  obj

      p1 <- ggsurvplot(f, main = "Kaplan-Meier Curve", 
                       palette = c("orange", "purple") #,#conf.int = TRUE,
                      # ggtheme = theme_bw() # Change ggplot2 theme
                       )
      ggplotly(p1[[1]])
      
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot2<-renderPlot({     
        
      f <- dat()$np  # Get the  data
  
      survdiffplot(f, col='darkgreen' )
      
     # survplot(f, conf='diffbands',col='purple',cex.aehaz=5,
      #         col.fill='blue'
       #                      , aehaz=TRUE, #times= c(5), 
        #      label.curves=list(method="arrow", cex=2), 
             #  label.curves=list(keys=1:2, cex=2),
         #      abbrev.label=TRUE, levels.only = FALSE)
    
        
    })

    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$mytable <- DT::renderDataTable({
      
      d=dat()$d
      
      d <- plyr::arrange(d, dt)
      
      DT::datatable(d, rownames=FALSE,
                    plugins = 'natural',
                   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
                              'Treatment'='trt'),
                    
                    options = list(
                   #  dom = 't',
                       columnDefs = list(list(type = 'natural', targets = c(1,2)))
                     )
                   ) %>%
        
        formatRound(
          columns= c("Time" ), 
          digits=4 )
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~maximum likelihood
    
    output$mytable2 <- DT::renderDataTable({
      
       d=dat()$d
       
       f <- dat()$f  # Get the  data
       y <- as.numeric(as.character(f$coefficients))
       guess=y  # we use the actual model hr estimate
      # sample <- random.sample()
      # hr=   sample$hr
   
       d <- plyr::arrange(d, dt)
      
      
       d$expB <- exp(guess*d$trt)
     
       d$part3 <- log(rev(cumsum(rev(d$expB))))
       d$part2 <- guess*d$trt
       d$likelihoodi <- d$e*(d$part2 - d$part3)
       d$LL <- sum(d$likelihoodi)
    
       
       
        datatable(d, rownames=FALSE,
                         plugins = 'natural',
                         colnames=c('Time' = 'dt', 
                                    'Event or censored' = 'e', 
                                    'Treatment'='trt',
                                    'Exponentiated HR'='expB',
                                   # 'A'='part1',
                                    'HR x trt'='part2',
                                    'Log(exp HR) of each risk set'='part3',
                                    'Individual likelihoods' ='likelihoodi',
                                   'Sum the Individual likelihoods to give log likelihood' ='LL'
                                    ),
                    
                    options = list(
                      #  dom = 't',
                      columnDefs = list(list(type = 'natural', targets = c(1,2)))
                    )
      ) %>%
         
         formatRound(
           columns= c("Time","Exponentiated HR", 
                      #"A",
                      "HR x trt","Log(exp HR) of each risk set",'Individual likelihoods'), 
           digits=4 ) %>%
        formatRound(
        columns= c("Sum the Individual likelihoods to give log likelihood"), 
        digits=1 )
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
}


shinyApp(ui, server)