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
      
      
      S <- Surv(d$dt, d$e)

      return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S))
      
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
                                    textInput(inputId='hr', label='Hazard ratio', width = '90%' , ".75"),
                                  ) 
                                  
                                ) 
                             ),
                                 
                             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            

                             menuItem("Survival analysis", tabName = "OVERVIEW",  icon = icon("bar-chart-o"), selected = FALSE),
                             
                             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             menuItem("More",  startExpanded = FALSE,  icon = icon("bar-chart-o"),
                                      
                                      menuSubItem("a Partial log likelihood",        tabName = "RESULTS1"),
                                      menuSubItem("b Diagnostics",                tabName = "RESULTS3"),
                                      menuSubItem("c Explanation",               tabName = "HELP")
                                      
                             ),
                           
                            
                           
                           
                           menuItem("Code", icon = icon("bar-chart-o"),
                                    menuSubItem("Shiny",  
                                                icon = icon("send",lib='glyphicon'), 
                                                href = "https://raw.githubusercontent.com/eamonn2014/PARTIAL-LIKELIHOOD-METHOD/master/app.R"),
                                    
                                    
                                    menuSubItem("R",  
                                                icon = icon("send",lib='glyphicon'), 
                                                href = "https://raw.githubusercontent.com/eamonn2014/PARTIAL-LIKELIHOOD-METHOD/master/cox%20calculations.R") 
                                    
                                    
                                    
                                    # menuSubItem("Click for bells and whistles main app.",  
                                    #             icon = icon("send",lib='glyphicon'), 
                                    #             href = "https://eamonn3.shinyapps.io/LoQs/")
                           )
                           
                           
                            )
                           
                           
                           ),
      
      
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
                        ,plotlyOutput("plot1", height = "720px")
                    )
                    
                    ,box(
                        title='Difference in two Kaplan-Meier estimates with approximate confidence bands for differences'
                        ,status = "primary"
                        ,solidHeader = TRUE 
                        ,collapsible = TRUE 
                        ,plotOutput("plot2", height = "720px")
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
               title = "Partial Likelihood by hand! Here we show how to calculate the log likelihood given the actual model HR. 
                In reality a starting HR is supplied and an iterative process used to search for the value that maximises the log partial likelihood function.
               So the Beta obtained is the value that maximizes log PL. Note the estimate depends only on the
ranks of the event times, not their numerical values!"
                     ,status = "primary"
                     ,solidHeader = TRUE 
                     ,collapsible = TRUE 
                     , DT::dataTableOutput("mytable2")
                            ))),
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   tabItem("HELP", 
           box("", 
               textOutput("help"),
               
               
               br(),
               textOutput("help2"),
               withMathJax(),
               p(strong("To do spieglehalter explanantion/ equations and partial likelihood example")) ,
               
               p("where $a=\\sqrt{b}$ is the dependent variable (operationalized as the hazard rate at time t for subject i), x1 to xk are k independent variables or covariates,
and 1 to k are the regression coefficients; h0(t) is a baseline hazard
function and is left unspecified. The baseline hazard function can be
thought of as the hazard function for an individual whose covariates all
have values of 0"),  
               
               
               p("$$\\begin{align}
                     h_{i}  {(t)} = h_{0} {(t)} {exp} ({\\beta_1}{x_{i1}} + \\cdots{\\beta_k}{x_{ik}})
                      \\end{align}$$"),
               
               p("$$\\begin{align}
                     \\log h_{i}  {(t)} ={\\alpha{(t)}} + {\\beta_1}{x_{i1}} + \\cdots{\\beta_k}{x_{ik}}
                      \\end{align}$$"),
               
               
               
               p("$$\\begin{align}
                      \\sigma^2 \\\\
                      \\end{align}$$"),
               p(""), 
               p(strong("using a pooled estimate of the variance in each group ",HTML(" <em>i</em>")," this is just the mean variance")),
               
               p("$$\\begin{align}
                      s_p^2 =\\frac{\\sum  s_i^2}{I} \\approx \\sigma^2 \\\\
                      \\end{align}$$"),
               p(""), 
               
               p(strong("Here is the trick, we have another way to estimate the population variance, if the group means do not differ the sample means are normally 
                distributed with variance:")),
               
               p("$$\\begin{align}
                      \\frac{\\sigma^2}{n} = s_\\bar{Y}^2 \\\\
                      \\end{align}$$"),
               p(""), 
               
       
               
               
               
               
               )
   ),
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   tabItem("RESULTS3",
           fluidRow(
             box(
               title = "Log-log survivor plot; log[-log S(t)] as the vertical axis, and
log time as the horizontal axis"
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
               ,plotOutput("plot3", height = "720px")
             )
             
             ,box(
               title='
Repeated Cox regression coefficients estimates and confidence limits within time intervals. 
               The log hazard ratios are plotted against the mean failure/censoring time within the interval'
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
               ,plotOutput("plot4", height = "720px")
             )))
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   
   
   )
   ))
   


# create the server functions for the dashboard  
server <- function(input, output) { 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # https://stackoverflow.com/questions/55043092/r-shinydashboard-display-sum-of-selected-input-in-a-valuebox
    output$value1 <- renderValueBox({
      
      valueBox(
        value =  tags$p(paste0(formatz0(setUpByName())," / ",formatz0(setUpByNamea()) ," / ",formatz00(setUpByNameb()) ," / ",formatz1(setUpByNamec()) ," / ",formatz2(setUpByNamea()/setUpByNameb()  )    )
                        ,style = "font-size: 100%;")
        ,subtitle = tags$p('N; Events (a); Exposure (b); Median surv.; Hazard (a/b)', style = "font-size: 150%;")
        ,icon = icon("stats",lib='glyphicon')
        ,color = "yellow" )
        
    })
    
    output$value2 <- renderValueBox({
        
        valueBox(
          value =  tags$p(paste0(formatz0(setUpByName2())," / ",formatz0(setUpByName2a()) ," / ",formatz00(setUpByName2b()) ," / ",formatz1(setUpByName2c()) ," / ",formatz2(setUpByName2a()/setUpByName2b()  )    )
          ,style = "font-size: 100%;")
          ,subtitle = tags$p('N; Events (a); Exposure (b); Median surv.; Hazard (a/b)', style = "font-size: 150%;")
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
  
        return(list(  d=res$d, f=res$f, f1=res$f1, sf=res$sf, np=res$np , LL1=res$LL1, LL0=res$LL0, S=res$S))
        
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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
  
      survdiffplot(f, col='darkgreen' , xlab="Time")
      
     # survplot(f, conf='diffbands',col='purple',cex.aehaz=5,
      #         col.fill='blue'
       #                      , aehaz=TRUE, #times= c(5), 
        #      label.curves=list(method="arrow", cex=2), 
             #  label.curves=list(keys=1:2, cex=2),
         #      abbrev.label=TRUE, levels.only = FALSE)
    
        
    })
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot3 <- renderPlot({
      
      d <- dat()$d  # Get the  obj
      S <- dat()$S
      f <- dat()$f
      
      #plot(survfit(S~ d$trt), col=c("purple", "orange"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)" , lwd=3)
      survplot(f,logt=TRUE, loglog=TRUE, col=c("orange", "purple"))   #Check for Weibull-ness (linearity)
      
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot4<-renderPlot({     
      
      d <- dat()$d  # Get the  obj
      S <- dat()$S

      hazard.ratio.plot(d$trt, S, e=20, legendloc='ll', xlab='Time', antilog=FALSE, col='blue', smooth=TRUE)
      
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
       d$part2 <- guess*d$trt
       d$part3 <- log(rev(cumsum(rev(d$expB))))
      
       d$likelihoodi <- d$e*(d$part2 - d$part3)
       d$LL <- sum(d$likelihoodi)
    
       
       
        datatable(d, rownames=FALSE,
                         plugins = 'natural',
                         colnames=c('Time' = 'dt', 
                                    'Event or censored' = 'e', 
                                    'Treatment'='trt',
                                    'HR'='expB',
                                   # 'A'='part1',
                                    'logHR x trt'='part2',
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
           columns= c("Time","HR", 
                      #"A",
                      "logHR x trt","Log(exp HR) of each risk set",'Individual likelihoods'), 
           digits=4 ) %>%
        formatRound(
        columns= c("Sum the Individual likelihoods to give log likelihood"), 
        digits=1 )
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    output$help <- renderText({
      HTML("xxxxxxxxxxxxxxxxxx'")
    })
    output$help2 <- renderText({
      HTML("xxxxxxxxxxxxxxx.")
    })
   
}


shinyApp(ui, server)