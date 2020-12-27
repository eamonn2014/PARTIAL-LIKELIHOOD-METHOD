#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load the required packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# https://thestatsgeek.com/2014/03/28/interpreting-changes-in-hazard-and-hazard-ratios/

  set.seed(333) # reproducible
  library(shiny)
  require(shinydashboard)
  library(ggplot2)
  library(dplyr)
  library(directlabels)
  library(Hmisc)
  library(ggplot2)
  library(tidyverse)
  library(plotly)
  library(survminer)
  library(rms)
  # library(scales) # For the trans_format function
  # library(shinyalert)
  library(DT)
  library(survival)
  options(max.print=1000000)    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # function to format 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  formatz4 <- function(x){
    sprintf(x, fmt = '%#.4f')  
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
#
coxdata <- function(n, allocation, hr, baseline) { 
  
  #n=1000; allocation =.5; hr=2; baseline=.4
  
  trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
  
  cens <- 15*runif(n)
  
  h <- baseline*exp(log(hr)*(trt==1))  # hazard function h(t)
  
  dt <- -log(runif(n))/h
  
  label(dt) <- 'Follow-up Time'
  
  e <- ifelse(dt <= cens,1,0)
  
  dt <- pmin(dt, cens)
  
  units(dt) <- "Year"
  
  d <<- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?
  
  dd <<- datadist(d)
  options(datadist='dd')
  
  foo <-d
  # S <- Surv(dt,e)
  f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
  f0 <- f$coefficients[[1]] #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  LL1 <- f$loglik[2]
  LL0 <- f$loglik[1]
  
  sf <- summary(f)
  
  f1 <- survfit(Surv(dt,e) ~ trt, data = d)
  
  np <- npsurv(Surv(dt,e) ~ trt,d)
  
  S <- Surv(d$dt, d$e)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  d <- plyr::arrange(d,dt)
  d$dt <- dd<- NULL
  d$dt <- sort(2*rexp(nrow(d)))# new times
  
  dx <<- datadist(d)
  options(datadist='dx')
  f0a <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
  f0a <- f0a$coefficients[[1]]
  f2 <- survfit(Surv(dt ,e)  ~ trt, data = d)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S,                  
              
         f0=f0,  f2=f2, f0a=f0a, foo=foo))
  
}

# dummy <- coxdata(n=1000, allocation =.5, hr=2, baseline=.4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function that can use the information from coxdata function to learn about the 
# behind the scenes cal in Cox PH, we allow a guess at HR and see the log likelihood
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## enter hr guess, data, dead var, trt var and treat var

loglike2 <- function(x, dat, dead, indep , time) {
  
  dd <- dat          # make another data object
  dd$dead <- dead    # take the key variables to run Cox PH
  dd$indep <- indep
  dd$time <- time
  
  ## run the analysis to get hr and log like
  ddd <<- datadist(dd)
  options(datadist='ddd')
  
  S <- Surv(time,dead)  # run Cox PH
  f <- cph(S ~  indep, x=TRUE, y=TRUE,dd)
  
  #~~~~~~~~~~extract hr and loglikelihood at null and maximised log likelihood
  dd$hr <- exp(f$coefficients)
  dd$lognull <- f$loglik[[1]]
  dd$lognmax <- f$loglik[[2]]
  
  #~~~~~~~~~using our guess x calculate log likelihood by jand
  dd$expB <- exp(x*dd$indep)
  dd$part1 <- dd$dead  
  dd$part2 <- x*dd$indep 
  
  dd <- plyr::arrange(dd,time)
  
  dd$part3 <- log(rev(cumsum(rev(dd$expB))))
  
  dd$guess <- exp(x)
  dd$likelihoodi <- dd$part1*(dd$part2 - dd$part3)
  dd$L <- sum(dd$likelihoodi)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dd <- as.data.frame(dd)
  dd$dead <- dd$indep <- dd$part1 <- dd$part2 <- dd$part3 <- NULL
   
  return(dd)
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Start app
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ui <- dashboardPage(  title="Survival Analysis",
  # Dashboard header carrying the title of the dashboard,
  
  dashboardHeader(title = h4(HTML("Cox proportional hazards<br/>and partial log likelihood"))),
  #Sidebar content of the dashboard
  sidebar <- dashboardSidebar(width=300,
                            
                               sidebarMenu(
                                
                                id = "tabs",
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               
                                br(),
                                tags$head(
                                  tags$style(HTML('#resample{background-color:palegreen}'))
                                ),
                                actionButton("resample"," Hit to sample another data set", icon = icon("th"),  width =250  ),
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
                                
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                menuItem("Analyses",  startExpanded = FALSE,  icon = icon("bar-chart-o"),
                                         #~~~~~~~~~~~~~~~~~~~~~~~~
                                         menuSubItem("Kaplan Meier (landing page)",    tabName = "OVERVIEW",  icon = icon("bar-chart-o"), selected = FALSE),
                                         menuSubItem("KM diagnostics",                 tabName = "RESULTS2",  icon = icon("bar-chart-o")),
                                         menuSubItem("Cox proportional hazards",       tabName = "RESULTS3",  icon = icon("bar-chart-o")),
                                         menuSubItem("Hazard ratio over time",         tabName = "RESULTS4",  icon = icon("bar-chart-o")),
                                         menuSubItem("Partial log likelihood",         tabName = "RESULTS1",  icon = icon("table")),
                                         #~~~~~~~~~~~~~~~~~~~~~~~~
                                         menuSubItem("Only ranks of event times needed!", tabName = "OVERVIEW2",  icon = icon("bar-chart-o"), selected = FALSE),
                                         #~~~~~~~~~~~~~~~~~~~~~~~~
                                         menuSubItem("Model assumptions", tabName = "OVERVIEW3",  icon = icon("bar-chart-o"), selected = FALSE),
                                         #~~~~~~~~~~~~~~~~~~~~~~~~
                                         menuSubItem("KM lifetable",                   tabName = "KMTABLE",  icon = icon("list-alt")),
                                         #~~~~~~~~~~~~~~~~~~~~~~~~
                                         menuItem("Partial likelihood exercise",  startExpanded = FALSE,    icon = icon("table"),
                                                 
                                                  tags$div(
                                                    textInput(inputId="guess", label='Enter a guess at Hazard Ratio (defaulted to null)', width = '90%' , value="1"),
                                                  ),
                                                  
                                                  menuSubItem("Hit to reveal Partial log likelihood",  tabName = "partial")
                                         )
                                         #~~~~~~~~~~~~~~~~~~~~~~~~
                                       
                                         
                                        
                               ),
                                #
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               menuItem("Change in hazard",  startExpanded = FALSE,    icon = icon("table")  ,
                                        
                                        tags$div(
                                          textInput(inputId="base", label='Enter a baseline hazard', width = '90%' , value="0.03"),
                                        ),
                                        tags$div(
                                          textInput(inputId="cens", label='Enter a censoring hazard', width = '90%' , value="0.02"),
                                        ),
                                        tags$div(
                                          textInput(inputId="hr2", label='Enter a hazard ratio', width = '90%' , value="1.2"),
                                        ),
                                        tags$div(
                                          textInput(inputId="per", label='Enter a survival probability',        width = '90%' , value="0.70"),
                                          textInput(inputId="per2", label='Enter another survival probability', width = '90%' , value="0.50"),
                                        ),
                                        
                                     menuSubItem("Hit to reveal change in hazard",  tabName = "Change"),
                                     menuSubItem("Hit to reveal hazard",  tabName = "Changeh")    ,
                                     menuSubItem("Hit to reveal KM & smoothed survival curve",  tabName = "Changeh2") 
                            ),
                               
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            menuItem("Power",  startExpanded = FALSE,    icon = icon("table")  ,
                                     
                                     
                                     textInput('ss', width = '90%' ,
                                               strong("Enter two survival probs & times"), "0.7, 0.5, 17, 23"),
                                     
                                      textInput('tt', width = '90%' ,
                                               strong("Enter ctrl, intervention sample size"), "500, 500"),
                                     
                                     textInput('af',  width = '90%' ,
                                               strong("Accrual time and follow up"), "36, 60"),
                                     # Here, let us accrue patients over three years, and
                                     #follow them for an additional seven years
                                   
                                     tags$div(
                                       textInput(inputId="t2", label='intervention arm non-compliance', width = '90%' , value="0.1"),
                                     ),
                                    
                                     
                                     ##
                                     splitLayout(
                                       tags$div(
                                         textInput(inputId="hrx", label='hazard ratio',   value="1.2"),
                                       ),
                                       tags$div(
                                         textInput(inputId="sim", label='No. of simulations',  value="500"),
                                       )
                                       
                                     ),
                                     
                                     ##
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     

                                   menuSubItem("Hit to reveal power",  tabName = "power")
                                  
                            ),
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               
                            menuItem("Parametric Survival Distributions",  startExpanded = FALSE,    icon = icon("table"),
                                     
                                     tags$div(
                                       textInput(inputId="shape", label='Weibull shape', width = '90%' , value="1"),
                                     ),
                                     
                                     tags$div(
                                       textInput(inputId="scale", label='Weibull scale', width = '90%' , value="0.03"),
                                     ),
                                     
                                     menuSubItem("Hit to reveal Survival hazard relationship",  tabName = "survhaz")
                            ),
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               
                               menuItem("Explanation",                    tabName = "HELP",icon = icon("bar-chart-o"), selected = FALSE),
                               menuItem("Wiki", tabName = "Wiki",  icon = icon("bar-chart-o"), selected = FALSE),
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               menuItem("Code", icon = icon("bar-chart-o"),
                                         menuSubItem("Shiny",  
                                                     icon = icon("send",lib='glyphicon'), 
                                                     href = "https://raw.githubusercontent.com/eamonn2014/PARTIAL-LIKELIHOOD-METHOD/master/app.R"),
                                         
                                         menuSubItem("R",  
                                                     icon = icon("send",lib='glyphicon'), 
                                                     href = "https://raw.githubusercontent.com/eamonn2014/PARTIAL-LIKELIHOOD-METHOD/master/cox%20calculations.R") 
                                  ),
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                         menuItem("References", icon = icon("bar-chart-o"),
                                         
                                         menuSubItem(h5(HTML( "Regression Models and Life-Tables")),  
                                                     icon = icon("send",lib='glyphicon'), 
                                                     href = "http://www.stat.cmu.edu/~ryantibs/journalclub/cox_1972.pdf"),
                                         
                                         menuSubItem(h5(HTML( "Individual survival time prediction <br/>using statistical models")),
                                                     icon = icon("send",lib='glyphicon'), 
                                                     href = "https://jme.bmj.com/content/medethics/31/12/703.full.pdf") ,
                                         #dashboardHeader(title = h4(HTML("This title<br/>is just way too long")))
                                         
                                         menuSubItem( h5(HTML("Can we say whether a drug would <br/>have enabled someone to <br/>live longer? Sadly not")),  
                                                      icon = icon("send",lib='glyphicon'), 
                                                      href = "https://understandinguncertainty.org/node/759"),
                                         
                                         menuSubItem( h5(HTML("Analysis of time-to-event for observational studies <br/>Guidance to the use of intensity models")),  
                                                      icon = icon("send",lib='glyphicon'), 
                                                      href = "https://github.com/eamonn2014/PARTIAL-LIKELIHOOD-METHOD/blob/master/Analysis%20of%20time-to-event%20for%20observational%20studies.pdf"),
                                         
                                         menuSubItem( h5(HTML("Cox's proportional hazards regression")),  
                                                      icon = icon("send",lib='glyphicon'), 
                                                      href = "https://influentialpoints.com/Training/coxs_proportional_hazards_regression_model-principles-properties-assumptions.htm"),
                                         
                                         menuSubItem( h5(HTML("Frank Harrell cph function")),  
                                                      icon = icon("send",lib='glyphicon'), 
                                                      href = "https://rdrr.io/cran/rms/man/cph.html")
                                         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                          
                                )
                             )
                              
    ),
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dashboardBody(
    # https://stackoverflow.com/questions/54876731/inline-latex-equations-in-shiny-app-with-mathjax
    tags$head(
      tags$link(rel="stylesheet", 
                href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css", 
                integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
                crossorigin="anonymous"),
      HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
      HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
      HTML('
    <script>
      document.addEventListener("DOMContentLoaded", function(){
        renderMathInElement(document.body, {
          delimiters: [{left: "$", right: "$", display: false}]
        });
      })
    </script>')
    ),
    
    fluidRow(
       valueBoxOutput("value1")
      ,valueBoxOutput("value2")
      ,valueBoxOutput("value3")
    ),
    
    tabItems(

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      tabItem("Wiki", 
              fluidRow(
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                box(  
                  title='Wiki'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  #textOutput("help"),
                  , p("Exploring the Cox proportional hazards model and the partial likelihood function"),
                )
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ,box(
                  title=' '
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE ,
                  # ,plotOutput("plot4", height = "720px")
                  p("xxxxxxxxxxxxxxxxxxxx"),
                ),  # box end
              )
      ),
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      tabItem("OVERVIEW",
              fluidRow(
                box(
                  title = "Kaplan-Meier Curve"
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotlyOutput("plot1", height = "720px"),
                  
                  #verbatimTextOutput("help2") ,
                  #p("zzzzzzzzzzzzzzz"),
                  h5(textOutput("Staff_name"))
                )
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ,box(
                  title='Difference in two Kaplan-Meier estimates with approximate confidence bands for differences'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("plot2", height = "720px")
                ))),               
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      tabItem("OVERVIEW2",
              fluidRow(
                box(
                  title =   "Kaplan-Meier curve"   # uiOutput('product'), 
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotlyOutput("plot5a", height = "720px")
                  ,h5(textOutput("info2"))
                )
                
                ,box(
                  title="KM, rank order of events preserved, sorted random values replace original numerical event times"
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotlyOutput("plot5b", height = "720px"),
                  h5(textOutput("info"))
                ))),   
      
   
   tabItem("OVERVIEW3",
           fluidRow(
             box(
               title =   "Altschuler-Nelson-Fleming-Harrington non parametric survival estimates and Cox-Breslow estimates"   # uiOutput('product'), 
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
               ,plotOutput("FH", height = "720px")
               #,h5(textOutput("info2"))
             )
             
             ,box(
               title=" "
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
               ,p("If the predicted survival curves from the fitted Cox model are in good agreement
               with the nonparametric estimates, this is evidence of verifying the PH assumption
for for these data")
              # ,plotlyOutput("plot5b", height = "720px"),
               #h5(textOutput("info"))
             ))),   
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      tabItem("RESULTS1",
              fluidRow(        
                #  box(
                #    title = "Data listing"
                #      ,status = "primary"
                #       ,solidHeader = TRUE 
                #       ,collapsible = TRUE 
                #    , DT::dataTableOutput("mytable")
                #)
                
                #,
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                box(width=8,
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
      # new tab
      tabItem("partial",
              fluidRow(
                box(
                  width=7,
                  # background = "green",
                  title = "Can you enter a HR that maximises the Log Likelihood? Hint: try the model's estimate of the population HR"
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  , DT::dataTableOutput("exercise")
                )
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ,box(
                  width=5,
                  title='Small sample of N=10 exemplifying partial likelihood calculations'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  , DT::dataTableOutput("exercise2")
                  ,p("")
                  ,p("
                In this small dataset we can look at the calculations in detail.
                Above time is sorted. The individual likelihood is overridden 
                if a patient is censored. As we are not working on the log scale 
                here the indivdual likelihoods column equals 1 (Log scale this would be 0).
                
                For a non censored patient the numerator is the (maximised) 
                HR x treatment indicator and the denominator is the sum of the 
                numerator columns for patients still in the risk set. 
                For example for the first patient (who is not censored) we sum all numerator 10 values. 
                The individual likelihoods are then Num./Den. Hence the censored individuals are included 
                in the summation over the risk sets at event times that occur before a censored time.
                 
                Because we have not logged the data, multiply all individual likelihoods 
                to give the log likelihood. This should be ok with 10 samples but will 
                cause numerical problems with more patients, so it is advisable to do the 
                calculations on the log scale.
                 
                Note the last person in the study always contributes 1 as they are the only person in their risk set. Therefore there is no 
                effect on the HR estimate whether
                they experience the event or are censored!
                
                Here we only log the final column.")
                ))),        
          
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      tabItem("HELP", 
              fluidRow(
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                box(  
                  title='Proportional hazards model'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE ,
                  textOutput("help"),
                  #      br(),
                  #textOutput("help2"),
                  withMathJax(),  # need this to be stated
                  #  p(strong("To do spieglehalter explanantion/ equations and partial likelihood example/cox paper link")) ,
                  
                  p("The proportional hazard model can be expressed as:"),
                  
                  p("$$\\begin{align}
                     h_{i}  {(t)} = h_{0} {(t)} {exp} ({\\beta_1}{x_{i1}} + \\cdots{\\beta_k}{x_{ik}})
                     \\end{align}$$"),
                  
                  p("where $h_{i}  {(t)}$ is the dependent variable (operationalized as the hazard rate at time t for subject i), ${x_{1}}$  to 
                ${x_{k}}$  are k independent variables or covariates,
                and $\\beta_{1}$ to $\\beta_{k}$ are the regression coefficients; $h_{0} {(t)}$ is a baseline hazard
                function and is left unspecified. The baseline hazard function can be
                thought of as the hazard function for an individual whose covariates all
                have values of 0."),  
                  
                  p("taking logs:"),   
                  
                  p("$$\\begin{align}
                     \\log h_{i}  {(t)} ={\\alpha{(t)}} + {\\beta_1}{x_{i1}} + \\cdots{\\beta_k}{x_{ik}}
                      \\end{align}$$"),
                  
                  p("Two features of the model are worth noting: (a) the model does not
              require assumptions about the underlying distribution of the survival
              times (i.e., no matter what the actual form of the survival distribution
              is - exponential, Weibull, Gompertz, standard gamma, generalized
              gamma, log-normal, or log-logistic - one can run the same Cox
              regression model for all) and (b) the model assumes a constant ratio of
              the hazards for any two individuals."),
                  
                  p("The second feature gives the model its name: proportional hazards
              model. Because there is no requirement for understanding the underlying survival distribution and because of the proportional hazards
              assumption, the model is also known as a semiparametric model. The second feature gives the model its name: proportional hazards
              model. Because there is no requirement for understanding the underlying survival distribution, and because of the proportional hazards
              assumption, the model is also known as a semiparametric model"),
                  
                  p("The proportional hazards assumption is that the hazard for any individual in a sample is a 
               fixed proportion of the hazard for any other individual,
              and the ratio of the two hazards is constant over time. Precisely, it means
              that in a log(hazard) plot, the log(hazard) curves for any two individuals
              should be strictly parallel. What is important here is that with this assumption, $h_{0} {(t)}$, 
              the baseline hazard function cancels out from the formula
              expressing a hazard ratio for any two individuals i and j, as follows (we can cross out the $h_{0} {(t)}$:"),
                  
                  p("$$\\begin{align}
                     \\frac{ h_{i} (t)}  {h_{j} (t)} = 
                       \\frac{  
                          h_{0} (t) {exp} ({\\beta_1}{x_{i1}} + \\cdots+{\\beta_k}{x_{ik}})}
                       {  h_{0} (t) {exp} ({\\beta_1}{x_{j1}} + \\cdots+{\\beta_k}{x_{jk}})}
                      \\end{align}$$"),
                  
                  p("rearranging:"),   
                  
                  p("$$\\begin{align}
                           {exp} ({\\beta_1}({x_{i1}}-{x_{j1}} +\\cdots+{\\beta_k}{x_{ik}}-{x_{jk}} ))
                     \\end{align}$$"),
                  
                 )
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ,box(
                  title='Partial likelihood'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE ,
                  # ,plotOutput("plot4", height = "720px")
                  p("Now we develop the partial likelihood function. First sort the data
                in an ascending order by the study time, so that the first subject in our
                 sample has the shortest study time or highest hazard rate $h_{1}$, the second
                 subject has the next shortest study time or second highest hazard rate $h_{2}$,
                 and so on, until the last or the $\\it{n}$th subject who has the longest study time
                 or lowest hazard rate $h_{n}$"),
                    
                  p("$$\\begin{align}
                           {h_{1} (t)} = {h_{0} (t)} exp({\\beta}{x_{1}})
                      \\end{align}$$"),
                  
                  p("where x1 denotes the value of x for subject 1. Likewise, we can have a
                   similar expression of hazard functions for all subjects. Taking the sum of
                   hazards over all sample subjects, we obtain a risk set as:"),
                            
                  p("$$\\begin{align}
                           {h_{1} (t)} + {h_{2} (t)} + {h_{3} (t)} + \\cdots+  {h_{n} (t)}
                   \\end{align}$$"),
                  
                  p("The likelihood for individual 1 to have the event at time t is simply the
                    ratio of hazard over the risk set, or is the hazard for subject 1 at time t
                    divided by the sum of the hazards for all subjects who are at risk of having
                    the event at time t. That is,"),
                  
                  p("$$\\begin{align}
                     {L_{1}} = \\frac{h_{1} (t)}    {  {h_{1} (t)} + {h_{2} (t)} + {h_{3} (t)} + \\cdots+  {h_{n} (t)} }
                   \\end{align}$$"),
                  
                  p("For the 2nd subject, note the first subject is no longer in the risk set"),
                  
                  p("$$\\begin{align}
                    {L_{2}} = \\frac{h_{2} (t)}    {  {h_{2} (t)} + {h_{3} (t)} + {h_{4} (t)} + \\cdots+  {h_{n} (t)} }
                   \\end{align}$$"),
                  
                  p("Substituting and cancelling all $h_{0} {(t)}$ for subject 1"),
                  
                  p("$$\\begin{align}
                   {L_{1}} = \\frac {exp({\\beta}{x_{1}})}    {  
                   exp({\\beta}{x_{1}}) + exp({\\beta}{x_{2}}) + exp({\\beta}{x_{3}}) + \\cdots+  exp({\\beta}{x_{n}}) 
                     }
                   \\end{align}$$"),
                  
                  p("Substituting and cancelling all $h_{0} {(t)}$ for subject 2, note the first subject is no longer in the risk set"),
                  
                  p("$$\\begin{align}
                   {L_{2}} =  \\frac {exp({\\beta}{x_{2}})}    {  
                  exp({\\beta}{x_{2}}) + exp({\\beta}{x_{3}}) + exp({\\beta}{x_{4}}) + \\cdots+  exp({\\beta}{x_{n}}) 
                     }
                 \\end{align}$$"),
                  
                  p("Three points to note (a) the
                    baseline hazard is canceled out; (b) as a result, the likelihood function
                    is solely expressed by $\\beta—the coefficient to be estimated and the
                    predictor; and (c) the model takes the information of censored
                    cases into account when building the likelihood function—censored
                    cases are not excluded, and their information (i.e., the hazard functions)
                    is built into the construction of the risk set"),
                  
                  p("Writing the partial likelihoods for each of the n subjects and multiplying all
                     these partial likelihoods together, we obtain the sample partial likelihood:"),
                  
                  p("$$\\begin{align}
                          {  \\it{PL} = \\prod_{i=1}^n {L_i} =  {L_1} \\times {L_2} \\times \\cdots \\times {L_n} }
                     \\end{align}$$"),

                  p("Each censored subject $\\it{j}$ contributes a likelihood of value 1, or $L_{j}^0=1$"),
                  
                  p("$$\\begin{align}
                 {\\it{PL} = \\prod_{i=1}^n }
                   \\left[ 
                    \\frac{
                 e^{\\beta x}
                       }{
                    \\sum_{j=1}^n {Y_{ij}}
                e^{\\beta x}
                   }
                  \\right]^{\\delta_i}
                  \\end{align}$$"),
                  
                  p("
                   where ${Y_{ij} =1}$ if ${t_{j}} \\gt {t_{i}}$; and ${Y_{ij} =0}$ if ${t_{j}} \\lt {t_{i}}$. 
                   Here ${Y_{ij}}$ serves as a switcher (i.e., on
                  and off), indicates that the study times are rank ordered, and signifies that
                  the estimating algorithm should not use the information for those whose
                  events occurred at a point earlier than the current ith subject in the calculation of the risk set 
                  (i.e., the formation of denominators). This makes sense
                  because those who already had the events have exited the set and are no
                  longer elements of the risk set.  $\\delta_i$ is the censoring code
                  and takes the value of 0 or 1. If $\\delta_i=0$ (i.e., the study time is censored), then
                  the whole partial likelihood for this subject equals value 1 otherwise a non one value. "),
                  
                  p("To avoid numerical problems we work on the log scale and so the partial log likelihood becomes: "),

                  p("$$\\begin{align}
                    {log\\it{PL} = 
                   \\sum_{i=1}^n }  {\\delta_i}
                   \\left[ 
                  {\\beta x_i} - log(
                  \\sum_{j=1}^n {Y_{ij}}
                   e^{\\beta x_j})
                  \\right]^{\\delta_i}
                  \\end{align}$$"),
                  
                  p("With this log partial likelihood function, we can search for the
                    best estimate of $\\beta$. Typically a starting value of $\\beta$ is plugged 
                    into the right-hand side of the equation to obtain a first 'guess’ of log
                    PL. Through several iterations, we should fing that further modification of $\\beta$ is no longer necessary because the difference between the
                    current log PL and the log PL obtained from the previous iteration is
                    less than a predetermined value called the convergence criterion,
                    typically a very small value such as .000001. We stop
                    searching, and the $\\beta$ so obtained is the best that maximizes log PL.
                    Using this $\\beta$, the likelihood of reproducing the sample data is maximum or optimal")
                                      
                ),  # box end
                
              )
              
      ),
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       tabItem("RESULTS3",
              fluidRow(
                
                box(width=6,
                  title = "Plot estimated survival curve from Cox proportional hazards"
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("plot2x", height = "720px")
                  ,h5(textOutput("Staff_name2"))
                ),
                
                box(width=6,
                  title = "Log-log survivor plot; log[-log S(t)] against
                                log time based on Cox PH model"
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("plot3", height = "720px")
                  ,p("
                  The function g(u) = log(-log(u)) is called the complementary log-log transformation, 
                  and has the effect of changing the range from (0,1) for u 
                  to (-$\\infty$ to $\\infty$) for g(u). A plot of g[$S_1$(t)] and g[$S_0$(t)] 
                  versus t or log(t) will yield two parallel curves separated by $\\beta$ if the 
                  proportional hazards assumption is correct.
                  
                  Since the survival functions are less than 1, their logarithms are negative.
                  Thus, we must negate them before we take a second logarithm.

                  ")
                )
 
                )),
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   tabItem("RESULTS4",
           fluidRow(
           
             box(width=6,
                  title='
                Repeated Cox regression coefficients estimates and confidence limits within time intervals'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("plot4", height = "720px")
                 ,p("The log hazard ratios are plotted against the mean failure/censoring time within the interval")
             )
             
             ,box(width=6,
                  title='Plot of how the log HR is estimated to vary over time, based on the Schoenfeld residuals'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("plot2y", height = "720px")
             ))),
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   tabItem("survhaz",
           fluidRow(
             
             box(width=6,
                 title='
               Survival probability'
                 ,status = "primary"
                 ,solidHeader = TRUE 
                 ,collapsible = TRUE 
                 ,plotOutput("survhaz", height = "720px")
                 ,p("The exponential distribution, the simplest survival distribution, has a constant hazard. The
exponential distribution is a special case of the Weibull distribution with shape =1. The shape of the hazard is monotone increasing
for shape > 1 and monotone decreasing for shape < 1.
                 The exponential distribution is easy to work with, but the constant hazard
assumption is not often appropriate for describing the lifetimes. The Weibull distribution offers more flexibility in modeling survival
data. Other
parametric families of survival distributions include the gamma, log-normal, log-logistic, Pareto, and many more.")
             )
             
             ,box(width=6,
                  title='Hazard function'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("survhaz2", height = "720px")
                  ,p("To plot the hazard function with shape and scale, as shown by the left curve
we define the hazard function as the p.d.f. divided by the survival function.")
             ))),
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   tabItem("Change",
           fluidRow(
             
             box(width=6,
                 title='
                KM curve based on user inputs, reference curve in blue'
                 ,status = "primary"
                 ,solidHeader = TRUE 
                 ,collapsible = TRUE 
                 ,plotOutput("plot1c", height = "720px")
                 ,p("KM curve based on user inputs, reference curve in blue")
             )
             
             ,box(width=6,
                  title='Exponential curve based on user inputs, reference curve in blue'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  ,plotOutput("plot1d", height = "720px")
                  ,h5(textOutput("info4"))
                  ,h5(textOutput("info5"))
                  # ,h5(textOutput("info4"))
                  # ,h5(textOutput("info5"))
                 # ,plotOutput("plot2y", height = "720px")
             ))),
   
   
   tabItem("Changeh",
           fluidRow(
             
             box(width=6,
                 title='
                xxxxxxxxxxxxxxxxxxxxxxxx'
                 ,status = "primary"
                 ,solidHeader = TRUE 
                 ,collapsible = TRUE 
                 ,plotOutput("ploth", height = "720px")
                 #,p("KM curve based on user inputs, reference curve in blue")
             )
             
             ,box(width=6,
                  title='Smoothed and step function estimates of the hazard function'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  #,plotOutput("plot1d", height = "720px")
                  #,h5(textOutput("info4"))
                  #,h5(textOutput("info5"))
                  # ,h5(textOutput("info4"))
                  # ,h5(textOutput("info5"))
                  ,plotOutput("ploth1", height = "720px")
             ))),
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   tabItem("Changeh2",
           fluidRow(
             
             box(width=6,
                 title=' Kaplan-Meier and smoothed survival curve estimate'
                 ,status = "primary"
                 ,solidHeader = TRUE 
                 ,collapsible = TRUE 
                 ,plotOutput("ploth2", height = "720px")
                 #,p("KM curve based on user inputs, reference curve in blue")
             )
             
             ,box(width=6,
                  title='xxxxxxxxxxxxxxxxxxxxxxxxxx'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  #,plotOutput("plot1d", height = "720px")
                  #,h5(textOutput("info4"))
                  #,h5(textOutput("info5"))
                  # ,h5(textOutput("info4"))
                  # ,h5(textOutput("info5"))
                #  ,plotOutput("ploth2", height = "720px")
             
             ))),
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   tabItem("power",
           fluidRow(
             
             box(width=6,
                 title='Survival curves based on a Weibull distribution for use in a power simulation'
                 ,status = "primary"
                 ,solidHeader = TRUE 
                 ,collapsible = TRUE 
                 ,plotOutput("powerp1", height = "720px")
                 ,p("The solid curve is for the control group and the dashed curve is
for the intervention group")
                 ,h5(verbatimTextOutput("powerp2"))
             )
             
             ,box(width=6,
                  title='Two-Group Event Time Comparison, one simulated realisation'
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  #,plotOutput("plot1d", height = "720px")
                 
                  #,h5(textOutput("info5"))
                  # ,h5(textOutput("info4"))
                  # ,h5(textOutput("info5"))
                   ,plotlyOutput("powerp3", height = "720px")
                  ,p("Random uniform censoring times. It is assumed that both treatment groups have the same censoring distribution.")
                  ,h5(verbatimTextOutput("powerp5"))
                  ,h5(verbatimTextOutput("powerp4"))
             ))),
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   tabItem("KMTABLE",
           fluidRow(
             box(
               width=6,
               # background = "green",
               title = "KM survival table"
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
               , DT::dataTableOutput("KM")
             )
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ,box(
               width=6,
               title='Cumulative hazard'
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
              , DT::dataTableOutput("CHAZ")
               ,p("")
               ,p("
                xxxxxxxxxxxxxxxxxxxx")
             ))),        
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   tabItem("RESULTS2",
           fluidRow(
             box(width=4,
               title = "Cumulative proportion" 
               ,status = "primary"
               ,solidHeader = TRUE 
               ,collapsible = TRUE 
               ,plotlyOutput("plot99a", height = "720px")
               ,p("This presentation is the the KM plot but with 1-survival probabilty on the y axis.")
             )
             
             ,box(width=4,
               title="Cumulative Hazard" 
               ,status = "primary"
               ,solidHeader = TRUE
               ,collapsible = TRUE
               ,plotlyOutput("plot99b", height = "720px")
              # ,p("We add a red dotted line the slope of which equals the true baseline cumulative hazard in treatment group 0 
               #   and a blue dashed line the slope of which equals the true cumulative hazard in treatment group 1...the slope equals true baseline x true hr")
             
               ,h5(textOutput("info3"))
               ,p("It is advisable to carry out the above graphical procedure before fitting a Cox regression model. 
               The interpretation of cumulative hazards $H(t)$ can be tricky, one way to think of cumulative hazards is the number of events  
               exepected by time $t$ if the event were a repeatable process.
               The Stata manual's interpretation of cumulative hazards has more: 'The cumulative hazard function, $H(t)$, has much more to offer than merely an
                  intermediate calculation to derive a 
survivor function from a hazard function. Hazards are rates, and in that respect they are not unlike the 
RPM-revolutions per minute-of an automobile engine. 
Cumulative hazards are the integral from zero to t of the hazard rates. 
Because an integral is really just a sum, a cumulative hazard is like the total number of revolutions
an automobile's engine makes over a given period. 
We could form the cumulative revolution function by integrating RPM over time. 
If we let a car engine run at a constant 2,000 RPM for 2 minutes, then the cumulative revolution function at time 2
minutes would be 4,000, meaning the engine would have revolved 4,000 times over that
period. Similarly, if a person faced a constant hazard rate of 2,000/minute (a big risk)
for 2 minutes, she would face a total hazard of 4,000. Going back to the car engine, if
we raced the engine at 3,000 RPM for 1 minute and then let it idle at 1,000 for another,
the total number of revolutions would still be 4,000. If our fictional risk taker faced a
hazard of 3,000/minute for 1 minute and then a hazard of 1,000/minute for another,
the total risk would still be 4,000.
Now let's stick with our fictional friend. Whatever the profile of risk, 
if the cumulative hazard is the same over a 2-minute period, then the probability of the event
(presumably death) occurring during that 2-minute period is the same.
Let's understand the units of this measurement of risk. In this, cumulative hazards
are more easily understood than the hazard rates themselves. Remember that 
$S(t) = \\exp( -H(t))$, so our fictional friend has a probability of surviving the 2-minute interval
of $\\exp(-4000) = 0$: our friend is going to die. One may similarly calculate the probability
of survival given other values for the cumulative hazard.
Probabilities, however, are not the best way to think about cumulative hazards.
Another interpretation of the cumulative hazard is that it records the number of times
we would expect (mathematically) to observe failures over a given period, if only the
failure event were repeatable. With our fictional friend, the cumulative hazard of 4,000
over the 2-minute period means that we would expect him to die 4,000 times if, as in
a video game, each time he died we could instantly resurrect him and let him continue
on his risky path.'")
               ) 
             
             ,box(width=4,
               title="Complementary log−log" 
               ,status = "primary"
               ,solidHeader = TRUE
               ,collapsible = TRUE
               ,plotlyOutput("plot99c", height = "720px")  
               ,p("
                  The function g(u) = log(-log(u)) is called the complementary log-log transformation, 
                  and has the effect of changing the range from (0,1) for u 
                  to (-$\\infty$ to $\\infty$) for g(u). A plot of g[$S_1$(t)] and g[$S_0$(t)] 
                  versus t or log(t) will yield two parallel curves separated by $\\beta$ if the 
                  proportional hazards assumption is correct.
                  
                  Since the survival functions are less than 1, their logarithms are negative.
                  Thus, we must negate them before we take a second logarithm.

                  ")
                ) 
             
             
             )  #fluidrow
   ) #tabitem
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
    )
  ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the server functions for the dashboard  
server <- function(input, output) { 
  output$res <- renderText({
    paste("You've selected:", input$tabs)
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # https://stackoverflow.com/questions/55043092/r-shinydashboard-display-sum-of-selected-input-in-a-valuebox
  output$value1 <- renderValueBox({
    
    valueBox(
      value =  tags$p(paste0(formatz0(setUpByName())," / ",formatz0(setUpByNamea()) ," / ",formatz00(setUpByNameb()) ," / ",formatz1(setUpByNamec()) ," / ",formatz2(setUpByNamea()/setUpByNameb()  )    )
                      ,style = "font-size: 100%;")
      ,subtitle = tags$p('N; Events (a); Exposure (b); Median surv.; Hazard (a/b)', style = "font-size: 150%;")
      ,icon = icon("stats",lib='glyphicon')
      ,color = "red" )
    
  })
  
  output$value2 <- renderValueBox({
    
    valueBox(
      value =  tags$p(paste0(formatz0(setUpByName2())," / ",formatz0(setUpByName2a()) ," / ",formatz00(setUpByName2b()) ," / ",formatz1(setUpByName2c()) ," / ",formatz2(setUpByName2a()/setUpByName2b()  )    )
                      ,style = "font-size: 100%;")
      ,subtitle = tags$p('N; Events (a); Exposure (b); Median surv.; Hazard (a/b)', style = "font-size: 150%;")
      ,icon = icon("stats",lib='glyphicon')
      ,color = "teal")
    
  })
  
  output$value3 <- renderValueBox({
    
    valueBox(
      value =  tags$p(paste0(formatz2(setUpByName4())," ( ",formatz2(setUpByName5()),", ",formatz2(setUpByName6())," ) " ," ; ",formatz1(setUpByNameLL()))
                      ,style = "font-size: 100%;")
      ,subtitle = tags$p(paste0("Hazard ratio trt 1 v 0 with 95% conf. ; log likelihood"), style = "font-size: 150%;")
      ,icon = icon("education",lib='glyphicon')
      ,color = "green")
    
  }) 
  
      
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This is where a new sample is instigated and inputs converted to numeric
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  random.sample <- reactive({
    
    foo <- input$resample
    
    n <- as.numeric(input$n)
    
    allocation <-as.numeric(input$allocation)
    
    hr <- as.numeric(input$hr)
    
    baseline <- as.numeric(input$baseline)
    
    ###############################################
    
    base <- as.numeric(input$base)
    cens <- as.numeric(input$cens)
    hr2 <-  as.numeric(input$hr2)
    per <-  as.numeric(input$per)
    per2 <-  as.numeric(input$per2)
    ###############################################
    
    return(list(  
      
      n=n,
      allocation =allocation,
      hr=hr,
      baseline=baseline,
      
      base =  base,
      cens =  cens,
      hr2  =  hr2,
      per=per,
      per2=per2
      
    ))
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  power <- reactive({
    
    foo <- input$resample
    
    ss <- as.numeric(unlist(strsplit(input$ss,",")))
    
    tt <- as.numeric(unlist(strsplit(input$tt,",")))
    
    # Here, let us accrue patients over three years, and
    #follow them for an additional seven years
    
    af <- as.numeric(unlist(strsplit(input$af,",")))  
    
    hrx <- as.numeric(unlist(strsplit(input$hrx,",")))
    
    nonc <- as.numeric(unlist(strsplit(input$t2,",")))
    
    sim <- as.numeric(unlist(strsplit(input$sim,",")))
    
    return(list(  
      
      ss1=ss[1],
      ss2=ss[2],
      prob1=ss[3],
      prob2=ss[4],
      
      nc=tt[1],
      ni=tt[2],
      
      AA=af[1],
      FF=af[2],
       
      hrx=hrx[1],
      
      nonc=nonc[1],
      
      sim=sim[1]
    
    ))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  power1 <- reactive({
    
    sample <- power()
    
    t1=sample$ss1
    t2=sample$ss2
    
    s1=sample$prob1
    s2=sample$prob2
    
    N1=sample$nc
    N2=sample$ni
    
    start=sample$AA
    fini=sample$FF
    
    hrx=sample$hrx
    drop=sample$nonc
    sim=sample$sim
    
    
    library(Hmisc)
    
    Weib.p <- Weibull2(c(s1,s2),c(t1,t2))
    
    rcens <- function(n) runif(n, start, fini)#
    
    ff.dropout <- Quantile2(Weib.p,hratio=function(x) hrx,
                            dropout=function(x) drop)
    
    #plot(ff.dropout)
    
    rcontrol <-      function(n) ff.dropout(n, what='control')
    rintervention <- function(n) ff.dropout(n, what='intervention')
    
    x<-spower(rcontrol, rintervention, rcens, pr=FALSE,
              nc=N1, 
              ni=N2,
              test=logrank, nsim=sim, alpha=0.025, cox=TRUE)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # get one simulation realisation and plot it for information
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ##lifted from spower function
    
    yc <- rcontrol(N1)
    yi <- rintervention(N2)
    cens <- rcens(N1+ N2)
    group <- c(rep(1, N1), rep(2, N2))
    y <- c(yc, yi)
    maxfail <- 0
    maxcens <- 0
    maxfail <- max(maxfail, max(y))
    maxcens <- max(maxcens, max(cens))
    S <- cbind(pmin(y, cens), 1 * (y <= cens))
    nexceed <- 0 + (logrank(S, group) > .025)
    fit <- coxph.fit(as.matrix(group), S, strata = NULL, 
                     offset = NULL, init = NULL, control = coxph.control(iter.max = 10, 
                                                                         eps = 0.0001), method = "efron", 
                     rownames = NULL)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    d <- cbind(S, group)
    d <- data.frame(d)
    names(d) <- c("dt","e", "trt")
    dd <<- datadist(d)
    options(datadist='dd')
    
    fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
    
    
    # survplot(fit, type="kaplan-meier", conf.int = TRUE, 
    #          col.fill=c("firebrick1","cyan2"), grid=TRUE, what='survival')
    
    
    f1 <- survfit(Surv(dt,e) ~ trt, data = d)
   
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return(list(x=x, f=ff.dropout , f1=f1, fit=fit)) 
    
 
  })
  
  
  
  # estimating hazard plot
  output$powerp1 <-renderPlot({     
    
    H=power1()$f
    plot(H)
    
  })
  
  
  output$powerp3 <-renderPlotly({     
    
    f1= power1()$f1
    fit=power1()$fit
    
    p1 <- ggsurvplot(f1, main = "Kaplan-Meier Curve",
                     
                  
                     legend.title = "Trt."
               
                     #,xlab=paste0("Time : HR=",round(exp(fit$coefficients),4))
                     ,xlab= "Time"
                      
    )
    ggplotly(p1[[1]] )
  })
  
 
  
  output$powerp2 <- renderPrint({  # renderText not so useful
    
    sample <- random.sample()
    
    x <- power1()$x
    
       return(print(x, digits=4))
    
  })
  
  
  output$powerp4 <- renderPrint({  # renderText not so useful
    
    sample <- random.sample()
    
    x <- power1()$fit
    
    return(print(summary(x), digits=4))
    
  })
  
  
  
  output$powerp5 <- renderPrint({  # renderText not so useful
    
    sample <- random.sample()
    
    x <- power1()$fit
    
    return(print((x), digits=4))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The change in  hazard tab, data generation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  datc <- reactive({
    
    sample <- random.sample()
    
    n =        sample$n
    lambdaT =  sample$base  # event hazard
    lambdaC =  sample$cens  # censoring haz
    beta1 =    sample$hr2   # hr
    per =      sample$per   # survival probability
    per2 =     sample$per2  # survival probability
    
 #  n=800; lambdaT=14.4; lambdaC=12; beta1=1.2; per=0.7 ;  # use this for on the fly
    
    beta1 <- log(as.numeric(beta1))  # log hr
     
    x1 = sample(0:1, n,replace=TRUE)  # assign trt randomly
    
    # create distributions
    T = rweibull(n, shape=1, scale=1/lambdaT)*exp(-beta1*x1)   # say if lambda is entered as 14 big hazard need to use 1/14 here?
    C = rweibull(n, shape=1, scale=1/lambdaC)                  #censoring time
    time = pmin(T,C)  # observed time is min of censored and true, pmin useful here
    event = time==T   # set to 1 if event is observed
    
    # run cox regression
    require(survival)
    f <- coxph(Surv(time, event)~ x1 , method="breslow")
    survfit <- survfit(Surv(time,event) ~ x1)
   # f
    #plot(survfit, ylab="Survival probability", xlab="Time", col=c('blue','red'))
    #run weibull regression
    w <- survreg(formula = Surv(time, event) ~ x1, dist = "w", control = list(maxiter = 90) )
    
    # grab the hr estimates from weibull
    hr <-  (c(w$coefficient[2],  confint(w)[2,]))
    hr <- exp(hr)      #exp(-coef(f))^exp(coef(f)["shape"])
    # grab the hr estimates from Cox
    hrc <- exp(c(f$coefficient,  confint(f) ) )
    
    # capture for later
    ss <- Surv(time, event)
    
    #  for practice run parametric regression using harrell package
    dd <<- datadist( x1=x1)
    options(datadist='dd')
    f.exp  <- psm(ss  ~ x1, dist ='exponential')
    fw    <-  psm(ss  ~ x1,  dist ='weibull')
 
    d <- table(x1,event)
    ev <- d[,"TRUE"]  # number of events
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # new, plotting the hazard, Dirk F moore page 38
    
    library(muhaz)
    
    # plotH
    result.simple <- muhaz(time, event, max.time=max(time),
                           bw.grid=2.25, bw.method="global", b.cor="none")
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PlotH1
    result.pe5 <- pehaz(time, event, width=5, max.time=max(time))  # in the text these are months, so will need to adjust? 
    result.pe1 <- pehaz(time, event, width=1, max.time=max(time))
    
    result.smooth <- muhaz(time, event, bw.smooth=max(time),
                           b.cor="left", max.time=max(time))
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #One use of smoothing the hazard function is to obtain a smooth estimate of the
    # survival function
    
    haz <- result.smooth$haz.est
    times <- result.smooth$est.grid
    survx <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
    
   # We may compare our
   # smoothed survival estimate to the Kaplan-Meier estimate as follows:
    result.km <- survfit(Surv(time, event) ~ 1,
                           conf.type="none")
    
    return(list(s=survfit, f=f ,w=w , hr=hr, hrc=hrc, ev=ev, d=d, result.simple=result.simple,
                result.pe5= result.pe5,result.pe1=result.pe1,result.smooth=result.smooth,  result.km= result.km,
               time=time, event=event , survx=survx, times=times)) 
    
  })
    
  # right plot....................., 
  output$plot1d<-renderPlot({     
    
    sample <- random.sample()
    
    n =        sample$n
    lambdaT =  sample$base  # event hazard
    lambdaC =  sample$cens  # censoring haz
    beta1 =    sample$hr2   # hr
    per =      sample$per   # survival probability
    per2 =     sample$per2  # survival probability
    
    s <- datc()$s  # survfit object
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    yo <- abs(100*((beta1/1)-1))
    
    wordd <- ifelse(beta1 < 1,"will reduce", 
                    ifelse(beta1 > 1, "will increase","will not change")) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot weibull density
    x <- curve(dweibull(x, shape=1, scale = 1/lambdaT), from=0, to=max(s$time))  
    x$y <- x$y/max(x$y)    #scale the weibull to max of 1 
    plot(x, type = "l", lty = 1, col='blue' , ylab="Survival probability", xlab="Time" , ylim=c(0,1))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.58, y = .99,                # Text with different color & size
         paste0(" For the blue reference curve:"),
         col = "#1b98e0",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.655, y = .97,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per)* 1/lambdaT), " months the survival probability is ",per* 100,"%"),
         col = "#1b98e0",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.66, y = .95,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per2)* 1/lambdaT), " months the survival probability is ",per2* 100,"%"),
         col = "#1b98e0",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    text(max(s$time)*.57, 0.92, expression( paste(
      "Using S(1)t = S(0)t"^{exp(beta.x)}
    )), cex = 1)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.7, y = .88,                # Text with different color & size
         paste0(" Postulating treatment " ,wordd," the hazard by ",yo,"%"),
         col = "red",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.69, y = .86,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per)*  1/lambdaT), " months the survival probability becomes ",formatz00(100*(per)^(beta1)) ,"%"),
         col = "red",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.695, y = .84,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per2)*  1/lambdaT), " months the survival probability becomes ",formatz00(100*(per2)^(beta1)) ,"%"),
         col = "red",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    # with the effect of trt and plot
    y <- x$y^(beta1)   # add another line based on S1(t) = S(0) ^exp(B)
    lines(y~x$x, col='red')     
 
    # add arrows to explain  
    arrows(                                  # x start
      -log(per2) *1/lambdaT,                   # time  
      per2 ,#                                # surv prob  
      -log(per2) *1/lambdaT ,  
      per2^beta1 ,
      col="black", lty=1 )       
    
    arrows(                                  # x start
      -log(per) *1/lambdaT,    
      per, # 
      -log(per) *1/lambdaT ,  
      per^beta1 ,
      col="black", lty=1 )    
    
   
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # left plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$plot1c<-renderPlot({     
    
    sample <- random.sample()
    
    n =        sample$n
    lambdaT =  sample$base  # event hazard
    lambdaC =  sample$cens  # censoring haz
    beta1 =    sample$hr2   # hr
    per =      sample$per   # survival probability
   # per2 =     sample$per2  # survival probability
    
    hr <- datc()$hr     # Weibull hr
    hrc <- datc()$hrc   # cox hr
    
    ev <- datc()$ev    # events
    
    N <- datc()$d      # sample size
    
    f <- datc()$f      # Get the data
    s <- datc()$s      # survfit object
    
    plot(s, ylab="Survival probability", xlab="Time", col=c('blue','red'))  # plot the survfit object
    
    text(x = max(s$time)*.58, y = .99,                # Text with different color & size
         paste0("HR from Weibull=",formatz2(1/hr[1])," 95%CI ( ",formatz2(1/hr[3]),", ",formatz2(1/hr[2])," )"),
         col = "#1b98e0",
         cex = 1)
    text(x = max(s$time)*.58, y = .97,                # Text with different color & size
         paste0("HR from Cox PH=",formatz2(hrc[1])," 95%CI ( ",formatz2(hrc[2]),", ",formatz2(hrc[3])," )"),
         col = "#1b98e0",
         cex = 1)
    text(x = max(s$time)*.54, y = .95,                # Text with different color & size
         paste0( "Actual N= ",sum(N), ", Actual events= ",sum(ev)),
         col = "#1b98e0",
         cex = 1)
    
    
    ## power 
    A <- 0.05 # alpha
    B <- 0.1  # beta
    
    s1=per
    s2=per^(beta1)
    
   # (-log(per)* lambdaT)  # time of s1 s2
    
    f <- ((qnorm(1-A/2))+qnorm(1-B))^2
    d1 <- (4*f)/ (log(beta1)^2)
     
    N2 <- d1/(1-(s1+s2)/2)

    
    text(x = max(s$time)*.667, y = .93,                # Text with different color & size
         paste0( "Alpha 2 sided ",A," Power ",1-B," events required= ",ceiling(d1), ", N= ",ceiling(N2)),
         col = "#1b98e0",
         cex = 1)
    
    
    
    
    
    # abline(h=.5,                lty=2, col='blue')
    # abline(v=-log(.5)*lambdaT,  lty=2, col='blue')
    # abline(h=per,               lty=2, col='red')
    # abline(v=-log(per)*lambdaT, lty=2, col='red')
    
    
    #lines(x=.5, y=-log(.5)*lambdaT , col="pink")
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
 #   segments(0,                0.5, -log(.5)*lambdaT , 0.5 , col="blue", lty=2 ) #h
    # segments( -log(.5)*lambdaT,  0, -log(.5)*lambdaT , 0.5 , col="blue", lty=2 ) #v
    # 
    # segments(0,                per , -log(per)*lambdaT , per , col="blue", lty=2 )  #h
    # segments( -log(per)*lambdaT,  0, -log(per)*lambdaT , per , col="blue", lty=2 )  #v
    # 
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # 
    # segments(0,                              # x start
    #          per^beta1 ,                     # x finish
    #          -log(per^beta1)*lambdaT ,       # y start 
    #          per^beta1  , col="red", lty=2 ) # y finsh
    # 
    # 
    # segments(0,                                   # x start
    #          0.5^beta1,                           # x finish
    #          -log(.5^beta1) *lambdaT ,            # y start
    #          0.5^beta1 ,  col="red", lty=2 )      # y finish
    
    
    # segments(                                  # x start
    #          -log(.5) *lambdaT,    
    #          0.5^beta1, # 
    #           -log(.5) *lambdaT ,  
    #          0.5 ,
    #           col="purple", lty=1 )       
    # 
    # segments(                                  # x start
    #   -log(per) *lambdaT,    
    #   per^beta1, # 
    #   -log(per) *lambdaT ,  
    #   per ,
    #   col="darkgreen", lty=1 )    
    # 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  # estimating hazard plot
  output$ploth <-renderPlot({     
 
    H=datc()$result.simple
    plot(H)
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  output$ploth1 <-renderPlot({     
    
    H=datc()$result.pe5
    H1=datc()$result.pe1
    HS=datc()$result.smooth
    
    
    plot(H,  col="green")
    lines(H1, col='blue')
    lines(HS, col='red')
    
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  output$ploth2 <-renderPlot({     
    
    H2=datc()$result.km
      
    times=datc()$times
    result.smooth=datc()$result.smooth
    
    haz <- result.smooth$haz.est
    times <- result.smooth$est.grid
    surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
    
    plot(H2, conf.int=F, mark="|", xlab="Time",  # confint true throws an error
         #xlim=c(0,30), 
         ylab="Survival probability")
    lines(surv ~ times[1:(length(times) - 1)], col='blue')
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  
  output$survhaz <-renderPlot({     
    
    alpha<- as.numeric(input$shape)
    
    lambda <- as.numeric(input$scale)
    
      weibSurv <- function(t, shape, scale) pweibull(t, shape=shape,
                                                   scale=scale, lower.tail=F)
    curve(weibSurv(x, shape=alpha, scale=1/lambda), from=0, to=2/lambda,
          ylim=c(0,1), ylab='Survival probability', xlab='Time')
    
    text(x = 2/lambda*.85, y = .95,                # Text with different color & size
         paste0(" Shape ", alpha,  ""),
         col = "blue",
         cex = 1.2)
    
    text(x = 2/lambda*.85, y = .9,                # Text with different color & size
         paste0("  Scale ", lambda ,""),
         col = "blue",
         cex = 1.2)
    
    
    text(x = 2/lambda*.85, y = .85,                # Text with different color & size
         paste0(" Mean ", formatz2(gamma(1 + 1/alpha)/lambda),  ""),
         col = "blue",
         cex = 1.2)
    
    
    
    text(x = 2/lambda*.85, y = .8,                # Text with different color & size
         paste0( " Median ",  formatz2((log(2)^(1/alpha))/lambda ),""),
         col = "blue",
         cex = 1.2)
    
    
    
    
    
      })
 
  output$survhaz2 <-renderPlot({     
    
    alpha<- as.numeric(input$shape)
    
    lambda <- as.numeric(input$scale)
    
    weibHaz <- function(x, shape, scale) 
      dweibull(x, shape=shape, scale=scale)/
      pweibull(x, shape=shape, scale=scale,lower.tail=F)
    
    curve(weibHaz(x, shape=alpha, scale=1/lambda), from=0, to=2/lambda,
          ylab='Hazard', xlab='Time', col="red")
    
  })
  
  
  
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text below right plot in changing hazard 
  output$info4 <- renderText({  
    
    sample <- random.sample()
    
    n = sample$n
    
    lambdaT =  sample$base
    lambdaC =  sample$cens
    beta1 =    sample$hr2
    per=       sample$per
    per2=       sample$per2
    
    #s <-  datC()$w
    #hr <- datC()$hr
    
    yo <- abs(100*((beta1/1)-1))
    
    wordd <- ifelse(beta1 < 1,"will reduce", 
                    ifelse(beta1 > 1, "will increase ","will not change")) 
     
    c(paste0(" 
    With a shape parameter of 1 and a scale parameter of ",lambdaT," a survival curve is created based on a Weibull distribution.
    Here we are measuring time in months.
    The time to reach median survival is equal to -log(0.5) x ",lambdaT,". 
             This equates to a median survival of ", formatz2(-log(.5)*lambdaT)," months. 
             
             Replacing 0.5 with a desired survival percentile will return the associated time. 
             Enter survival probabilities in the two boxes on the left. 
             
             The time at which the survival probability is ",per* 100,"% is ",  formatz2(-log(per)* lambdaT), " months.
             The time at which the survival probability is ",per2*100,"% is ",  formatz2(-log(per2)*lambdaT), " months"
      ))
               
     
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  # another piece of text below right plot in changing hazard 
  
  output$info5 <- renderText({  
    
    sample <- random.sample()
    
    n=          sample$n
    
    lambdaT =  sample$base
    lambdaC =  sample$cens
    beta1 =   sample$hr2
    per = sample$per
    per2 = sample$per2
    
    yo <- abs(100*((beta1/1)-1))
    
    wordd <- ifelse(beta1 < 1,"will reduce", 
                    ifelse(beta1 > 1, "will increase ","will not change")) 
     
    
    c(paste0("Now if we are postulating that a new treatment " ,wordd," the hazard by ",yo,"% we can use 
               the fact S1(t) = S0(t)^exp Bx, with exp Bx=", (beta1)," the 
              probability of survival at ",  formatz2(-log(per)*  lambdaT), " months now becomes ",formatz2((per)^(beta1)) ,".
          The probability of survival at ",  formatz2(-log(per2)* lambdaT), " months now becomes ",formatz2((per2)^(beta1)) ,"
          . See the arrows in the plot showing the changes in survival. 
             " ) )
    
    })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
     
   
    # theoretical median, self learning txt p301
 #   -log(.5)*lambdaT   # qweibull(.5,1,15.6)           ## median
  #  -log(.7)*lambdaT                                   ## survival prob 70%
    
    ## surv rates 5 year  survival 70%  
    ## surv rates 10 year survival 50%
    
    ## detect 20% increase in haz rate, what does this inply for 5 n 10 year surv rates
    ## S1(t) = So(t)expBx =0.7^1.2=65.2 
    ###                    0.5^1.2=43.5%
    ## so 5 yr drops to 65.2 and 10 yr drops to 43.5
    
    
    #mfit$surv
  #  survfit <- survfit(Surv(time,event) ~ x1)
  #  plot(survfit, ylab="Survival probability", xlab="Time")
    #survfit
    #summary(survfit)
   # quantile(survfit, probs = c(.3,.348,0.5, .565,0.7), conf.int = FALSE)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # summary(survfit, times=5)
  #  summary(survfit, times=10)
    
    #return(list(  d=res$d, f=res$f, f1=res$f1, sf=res$sf, np=res$np , LL1=res$LL1, LL0=res$LL0, S=res$S, res=res
                  
                     
  #  ))
    
  #})
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
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
    
    return(list(  d=res$d, f=res$f, f1=res$f1, sf=res$sf, np=res$np , LL1=res$LL1, LL0=res$LL0, S=res$S, res=res
                  
                  
                  ,f0a=res$f0a, f0=res$f0, f2=res$f2 #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
                  
                  ))
    
  })
  
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
  output$plot5a <- output$plot1 <- renderPlotly({

    f <- dat()$f1  # Get the  obj

    sf  <- dat()$sf
    X <- as.numeric(as.character(sf[2,c("Effect")]))
    Y <- as.numeric(as.character(sf[2,c("Lower 0.95")]))
    Z <- as.numeric(as.character(sf[2,c("Upper 0.95")]))
    f0 <- dat()$f0

    p1 <- ggsurvplot(f, main = "Kaplan-Meier Curve",
                     
                   #  text = paste("Time: ", time),
                     
                     #palette = c("orange", "purple")
                     legend.title = "Trt."
                    # , xlab= paste0("Time" )
                     ,xlab=paste0("Time : HR=",formatz4(exp(f0)) )
                     # , xlab= paste0("Time: HR = ", formatz2(X),", 95%CI( ",formatz2(Y),", ",formatz2(Z)," )" )
                     #,#conf.int = TRUE,
                     # ggtheme = theme_bw() # Change ggplot2
                  #   ,text = paste0("wt: ", round(wt), "</br></br>qsec: ", round(qsec))
    )
    ggplotly(p1[[1]] )

  })
  
  # p <- ggplot(mtcars, 
  #             aes(x = wt, y = qsec, 
  #                 text = paste0("wt: ", round(wt), "</br></br>qsec: ", round(qsec)))) +
  #   geom_point()
  # 
  # ggplotly(p, tooltip = "text")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$plot99a <- renderPlotly({
    
    f <- dat()$f1  # Get the survfit obj
 
  p1 <-  ggsurvplot(f, fun = "event",   main = "Cumulative proportion",  
                    legend.title = "Trt.")#
                  #  palette = c("orange", "purple"))
  ggplotly(p1[[1]])
   
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot99b <- renderPlotly({
    
    f <- dat()$f1  # Get the survfit obj
   
    sample <- random.sample()
    
    hr=        sample$hr
    baseline=  sample$baseline
    
    p1 <- ggsurvplot(f, fun = "cumhaz",  main = "Cumulative Hazard"  ,
                    legend.title = "Trt.")  
    
    p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline,    linetype="dotted", col='red') 
    p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline*hr, linetype="dotted", col='blue')

   # p1$plot <- p1$plot + scale_y_continuous(trans = 'log')
    
  
    
    ggplotly(p1[[1]])
    
  })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot99c <- renderPlotly({
      
      f <- dat()$f1  # Get the survfit obj
      
      p1 <- ggsurvplot(f, fun = "cloglog",  
                       main = "Complementary log−log" ,
                       legend.title = "Trt.")
                      # legend.labs = c("0", "1")) 
                      #palette = c("jco"))
      ggplotly(p1[[1]])
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$plot2<-renderPlot({     
    
    f <- dat()$np  # Get the  data
    survdiffplot(f, col='darkgreen' , xlab= "Time")
    
    # survplot(f, conf='diffbands',col='purple',cex.aehaz=5,
    #         col.fill='blue'
    #                      , aehaz=TRUE, #times= c(5), 
    #      label.curves=list(method="arrow", cex=2), 
    #  label.curves=list(keys=1:2, cex=2),
    #      abbrev.label=TRUE, levels.only = FALSE)
    
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   output$plot2x <-renderPlot({     # Cox
    
    d <- dat()$d
    
    dd <<- datadist(d)
    options(datadist='dd')
    
    fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
    
    survplot(fit, type="kaplan-meier", conf.int = TRUE, 
             col.fill=c("firebrick1","cyan2"), grid=TRUE, what='survival')
    
  })
  
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   output$plot4<-renderPlot({     
     
     d <- dat()$d  # Get the  obj
     #S <- dat()$S
     
     S <- Surv(d$dt, d$e)
     
     hazard.ratio.plot(d$trt, S, e=20, legendloc='ll', xlab='Time', antilog=FALSE, col='blue', smooth=TRUE,
                       ylim=c(-4,4), ylab=c("Log hazard ratio"))
     
   })

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   output$plot2y <-renderPlot({     # Cox
     
     d <- dat()$d
     
     fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
     
     plot(cox.zph(fit, transform="identity" ), ylim=c(-4,4), ylab=c("Log hazard ratio"))
     
   })
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   output$plot3 <- renderPlot({
    
    d <- dat()$d  # Get the  obj
    S <- dat()$S
    f <- dat()$f
    
    #plot(survfit(S~ d$trt), col=c("purple", "orange"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)" , lwd=3)
    survplot(f,logt=TRUE, loglog=TRUE, 
             col=c("red", "lightblue")
             
           # col=c("orange", "purple")
             )   #Check for Weibull-ness (linearity)
    
  })

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$plot5b <- renderPlotly({
    
    fx <-  dat()$f2 # Get the  obj #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    f0a <- dat()$f0a              #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    p1 <- ggsurvplot(fx, main = "Kaplan-Meier Curve", legend.title = "Trt.",
                    # palette = c("orange", "purple")  ,
                     xlab=paste0("Time : HR=", formatz4(exp(f0a)))
                     # ggtheme = theme_bw() # Change ggplot2 theme
    )
    ggplotly(p1[[1]])
    
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
  
  #~~~~~~~~~~~~~~~~~~~~~~KM table
  output$KM <- DT::renderDataTable({
    
    library(survminer)
    
    require(ggfortify)
    d=dat()$d
    
    KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
    
    KM <- fortify(KM_fit) # fortify(KM_fit, fun = "cumhaz")'
    
    DT::datatable(KM, rownames=FALSE,
                  plugins = 'natural',
               #   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
               #              'Treatment'='trt'),
                  
                  options = list(
                    #  dom = 't',
                    columnDefs = list(list(type = 'natural', targets = c(1,2)))
                  )
    ) %>%
      
      formatRound(
        columns= c("time" ,"surv","std.err","upper","lower"), 
        digits=4 )
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$CHAZ <- DT::renderDataTable({
    
    library(survminer)
    
    require(ggfortify)
    d=dat()$d
    
    KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
    
    KM <-  fortify(KM_fit, fun = "cumhaz")
    
    DT::datatable(KM, rownames=FALSE,
                  plugins = 'natural',
                  #   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
                  #              'Treatment'='trt'),
                  
                  options = list(
                    #  dom = 't',
                    columnDefs = list(list(type = 'natural', targets = c(1,2)))
                  )
    ) %>%
      
      formatRound(
        columns= c("time" ,"surv","std.err","upper","lower"), 
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
                         'Individual likelihoods' ='likelihoodi',
                         'logHR x trt'='part2',
                         'Log(exp HR) of each risk set'='part3',
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
  
  output$exercise <- DT::renderDataTable({
    
    dummy = dat()$res
    
    g <- log(as.numeric(input$guess))
    
    foo <- loglike2(g, dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
    foo$time =NULL
    foo$expB = NULL
    
    datatable(foo, rownames=FALSE,
              plugins = 'natural',
              colnames=c('Time' = 'dt', 
                         'Event or censored' = 'e', 
                         'Treat.'='trt',
                         'Model Hazard Ratio'='hr',
                         'Null Log Likelihood'='lognull',
                         'Maximised Log Likelihood'='lognmax',
                         'HR guess' ='guess',
                         'Individual likelihoods' ='likelihoodi',
                         'Sum the Individual likelihoods to give log likelihood based on guess' ='L'
              ),
              
              options = list(
                #  dom = 't',
                columnDefs = list(list(type = 'natural', targets = c(1,2)))
              )
    ) %>%
      
      formatRound(
        columns= c("Time","Model Hazard Ratio",  
                   'Individual likelihoods'
        ),
        digits=4 ) %>%
      formatRound(
        columns= c( 'Null Log Likelihood',
                    'Maximised Log Likelihood',
                    'Sum the Individual likelihoods to give log likelihood based on guess'), 
        digits=0 )
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$exercise2 <- DT::renderDataTable({
    
    sample <- random.sample()

    allocation=sample$allocation
    hr=        sample$hr
    baseline=  sample$baseline
    
    res <- coxdata(n=10, allocation=allocation, hr=hr, baseline =baseline)
    
    d <- res$d
    d <- plyr::arrange(d, dt)  # sort by time
    
    # Calculate Li for everyone
    d$Numerator   <- exp(res$f$coefficients[[1]] * d$trt)
    d$Denominator <- (rev(cumsum(rev(d$Numerator))))
    d$Li          <- d$Numerator/d$Denominator
    
    # all censored contribute 1 (on multiplicative scale)
    d$Li<- ifelse(d$e %in% 1,d$Li,1)
    
    # get the product of all and log answer
    d$LL <- log(prod(d$Li))  
    d
    # model LL, prove we have ecalc correctly
    res$f$loglik
    
    
    datatable(d, rownames=FALSE,
              plugins = 'natural',
              colnames=c('Time' = 'dt', 
                         'Event or censored' = 'e', 
                         'Treat.'='trt',
                         'Num.'='Numerator',
                         'Den.'='Denominator',
                         'Individual likelihoods'='Li',
                         'log of product of Individual likelihoods to give maximizes log likelihood' ='LL'

              ),
              
              options = list(
                #  dom = 't',
                columnDefs = list(list(type = 'natural', targets = c(1,2)))
              )
    ) %>%
      
      formatRound(
        columns= c("Time","Num.","Den.", 
                   'Individual likelihoods', 'log of product of Individual likelihoods to give maximizes log likelihood'
        ),
        digits=4 ) %>%
         formatRound(
         columns= c( 'Event or censored',
                    'Treat.'), 
        digits=0 )
  })
  
  output$help <- renderText({
    HTML(" ")
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$Staff_name2 <- output$Staff_name <- renderText({  
    
    sf  <- dat()$sf
    X <- as.numeric(as.character(sf[2,c("Effect")]))
    Y <- as.numeric(as.character(sf[2,c("Lower 0.95")]))
    Z <- as.numeric(as.character(sf[2,c("Upper 0.95")]))
    
    Xp  <- X/(X+1)
    Yp  <- Y/(Y+1)
    Zp  <- Z/(Z+1)
    
    wordup <- ifelse(X>1,"higher", "")
    
    wordup2 <- ifelse(X>1,"increase", "reduction")
    
    paste0( "From the Cox model the estimated hazard ratio is "
            , formatz2(X),", 95%CI ( ",formatz2(Y),", ",formatz2(Z),
            " ) comparing treatment 1 to 0. 
            
             A hazard ratio of  ", formatz2(X)," means that, in each unit of time, someone 
            treated in group 1 has ", formatz00(abs(X/1-1)*100),"% ", wordup ," of the chance of experiencing the event of interest
            in the following unit of time as they would were they taking treatment 0.
            
            There is an estimated ", formatz00(abs(X/1-1)*100),"% ", wordup2 ," in the hazard of the outcome by the Cox model. 
            
            Equivalently, the hazard ratio is equal to the odds that a patient in treatment group 1 experiences the event of interest before a
            a patient in treatment group 0.
            
           Therefore we can reformulate the hazard ratio, possibly more intuitively, as
            the probability that a patient in treatment 
            group 1 experiences the event before a patient in treatment group 0, which is: "
                  , formatz2(Xp),", 95%CI ( ",formatz2(Yp),", ",formatz2(Zp),").        ")
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # frank Harrell rms page 479
   output$FH <- renderPlot({
  
    sample <- random.sample()
     
    d <- dat()$d

    trt <-  d$trt
    e   <-  d$e
    dt  <-  d$dt
    d   <-  d$d
    
    limx <- quantile(dt, prob=.99)
    limx <- max(dt)*.8
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    S <- Surv(dt,e)
    f <- npsurv(S ~ trt)
    
    for (meth in c('exact','breslow','efron')) {
      
      g <- cph(S  ~ trt, method=meth, surv=TRUE, x=TRUE, y=TRUE)
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.exp  <- psm(S  ~ trt, dist ='exponential')
    fw    <-  psm(S  ~ trt,  dist ='weibull')
    phform <- pphsm(fw)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    co <- gray(c(0,.8))
    co <- c("red", "lightblue")
    survplot(f, lty=c(1,1)   , lwd=c(1,3), col=co,           label.curves=FALSE, conf='none')
    survplot(g, lty=c(3,3)   , lwd=c(1,3), col=co, add=TRUE, label.curves=FALSE, conf.type='none')
    
    legend(c(limx,160),c(.38,.99),
      c('Nonparametric estimates', 'Cox-Breslow estimates'),
           lty=c(1,3), bty='n',    cex=1.0) # col=co 
    
    legend(c(limx,160),c(.18,.89), 
           c('Trt 0','Trt 1'), lwd=c(1,3), col=co, bty='n',  cex=1.0)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  })
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$info <- renderText({  
    
    c("Because the model
      depends only on ranks, any transformation of the event
      times that preserves the order will leave the coefficient
      estimates unchanged as seen above.")
    
  })
  
  output$info2 <- renderText({  
    
    c("The regression coefficients of the proportional hazards 
      model are estimated without having to specify the 
      baseline hazard function (distribution-free approach), 
      and the estimates depend only on the
      ranks of the event times, not their numerical values.")
    
  })
  
  output$info3 <- renderText({  
    
    sample <- random.sample()
    
    hr=        sample$hr
    baseline=  sample$baseline
    
    c(paste0("For treatment group 0 we add a red dotted line, 
             the slope of which equals the true baseline cumulative hazard ",baseline
             ," and a blue dotted line, the slope of which equals the true
             cumulative hazard in treatment group 1...the slope being the true baseline x true hr ", 
             baseline*hr,""))
  
  })
  
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

shinyApp(ui, server)