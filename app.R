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
                                         menuItem("Kaplan Meier (landing page)",                      tabName = "OVERVIEW",  icon = icon("bar-chart-o"), selected = FALSE),
                                         menuSubItem("KM diagnostics",                 tabName = "RESULTS2",  icon = icon("bar-chart-o")),
                                         menuSubItem("Cox proportional hazards",       tabName = "RESULTS3",  icon = icon("bar-chart-o")),
                                         menuSubItem("Hazard ratio over time",         tabName = "RESULTS4",  icon = icon("bar-chart-o")),
                                         menuSubItem("Partial log likelihood",         tabName = "RESULTS1",  icon = icon("table")),
                                         menuItem("Only ranks of event times needed!", tabName = "OVERVIEW2",  icon = icon("bar-chart-o"), selected = FALSE),
                                         menuItem("Model assumptions", tabName = "OVERVIEW3",  icon = icon("bar-chart-o"), selected = FALSE),
                                         menuSubItem("KM lifetable",                   tabName = "KMTABLE",  icon = icon("list-alt")),
                                         menuItem("Partial likelihood exercise",  startExpanded = FALSE,    icon = icon("table"),
                                                 
                                                  tags$div(
                                                    textInput(inputId="guess", label='Enter a guess at Hazard Ratio (defaulted to null)', width = '90%' , value="1"),
                                                  ),
                                                  
                                                  menuSubItem("Partial log likelihood reveal",  tabName = "partial")
                                         ),
                                         menuSubItem("Explanation",                    tabName = "HELP")
                                ),
                                
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
                  title = "Can you enter a HR that maximises the Log Likelihood? Hint: try thr model HR"
                  ,status = "primary"
                  ,solidHeader = TRUE 
                  ,collapsible = TRUE 
                  , DT::dataTableOutput("exercise")
                )
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ,box(
                  width=5,
                  title='Small example of N=10 exemplifying partial likelihood calculations'
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
    
    p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline, linetype="dotted", col='red') 
    p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline*hr, linetype="dotted", col='blue')

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
    
    paste0( "The estimated hazard ratio is "
            , formatz2(X),", 95%CI ( ",formatz2(Y),", ",formatz2(Z),
            " ) comparing treatment 1 to 0. 
            
             A hazard ratio of  ", formatz2(X)," means that, in each unit of time, someone 
            treated in group 1 has ", formatz00(abs(X/1-1)*100),"% ", wordup ," of the chance of experiencing the event of interest
            in the following unit of time as they would were they taking treatment 0.
            
            There is an estimated ", formatz00(abs(X/1-1)*100),"% ", wordup2 ," in the hazard of the outcome. 
            
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