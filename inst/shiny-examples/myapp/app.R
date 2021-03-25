library(grid)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(gghighlight)
library(plotly)
library(matlab)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(MSnbase)
library(pRoloc)

#this app still needs some improvement : some warnings, use more synthetic code, espacially in the map part
#I still need to put my function on this part (implement plotly in Alldatavisu or in another function)
#also maybe print a loading on the app when clustering is running
#it's on my to do list


setStockcol(NULL)                       #set up of color palette for better visualization
setStockcol(paste0(getStockcol(), 70))

#save the reduction and clustering method available from pRoloc package
reducmeth <- c("PCA", "MDS", "kpca", "nipals", "t-SNE", "lda", "umap")
Clusmeth <- c("knn", "svm", "ksvm", "naiveBayes", "rf", "nnet", "perTurbo", "tagm.mcmc.allocation")

library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinybusy)

#set the aspect of the spinner when a graph is loading
options(spinner.color = "#60c73e", spinner.color.background = "000000", spinner.size = 2)

ui <- dashboardPage(
  skin = "green",

  dashboardHeader(title = "Cellular localisation of proteins", titleWidth = 350),

  #the side bar
  dashboardSidebar(
    width = 350,
    fluidRow(fileInput("file", label= h3("Select a File for your analysis (not yet implemented"), accept=NULL)),
    menuItem("Data", tabName = "data", icon = icon("database"),
             selectInput("datapack", label = h3("Select data from the package pRolocdata to visualise"),
                         choices = alldatapRoloc,
                         selected = "tan2009r1")
             ),
    menuItem("Reduction method", tabName = "Rmeth", icon = icon("object-group"),
             selectInput("Redmet", label = h3("Choose your Reduction method to visualise the cellular map"),
                         choices = reducmeth,
                         selected = "PCA"),
             conditionalPanel(condition = "input.Redmet == 'PCA' ",
                              checkboxInput("paret", label = "Visualise the pareto diagramm from the PCA", value = FALSE)
                              ),
             conditionalPanel(condition = "input.Redmet == 't-SNE' | input.Redmet == 'umap' ",
                              numericInput("yseed", label = "Choose the seed for the stochastic calculation", value = 500)
                              ),

             #the axes on which we can see the plot depends on the reduction method
             conditionalPanel(condition = "input.Redmet != 'MDS' & input.Redmet != 'umap' ",
                              numericInput("axe1", label = h4("Choose the first axe to visualise the plot"),
                                           value = 1, min = 1),
                              numericInput("axe2", label = h4("Choose the second axe to visualise the plot"),
                                           value = 2, min = 1))),

    menuItem("Clustering method", tabName = "Cmeth", icon = icon("object-group"),
             selectInput("Cmet", label = h3("Choose your clustering method"),
                         choices = Clusmeth,
                         selected = "kNN")
             )
    ),

  dashboardBody(
    tabsetPanel(type = "tabs",
                #the first panel; for the protein profiles
                tabPanel("Markers Intensity",
                         fluidRow(
                           box(title="Graphic", width=12, status="success", solidHeader=TRUE,
                               actionButton("profbutton", "See protein profile"),
                               checkboxInput("inter_prof", label = "Interactive profile", value = FALSE),
                               conditionalPanel(condition = "input.inter_prof == false",
                                                withSpinner(plotOutput(outputId="int_a"), type =6)),
                               conditionalPanel(condition = "input.inter_prof == true",
                                                withSpinner(plotlyOutput(outputId="int_prof_a"), type =6)),

                               conditionalPanel(condition = "input.One_pr_pf == false",  uiOutput("organ")),
                               checkboxInput("clus_prof",
                                             label = "Visualise the proteins profiles with their new assignation", value = FALSE),
                               checkboxInput("One_pr_pf", label = "Visualise profile from spefics proteins", value = FALSE),
                               conditionalPanel(condition = "input.One_pr_pf == true", uiOutput("prot_for_pf")),

                               #including specific options fro data from Borner and al.
                               conditionalPanel(
                                 condition = "input.datapack == 'alldyn' | input.datapack == 'alldyn_two' | input.datapack == 'alldyn_mean' ",
                                 checkboxInput("dyn_prof",
                                               label = "Visualise profile from one prot between the two conditions", value = FALSE),
                                 conditionalPanel(condition = "input.dyn_prof == true",
                                                  plotOutput("Dyn_profile"), uiOutput("prot_profi")
                                                  )
                                 ),
                               )
                           ),

                         ),

                #second panel, the cellular map
                tabPanel("Cellular map",
                         fluidRow(column(6, textInput("text2",label=h1("Enter a title for the graphic"), value="Cellular map"))),
                         tags$hr(),
                         fluidRow(box(title="Graphic", width=12, status="success", solidHeader=TRUE,
                                      actionButton("mapbutton", "See map"),
                                      checkboxInput("Int", label = "Intercativ graphic", value = TRUE),
                                      checkboxInput("Mean", label = "See mean points", value = FALSE),
                                      tags$hr(),

                                      conditionalPanel(condition = "input.Int == true",
                                                       uiOutput("sel_pr"),
                                                       tags$hr(),
                                                       withSpinner(plotlyOutput(outputId="map"), type = 6)),

                                      conditionalPanel(condition = "input.Int == false",
                                                       withSpinner(plotOutput("mapin"), type = 6)),

                                      #specific conditions for data from Borner and al.
                                      conditionalPanel(
                                        condition = "input.datapack == 'alldyn' | input.datapack == 'alldyn_two' | input.datapack == 'alldyn_mean' ",
                                        #choose the proteins to select when plotly
                                        selectInput("cond_s", label = "Choose the condition when you select proteins on the plots",
                                                    choices = list("C1", "C2", "C3", "E1", "E2", "E3"), selected = "C1"),

                                        radioButtons(inputId = "hili_all", label = "Highlight specifics proteins",
                                                     choices = c("No highlight" = 1,
                                                                 "Highlight all prot from one condition" = 2,
                                                                 "Highlight one prot from all condtions" = 3,
                                                                 "Highlight one prot from one condition" = 4),
                                                     selected = 1)),

                                      #choose the proteins to highlight
                                      conditionalPanel(condition = "input.hili_all == 2",
                                                       selectInput("cond", label = "Choose a condition to highlight",
                                                                   choices = list("C1"=1, "C2"=2, "C3"=3, "E1"=4, "E2"=5, "E3"=6),
                                                                   selected = 1)
                                                       ),
                                      conditionalPanel(condition = "input.hili_all == 3",
                                                       selectizeInput("uniprot",
                                                                   label = "Choose a protein from all condtions to highlight",
                                                                   choices = unique(fData(alldyn)$gene_name), multiple = TRUE,
                                                                   selected = unique(fData(alldyn)$gene_name)[1])
                                                       ),

                                      conditionalPanel(
                                        condition = "input.datapack != 'alldyn' & input.datapack != 'alldyn_two' & input.datapack != 'alldyn_mean' ",
                                        checkboxInput("hili", label = "Highlight a protein", value = FALSE)
                                        ),

                                      conditionalPanel(condition = "input.hili == true | input.hili_all == 4",
                                                       uiOutput("prot")),
                                      tags$hr(),

                                      #the second map : all proteins are clustered
                                      actionButton("clusbutton", "See clustering"),
                                      tags$hr(),

                                      conditionalPanel(condition = "input.Int == true",
                                                       uiOutput("sel_pr2"),
                                                       tags$hr(),
                                                       withSpinner(plotlyOutput("clusmap"), type = 6)),

                                      conditionalPanel(condition = "input.Int == false",
                                                       withSpinner(plotOutput("clusmapin"), type = 6)),

                                      #the pareto diagram from the PCA (if PCA selected)
                                      conditionalPanel(condition = "input.Redmet == 'PCA' & input.paret == true",
                                                       plotOutput("pareto")
                                                       )
                                      )
                                  )
                         )
                )
    )
  )

server <- function(input, output, session){
  #the data
  data_marker <- reactive({get(input$datapack)})

  #update the ui depending on the data that are selected

  observe({
    if (input$datapack == "alldyn" | input$datapack == "alldyn_two" | input$datapack == "alldyn_mean"){
      updateCheckboxInput(session, "hili", value = FALSE)
    }
    else
      updateRadioButtons(session, "hili_all", selected = 1)
  })

  observe({
    if (input$datapack == "alldyn_two" | input$datapack == "alldyn_mean"){
      updateSelectInput(session, "cond",
                        choices = list("C"=1, "E"=2), selected = 1)
    }
    if (input$datapack == "alldyn"){
      updateSelectInput(session, "cond",
                        choices = list("C1"=1, "C2"=2, "C3"=3, "E1"=4, "E2"=5, "E3"=6), selected = 1)
    }

  })

  observe({
    if (input$datapack == "alldyn_two" | input$datapack == "alldyn_mean"){
      updateSelectInput(session, "cond_s",
                        choices = list("C", "E"), selected = "C")
    }
    if (input$datapack == "alldyn"){
      updateSelectInput(session, "cond_s",
                        choices = list("C1", "C2", "C3", "E1", "E2", "E3"), selected = "C1")
    }

  })


  output$organ <- renderUI({
      selectizeInput("organe", label = "Choose an organelle to visalise the proteins profiles",
                     choices = as.character(unique(fData(data_marker())$markers))[order(as.character(unique(fData(data_marker())$markers)))],
                     options = list(maxOptions = 15000),
                     multiple = TRUE,
                     selected = as.character(unique(fData(data_marker())$markers))[order(as.character(unique(fData(data_marker())$markers)))][1])

  })

  output$prot_for_pf <- renderUI({
    selectizeInput("Prot_for_pf", label = "Choose a protein to visalise",
                   choices = rownames(fData(data_marker())),
                   options = list(maxOptions = 15000),
                   multiple = TRUE,
                   selected = rownames(fData(data_marker()))[1])

  })



  output$prot <- renderUI({
    if (input$hili == TRUE | input$hili_all == 4){
      selectizeInput("prot_1", label = "Choose a protein to highlight",
                     choices = rownames(fData(data_marker())),
                     selected = rownames(fData(data_marker()))[1],
                     multiple = TRUE,
                     options = list(maxOptions = 15000))
    }
    })

  output$prot_profi <- renderUI({
    if (input$datapack == "alldyn" | input$datapack == "alldyn_two" | input$datapack == "alldyn_mean"){
      selectizeInput("Prot_Profi", label = "Choose a protein",
                     choices = unique(fData(data_marker())$gene_name),
                     options = list(maxOptions = 15000),
                     selected = unique(fData(data_marker())$gene_name))
    }
    })

  #Plot the protein profile, interactive or not, with the ProfileInteract function
  int <- reactive({
    if (!input$inter_prof){
        if (input$clus_prof){
          if (input$One_pr_pf){
            ProfileInteract(data_markerfc(), mrk = input$Cmet, Organelle = input$organe, Interact = FALSE,
                            Clust = TRUE, one_pr = TRUE, protein = input$Prot_for_pf, mytitle = TRUE, TITLE = input$datapack)
          }

          else {
            ProfileInteract(data_markerfc(), mrk = input$Cmet, Organelle = input$organe, Interact = FALSE,
                            Clust = TRUE, mytitle = TRUE, TITLE = input$datapack)
          }

        }
        else {
          if (input$One_pr_pf){
            ProfileInteract(data_marker(), Organelle = input$organe, Interact = FALSE,
                            one_pr = TRUE, protein = input$Prot_for_pf, mytitle = TRUE, TITLE = input$datapack)
          }
          else {
            ProfileInteract(data_marker(), Organelle = input$organe, Interact = FALSE,
                            mytitle = TRUE, TITLE = input$datapack)
          }
        }

      }


    else
      NULL

    })


  int_prof <- reactive({

    if (input$inter_prof){
      if (input$clus_prof){
        if (input$One_pr_pf){
          ProfileInteract(data_markerfc(), mrk = input$Cmet, Organelle = input$organe, Clust = TRUE,
                          one_pr = TRUE, protein = input$Prot_for_pf, mytitle = TRUE, TITLE = input$datapack)
        }

        else {
          ProfileInteract(data_markerfc(), mrk = input$Cmet, Organelle = input$organe, Clust = TRUE,
                          mytitle = TRUE, TITLE = input$datapack)
        }
      }
      else {
        if (input$One_pr_pf){
          ProfileInteract(data_marker(), Organelle = input$organe,
                          one_pr = TRUE, protein = input$Prot_for_pf, mytitle = TRUE, TITLE = input$datapack)
        }

        else {
          ProfileInteract(data_marker(), Organelle = input$organe, mytitle = TRUE, TITLE = input$datapack)
        }

      }
    }

    else
      NULL
    })


  #this is for using the button
  prof1 <- reactiveValues(
    ch = NULL
  )
  prof2 <- reactiveValues(
    ch = NULL
  )

  observeEvent(input$profbutton, {
    if (input$inter_prof == FALSE){
      prof1$ch <- int()
    }
    else {
      prof1$ch <- NULL
    }
  }, ignoreInit = TRUE, ignoreNULL = FALSE)

  observeEvent(input$profbutton, {
    if (input$inter_prof == TRUE){
      prof2$ch <- int_prof()
    }
    else {
      prof2$ch <- NULL
    }
  }, ignoreInit = TRUE, ignoreNULL = FALSE)

  output$int_a <- renderPlot({
    if (input$inter_prof == FALSE){
      prof1$ch
    }
    else
      NULL
  })
  output$int_prof_a <- renderPlotly({
    if (input$inter_prof == TRUE){
      prof2$ch
    }
    else
      NULL
  })


  #plot the box plot when data from Borner and al. are selected, with BoxProfile function
  output$Dyn_profile <- renderPlot({
    if (input$dyn_prof == TRUE){
      BoxProfile(data_marker(), prot = input$Prot_Profi)$graph
    }
    else
      NULL
  })

  #update the ui depending on the reduction method that is selected
  observe({
      if (input$Redmet != "t-SNE"){
        updateNumericInput(session, "axe1", max = ncol(exprs(data_marker())))
        updateNumericInput(session, "axe2", max = ncol(exprs(data_marker())))
        }
      if (input$Redmet == "t-SNE"){
        updateNumericInput(session, "axe1", max = 3)
        updateNumericInput(session, "axe2", max = 3)
        }
      })


  #clustered the data, depending on the clustering method that is selected; method from pRoloc package
  param <- reactive({
    if (input$Cmet == "svm")
      param <- svmOptimisation(data_marker(), times = 3)
    if (input$Cmet == "knn")
      param <- knnOptimisation(data_marker(), k = seq(11, 30, 2), times = 3)
    if (input$Cmet == "ksvm")
      param <- ksvmOptimisation(data_marker(),  times = 3)
    if (input$Cmet == "naiveBayes")
      param <- nbOptimisation(data_marker(), times = 3)
    if (input$Cmet == "nnet")
      param <- nnetOptimisation(data_marker(), times = 3)
    if (input$Cmet == "rf")
      param <- rfOptimisation(data_marker(), times = 3)
    if (input$Cmet == "tagm.mcmc.allocation")
      param <- tagmMcmcTrain(data_marker(), numIter = 2000, burnin = 200,
                             method = "MCMC", thin = 10, numChains = 4)
    if (input$Cmet == "perTurbo")
      param <- perTurboOptimisation(data_marker(), pRegul = 2^seq(-2,2,2), sigma = 10^seq(-1,1,1),
                                    inv = "Inversion Cholesky",
                                    reg = "tikhonov", times = 3)
    param
    })


  data_markerfc <- reactive({
    if (input$Cmet == "svm")
      data_markerfc <- svmClassification(data_marker(), param())
    if (input$Cmet == "knn")
      data_markerfc <- knnClassification(data_marker(), param())
    if (input$Cmet == "ksvm")
      data_markerfc <- ksvmClassification(data_marker(), param())
    if (input$Cmet == "naiveBayes")
      data_markerfc <- nbClassification(data_marker(), param())
    if (input$Cmet == "nnet")
      data_markerfc <- nnetClassification(data_marker(), param())
    if (input$Cmet == "rf")
      data_markerfc <- rfClassification(data_marker(), param(), mtry = c(2,5,10))
    if (input$Cmet == "tagm.mcmc.allocation")
      data_markerfc <- tagmMcmcPredict(data_marker(), params = tagmMcmcProcess(param()))
    if (input$Cmet == "perTurbo")
      data_markerfc <- perTurboClassification(data_marker(), param())

    data_markerfc
  })

  #saving the clustering score
  ptsze <- reactive({

    if (input$Cmet == "svm")
      ptsze <- exp(fData(data_markerfc())$svm.scores) - 1
    if (input$Cmet == "knn")
      ptsze <- exp(fData(data_markerfc())$knn.scores) - 1
    if (input$Cmet == "ksvm")
      ptsze <- exp(fData(data_markerfc())$ksvm.scores) - 1
    if (input$Cmet == "naiveBayes")
      ptsze <- exp(fData(data_markerfc())$naiveBayes.scores) - 1
    if (input$Cmet == "nnet")
      ptsze <- exp(fData(data_markerfc())$nnet.scores) - 1
    if (input$Cmet == "rf")
      ptsze <- exp(fData(data_markerfc())$rf.scores) - 1
    if (input$Cmet == "tagm.mcmc.allocation")
      ptsze <- exp(fData(data_markerfc())$tagm.mcmc.probability) - 1
    if (input$Cmet == "perTurbo")
      ptsze <- exp(fData(data_markerfc())$perTurbo.scores) - 1


    #t-SNE only take unique features -> update the data
    if (input$Redmet == "t-SNE"){
      if (length(unique(data_marker()@assayData$exprs)[,1]) == length(fData(data_marker())$markers))
        ptsze <- ptsze
      else
        ptsze <- ptsze[-which(duplicated(data_marker()@assayData$exprs) == TRUE)]

      }

    ptsze
    })

  #saving the proteins assignations
  dataclus <- reactive({
    if (input$Cmet == "svm")
      dataclus <- fData(data_markerfc())$svm
    if (input$Cmet == "knn")
      dataclus <- fData(data_markerfc())$knn
    if (input$Cmet == "ksvm")
      dataclus <- fData(data_markerfc())$ksvm
    if (input$Cmet == "naiveBayes")
      dataclus <- fData(data_markerfc())$naiveBayes
    if (input$Cmet == "nnet")
      dataclus <- fData(data_markerfc())$nnet
    if (input$Cmet == "rf")
      dataclus <- fData(data_markerfc())$rf
    if (input$Cmet == "tagm.mcmc.allocation")
      dataclus <- fData(data_markerfc())$tagm.mcmc.allocation
    if (input$Cmet == "perTurbo")
      dataclus <- fData(data_markerfc())$perTurbo

    dataclus
  })

  #saving the mean clustering score
  loss <- reactive({
    if (input$Cmet == "svm")
      loss <- mean(fData(data_markerfc())$svm.scores)
    if (input$Cmet == "knn")
      loss <- mean(fData(data_markerfc())$knn.scores)
    if (input$Cmet == "ksvm")
      loss <- mean(fData(data_markerfc())$ksvm.scores)
    if (input$Cmet == "naiveBayes")
      loss <- mean(fData(data_markerfc())$naiveBayes.scores)
    if (input$Cmet == "nnet")
      loss <- mean(fData(data_markerfc())$nnet.scores)
    if (input$Cmet == "rf")
      loss <- mean(fData(data_markerfc())$rf.scores)
    if (input$Cmet == "tagm.mcmc.allocation")
      loss <- mean(fData(data_markerfc())$tagm.mcmc.probability)
    if (input$Cmet == "perTurbo")
      loss <- mean(fData(data_markerfc())$perTurbo.scores)

    loss
  })

  #the plot, interactive or not, clustered or not (reactiveValues in order to use an action button)
  gra1 <- reactiveValues(
    ch = NULL
  )

  #same as Alldatavisu function but now include the plotly option
  map_Int <- reactive({
    if (input$Int){
      prt <- NULL
      if (input$hili_all == 2){
        cd <- TRUE
      }
      else{
        cd <- FALSE
      }
      if (input$hili_all == 3){
        vc <- TRUE
        prt <- input$uniprot
      }
      else{
        vc <- FALSE
      }
      if (input$hili | input$hili_all == 4){
        pr <- TRUE
        prt <- input$prot_1
      }
      else{
        pr <- FALSE
      }


      g <- AllDatavisuInt(data_marker(), redmet = input$Redmet, cmet = "markers", unknow = TRUE, Interact = TRUE,
                     highpr = pr, highcond = cd, vect = vc, proteins = prt, Mean_point = input$Mean,
                     Title = input$text2, yourseed = input$yseed, ax = c(input$axe1, input$axe2),
                     mysubtitle = TRUE, subtitle = input$datapack)
    }
    else {
       g <- NULL
    }

    g


   })

  gra2 <- reactiveValues(
    ch = NULL
  )


  #the non interactive graph
  map_nIn <- reactive({
    if (input$Int){
      g <- NULL
    }

    else {
      prt <- NULL
      if (input$hili_all == 2){
        cd <- TRUE
      }
      else{
        cd <- FALSE
      }
      if (input$hili_all == 3){
        vc <- TRUE
        prt <- input$uniprot
      }
      else{
        vc <- FALSE
      }
      if (input$hili | input$hili_all == 4){
        pr <- TRUE
        prt <- input$prot_1
      }
      else{
        pr <- FALSE
      }


      g <- AllDatavisuInt(data_marker(), redmet = input$Redmet, cmet = "markers", unknow = TRUE, Interact = FALSE,
                          highpr = pr, highcond = cd, vect = vc, proteins = prt, Mean_point = input$Mean,
                          Title = input$text2, yourseed = input$yseed, ax = c(input$axe1, input$axe2),
                          mysubtitle = TRUE, subtitle = input$datapack)

    }

    g

  })

  #update only when action button is pressed
  observeEvent(input$mapbutton, {
    if (input$Int == TRUE){
      gra1$ch <- map_Int()
    }
    else {
      gra1$ch <- NULL
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$mapbutton, {
    if (input$Int == FALSE){
      gra2$ch <- map_nIn()
    }
    else {
      gra2$ch <- NULL
    }
  }, ignoreInit = TRUE, ignoreNULL = FALSE)

  #same thing as before but now using the clustered data
  output$map <- renderPlotly({
    if (input$Int == TRUE){
      gra1$ch
    }
    else
      NULL
    })
  output$mapin <- renderPlot({
    if (input$Int == FALSE){
      gra2$ch
    }
    else
      NULL
    })


  newgr <- reactiveValues(
    ch = NULL
  )

  truc <- reactive({

    if (input$Int == FALSE){
      g1 <- NULL
    }

    else {
      prt <- NULL
      if (input$hili_all == 2){
        cd <- TRUE
      }
      else{
        cd <- FALSE
      }
      if (input$hili_all == 3){
        vc <- TRUE
        prt <- input$uniprot
      }
      else{
        vc <- FALSE
      }
      if (input$hili | input$hili_all == 4){
        pr <- TRUE
        prt <- input$prot_1
      }
      else{
        pr <- FALSE
      }


      g1 <- AllDatavisuInt(data_markerfc(), redmet = input$Redmet, cmet = input$Cmet, unknow = TRUE, Interact = TRUE,
                          highpr = pr, highcond = cd, vect = vc, proteins = prt, Mean_point = input$Mean,
                          Title = input$text2, yourseed = input$yseed, ax = c(input$axe1, input$axe2), Source = "BB",
                          mysubtitle = TRUE, subtitle = input$datapack)
    }

    g1

})

  newgrin <- reactiveValues(
    ch = NULL
  )

  #the non interactive graph
  trucin <- reactive({

    if (input$Int == TRUE){
      g1 <- NULL
    }

    else {
      prt <- NULL
    if (input$hili_all == 2){
      cd <- TRUE
    }
    else{
      cd <- FALSE
    }
    if (input$hili_all == 3){
      vc <- TRUE
      prt <- input$uniprot
    }
    else{
      vc <- FALSE
    }
    if (input$hili | input$hili_all == 4){
      pr <- TRUE
      prt <- input$prot_1
    }
    else{
      pr <- FALSE
    }


    g1 <- AllDatavisuInt(data_markerfc(), redmet = input$Redmet, cmet = input$Cmet, unknow = TRUE, Interact = FALSE,
                        highpr = pr, highcond = cd, vect = vc, proteins = prt, Mean_point = input$Mean,
                        Title = input$text2, yourseed = input$yseed, ax = c(input$axe1, input$axe2), Source = "BB",
                        mysubtitle = TRUE, subtitle = input$datapack)
    }

    g1

  })


  observeEvent(input$clusbutton, {
    newgr$ch <- truc()
  }, ignoreInit = TRUE)

  observeEvent(input$clusbutton, {
    newgrin$ch <- trucin()
  }, ignoreInit = TRUE)


output$clusmap <- renderPlotly({newgr$ch})
output$clusmapin <- renderPlot({newgrin$ch})


  #the pareto diagram from the PCA (plot2D is from pRoloc package)
  output$pareto <- renderPlot({
    if (input$Redmet == "PCA" & input$paret == TRUE){
      plot2D(data_marker(), method = "scree")
    }
    else
      NULL
  })

  #print the proteins selected when plotly is used, depending on certain conditions
  #the first plot (not clustered)
  output$sel_pr <- renderUI({
     select_data <- event_data("plotly_selected", source = "AA")

     if (is.null(select_data)) "Brushed points appear here (double-click to clear)"
     else {
       if (input$datapack == "alldyn" | input$datapack == "alldyn_two" | input$datapack == "alldyn_mean"){
         Sel_data <- select_data$key[which(!is.na(str_extract(select_data$key, paste0("_", input$cond_s))))]
         Sel_data <- str_remove(Sel_data, paste0("_", input$cond_s))
       }
       else {
         Sel_data <- select_data$key
       }

       paste("The proteins slected are :", paste(unique(Sel_data)[order(unique(Sel_data))], collapse = ", "))
     }
     })

  #the second plot (the clustered one)
  output$sel_pr2 <- renderUI({
    select_data <- event_data("plotly_selected", source = "BB")

    if (is.null(select_data)) "Brushed points appear here (double-click to clear)"

    else {
      if (input$datapack == "alldyn" | input$datapack == "alldyn_two" | input$datapack == "alldyn_mean"){
        Sel_data <- select_data$key[which(!is.na(str_extract(select_data$key, paste0("_", input$cond_s))))]
        Sel_data <- str_remove(Sel_data, paste0("_", input$cond_s))
      }
      else {
        Sel_data <- select_data$key
      }
      paste("The proteins slected are :", paste(unique(Sel_data)[order(unique(Sel_data))], collapse = ", "))
    }

  })
}



shinyApp(ui, server)



