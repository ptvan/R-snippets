library(shiny)
library(DT)
library(igraph)
library(intergraph)
library(ggnetwork)

shinyUI(fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      selectInput("layout_algo", label="layout algorithm",
                  choices = list("fruchtermanreingold",
                                 "kamadakawai",
                                 "spring",
                                 "circle",
                                 "princoord",
                                 "random"
                  )
      ),
      selectInput("score_type", label="score type",
                  choices = list("combined_score",
                                 "textmining",
                                 "experiments",
                                 "homology",
                                 "neighborhood"
                  )
      ),
      
      # conditionalPanel(
      #   condition = "input.score_type == 'combined_score'",
      #   sliderInput("score_threshold", label = "score threshold",
      #             min = 0, max = combinedMax, value = combinedMax/2, step = 50)
      # ),
      # conditionalPanel(
      #   condition = "input.score_type == 'textmining'",
      #   sliderInput("score_threshold", label = "score threshold",
      #             min = 0, max = textMax, value = textMax/2, step = 50)
      # ),
      # conditionalPanel(
      #   condition = "input.score_type == 'experiments'",
      #   sliderInput("score_threshold", label = "score threshold",
      #             min = 0, max = expMax, value = expMax/2, step = 50)
      # ),
      # 
      # conditionalPanel(
      #   condition = "input.score_type == 'homology'",
      #   sliderInput("score_threshold", label = "score threshold",
      #             min = 0, max = hoMax, value = hoMax/2, step = 50)
      # ),
      # conditionalPanel(
      #   condition = "input.score_type == 'neighborhood'",
      #   sliderInput("score_threshold", label = "score threshold",
      #             min = 0, max = neighMax, value = neighMax/2, step = 50)
      # ),
      
      sliderInput("score_threshold", label = "score threshold",
                  min = 0, max = 900, value = 450, step = 50),
      
      downloadButton("downloadPlot", label="save network plot"),
      br(),
      br(),
      downloadButton("downloadGenes", label="save gene list"),
      
      tags$div(id ='placeholder')
      , width =3
    ), # end sidebarPanel
    mainPanel(  
      plotOutput("graphPlot"),
      DTOutput("table")
    ) # end mainPanel
  ) # end sidebarLayout
  
))
