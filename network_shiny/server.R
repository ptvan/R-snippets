library(shiny)
library(DT)
library(igraph)
library(intergraph)
library(ggnetwork)

allGraph <- readRDS("data/saved_igraph_obj.Rds")
pu <- colorRampPalette(c("purple","mediumorchid","purple4"))(100)                      
ye <- colorRampPalette(c("lemonchiffon", "khaki","yellow"))(100)

shinyServer(function(input, output) {
   
  output$downloadPlot <- downloadHandler(
    filename = function() { paste('My_network_', input$score_type, '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = plotInput(), device = "png")
    }
  ) 
  
  output$downloadGenes <- downloadHandler(
    filename = function() {
      paste(input$score_type, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(buildTable(), file, row.names = FALSE)
    }
  )
  
  
  buildTable <- reactive({
    grph <- igraph::delete.edges(allGraph, which(edge_attr(allGraph, input$score_type) < input$score_threshold))
    grph <- induced_subgraph(grph, v=which(igraph::degree(g=grph, v=V(grph))>1) )
    if (length(V(grph)) > 0 && length(E(grph) > 2)){
      tab <- data.frame(cbind(vertex_attr(grph)$geneName, vertex_attr(grph)$logFC))
      colnames(tab) <- c("gene", "logFC")
      tab$logFC <- as.numeric(levels(tab$logFC))[tab$logFC]
      tab
      
    }
    
  })
  
  output$table <- renderDT(
    # formatRound(buildTable(), c("logFC"), 3)
    buildTable()
  )
  
  plotInput <- reactive({
    
    
    set.seed(100)
    grph <- igraph::delete.edges(allGraph, which(edge_attr(allGraph, input$score_type) < input$score_threshold))
    grph <- induced_subgraph(grph, v=which(igraph::degree(g=grph, v=V(grph))>1) )
    if (length(V(grph)) > 0 && length(E(grph) > 2)){
      
      out <- ggplot(ggnetwork(grph, layout=eval(input$layout_algo)), aes(x=x,y=y,xend=xend,yend=yend)) +
        geom_edges(aes_string(size = eval(input$score_type)), color = "grey80") +
        geom_nodes(aes(color=logFC, fill=NULL), size = 12) +
        scale_color_gradientn(colours=c(pu,"black", ye), na.value = "grey98", limits = c(-2, 2)) +
        geom_nodetext(aes(label=geneName), size=3) +
        labs(title = paste0(length(V(grph)), " genes in network at ", input$score_type, ">", input$score_threshold)) +
        theme_blank()
      
    } else {
      showModal(modalDialog(title="Empty network", "Your filtering threshold has excluded all genes!"))
    }
    
  })
  
  output$graphPlot <- renderPlot({
    print(plotInput())
    
  }
  ) # end renderPlot
  
  
})
