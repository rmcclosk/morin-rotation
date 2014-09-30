library(shiny)

d <- read.table("freqs.dat", header=T)

shinyServer(function(input, output) {

  output$freqPlot <- renderPlot({
    plot(1:100, 1:100)
  })
})
