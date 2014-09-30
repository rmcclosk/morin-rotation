library(shiny)

shinyUI(fluidPage(
  titlePanel("Variant allele fractions through time"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient",
                  "Patient ID:",
                  levels(d$patient))
    ),

    mainPanel(
      plotOutput("freqPlot")
    )
  )
))
