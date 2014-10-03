library(shiny)

shinyUI(fluidPage(
  titlePanel("Variant allele fractions through time"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient", "Patient ID:", levels(d$patient), multi=T),
      selectInput("chrom", "Chromosome:", levels(d$chrom), multi=T),
      sliderInput("depth", "Minimum coverage depth:", min=0, max=100, step=1, value=0),
      selectInput("color", "Color by:", c("None", "Patient", "Chromosome")),
      sliderInput("n", "Number to show:", min=0, max=30, step=1, value=30),
      selectInput("order", "Order by:", c("Highest fraction", "Total change"))
    ),

    mainPanel(
      ggvisOutput("freqPlot")
    )
  )
))
