library(shiny)

shinyUI(fluidPage(
  titlePanel("Variant allele fractions through time"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient",
                  "Patient ID:",
                  levels(d$patient),
                  multi=F),
      selectInput("chrom",
                  "Chromosome:",
                  levels(d$chrom),
                  multi=T)
    ),

    mainPanel(
      plotOutput("freqPlot")
    )
  )
))
