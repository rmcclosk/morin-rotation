library(shiny)

shinyUI(fluidPage(
  titlePanel("Variant allele fractions through time"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient",
                  "Patient ID:",
                  levels(d$patient),
                  multi=T),
      selectInput("chrom",
                  "Chromosome:",
                  levels(d$chrom),
                  multi=T)
    ),

    mainPanel(
      ggvisOutput("freqPlot")
    )
  )
))
