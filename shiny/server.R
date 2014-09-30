library(shiny)

shinyServer(function(input, output) {

  output$freqPlot <- renderPlot({
    if (length(input$patient) > 0 & length(input$chrom) > 0) {
        pat.d <- subset(d, patient == input$patient & chrom %in% input$chrom)
        pat.d$sample <- factor(pat.d$sample, levels=unique(pat.d$sample))
        nsamples <- length(unique(pat.d$sample))

        plot(as.integer(pat.d$sample), pat.d$vaf, type="p", 
             xlim=c(0.5, length(unique(pat.d$sample))+0.5),
             xaxt="n", xlab="sample", ylab="variant allele fraction")
        y <- sapply(unique(pat.d$sample), function (s) pat.d[pat.d$sample==s, "vaf"])
        axis(side=1, at=1:nsamples, labels=levels(pat.d$sample))
        sapply(1:(ncol(y)-1), function (i) {
           segments(rep(i, nrow(y)), y[,i], rep(i+1, nrow(y)), y[,i+1])
        })
    }
  })

})
