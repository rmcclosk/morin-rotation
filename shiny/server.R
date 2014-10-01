library(shiny)

shinyServer(function(input, output) {

    plot.data <- reactive({
        pd <- subset(d, patient == input$patient & chrom == input$chrom)
        pd$sample <- factor(pd$sample, levels=unique(pd$sample))
        pd$vaf <- as.numeric(pd$vaf)
        pd
    })

    output$freqPlot <- renderPlot({
        pd <- plot.data()
        plot(as.integer(pd$sample), pd$vaf, type="p", 
             ylab="variant allele fraction",
             xlab="sample",
             xaxt="n",
             xlim=c(0.5, length(levels(pd$sample))+0.5))
        axis(1, at=1:length(levels(pd$sample)), labels=levels(pd$sample))
    })

    output$chromPlot <- renderPlot({
        plot.new()
    })

#    output$freqPlot <- renderPlot({
#        nsamples <- length(unique(pat.d$sample))
#        do.plot <- (length(input$patient) > 0 & length(input$chrom) > 0)
#
#        if (do.plot) {
#            plot(as.integer(pat.d$sample), pat.d$vaf, type="p", 
#                 xlim=c(0.5, length(unique(pat.d$sample))+0.5),
#                 xaxt="n", xlab="sample", ylab="variant allele fraction")
#            y <- sapply(unique(pat.d$sample), function (s) as.numeric(pat.d[pat.d$sample==s, "vaf"]))
#            print(y)
#            axis(side=1, at=1:nsamples, labels=levels(pat.d$sample))
#            sapply(1:(ncol(y)-1), function (i) {
#                segments(rep(i, nrow(y)), y[,i], rep(i+1, nrow(y)), y[,i+1])
#            })
#        }
#    })
#
#    output$chromPlot <- renderPlot({
#        pat.d <- subset(d, patient == input$patient & chrom %in% input$chrom)
#        pat.d$sample <- factor(pat.d$sample, levels=unique(pat.d$sample))
#        nsamples <- length(unique(pat.d$sample))
#        do.plot <- (length(input$patient) > 0 & length(input$chrom) > 0)
#
#        if (do.plot) {
#            x <- sapply(unique(pat.d$sample), function (s) as.numeric(pat.d[pat.d$sample==s, "pos"]))
#            y <- sapply(unique(pat.d$sample), function (s) as.numeric(pat.d[pat.d$sample==s, "vaf"]))
#            sapply(1:nsamples, function (i) {
#                plot(x[,i], y[,i], type="p", col=i)
#            })
#        }
#    })

})
