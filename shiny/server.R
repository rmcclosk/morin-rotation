library(shiny)
library(ggvis)

shinyServer(function(input, output) {

    plot.data <- reactive({
        pd <- subset(d, patient %in% input$patient & chrom %in% input$chrom)

        # can't have an empty plot...
        if (nrow(pd) == 0) {
            names <- colnames(pd)
            pd <- rbind(pd, rep(0, ncol(pd)))
            colnames(pd) <- names
        }
        pd[order(pd$pos, pd$sample.num),]
    })

    tooltip <- function (x) {
        pd <- isolate(plot.data())
        point <- pd[pd$key == x$key,]
        html <- paste0("<b>Patient: </b>", point$patient, "</b><br />")
        html <- paste0(html, "<b>Sample: </b>", point$sample, "</b><br />")
        html <- paste0(html, "<b>Chromosome: </b>", point$chrom, "</b><br />")
        html <- paste0(html, "<b>Position: </b>", point$pos, "</b><br />")
        html <- paste0(html, "<b>Reference: </b>", point$ref, "</b><br />")
        html <- paste0(html, "<b>Variant: </b>", point$alt, "</b><br />")
    }

    freqPlot.vis <- reactive({
        plot.data %>%
            ggvis(x=~factor(sample.num), y=~vaf) %>%
            layer_points(key := ~key) %>%
            group_by(pos, patient) %>%
            layer_paths() %>%
            add_tooltip(tooltip, "hover") %>%
            add_axis("x", title="sample") %>%
            add_axis("y", title="variant allele fraction") 
    })

    freqPlot.vis %>% bind_shiny("freqPlot")

})
