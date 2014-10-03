library(shiny)
library(ggvis)

shinyServer(function(input, output) {

    all.plot.data <- reactive({
        pd <- subset(d, patient %in% input$patient & 
                        chrom %in% input$chrom & 
                        depth >= input$depth)

        if (nrow(pd) == 0) {
            names <- colnames(pd)
            pd <- rbind(pd, rep(0, ncol(pd)))
            colnames(pd) <- names
        }
        pd[order(pd$pos, pd$sample.num),]
    })

    plot.data <- reactive({
        pd <- all.plot.data()

        # can't have an empty plot...
        if (nrow(pd) > 1) {
	        if (input$order == "Highest fraction") {
	            aggfun <- sum
	        } else if (input$order == "Total change") {
	            aggfun <- function (x) max(abs(diff(x)))
	        }
	        agg <- aggregate(vaf~patient+pos, pd, aggfun)
	        agg <- agg[order(-agg$vaf),]
	        agg <- head(agg, input$n)
	        pd <- merge(pd, agg, by=c("patient", "pos"), suffixes=c("", ".agg"))
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
        vis <- plot.data %>% 
            ggvis(x=~factor(sample.num), y=~vaf)

        if (input$color == "None") {
            vis <- vis %>% 
                layer_points(key := ~key) %>%
                group_by(pos, patient)  %>%
                layer_paths()
        } else if (input$color == "Patient") {
            vis <- vis %>% 
                layer_points(key := ~key, fill = ~factor(patient)) %>%
                group_by(pos, patient)  %>%
                layer_paths(stroke = ~factor(patient))
        } else if (input$color == "Chromosome") {
            vis <- vis %>% 
                layer_points(key := ~key, fill=~factor(chrom)) %>%
                group_by(pos, patient)  %>%
                layer_paths(stroke = ~factor(chrom))
        }

        vis %>%
            add_tooltip(tooltip, "hover") %>%
            add_axis("x", title="sample") %>%
            add_axis("y", title="variant allele fraction") 
    })

    freqPlot.vis %>% bind_shiny("freqPlot")
})
