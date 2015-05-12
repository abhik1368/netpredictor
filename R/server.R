require(shiny)
require(rCharts)
require(datasets)
require(tsne)

shinyServer(function(input, output, session) {
    data <-  reactive({
        inFile <- input$file1
        if(is.null(inFile)) {
            dataframe <- mtcars          
        } else {
            
            dataframe <- read.csv(
                inFile$datapath, 
                sep=",",
                quote='"',
                stringsAsFactors=FALSE
            )}
    })
    
    output$myChart<-renderChart({
        p1<-rPlot(input$x,input$y, data=mtcars,type="point",color=input$color,facet=input$facet)
        p1$addParams(dom="myChart")
        return(p1)
    })
}