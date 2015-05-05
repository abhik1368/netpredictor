library(shiny)

shinyUI(navbarPage("NetPredictor",
                   tabPanel("Run!!",
                            sidebarLayout(
                                sidebarPanel(
                                    fileInput('file1', 'Drug Target Matrix',
                                              accept=c('text/csv', 
                                                       'text/comma-separated-values,text/plain', 
                                                       '.csv')),
                                    fileInput('file2', 'Target Similarity Matrix',
                                              accept=c('text/csv', 
                                                       'text/comma-separated-values,text/plain', 
                                                       '.csv')),
                                    fileInput('file3', 'Chemical Similarity',
                                              accept=c('text/csv', 
                                                       'text/comma-separated-values,text/plain', 
                                                       '.csv')),
                                    helpText(HTML("<h4>Algorithm</h4>")),
                                    radioButtons("algo", "",
                                                 c("RWR"="rwr",
                                                   "NBI"="nbi",
                                                   "HGVI"="hgvi",
                                                   "NetCombo"="combo"
                                                   )),
                                    HTML("<hr />"),
#                                     conditionalPanel(
#                                         condition = "input.algo == 'rwr'",
#                                         sliderInput("Eta", "eta",
#                                                     min = 0, max = 1, value = 0.99, step = 0.01),
#                                         sliderInput("restart", "restart",
#                                                     min = 0, max = 1, value = 0.7, step = 0.1)
#                                     ),
#                                     conditionalPanel(
#                                         condition = "input.algo == 'nbi'",            
#                                         sliderInput("alpha", "alpha",
#                                                     min = 0, max = 1, value = 0.5, step = 0.1),
#                                         sliderInput("lamda", "lamda",
#                                                     min = 0, max = 1, value = 0.5, step = 0.1)
#                                     ),
#                                     conditionalPanel(
#                                         condition = "input.algo == 'hgvi'",
#                                         helpText(HTML("<div style=\"text-indent: 25px\">This sample dataset contains every four-seam fastball and cutting fastball thrown by Mariano Rivera and Phil Hughes over the 2011 season.</div>"))
#                                        # sliderInput("alpha", "alpha",
#                                     #                min = 0, max = 1, value = 0.4, step = 0.1),
#                                      #   sliderInput("steps", "steps",
#                                      #               min = 0, max = 100, value = 20, step = 10)
#                                         
#                                     ),
#                                     HTML("<hr />"),
                                    submitButton("Run Prediction !")
                                ), mainPanel(
                                    tableOutput('contents')
                                ))

                   ),
                   tabPanel("Summary",
                            verbatimTextOutput("summary")
                   ),
                   tabPanel("About",
                            verbatimTextOutput("About")

                              ) )
                   )