shinyUI(
  fluidPage(
    headerPanel("TCGA Analysis"),
    sidebarLayout(
      sidebarPanel(
        selectInput("operation","Select Plot Type",choices = c("Boxplot","Survival")),
        conditionalPanel(condition = "input.operation=='Boxplot'",
                         selectInput("cancerType","Select one or more cancer types",
                                     choices = substr(list.files(data.dir,"Expression/",pattern=".txt"),1,
                                                      nchar(list.files(data.dir,"Expression/",pattern=".txt"))-4),
                                     multiple = TRUE),
                         textInput("gene","Please input one gene of interest",placeholder = "Gene name ie TRIM21")
        ),
        conditionalPanel(condition="input.operation=='Survival'",
                         selectInput("cancerType2","Select a cancer type",
                                     choices = substr(list.files(data.dir,"Expression/",pattern=".txt"),1,
                                                      nchar(list.files(data.dir,"Expression/",pattern=".txt"))-4)),
                         textInput("gene2","Please input one gene of interest",placeholder = "Gene name ie BRCA1"),
                         selectInput("operation2","Select automatic or manual cutoff", choices = c("Automatic", "Manual"))
                         
        ),
        conditionalPanel(condition="input.operation2=='Manual' & input.operation=='Survival'",
                         sliderInput("cut.off","Select a cutoff value",min = .10,max=.50,step = .05,value = .25)
        ),
        actionButton("go","GO")
        
      ),
      mainPanel(
        plotOutput("myPlot"),
        downloadButton('downloadplot',"Download Figure")
      )
      
    )
  )
)

