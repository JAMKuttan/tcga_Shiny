shinyUI(
  pageWithSidebar(
    headerPanel("TGCA Analysis"),
    sidebarPanel(
      selectInput("operation","Select Plot Type",choices = c("Pancan","Survival")),
      conditionalPanel(condition = "input.operation=='Pancan'",
                       selectInput("cancerType","Select one or more cancer types",
                                   choices = substr(list.files("../Expression/",pattern=".txt"),1,
                                                    nchar(list.files("../Expression/",pattern=".txt"))-4),
                                   multiple = TRUE),
                       textInput("gene","Please input one gene of interest",placeholder = "Gene name ie 'TRIM21'")
                       ),
      conditionalPanel(condition="input.operation=='Survival'",
                       selectInput("cancerType2","Select a cancer type",
                                   choices = substr(list.files("../Expression/",pattern=".txt"),1,
                                                    nchar(list.files("../Expression/",pattern=".txt"))-4)),
                       textInput("gene2","Please input one gene of interest",placeholder = "Gene name ie 'TRIM21'",value="TRIM21"),
                       sliderInput("cut.off","Select a cutoff value",min = .10,max=.50,step = .05,value = .25)
      ),
      actionButton("go","GO")
      #conditionalPanel(condition = "input.operation=='Pancan'"
                       
      #                 )
      #selectInput("Distribution", "Select dist type", 
      #            choices=c("normal","expo")),
      #sliderInput("sampleSize","Select sample size: ", 
      #            min=100,max=5000,value=1000,step=100),
      #conditionalPanel(condition = "input.Distribution=='normal'",
      #                 textInput("Mean","Plase Select the mean",10),
      #                 textInput("sd","Please Select Standard Deviation",3)),
      #conditionalPanel(condition="input.Distribution == 'expo'",
      #                 textInput("lambda","Pleas select Expo Lambda",1)
      #)
    ),
    mainPanel(
      plotOutput("myPlot")
    )
    
  )
)

