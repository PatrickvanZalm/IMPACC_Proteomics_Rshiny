

library(shiny)
library(shinyWidgets)
library(shinyjs)
library(bslib)
library(markdown)

ui <- shinyUI(navbarPage(theme = bs_theme(version = 4, bootswatch = "flatly"), #Set a theme. More can be selected at bootswatch.com
                         
                         "IMPACC Proteomics Analysis", #Set the Title
                         
                         #make first Tab. Select the dataset (DDA, DIA, Targeted)
                         tabPanel("Select-Dataset", #####
                                  
                                  #Sidebar Items
                                  sidebarLayout(
                                    sidebarPanel(
                                      
                                      #Select the Dataset
                                      radioButtons("DataSelect", "Select Dataset",
                                                   c("DDA"= "DDA_Data",
                                                     "DIA"= "DIA_Data",
                                                     "Targeted" = "Targeted_Data")),
                                      
                                      #Select If URS samples should be removed
                                      radioButtons("URSFilter", "Remove URS Samples?",
                                                   c("Yes" = TRUE,
                                                     "No" = FALSE)
                                      ),
                                      
                                      #select for Emory
                                      radioButtons("EmoryFilter", "Remove Emory Samples?",
                                                   c("Yes" = TRUE,
                                                     "No" = FALSE)                                
                                      ),
                                      
                                      #select for Non COVID samples
                                      radioButtons("NonPositive", "Remove SARS-CoV-2 Negative?",
                                                   c("Yes" = TRUE,
                                                     "No" = FALSE)
                                      ),
                                      
                                      #select normalisation (only for DDA its nice)
                                      selectInput("Normalization", "Choose Normalisation (DDA)",
                                                  choices = list("None" = 1, "Median" = 2, "VSN" = 3),
                                                  selected = 1),
                                      
                                      #cut off filter (50% is standard)
                                      sliderInput("CutOffFilter", "Cut Off Filter (% NA)",
                                                  min = 2, max = 99, value = 50),
                                      
                                      
                                    ),
                                    
                                    #Main Panel
                                    mainPanel(
                                      tabsetPanel(type = "tabs",
                                                  
                                                  #Panel 1, Data selection info
                                                  tabPanel("Dataset Info",
                                                           plotOutput("DataSelectionPlot")
                                                  ),
                                                  
                                                  # #Panel 2, Plot XXXXXXX
                                                  # tabPanel("XXXXX",
                                                  #          plotOutput("SampleBarPlot")
                                                  # ),
                                                  
                                                  tabPanel("Venny (DDA vs. DIA vs. Targeted)",
                                                           plotOutput("VennyThreeMethods")
                                                  )
                                                  
                                      )
                                      
                                    )
                                  )
                                  
                         ),
                         
                         
                         
                         tabPanel("Clinical Info", #####
                                  
                                  #Add a sidebar
                                  sidebarLayout(
                                    sidebarPanel(
                                      sliderInput("HistBins", "Number of Bins for Histogram",
                                                  min = 5, max = 50, value = 25)
                                    ),
                                    
                                    mainPanel(
                                      
                                      
                                      
                                      tabsetPanel(type = "tabs",
                                                  tabPanel("Table Info",
                                                           htmlOutput("BasicClinicalInfo"),
                                                           
                                                  ),
                                                  
                                                  tabPanel("Clinical group: day 0 & 14",
                                                           plotOutput("dayonetofourteen", height = "600px")
                                                  ),
                                                  
                                                  
                                                  tabPanel("day 0 flow to day 14",
                                                           plotOutput("alluvianPlot")
                                                  ),
                                                  
                                                  tabPanel("Number of Samples per visit",
                                                           plotOutput("ClinicalPlot_Bar")
                                                  ),
                                                  
                                                  tabPanel("Age Distribution (Histogram)",
                                                           plotOutput("AgeHistogram")
                                                  ),
                                                  tabPanel("Sex",
                                                           plotOutput("ClinicalSex")
                                                  ),
                                                  tabPanel("Enrollment Sites",
                                                           plotOutput("SampleLocation")
                                                  ),
                                                  tabPanel("Ethnicity",
                                                           plotOutput("ClinicalRace"))
                                                  
                                                  
                                                  
                                      )
                                      
                                    )
                                    
                                    
                                  )
                                  
                                  
                                  
                                  
                                  
                                  
                         ),
                         
                         tabPanel("Clustering", #####
                                  
                                  sidebarLayout(
                                    sidebarPanel(
                                      
                                      selectInput("Scaler", "Should Data be Scaled?",
                                                  choices = list("Yes" = TRUE, "No" = FALSE),
                                                  selected = TRUE),
                                      
                                      selectInput("Imputation", "Choose an imputation Method",
                                                  choices = list("Half Min Value" = 1, "Replace with 0" = 2),
                                                  selected = 1),
                                      
                                      selectInput("Colouring", "What Clinical Factor Should be Coloured?",
                                                  choices = list("Visit Numbers" = "event_type",
                                                                 "Respiratory Status (day 0)" = "respiratory_status",
                                                                 "Respiratory Status (day 14)" = "respiratory_status_day14",
                                                                 "Enrollment Site" = "enrollment_site",
                                                                 "Age" = "admit_age",
                                                                 "Sex" = "sex",
                                                                 "Race" = "race",
                                                                 "Clinical Group (day 0)" = "EndpointDay0",
                                                                 "Clinical Group (day 14)" = "EndpointDay14",
                                                                 "Death" = "death")),
                                      
                                      selectInput("ClusterElipse", "Add Elipse on Plot?",
                                                  choices = list("Yes" = TRUE, "No" = FALSE),
                                                  selected = "No"),
                                      
                                      sliderInput("Perplexity", "t-SNE Perplexity",
                                                  min = 1, max = 200, value = 15),
                                      
                                      
                                      pickerInput("pvcavariables","Select Variables for PVCA",
                                                  choices= colnames(clinical_data %>% dplyr::select(- samples_collected,
                                                                                             - sample_id, 
                                                                                             - event_type,
                                                                                             - event_id,
                                                                                             - event_location,
                                                                                             - participant_type,
                                                                                             -sample_type
                                                  )),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = colnames(clinical_data %>% dplyr::select(- samples_collected,
                                                                                               - sample_id, 
                                                                                               - event_type,
                                                                                               - event_id,
                                                                                               - event_location,
                                                                                               - participant_type,
                                                                                               - sample_type,
                                                                                               - sex,
                                                                                               - pcr_date,
                                                                                               - EndpointDay0,
                                                                                               - event_date,
                                                                                               - ethnicity,
                                                                                               - symptom_date,
                                                                                               - EndpointDay14
                                                                                               
                                                                                               
                                                  )))
                                      
                                      
                                      
                                    ),
                                    
                                    mainPanel(
                                      
                                      tabsetPanel(type = "tabs",
                                                  
                                                  tabPanel("PCA",
                                                           plotOutput("PCAPlot")),
                                                  tabPanel("t-SNE",
                                                           plotOutput("tSNE")),
                                                  tabPanel("UMAP",
                                                           plotOutput("UMAP")),
                                                  tabPanel("PVCA",
                                                           plotOutput("PVCA")),
                                                  useShinyjs()
                                                  
                                                  
                                                  
                                      )
                                    )
                                    
                                    
                                  )
                                  
                                  
                         ),
                         
                         tabPanel("Longitudinal Projectories", #####
                                  
                                  
                                  
                                  sidebarLayout(
                                    sidebarPanel(
                                      pickerInput("LongFilter","Select visits to include in Analysis",
                                                  choices= c("Visit 1", "Visit 2", "Visit 3", "Visit 4", "Visit 5", "Visit 6"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected =c("Visit 1", "Visit 2", "Visit 3")),
                                      
                                      radioButtons("AnovaType", "Type of ANOVA",
                                                   c("ANOVA" = 1,
                                                     "within-subjects ANOVA (only parametric)" = 2)),
                                      
                                      radioButtons("AnovaParametric", "Parametric / non-parametric",
                                                   c("parametric" = 1,
                                                     "Non-Parametric" = 2)),
                                      pickerInput("SpaghettiFacet","Select Variable for Split Spaghetti Plot",
                                                  choices= colnames(clinical_data %>% dplyr::select(- samples_collected,
                                                                                             - sample_id, 
                                                                                             - event_type,
                                                                                             - event_id,
                                                                                             - event_location,
                                                                                             - participant_type,
                                                                                             -sample_type,
                                                                                             -pcr_date,
                                                                                             -admit_age,
                                                                                             -event_date,
                                                                                             -symptom_date,
                                                                                             -participant_id,
                                                                                             -enrollment_site,
                                                                                             -race
                                                  )), options = list(`actions-box` = TRUE),multiple = F, selected = "sex"),
                                      pickerInput("BoxPlotProtein", "Select Protein",
                                                  choices = c(0))
                                      
                                      
                                      #NOTE TO SELF. Find a way to easily search all from a reactive(dynamic DF).... Based on that
                                      #We will build boxplots with post hoc tests.
                                      
                                    ),
                                    
                                    mainPanel(
                                      tabsetPanel(type = "tabs",
                                                  tabPanel("Number of Samples per Visit",
                                                           plotOutput("Longit_SamplePerVisit")
                                                  ),
                                                  
                                                  tabPanel("Projectories",
                                                           plotlyOutput("SpaghettiLong_all"),
                                                           plotlyOutput("Longit_volcano")
                                                  ),
                                                  
                                                  tabPanel("BoxPlot & Statistics",
                                                           plotOutput("BoxplotComparisons"),
                                                           DT::dataTableOutput("Longitudinal_Table")
                                                  ),
                                                  
                                                  tabPanel("Projectories (Split)",
                                                           plotlyOutput("SpaghettiLong_Facet", height = "750px")),
                                                  
                                                  tabPanel("Enrichment Analysis",
                                                           plotOutput("EnrichmentPlot"),
                                                           DT::dataTableOutput("EnrichmentLongitudinal")),
                                                  
                                                  tabPanel("GO Similarity",
                                                           plotOutput("GoPlotEnrich" , height = "600px"))
                                                           
                                                           
                                                  
                                                  
                                                  
                                      )
                                    )
                                    
                                  )
                         ),
                         
                         tabPanel("VolcanoPlotter", ######
                                  
                                  
                                  
                                  sidebarLayout(
                                    sidebarPanel(
                                      sliderInput("Volcano_eventDate", "Filter for Event date (days)",
                                                  min = min(clinical_data$event_date), max = max(clinical_data$event_date),
                                                  value = c(min(clinical_data$event_date),max(clinical_data$event_date))),
                                      
                                      pickerInput("Volcano_EventType","Select visits to include in Volcanoplot",
                                                  choices= c("Escalation 1", "Visit 1", "Visit 2", "Visit 3", "Visit 4", "Visit 5", "Visit 6"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Visit 1", "Visit 2", "Visit 3")),
                                      
                                      pickerInput("Volcano_Endpoint0","Select Endpoint day 0 to include in Volcanoplot",
                                                  choices= c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)")),
                                      
                                      pickerInput("Volcano_Endpoint14","Select Endpoint day 14 to include in Volcanoplot",
                                                  choices= c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)")),
                                      
                                      pickerInput("Volcano_Sex","Select sex to include in Volcanoplot",
                                                  choices= c("Female", "Male"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Female", "Male")),
                                      
                                      pickerInput("Volcano_death","Select sex to include in Volcanoplot",
                                                  choices= c("Alive" = FALSE, "Dead" = TRUE),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Alive" = FALSE, "Dead" = TRUE)),
                                      
                                      sliderInput("Volcano_Age", "Filter for age of Participants in Volcanoplot",
                                                  min = min(clinical_data$admit_age), max = max(clinical_data$admit_age),
                                                  value = c(min(clinical_data$admit_age),max(clinical_data$admit_age))),
                             
                                      pickerInput("GroupCompare","Which variable do you want to Compare?",
                                                  choices= c("Visit Number" = "event_type",
                                                             "Endpoint day 0" = "EndpointDay0", "Endpoint day 14" = "EndpointDay14",
                                                             "Sex" = "sex", "Death" = "death", "Age" = "age"),
                                                  options = list(`actions-box` = TRUE),multiple = F),
                                      
                                      pickerInput("GroupA","Select which for group A (DONT OVERLAP WITH B)",
                                                  choices= c(0),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected =c(0)),
                                      
                                      pickerInput("GroupB","Select which for group B (DONT OVERLAP WITH A)",
                                                  choices= c(0),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected =c(0)),
                                      
                                      sliderInput("Split_Age", "Filter for age of Participants in Volcanoplot if age is selected. Group A will include selected number",
                                                  min = min(clinical_data$admit_age), max = max(clinical_data$admit_age),
                                                  value = 65),
                                      
                                      
                                      
                                    ),
                                    
                                    mainPanel(
                                      tabsetPanel(type = "tabs",
                                                  tabPanel("Selected Number of Samples",
                                                           plotOutput("SelectedSamples_Volcano"),
                                                  ),
                                                  tabPanel("VolcanoPlot",
                                                           plotlyOutput("VolcanoPlotly"),
                                                           DT::dataTableOutput("VolcanoTable")),
                                                  
                                                  tabPanel("Enrichment Analysis",
                                                           plotOutput("EnrichmentPlot_Volcano"),
                                                           DT::dataTableOutput("EnrichmentVolcano")),
                                                  
                                                  tabPanel("GO Similarity",
                                                           plotOutput("GoPlotEnrich_Volcano" , height = "600px"))
                                                  
                                                  )
                                    )
                                    
                                  )
                         ),
                         
                         
                         
                         tabPanel("Heatmap", ######
                                  
                                  
                                  
                                  sidebarLayout(
                                    sidebarPanel(
                                      sliderInput("Heatmap_eventDate", "Filter for Event date (days)",
                                                  min = min(clinical_data$event_date), max = max(clinical_data$event_date),
                                                  value = c(min(clinical_data$event_date),max(clinical_data$event_date))),
                                      
                                      pickerInput("Heatmap_EventType","Select visits to include in Heatmap",
                                                  choices= c("Escalation 1", "Visit 1", "Visit 2", "Visit 3", "Visit 4", "Visit 5", "Visit 6"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Visit 1", "Visit 2", "Visit 3")),
                                      
                                      pickerInput("Heatmap_Endpoint0","Select Endpoint day 0 to include in Heatmap",
                                                  choices= c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)")),
                                      
                                      pickerInput("Heatmap_Endpoint14","Select Endpoint day 14 to include in Heatmap",
                                                  choices= c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Discharged(1-2)", "Hospitalized(3-4)", "Hospitalized(5-6)", "Dead(7)")),
                                      
                                      pickerInput("Heatmap_Sex","Select sex to include in Heatmap",
                                                  choices= c("Female", "Male"),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Female", "Male")),
                                      
                                      pickerInput("Heatmap_death","Select sex to include in Heatmap",
                                                  choices= c("Alive" = FALSE, "Dead" = TRUE),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected = c("Alive" = FALSE, "Dead" = TRUE)),
                                      
                                      sliderInput("Heatmap_Age", "Filter for age of Participants in Heatmap",
                                                  min = min(clinical_data$admit_age), max = max(clinical_data$admit_age),
                                                  value = c(min(clinical_data$admit_age),max(clinical_data$admit_age))),
                                      
                                      pickerInput("Heatmap_Grouping","Do you want to group (mean) some factors?",
                                                  choices= c("Yes" = T, "No" = F),
                                                  options = list(`actions-box` = TRUE),multiple = F,
                                                  selected = c("Yes" = TRUE)),
                                      
                                      pickerInput("Heatmap_groupby","Which Grouping?",
                                                  choices= c("Visit Number" = "event_type",
                                                             "Endpoint day 0" = "EndpointDay0", "Endpoint day 14" = "EndpointDay14",
                                                             "Sex" = "sex", "Death" = "death"),
                                                  options = list(`actions-box` = TRUE),multiple = T)
                                                  ,
                               
                                      
                                      pickerInput("Heatmap_Protein","Select Proteins to include in Heatmap (significant selected automatically)",
                                                  choices= c(0),
                                                  options = list(`actions-box` = TRUE),multiple = T,
                                                  selected =c(0)),
                                      
                                    ),
                                    
                                    mainPanel(
                                      tabsetPanel(type = "tabs",
                                                  tabPanel("Selected Number of Samples",
                                                           plotOutput("SelectedSamples")
                                                  ),
                                                  
                                                  tabPanel("HeatMap",
                                                           plotOutput("HMT", height = "750px")),
                                                  
                                                  tabPanel("ValueTable",
                                                           #DT::dataTableOutput("HeatMapSuper")
                                                  )
      
                                      )
                                    )
                                    
                                  )
                         ),
                         
            
                         tabPanel("", ######
                                  
                                  
                                  
                                  sidebarLayout(
                                    sidebarPanel(),
                                    
                                    mainPanel(
                                      tabsetPanel()
                                    )
                                    
                                  )
                         )
                         
                         
) #navbarpage
)#shinyUI