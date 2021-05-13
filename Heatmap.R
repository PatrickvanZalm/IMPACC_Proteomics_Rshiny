#Filters on the visit numbers you selected. If you choose within subject ANOVA, it also only keeps participants that has all the selected visits.
Heatmap_filterData <- reactive({
  #Apply all Filters
  dataset <- selectedDataset_Clinical()
  
  dataset <- dataset %>% filter(event_date >= input$Heatmap_eventDate[1] & event_date <= input$Heatmap_eventDate[2])
  
  dataset <- dataset %>% filter(event_type %in% input$Heatmap_EventType)
  
  dataset <- dataset %>% filter(EndpointDay0 %in% input$Heatmap_Endpoint0)
  
  dataset <- dataset %>% filter(EndpointDay14 %in% input$Heatmap_Endpoint14)
  
  dataset <- dataset %>% filter(sex %in% input$Heatmap_Sex)
  
  dataset <- dataset %>% filter(death %in% input$Heatmap_death)
  
  dataset <- dataset %>% filter(admit_age >= input$Heatmap_Age[1] & admit_age <= input$Heatmap_Age[2])
  
  return(dataset)
  
})

observeEvent(StatResults(), #Observe changes in the reactive DF.
             updatePickerInput(session = session, inputId = "Heatmap_Protein",
                               choices = StatResults()$Protein,
                               selected = (StatResults() %>% filter(StatResults()$Padjusted < 0.05))$Protein)
)

output$SelectedSamples <- renderPlot({
  
  Heatmap_filterData() %>% tally() %>%
    ggplot(aes(1, y = n))+
    geom_bar(stat = "identity") +
    theme_light() + 
    geom_text(aes(label= n), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(legend.position = "none") + 
    ggtitle("Number of Samples Selected") +
    xlab("") + ylab("") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
})


#Plot the heatmap.
output$HeatMapSuper <- DT::renderDataTable({
  
  if (input$Heatmap_Grouping ==T){
    
    #copy data
    grouped_data <- Heatmap_filterData()
    
    grouped_data[,22:ncol(grouped_data)] <- scale(grouped_data[,22:ncol(grouped_data)])
    
    GroupedData <- grouped_data %>% dplyr::group_by(across(all_of(input$Heatmap_groupby))) %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>% ungroup()
    
    clinical <- t(as.matrix(Heatmap_filterData() %>% dplyr::select(colnames(clinical_data))))
    proteins <- t(as.matrix(scale(Heatmap_filterData() %>% dplyr::select(!!!input$Heatmap_Protein))))
                  
    
    return(proteins)
  }
  
  else if (input$Heatmap_Grouping == F){
    Heatmap_filterData()
  }
  

})


output$HMT <- renderPlot({
  
  if (input$Heatmap_Grouping ==T){
    
    #copy data
    grouped_data <- Heatmap_filterData()
    #Select only numeric data and scale it
    grouped_data[,22:ncol(grouped_data)] <- scale(grouped_data[,22:ncol(grouped_data)])
    #Group based on input. Calculate mean. Ungroup
    GroupedData <- grouped_data %>% dplyr::group_by(across(all_of(input$Heatmap_groupby))) %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>% ungroup()
    
    #Split the two. We select proteins based on input and the group_by factors --> these should be of interest here
    clinical <- t(as.matrix(GroupedData %>% dplyr::select(!!!input$Heatmap_groupby)))
    proteins <- t(as.matrix(GroupedData %>% dplyr::select(!!!input$Heatmap_Protein)))
    
    #how to make an x amount of Annotations?!
    #Try to make each of them.
    #Make a list with all of them and if not an item --> NULL.
    Annot_EndpointDay0 <- try(HeatmapAnnotation(EndpointDay0 = clinical["EndpointDay0",],
                                                col = list(EndpointDay0 = c("Discharged(1-2)" = "#A5D6A7",
                                                                            "Hospitalized(3-4)" = "#FFF59D",
                                                                            "Hospitalized(5-6)" = "#FFAB91"))))
    
    Annot_EndpointDay14 <- try(HeatmapAnnotation(EndpointDay14 = clinical["EndpointDay14",],
                                                col = list(EndpointDay14 = c("Discharged(1-2)" = "#4CAF50",
                                                                             "Hospitalized(3-4)" = "#FFEB3B",
                                                                             "Hospitalized(5-6)" = "#FF5722",
                                                                             "Dead(7)" = "#424242"))))
    Annot_sex <- try(HeatmapAnnotation(sex = clinical["sex",],
                                       col = list(sex = c("Female" = "#DC2543", "Male" = "#1E94A0"))))
    
    Annot_death <- try(HeatmapAnnotation(death = clinical["death",],
                                         col = list(death = c("TRUE" = "black", "FALSE" = "white"))))
    
    Annot_eventType <- try(HeatmapAnnotation(VisitNumber = clinical["event_type",]))
             
    #Make empty list      
    AnnotList <-  HeatmapAnnotation(foo = anno_empty(border = TRUE))
    
    #Repeat this 5 times
    if (exists("Annot_EndpointDay0")){
      try(AnnotList <- c(AnnotList, Annot_EndpointDay0))
    }
    
    
    if (exists("Annot_EndpointDay14")){
      try(AnnotList <- c(AnnotList, Annot_EndpointDay14))
    }
    
    if (exists("Annot_sex")){
      try(AnnotList <- c(AnnotList, Annot_sex))
    }
    
    if (exists("Annot_death")){
      try(AnnotList <- c(AnnotList, Annot_death))
    }
    
    if (exists("Annot_eventType")){
      try(AnnotList <- c(AnnotList, Annot_eventType))
    }
    
    #Add thingy
    heatmapPlot <- Heatmap(proteins, top_annotation = AnnotList)
    
    
    draw(heatmapPlot)
    
    
  }
  
  else if (input$Heatmap_Grouping == F){
    
    clinical <- t(as.matrix(Heatmap_filterData() %>% dplyr::select(colnames(clinical_data))))
    proteins <- t(as.matrix(scale(Heatmap_filterData() %>% dplyr::select(!!!input$Heatmap_Protein))))
    
    Annot_EndpointDay0 <- HeatmapAnnotation(EndpointDay0 = clinical["EndpointDay0",],
                                                col = list(EndpointDay0 = c("Discharged(1-2)" = "#A5D6A7",
                                                                            "Hospitalized(3-4)" = "#FFF59D",
                                                                            "Hospitalized(5-6)" = "#FFAB91")))
    
    Annot_EndpointDay14 <- HeatmapAnnotation(EndpointDay14 = clinical["EndpointDay14",],
                                                 col = list(EndpointDay14 = c("Discharged(1-2)" = "#4CAF50",
                                                                              "Hospitalized(3-4)" = "#FFEB3B",
                                                                              "Hospitalized(5-6)" = "#FF5722",
                                                                              "Dead(7)" = "#424242")))
    
    Annot_sex <- HeatmapAnnotation(sex = clinical["sex",], col = list(sex = c("Female" = "#DC2543", "Male" = "#1E94A0")))
    
                     
    
    Annot_death <- HeatmapAnnotation(death = clinical["death",],
                                         col = list(death = c("TRUE" = "black", "FALSE" = "white")))
    
    Annot_eventType <- HeatmapAnnotation(VisitNumber = clinical["event_type",])
    
    Annot_Age <- HeatmapAnnotation(age = anno_points(as.numeric(clinical["admit_age",])))
    
    Annotlist <- c(Annot_death, Annot_Age, Annot_sex, Annot_eventType, Annot_EndpointDay0, Annot_EndpointDay14)
    
    heatmapPlot <- Heatmap(proteins, top_annotation = Annotlist)
    
    draw(heatmapPlot)
    
  }
  
  
})







