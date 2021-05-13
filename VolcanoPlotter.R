#Filters on the visit numbers you selected. If you choose within subject ANOVA, it also only keeps participants that has all the selected visits.
Volcano_filterData <- reactive({
  #Apply all Filters
  dataset <- selectedDataset_Clinical()
  
  dataset <- dataset %>% filter(event_date >= input$Volcano_eventDate[1] & event_date <= input$Volcano_eventDate[2])
  
  dataset <- dataset %>% filter(event_type %in% input$Volcano_EventType)
  
  dataset <- dataset %>% filter(EndpointDay0 %in% input$Volcano_Endpoint0)
  
  dataset <- dataset %>% filter(EndpointDay14 %in% input$Volcano_Endpoint14)
  
  dataset <- dataset %>% filter(sex %in% input$Volcano_Sex)
  
  dataset <- dataset %>% filter(death %in% input$Volcano_death)
  
  dataset <- dataset %>% filter(admit_age >= input$Volcano_Age[1] & admit_age <= input$Volcano_Age[2])
  
  return(dataset)
  
})

#Observe which variable is selected.
observeEvent(input$GroupCompare, #Observe changes in the reactive DF.
             
             if (input$GroupCompare == "EndpointDay0"){
               try(updatePickerInput(session = session, inputId = "GroupA",
                                     choices = unique(unfactor(clinical_data$EndpointDay0)))
                   )
             }
             
             else if (input$GroupCompare == "EndpointDay14"){
               try(updatePickerInput(session = session, inputId = "GroupA",
                                     choices = unique(unfactor(clinical_data$EndpointDay14)))
               )
             }
             
             else{
               try(updatePickerInput(session = session, inputId = "GroupA",
                                     choices = unique(Volcano_filterData() %>% dplyr::select(!!input$GroupCompare))))
             }
             
)

#Same for GroupB
observeEvent(input$GroupCompare, #Observe changes in the reactive DF.
             
             if (input$GroupCompare == "EndpointDay0"){
               try(updatePickerInput(session = session, inputId = "GroupB",
                                     choices = unique(unfactor(clinical_data$EndpointDay0)))
               )
             }
             else if (input$GroupCompare == "EndpointDay14"){
               try(updatePickerInput(session = session, inputId = "GroupB",
                                     choices = unique(unfactor(clinical_data$EndpointDay14)))
               )
             }
             
             else{
               try(updatePickerInput(session = session, inputId = "GroupB",
                                     choices = unique(Volcano_filterData() %>% dplyr::select(!!input$GroupCompare))))
             }
             
)
Volcano_SelectedGroup <- reactive({
  
  ProteinNamesTtest <- colnames(Volcano_filterData())[22:ncol(Volcano_filterData())]
  
  if (input$GroupCompare == "age"){
    #SubGroup them. Add column with which Grouping Info.
    GroupA <- Volcano_filterData()[Volcano_filterData()$admit_age <= input$Split_Age, ]
    GroupB <- Volcano_filterData()[Volcano_filterData()$admit_age > input$Split_Age, ]
    
  }
  
  else{
    #SubGroup them. Add column with which Grouping Info.
    GroupA <- Volcano_filterData()[Volcano_filterData()[[input$GroupCompare]] %in% input$GroupA, ]
    GroupB <- Volcano_filterData()[Volcano_filterData()[[input$GroupCompare]] %in% input$GroupB, ]
  }
  
  GroupA$Group <- "A"
  GroupB$Group <- "B"
  
  return(rbind(GroupA, GroupB))
})

Volcano_Significant <- reactive({
  
  #ID proteins
  ProteinNamesTtest <- colnames(Volcano_filterData())[22:ncol(Volcano_filterData())]
  
  #longFormat it
  SelectedGroups <- Volcano_SelectedGroup() %>% dplyr::select(Group, ProteinNamesTtest) %>% 
    gather(ProteinNamesTtest, key = "Protein", value = "Intensity")
  
  #run T.test
  ttestResults <- result <- group_by(SelectedGroups, Protein)
  ttestResults <- do(result, tidy(t.test(.$Intensity ~ .$Group)))
  
  #Add Fold Change and P.adjusted
  ttestResults$FC <- ttestResults$estimate1 - ttestResults$estimate2
  ttestResults$Padjusted <- p.adjust(ttestResults$p.value, method = "BH")
  
  #select interesting columns
  ttestResults <- ttestResults %>% dplyr::select(Protein, FC, p.value, Padjusted)
  
  return(ttestResults)
})


output$SelectedSamples_Volcano <- renderPlot({
  
  Volcano_SelectedGroup() %>%
    group_by(Group) %>%
    tally() %>%
    ggplot(aes(x = Group, y = n)) +
    geom_bar(sta = "identity") +
    theme_light() + 
    geom_text(aes(label= n), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(legend.position = "none") + 
    ggtitle("Number of Samples Selected per Group") +
    xlab("") + ylab("")
  
})

output$VolcanoPlotly <- renderPlotly({
  
  #Copy, set true false for significant
  dataset <- Volcano_Significant()
  dataset$Significant <- dataset$Padjusted < 0.05
  
  ggplotly(
    ggplot(dataset, aes(x=FC, y= -log10(Padjusted), label=Protein)) +
      geom_point(aes(fill=Significant), size = 2.3, alpha = .8) +
      scale_fill_manual(values=c("Grey", "orange")) + 
      theme_light() +
      xlab("Fold Change (Group A - Group B)") + ylab("-log10(adjusted p-value")
  )
  
})


output$VolcanoTable <- DT::renderDataTable({
  Volcano_Significant()
})


output$EnrichmentPlot_Volcano <- renderPlot({
  #Map Significant Proteins
  Significant <- Volcano_Significant()[Volcano_Significant()$Padjusted < 0.05, 1]
  
  #run our proteins through the DB. Make STRING thingy
  string_db$plot_network(Significant)
  
})

output$EnrichmentVolcano <- DT::renderDataTable({
  
  #Map Significant Proteins
  Significant <- Volcano_Significant()[Volcano_Significant()$Padjusted < 0.05, 1]
  
  #run our proteins through the DB
  enrichmentResults_Long <- string_db$get_enrichment(Significant)
  #filter only GO process, KEGG and RCTM
  enrichmentResults_Long <- enrichmentResults_Long %>%
    filter(category == "Process" | category ==  "RCTM" | category ==  "KEGG") %>%
    dplyr::select(-ncbiTaxonId, -number_of_genes, -number_of_genes_in_background, -preferredNames) %>%
    relocate(inputGenes, .after = description)
  
  DT::datatable(enrichmentResults_Long) %>% DT::formatStyle(columns = c(0:4), fontSize = '80%')
  
})

output$GoPlotEnrich_Volcano <- renderPlot({
  
  #Map Significant Proteins
  Significant <- Volcano_Significant()[Volcano_Significant()$Padjusted < 0.05, 1]
  
  #run our proteins through the DB
  enrichmentResults_Long <- string_db$get_enrichment(Significant)
  #filter only GO process, KEGG and RCTM
  enrichmentResults_Long <- enrichmentResults_Long %>%
    filter(category == "Process") 
  
  enrichmentResults_Long$term <- gsub("\\.", ':', enrichmentResults_Long$term)
  
  mat <- GO_similarity(enrichmentResults_Long$term)
  
  cl <- cluster_terms(mat)
  
  plotSimi <- ht_clusters(mat, cl, fontsize = runif(30, min = 15, max = 25))
  
  draw(plotSimi)
  
})


