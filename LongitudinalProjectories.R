#Filters on the visit numbers you selected. If you choose within subject ANOVA, it also only keeps participants that has all the selected visits.
dataset_longitudinal <- reactive({
  
  #initial Copy of Data
  dataset <- selectedDataset_Clinical()
  
  dataset_eventFilter <- dataset %>% filter(event_type %in% input$LongFilter)
  
  #if within we only want to keep participants we have data for across the selected visits.
  if (input$AnovaType == 2){
    dataset_partFilter <- dataset_eventFilter[dataset_eventFilter$participant_id %in% names(which(table(dataset_eventFilter$participant_id) == length(unique(dataset_eventFilter$event_type)))), ]
  }
  
  else if (input$AnovaType == 1){
    dataset_partFilter <- dataset_eventFilter
  }
  
  return(dataset_partFilter)
})

StatResults <- reactive({
  
  #grab the dataset
  dataset <- dataset_longitudinal()
  
  #Set Protein Names 
  Proteinnames_anova <- colnames(dataset)[22:ncol(dataset)]
  
  #LongFormat the data
  dataset_long <- dataset %>% dplyr::select(event_type, participant_id, Proteinnames_anova) %>%
    gather(Proteinnames_anova, key = "Protein", value = "Intensity") 
  
  #We do a split of the dataframe based on unique ProteinID. Then we run the Anova.
  #Based on some input values the formulas change slightly
  
  
  #Normal ANOVA, if parametric is choosen
  if (input$AnovaType == 1 && input$AnovaParametric == 1){
    dataset_anovaResults <- lapply(split(dataset_long, dataset_long$Protein), aov, formula = Intensity ~ event_type)
    
    anovaResults <- cbind(
      data.frame(sapply(dataset_anovaResults, function(x) summary(x)[[1]][["F value"]][[1]])),
      data.frame(sapply(dataset_anovaResults, function(x) summary(x)[[1]][["Pr(>F)"]][[1]]))) %>%
      rownames_to_column()
    
  }
  
  #Within Subject, if parametric is choosen
  else if (input$AnovaType == 2 && input$AnovaParametric == 1){
    
    dataset_anovaResults <- lapply(split(dataset_long, dataset_long$Protein), aov, formula = Intensity ~ event_type + Error(participant_id/event_type))
    
    anovaResults <- cbind(
      data.frame(sapply(dataset_anovaResults, function(x) summary(x)$`Error: participant_id:event_type`[[1]]$'F value'[1])),
      data.frame(sapply(dataset_anovaResults, function(x) summary(x)$`Error: participant_id:event_type`[[1]]$'Pr(>F)'[1]))) %>%
      rownames_to_column()
    
  }
  
  #Kruskal Wallis (Non Para ANOVA)
  else if (input$AnovaType == 1 && input$AnovaParametric == 2){
    dataset_anovaResults <- lapply(split(dataset_long, dataset_long$Protein), function(x){kruskal.test(Intensity ~ event_type, data = x)})
    
    anovaResults <- cbind(
      data.frame(sapply(dataset_anovaResults, function(x) x$statistic )),
      data.frame(sapply(dataset_anovaResults, function(x) x$p.value))) %>%
      rownames_to_column()
    
    anovaResults$rowname <- substr(anovaResults$rowname, 1, nchar(anovaResults$rowname)-27)
    
  }
  
  #return empty data
  else if (input$AnovaType == 2 && input$AnovaParametric == 2){
    
    anovaResults <-  data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Dingetjes",
                                                                            "Dangetjes", "Dingen"))))
    
  }
  
  colnames(anovaResults) <- c("Protein", "Fvalue", "Pvalue")
  
  try(anovaResults$Padjusted <- p.adjust(anovaResults$Pvalue, method = "BH"))
  
  return(anovaResults)
})

############ OUTPUT ###############

output$Longit_SamplePerVisit <- renderPlot({
  dataset_longitudinal() %>% group_by(event_type) %>% tally() %>%
    ggplot(aes(x = event_type, y = n, fill = event_type)) +
    geom_bar(sta = "identity") +
    theme_light() + 
    geom_text(aes(label= n), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(legend.position = "none") + 
    ylab("Number of Samples") +
    xlab("") 
})

output$Longit_volcano <- renderPlotly({
  
  #Datasetcopy
  dataset <- StatResults()
  
  #SignificantTrueNo
  dataset$significant <- dataset$Padjusted < 0.05
  
  #Change Column Names. Plot in both cases.
  if (input$AnovaParametric == 1){
    colnames(dataset) <- c("Protein", "Fvalue", "Pvalue", "Padjusted", "significant")
    
    ggplotly(
      dataset %>% 
        ggplot(aes(x = Fvalue, y = -log10(Padjusted), label = Protein)) +
        geom_point(aes(fill = significant), size = 2.3, alpha = 0.6) +
        scale_fill_manual(values=c("Grey", "orange")) + 
        theme_light() +
        xlab("Anova F value") + ylab("-log10(adjusted p-value")
    )
    
  }
  
  else if (input$AnovaParametric == 2){
    colnames(dataset) <- c("Protein", "ChiSquared", "Pvalue", "Padjusted", "significant")
    
    ggplotly(
      dataset %>% 
        ggplot(aes(x = ChiSquared, y = -log10(Padjusted), label = Protein)) +
        geom_point(aes(fill = significant), size = 2.3, alpha = 0.6) +
        scale_fill_manual(values=c("Grey", "orange")) + 
        theme_light() +
        xlab("Kruskal Wallis Chi Squared") + ylab("-log10(adjusted p-value")
    )
    
  }
  
  
})

output$Longitudinal_Table <- DT::renderDataTable({
  StatResults()
})



output$SpaghettiLong_all <- renderPlotly({
  
  #Copy
  dataset_scaled <- dataset_longitudinal()
  
  dataset_scaled[,22:ncol(dataset_scaled)] <- scale(dataset_scaled[,22:ncol(dataset_scaled)])
  
  SignificantProteins <- (StatResults() %>% filter(StatResults()$Padjusted < 0.05))$Protein
  
  #Only keep columns that were significant
  significant <- dataset_scaled %>%
   dplyr::select(colnames(dataset_scaled)[1:21], SignificantProteins) %>% #filter for significant
    group_by(event_type) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::select(event_type, SignificantProteins) %>%
    gather(key = "ProteinID", value = "Average", SignificantProteins)
  
  ggplotly(significant %>% ggplot(aes(x = event_type, y = Average, group = ProteinID, color = "red")) +
             geom_line(alpha = 0.5)  +
             theme_light() +
             theme(legend.position = "none") +
             xlab(""))
  
})

output$SpaghettiLong_Facet <- renderPlotly({
  
  #Copy
  dataset_scaled <- dataset_longitudinal()
  
  dataset_scaled[,22:ncol(dataset_scaled)] <- scale(dataset_scaled[,22:ncol(dataset_scaled)])
  
  SignificantProteins <- (StatResults() %>% filter(StatResults()$Padjusted < 0.05))$Protein
  
  #Only keep columns that were significant
  significant <- dataset_scaled %>% dplyr::select(event_type, !!!input$SpaghettiFacet, SignificantProteins) %>%
    group_by_("event_type", input$SpaghettiFacet) %>%
    summarise_at(SignificantProteins, mean, na.rm = T) %>%
    gather(key = "ProteinID", value = "Average", SignificantProteins)
  
  plotly::ggplotly(significant %>% ggplot(aes(x = event_type, y = Average, group = ProteinID, colour = input$SpaghettiFacet)) +
                     geom_line(alpha = 0.5)  +
                     theme_light() +
                     theme(legend.position = "none") +
                     xlab("") +
                     facet_grid(input$SpaghettiFacet))
  
})

observeEvent(StatResults(), #Observe changes in the reactive DF.
             updatePickerInput(session = session, inputId = "BoxPlotProtein",
                               choices = (StatResults() %>% filter(StatResults()$Padjusted < 0.05))$Protein)
             
)


output$BoxplotComparisons <- renderPlot({
  
  Comparisons <- data.frame(t(data.frame(combn(input$LongFilter, 2)))) %>% rownames_to_column()
  
  Comparisons <- by(Comparisons, Comparisons$rowname, function(y) unlist(y[y != ''][-1]))
  
  #Based on the type of ANOVA you choose we have to do t.test, paired t.test or Wilcoxon Rank Sum.
  
  #Normal Anova, Parametric assumptions --> t.test
  if (input$AnovaType == 1 && input$AnovaParametric == 1){
    dataset_longitudinal() %>% 
      ggplot(aes(x = event_type, y = !!as.symbol(input$BoxPlotProtein), fill = event_type)) +
      geom_boxplot() +
      stat_compare_means(comparisons = Comparisons, method = "t.test") + # Add pairwise comparisons p-value
      theme(legend.position="right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1)) +
      xlab("")
  }
  
  #Paired Anova, Parametric assumptions --> paired t.test
  else if (input$AnovaType == 2 && input$AnovaParametric == 1){
    dataset_longitudinal() %>% 
      ggplot(aes(x = event_type, y = !!as.symbol(input$BoxPlotProtein), fill = event_type)) +
      geom_boxplot() +
      stat_compare_means(comparisons = Comparisons, paired = T, method = "t.test") + # Add pairwise comparisons p-value
      theme(legend.position="right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1)) +
      xlab("")
  }
  
  #Kruskal, non-Parametric assumptions --> Wilcoxon.
  else if (input$AnovaType == 1 && input$AnovaParametric == 2){
    dataset_longitudinal() %>% 
      ggplot(aes(x = event_type, y = !!as.symbol(input$BoxPlotProtein), fill = event_type)) +
      geom_boxplot() +
      stat_compare_means(comparisons = Comparisons,  method = "wilcox.test") + # Add pairwise comparisons p-value
      theme(legend.position="right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1)) +
      xlab("")
  }
})

output$EnrichmentPlot <- renderPlot({
  #Map Significant Proteins
  Significant <- StatResults()[StatResults()$Padjusted < 0.05, 1]
  
  #run our proteins through the DB. Make STRING thingy
  string_db$plot_network(Significant)
  
})

output$EnrichmentLongitudinal <- DT::renderDataTable({
  
  #Map Significant Proteins
  Significant <- StatResults()[StatResults()$Padjusted < 0.05, 1]
  
  #run our proteins through the DB
  enrichmentResults_Long <- string_db$get_enrichment(Significant)
  #filter only GO process, KEGG and RCTM
  enrichmentResults_Long <- enrichmentResults_Long %>%
    filter(category == "Process" | category ==  "RCTM" | category ==  "KEGG") %>%
    dplyr::select(-ncbiTaxonId, -number_of_genes, -number_of_genes_in_background, -preferredNames) %>%
    relocate(inputGenes, .after = description)
  
  DT::datatable(enrichmentResults_Long) %>% DT::formatStyle(columns = c(0:4), fontSize = '80%')
  
})
output$GoPlotEnrich <- renderPlot({
  
  #Map Significant Proteins
  Significant <- StatResults()[StatResults()$Padjusted < 0.05, 1]
  
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



