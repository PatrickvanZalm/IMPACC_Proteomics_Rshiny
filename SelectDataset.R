#### Functions ####

#Function for Normalization
PatrickNormalizer <- function(dataset){
  #We assume log2 transformed data
  dataset <- 2^dataset
  
  #Calc normalization Factor
  NormalisationFactor <- median(rowSums(dataset, na.rm = TRUE))
  
  #make copy
  dataset_normalized <- dataset
  
  #run for loop over each over the rows (samples)
  for(i in 1:nrow(dataset)) {
    
    #get sample i, apply normalisation factor on it. Change sample i in the normalized df.
    dataset_normalized[i,] <- dataset[i,] * (NormalisationFactor / sum(dataset[i,], na.rm = TRUE))
  }
  #return
  return(dataset_normalized)
  
}

# Function for Plotting CV  ####
CV_Plotter <- function(dataframe){ #function for making the plots
  #apple CV function
  dataframe_results <- data.frame(apply(dataframe, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100))
  #change colnames. Just to make life easier
  colnames(dataframe_results) <- c("Proteins")
  #Plot
  dataframe_results %>% ggplot(aes(x = rownames(dataframe_results), y = Proteins)) + geom_point() +
    theme_light() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
} 



#### Actual Shiny Code ############

#### INPUTS ####

#set up the Dataset
selectedDataset <- reactive({
  
  #Choose the dataset
  if (input$DataSelect == "DDA_Data"){
    dataset_select <- DDA_Data
  }
  else if (input$DataSelect == "DIA_Data"){
    dataset_select <- DIA_Data
  }
  else if (input$DataSelect == "Targeted_Data"){
    dataset_select <- Targeted_Data
  }
  
  #change to DF format
  dataset_select <- as.data.frame(dataset_select)
  #Set rownames and remove first column
  dataset_select_rownames <- dataset_select[,-1]
  rownames(dataset_select_rownames) <- dataset_select[,1]
  
  #Based on yes no for sample removal, we remove URS samples
  if (input$URSFilter == TRUE){
    dataset_select_rownames <- dataset_select_rownames[!grepl("URS", rownames(dataset_select_rownames)),]
  }
  
  #Based on yes no for sample removal, we remove EMT samples
  if (input$EmoryFilter == TRUE){
    dataset_select_rownames <- dataset_select_rownames[!grepl("EDT-PRT", rownames(dataset_select_rownames)),]
  }
  
  #Based on yes no for sample removal, we remove COVID19 Negative samples
  if (input$NonPositive == TRUE){
    try(dataset_select_rownames <- dataset_select_rownames %>% filter(!rownames(dataset_select_rownames) %in% COVID19Negative$sample_id))
  }
  
  return(dataset_select_rownames)
})

#Normalisation
selectedDataset_Normalized <- reactive({
  if (input$Normalization == 2){
    dataset_select_normalized <- PatrickNormalizer(selectedDataset())
  }
  if (input$Normalization == 3){
    dataset_select_normalized <- as.data.frame(NormalyzerDE::performVSNNormalization(as.matrix(2^selectedDataset())))
  }
  if (input$Normalization == 1){
    dataset_select_normalized <- selectedDataset()
  }
  
  return(dataset_select_normalized)
  
})

#Cut off filter
selectedDataset_Filtered <- reactive({
  dataset_select_filtered <- selectedDataset_Normalized() %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >= (100 -input$CutOffFilter))
  
  return(dataset_select_filtered)
  
})

#Connect with the clinical Data
selectedDataset_Clinical <- reactive({
  #put the rownames back as column again.
  selectedDataset_Filtered_rownametocolumn <- selectedDataset_Filtered() %>% rownames_to_column()
  
  #We only want to keep sample & proteomics info on the ones we have of BOTH of. Thus we use inner_join.
  dataset_select_clinical <- inner_join(clinical_data, selectedDataset_Filtered_rownametocolumn,  by = c("sample_id" = "rowname"))
  
  return(dataset_select_clinical)
  
})



#### OUTPUTS ####

#Venny Diagram of DDA, DIA, Targeted
output$VennyThreeMethods <- renderPlot({
  
  #Filter based on Input
  DDA_Filtered <- DDA_Data %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >= (100 -input$CutOffFilter))
  DIA_Filtered <- DIA_Data %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >= (100 -input$CutOffFilter))
  Targeted_Filtered <- Targeted_Data %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >= (100 -input$CutOffFilter))
  
  #Keep column names --> proteinID
  DDA_Prots <- colnames(DDA_Filtered[,2:ncol(DDA_Filtered)])
  DIA_Prots <- colnames(DIA_Filtered[,2:ncol(DIA_Filtered)])
  Targeted_Prots <- colnames(Targeted_Filtered[,2:ncol(Targeted_Filtered)])
  
  #List the three
  ProteinIDList <- list("DDA" = DDA_Prots, "DIA" = DIA_Prots, "Targeted" = Targeted_Prots)
  
  ProteinIDListLong <- data.frame(list("DDA" = ncol(DDA_Filtered[,2:ncol(DDA_Filtered)]),
                            "DIA" = ncol(DIA_Filtered[,2:ncol(DIA_Filtered)]),
                            "Targeted" = ncol(Targeted_Filtered[,2:ncol(Targeted_Filtered)]))) %>%
    gather(key = "Method", value = "Protein_Numbers")
  
  barplot <- try(ProteinIDListLong %>% ggplot(aes(x = Method, y = Protein_Numbers)) +
    geom_bar(stat = "identity") +
    theme_light())
  
  
  #Plot
  venny <- plot(euler(ProteinIDList, shape = "ellipse"), quantities = TRUE)
  
  ggarrange(barplot, venny, ncol = 2, nrow = 1)
  
})


output$DataSelection = renderText({
  print(paste("There is",
              dim(selectedDataset_Filtered())[1],
              "Samples Selected and",
              dim(selectedDataset_Filtered())[2],
              "Protein Selected."))
})

output$DataSelectionPlot <- renderPlot({
  items <- c("Samples", "Proteins")
  numbers <- c(dim(selectedDataset_Filtered())[1], dim(selectedDataset_Filtered())[2])
  
  DF <- data.frame(items, numbers)
  
  DF %>% ggplot(aes(x = items, y = numbers, fill = items)) + 
    geom_bar(stat ="identity") +
    theme_light() + 
    geom_text(aes(label= numbers), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(legend.position = "none") + 
    ggtitle("Number of Samples & Proteins") +
    xlab("") + ylab("") +
    scale_fill_manual(values = c("darkred", "darkblue"))
  
})

output$SampleBarPlot <- renderPlot({
  CV_Plotter(selectedDataset_Filtered())
  
})



