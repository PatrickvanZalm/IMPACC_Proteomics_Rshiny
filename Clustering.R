####
# START DATASET
# 

#PVCA function. Gotten from: https://github.com/dleelab/pvca

PVCA <- function(counts, meta, threshold, inter){
  
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)
  
  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}
  
  pred.list <- colnames(meta)
  meta <- droplevels(meta)
  
  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
  
  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit)
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}

dataset_imputed <- reactive({
  
  #initial Copy of Data
  dataset <- selectedDataset_Clinical()
  
  #subselect only protein column
  dataset_Proteins <- dataset[,22:ncol(dataset)]
  
  if (input$Imputation == 1){ #IMPUTE WITH HALF MIN
    
    #Set up Columns. Have to be ALL NUMERIC
    cols<- 1:ncol(dataset_Proteins)
    
    #Replace NA with 0
    datasetImputation <- dataset_Proteins %>% replace(is.na(.), 0)
    
    #Replace 0 with half minimum value per column
    datasetImputation[cols] <- lapply(datasetImputation[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/2))
    
    
  }
  
  else if (input$Imputation == 2){ #IMPUTE WITH 0
    
    #Set up Columns. Have to be ALL NUMERIC
    cols<- 1:ncol(dataset_Proteins)
    
    #Replace NA with 0
    datasetImputation <- dataset_Proteins %>% replace(is.na(.), 0)
    
  }
  
  if (input$Scaler == TRUE){
    
    dataset_scaled <- data.frame(t(scale(t(datasetImputation))))
    
  }
  
  else if(input$Scaler == FALSE){
    
    dataset_scaled <- datasetImputation
    
  }
  
  #Connect With Clinical Data
  
  
  return(dataset_scaled)
  
})



output$textys = renderText({
  
  print(paste("There is",
              dim(dataset_imputed())[1],
              "Samples Selected and",
              dim(dataset_imputed())[2],
              "Protein Selected."))
  
})

output$PCAPlot <- renderPlot({
  
  pca1 <- prcomp(dataset_imputed()[2:ncol(dataset_imputed())]) #transpose?
  scores <- data.frame(pca1$x) #Grab PC scores
  
  #calculate variance Percentages
  eigs <- pca1$sdev^2
  variance_percentage <- (eigs / sum(eigs))*100
  pc1var <- round(variance_percentage[1],digits=0)
  pc2var <- round(variance_percentage[2],digits=0)
  pc3var <- round(variance_percentage[3],digits=0)
  
  #grab clinical data
  dataset_clinical <- selectedDataset_Clinical()[,1:21]
  dataset_clinical$respiratory_status <- as.factor(dataset_clinical$respiratory_status)
  dataset_clinical$respiratory_status_day14 <- as.factor(dataset_clinical$respiratory_status_day14)
  
  #combine the two
  scores_clinical <- cbind(dataset_clinical, scores)
  
  pca_plot <- ggplot(data = scores_clinical, aes_string(x = "PC1", y = "PC2", color = input$Colouring)) +
    theme_light()  +
    theme(legend.position="bottom") +
    xlab(paste('PC1',' (',pc1var,'% variance)',sep='')) +
    ylab(paste('PC2',' (',pc2var,'% variance)',sep=''))
  
  if (input$ClusterElipse == TRUE){
    pca_plot +  stat_ellipse(geom = "polygon", alpha = 0.05, aes_string(fill = input$Colouring)) +
      geom_point(size = 4)
  }
  else if (input$ClusterElipse == FALSE){
    pca_plot + geom_point(size = 4)
  }
  
})


output$tSNE <- renderPlot({
  
  
  tsne <- Rtsne::Rtsne(as.matrix(dataset_imputed()), dims = 2, perplexity= input$Perplexity, verbose=FALSE, max_iter = 500)
  
  scores <- data.frame(tsne$Y) #Grab tSNE
  
  
  #grab clinical data
  dataset_clinical <- selectedDataset_Clinical()[,1:21]
  dataset_clinical$respiratory_status <- as.factor(dataset_clinical$respiratory_status)
  dataset_clinical$respiratory_status_day14 <- as.factor(dataset_clinical$respiratory_status_day14)
  
  #combine the two
  scores_clinical <- cbind(dataset_clinical, scores)
  
  tsne_plot <- ggplot(data = scores_clinical, aes_string(x = "X1", y = "X2", color = input$Colouring)) +
    theme_light()  +
    theme(legend.position="bottom")
  
  if (input$ClusterElipse == TRUE){
    tsne_plot +  stat_ellipse(geom = "polygon", alpha = 0.05, aes_string(fill = input$Colouring)) +
      geom_point(size = 4)
  }
  else if (input$ClusterElipse == FALSE){
    tsne_plot + geom_point(size = 4)
  }
  
  
})


output$UMAP <- renderPlot({
  
  
  umap <- umap::umap(as.matrix(dataset_imputed()))
  
  scores <- data.frame(umap$layout) #Grab UMAP
  
  
  #grab clinical data
  dataset_clinical <- selectedDataset_Clinical()[,1:21]
  dataset_clinical$respiratory_status <- as.factor(dataset_clinical$respiratory_status)
  dataset_clinical$respiratory_status_day14 <- as.factor(dataset_clinical$respiratory_status_day14)
  
  #combine the two
  scores_clinical <- cbind(dataset_clinical, scores)
  
  umapPlot <- ggplot(data = scores_clinical, aes_string(x = "X1", y = "X2", color = input$Colouring)) +
    theme_light()  +
    theme(legend.position="bottom") 
  
  if (input$ClusterElipse == TRUE){
    umapPlot +  stat_ellipse(geom = "polygon", alpha = 0.05, aes_string(fill = input$Colouring)) +
      geom_point(size = 4)
  }
  else if (input$ClusterElipse == FALSE){
    umapPlot + geom_point(size = 4)
  }
})

output$PVCA <- renderPlot({
  
  
  #SELECT THE CLINICAL PARAMETERS YOU WANNA TEST. SET SAMPLE ID AS ROWNAMES
  pvca_select <- data.frame(selectedDataset_Clinical() %>% dplyr::select(!!!input$pvcavariables))
  rownames(pvca_select) <- selectedDataset_Clinical()$sample_id
  
  #SELECT THE QUANT.DATA. Transpose it. COLNAMES SHOULD BE SAME AS ROWNAMES OF THE ONE ABOVE.
  pvca_quant <-data.frame(t(dataset_imputed()))
  colnames(pvca_quant) <- rownames(pvca_select)
  
  #OUTPUT IS a few numbers. Based on the number of items you did select above.
  pvca_results <- data.frame(PVCA(counts = pvca_quant,
                                  meta = pvca_select, inter = F, threshold = 0.6))
  
  pvca_results <- pvca_results %>% rownames_to_column()
  colnames(pvca_results) <- c("Clinical_Variable", "PVCA")
  
  pvca_results %>% ggplot(aes(x = reorder(Clinical_Variable, PVCA), y = PVCA)) +
    geom_bar(stat = "identity") +
    theme_light() +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1))
  
  
})