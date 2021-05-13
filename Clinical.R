
output$BasicClinicalInfo <- renderTable({
  
  
  table1_produced <- table1(~ admit_age + sex+ race + ethnicity + enrollment_site + symptom_date + death +
                              event_type | EndpointDay0, data= selectedDataset_Clinical(), topclass="Rtable1-zebra")
  
  table1_produced
}, sanitize.text.function = function(x) x)





output$dayonetofourteen <- renderPlot({
  
  clin_group_0 <- selectedDataset_Clinical() %>% filter(event_type == "Visit 1") %>% group_by(EndpointDay0) %>% tally() %>% drop_na() %>%
    ggplot( aes(x= EndpointDay0, y=n, fill= EndpointDay0) ) + 
    geom_bar(stat="identity") + 
    theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency") +
    scale_fill_manual(values = c("#00B81F","orange",  "darkred", "black")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1)) +
    ggtitle("Clinical Grouping: day 0")
  
  resp_group_0 <- selectedDataset_Clinical() %>% filter(event_type == "Visit 1") %>% group_by(respiratory_status) %>% tally() %>% drop_na() %>%
    ggplot( aes(x=respiratory_status, y=n, fill=as.factor(respiratory_status)) ) + 
    geom_bar(stat="identity") + 
    theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency") +
    scale_fill_manual(values = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026", "black"))+
    ggtitle("Respiratory Status: day 0")
  
  clin_group_14 <- selectedDataset_Clinical() %>% filter(event_type == "Visit 1") %>% group_by(EndpointDay14) %>% tally() %>% drop_na() %>%
    ggplot( aes(x=EndpointDay14, y=n, fill=EndpointDay14) ) + 
    geom_bar(stat="identity") + 
    theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency") +
    scale_fill_manual(values = c("#00B81F","orange",  "darkred", "black")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1))+
    ggtitle("Clinical Grouping: day 14")
  
  resp_group_14 <- selectedDataset_Clinical() %>% filter(event_type == "Visit 1") %>% group_by(respiratory_status_day14) %>% tally() %>% drop_na() %>%
    ggplot( aes(x=respiratory_status_day14, y=n, fill=as.factor(respiratory_status_day14)) ) + 
    geom_bar(stat="identity") + 
    theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency") +
    scale_fill_manual(values = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026", "black"))+
    ggtitle("Respiratory Status: day 14")
  
  ggarrange(resp_group_0,clin_group_0,
            resp_group_14, clin_group_14,
            ncol = 2,
            nrow = 2,
            align = "hv")
  
})

output$alluvianPlot <- renderPlot({
  
  selectedDataset_Clinical() %>% filter(event_type == "Visit 1") %>% group_by(EndpointDay0, EndpointDay14) %>% tally() %>% drop_na() %>%
    ggplot(aes(y = n, axis1 = EndpointDay0, axis2 = EndpointDay14)) +
    geom_alluvium(aes(fill = EndpointDay14), width = 1/12) +
    geom_stratum() +
    scale_x_discrete(limits = c("Day 0", "Day 14"), expand = c(.05, .05)) +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
    theme_light() +
    scale_fill_manual(values = c("#00B81F","orange",  "darkred", "black")) +
    ylab("Number of Patients") +
    ggtitle("Patient Clinical Grouping Development") +
    theme(legend.position = "none")
  
})

output$ClinicalPlot_Bar <- renderPlot({
  selectedDataset_Clinical() %>% group_by(event_type) %>% tally() %>%
    ggplot(aes(x = event_type, y = n, fill = event_type)) +
    geom_bar(sta = "identity") +
    theme_light() + 
    geom_text(aes(label= n), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(legend.position = "none") + 
    ylab("Number of Samples") +
    xlab("") 
})

output$AgeHistogram <- renderPlot({
  
  hist_unsplit <- selectedDataset_Clinical() %>% ggplot(aes(x = admit_age)) +
    geom_histogram(aes(y=..density..), bins = input$HistBins, fill = "darkorange1") +
    geom_density(alpha = .2) +
    theme_light() + xlab("")
  
  hist_split <- selectedDataset_Clinical() %>% ggplot(aes(x = admit_age, fill = sex)) +
    geom_histogram(aes(y=..density..), bins = input$HistBins) +
    geom_density(alpha = .2) +
    theme_light() +
    scale_fill_manual( values = color_sex ) + 
    facet_grid(cols = vars(sex)) +
    theme(legend.position = "none")
  
  ggarrange(hist_unsplit, hist_split,
            ncol = 1,
            nrow = 2)
})

output$ClinicalSex <- renderPlot({
  selectedDataset_Clinical() %>% count(sex) %>% 
    ggplot( aes(x=sex, y=n, fill=sex) ) + 
    geom_bar(stat="identity") + 
    scale_fill_manual( values = color_sex ) + theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency")
  
})

output$SampleLocation <- renderPlot({
  selectedDataset_Clinical() %>% count(enrollment_site) %>% 
    ggplot( aes(x=enrollment_site, y=n, fill=enrollment_site) ) + 
    geom_bar(stat="identity") + theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency") +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1))
  
})

output$ClinicalRace <- renderPlot({
  selectedDataset_Clinical() %>% count(race) %>% 
    ggplot( aes(x=race, y=n, fill=race) ) + 
    geom_bar(stat="identity") + 
    scale_fill_manual( values = color_race) + theme_light() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.999, hjust=1))
  
})
