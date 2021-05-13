#Standard package
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(lme4)
library(broom)

#packages for Venny Diagram
library(eulerr)

#library for Alluvian Plot
library(ggalluvial)

#library for tables
library(DT)
library(table1)

#library for STRING analysis
library(STRINGdb)

#library for heatmaps
library(ComplexHeatmap)
library(simplifyEnrichment)
#library(InteractiveComplexHeatmap)
library(magick)

####Import the DataFrames #####
DDA_Data <- readr::read_csv(file = "B:\\IMPACC_Shiny\\IMPACC_Shiny\\Data\\DDA_Data.csv")

DIA_Data <- readr::read_csv("B:\\IMPACC_Shiny\\IMPACC_Shiny\\Data\\DIA_Data.csv")

Targeted_Data <- readr::read_csv("B:\\IMPACC_Shiny\\IMPACC_Shiny\\Data\\Targeted_Data.csv")


####Import the Clinical Data ####
clinical_event_data <- readr::read_csv("B:\\IMPACC_Shiny\\IMPACC_Shiny\\Data\\clinical_event_data.csv")
clinical_individ_data <- readr::read_csv("B:\\IMPACC_Shiny\\IMPACC_Shiny\\Data\\clinical_individ_data.csv")
clinical_sample_data <- readr::read_csv("B:\\IMPACC_Shiny\\IMPACC_Shiny\\Data\\clinical_sample_data.csv")


#Import StringDB
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=200, input_directory="")


#### combine them. #####
clinical_data <-inner_join(clinical_sample_data, clinical_event_data, by = c("event_id" = "event_id")) %>%    # Join the clinical_sample_data and clinical_event_data tables using "event_id" column
  mutate(
    participant_id = participant_id.x,
    participant_id.x = NULL,
    participant_id.y = NULL
  ) %>%   # Remove redundant "participant_id" columns
  inner_join(clinical_individ_data,
             by = c("participant_id" = "participant_id"))   # Append the clinical_individ_data table

#Add some info to clinical data on Endpoint Grouping.
#add column with specific info on initial endpoint day 0 info (factor it)
clinical_data <- add_column(clinical_data, 'EndpointDay0' =
                              clinical_data$respiratory_status,
                            .after = "respiratory_status")

#Endpoint to string grouping
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 1] <- "Discharged(1-2)"
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 2] <- "Discharged(1-2)"
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 3] <- "Hospitalized(3-4)"
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 4] <- "Hospitalized(3-4)"
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 5] <- "Hospitalized(5-6)"
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 6] <- "Hospitalized(5-6)"
clinical_data$EndpointDay0[clinical_data$EndpointDay0 == 7] <- "Dead(7)"

#order factor wise
clinical_data$EndpointDay0 <- factor(clinical_data$EndpointDay0,
                                     levels = c("Discharged(1-2)",
                                                "Hospitalized(3-4)",
                                                "Hospitalized(5-6)",
                                                "Dead(7)"))

#Now for Day 14

#add column with specific info on initial endpoint day 14 info (factor it)
clinical_data <- add_column(clinical_data, 'EndpointDay14' =
                              clinical_data$respiratory_status_day14,
                            .after = "respiratory_status_day14")

#Endpoint to string grouping
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 1] <- "Discharged(1-2)"
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 2] <- "Discharged(1-2)"
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 3] <- "Hospitalized(3-4)"
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 4] <- "Hospitalized(3-4)"
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 5] <- "Hospitalized(5-6)"
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 6] <- "Hospitalized(5-6)"
clinical_data$EndpointDay14[clinical_data$EndpointDay14 == 7] <- "Dead(7)"

#order factor wise
clinical_data$EndpointDay14 <- factor(clinical_data$EndpointDay14,
                                      levels = c("Discharged(1-2)",
                                                 "Hospitalized(3-4)",
                                                 "Hospitalized(5-6)",
                                                 "Dead(7)"))

COVID19Negative <- clinical_data %>% filter(sample_type == "PROTEOMICS/MET") %>% filter(participant_type == "COVID-19 Negative")
#### ####

#Setup some color schemes
color_sex <- c("Female" = "#DC2543", "Male" = "#1E94A0")

color_event_type <- c("Escalation 1" = "#FFFFCC",
                      "Visit 1" = "#FEE692",
                      "Visit 2" = "#FEBF5A",
                      "Visit 3" = "#FD8D3C",
                      "Visit 4" = "#F33C25",
                      "Visit 5" = "#C90822",
                      "Visit 6" = "#800026"
                      )

color_respiratory_status <- c("2" = "#FFFFCC",
                              "3" = "#FED976",
                              "4" = "#FD8D3C",
                              "5" = "#E31A1C",
                              "6" = "#800026",
                              "7" = "black")

color_race <- c("American Indian / Alaska Native" = "#F0A0FF",
                "Asian" = "#0075DC",
                "Black / African American" = "#993F00",
                "Multiple" = "#4C005C",
                "Native Hawaiian / Pacific Islander" = "#191919",
                "Other / Declined" = "#005C31",
                "Unknown / Unavailable" = "#2BCE48",
                "White" = "#FFCC99"
                )


# Define server logic ----
server <- function(input, output, session) {
  
  source("SelectDataSet.R", local = TRUE)
  source("Clinical.R", local = TRUE)
  source("Clustering.R", local = TRUE)
  source("data_analysis_template_codebase.R", local = TRUE)
  source("LongitudinalProjectories.R", local = TRUE)
  source("Heatmap.R", local = TRUE)
  source("VolcanoPlotter.R", local = TRUE)
  
  
}