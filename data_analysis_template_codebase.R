##### Codebase for the Data analysis template
##### Authors: Ravi K. Patel, Jeremy Gygi



# A function to fetch "current" or "legacy" version of the data. The function will return data files that were current on the indicated date.
fetch_file_names <- function(DATA_VERSION, dir, pattern = ".Counts.csv") {
  # Fetch all files within the dir
  files <- list.files(dir, recursive = T)
  
  # Select the files associated with the selection version of the data.
  if (DATA_VERSION == "current") {
    files <- files[grep("current", files)]
    
  } else {
    # Will store the date to use in legacy_date
    legacy_date <- ""
    # Today's date
    today_date <- Sys.Date()
    
    files_legacy <- files[grep("legacy", files)]
    dates <-
      as.Date(str_match(files_legacy, "/(\\d{4}-\\d{2}-\\d{2})/")[, 2])
    unique_reverse_sorted_dates <-
      sort(unique(dates), decreasing = T)
    
    # If there is no data in legacy folder, use today's date as legacy_date
    if (length(unique_reverse_sorted_dates) == 0) {
      legacy_date = today_date
      
    } else {
      # Else, find an appropriate version date.
      
      user_version_date <- as.Date(DATA_VERSION)
      # If the user's version date is same as today's date or is more recent than the most recent folder in "legacy" folder, use the today's date as the legacy date.
      if (today_date == user_version_date |
          unique_reverse_sorted_dates[1] < user_version_date) {
        legacy_date <- today_date
        
      } else {
        # Else, find the data version from legacy folder that was current on the user's version date.
        for (i in 1:length(unique_reverse_sorted_dates)) {
          if (user_version_date <= unique_reverse_sorted_dates[i]) {
            legacy_date <- unique_reverse_sorted_dates[i]
            
          }
        }
      }
    }
    
    
    if (legacy_date == today_date) {
      files <- files[grep("current", files)]
      
    } else {
      files <- files_legacy[dates == legacy_date]
      
    }
    
    
  }
  
  # Select a file that matches the pattern.
  matched_file_name <- files[grepl(pattern, files)]
  if (length(matched_file_name) > 1) {
    warning("More than one files match \"",
            pattern  ,
            "\" in \"",
            dir ,
            "\". Skipping this dataset.\n")
    return()
  } else if (length(matched_file_name) == 0) {
    #cat("No files match \"", pattern  ,"\" in \"", dir , "\".\n", sep = "")
    return()
  }
  
  return(file.path(dir, matched_file_name))
}


# A function to generate a named vector of colors for a categorical variable using the levels_vector (a vector of unique categories) provided by the user.
get_categorical_color_vector <- function(levels_vector) {
  return(alphabet(length(levels_vector)) %>% setNames(levels_vector))
}


# A function to generate a named vector of colors for an ordinal variable using the ordered_levels_vector (a vector of unique categories arranged in an increasing order) provided by the user. The user can provide name of RColorBrewer sequential palette or just assign a value between 1 and 18 to palette_index to choose an RColorBrewer sequential palette at that index.
get_ordinal_color_vector <-
  function(ordered_levels_vector,
           rcolorbrewer_palette_name = NULL,
           palette_index = 1,
           bidirectional = FALSE) {
    # A list of RColorBrewer sequential palettes
    rcolorbrewer_palettes <-
      brewer.pal.info %>% filter(category == "seq") %>% row.names() %>% rev()
    
    # Select a palette that matches user's selection
    if (is.null(rcolorbrewer_palette_name) & palette_index == 1) {
      rcolorbrewer_palette_name <- rcolorbrewer_palettes[1]
    }
    if (palette_index != 1) {
      rcolorbrewer_palette_name <- rcolorbrewer_palettes[palette_index]
    }
    
    # Create a color vector
    palette <-
      colorRampPalette(brewer.pal(9, rcolorbrewer_palette_name))
    color_vector <-
      palette(length(ordered_levels_vector)) %>% setNames(ordered_levels_vector)
    
    # If one of the values is "NA", use grey color for the category.
    # color_vector["NA"] <- "grey"
    return(color_vector)
  }


# Generate and return a color-key (legend) displaying the provided color scheme.
display_color_scheme <- function(color_scheme_vector, title = "") {
  p <- data.frame(color_scheme_vector = color_scheme_vector) %>%
    ggplot(aes(
      x = row.names(.),
      y = 1:nrow(.),
      fill = row.names(.)
    )) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_scheme_vector) +
    guides(fill = guide_legend(title = title))
  return(get_legend(p))
}


read_data_files <- function(file, row.names = NULL, header = FALSE) {
  if(file.access(file, mode = 4) == -1) {
    warning(file, " does not have a read permission.")
    return(data.frame())
  }
  if (grepl(".csv$", file)) {
    x <- tryCatch({
      data <- read.csv(file, row.names = row.names, check.names = F, header=header)
    },
    error = function(e) {
      warning("Cannot read ", file,". Please make sure it is a comma-separated file.\n")
      return(data.frame())
    })
  } else if (grepl(".tsv$", file)) {
    x <- tryCatch({
      data <- read.table(file, row.names = row.names, check.names = F, header=header)
    },
    error = function(e) {
      warning("Cannot read ", file,". Please make sure it is a tab-separated file.\n")
      return(data.frame())
    })
  } else {
    warning("Unrecognized file extension for ", file, ". Currently support only .csv and .tsv files.\n")
    return(data.frame())
  }
  return(data)
}


### Assign elements of a list to variables
assign_env_vars <- function(env) {
  lnames <- names(list)
  if( ! is.null( lnames ) ) {
    for( n in lnames ) {
      assign( n, list[[n]] )
    }
  }
}


#### Prepare clinical data
prepare_clinical_data <- function(data_env, DATA_VERSION, KEEP_COVID19_POS) {
  
  clinical_dir <- data_env[[ "data_dirs" ]][ "clinical_dir" ]
  
  ### Fetch clinical data files
  clinical_sample_file <- fetch_file_names( DATA_VERSION, clinical_dir, "sample.csv$")
  clinical_event_file <- fetch_file_names( DATA_VERSION, clinical_dir, "event.csv$")
  clinical_individ_file <- fetch_file_names( DATA_VERSION, clinical_dir, "individ.csv$")
  
  if (is.null(clinical_sample_file) |
      is.null(clinical_event_file) | is.null(clinical_individ_file)) {
    stop(
      "Error: Either you don't have access to the clinical data or the required clinical data files don't exist for the version of data you have selected. Please make sure that *sample.csv, *event.csv and *.individ.csv files exist in \"",
      clinical_dir,
      "\".\nExiting...\n",
      sep = ""
    )
  }
  
  message("Using the following clinical data files:\n")
  message("clinical_sample_file <- ", clinical_sample_file, "\n")
  message("clinical_event_file <- ", clinical_event_file, "\n")
  message("clinical_individ_file <- ", clinical_individ_file, "\n")
  
  ### Load clinical data
  clinical_sample_data <- read_data_files( clinical_sample_file, header = TRUE )
  clinical_event_data <- read_data_files( clinical_event_file, header = TRUE )
  clinical_individ_data <- read_data_files( clinical_individ_file, header = TRUE )
  
  ### Combine clinical data
  clinical_data <-
    inner_join(clinical_sample_data,
               clinical_event_data,
               by = c("event_id" = "event_id")) %>%    # Join the clinical_sample_data and clinical_event_data tables using "event_id" column
    mutate(
      participant_id = participant_id.x,
      participant_id.x = NULL,
      participant_id.y = NULL
    ) %>%   # Remove redundant "participant_id" columns
    inner_join(clinical_individ_data,
               by = c("participant_id" = "participant_id"))   # Append the clinical_individ_data table
  
  ### Keep data from COVID-19 Positive individuals
  if(KEEP_COVID19_POS)
    clinical_data <- clinical_data %>% filter(grepl("COVID-19 Positive", participant_type))
  
  data_env[[ "clinical_data" ]] <- clinical_data
  data_env[[ "clinical_sample_file" ]] <- clinical_sample_file
  data_env[[ "clinical_event_file" ]] <- clinical_event_file
  data_env[[ "clinical_individ_file" ]] <- clinical_individ_file
  
  return(data_env)
}


#### Prepare info for the assay data files
prepare_assay_files_info <- function(data_env) {
  ## Individual assay data files are fetched for the selected version of data and stored in a list for easy handling of these large number of data files.
  
  data_dirs <- data_env[[ "data_dirs" ]]
  for( n in names(data_dirs) ) {
    assign( n, data_dirs[ n ] )
  }
  
  ### A list containing, for individual assay files, file-name-pattern, directory, variable name to load the data to and description of the data file/type.
  assay_files_info <- list(
    "plasma_proteomics_targeted" = list(
      dir = plasma_proteomics_targeted_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "plasma_proteomics_targeted_counts",
      count_desc = "Count data from targeted proteomics of plasma",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "plasma_proteomics_targeted_rowfeature",
      rowfeature_desc = "RowFeature data from targeted proteomics of plasma"
    ),
    "plasma_proteomics_global_dda" = list(
      dir = plasma_proteomics_global_dda_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "plasma_proteomics_global_dda_counts",
      count_desc = "Count data from global proteomics (DDA) of plasma",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "plasma_proteomics_global_dda_rowfeature",
      rowfeature_desc = "RowFeature data from global proteomics (DDA) of plasma"
    ),
    "plasma_proteomics_global_dia" = list(
      dir = plasma_proteomics_global_dia_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "plasma_proteomics_global_dia_counts",
      count_desc = "Count data from global proteomics (DIA) of plasma",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "plasma_proteomics_global_dia_rowfeature",
      rowfeature_desc = "RowFeature data from global proteomics (DIA) of plasma"
    ),
    "serum_proteomics_global" = list(
      dir = serum_proteomics_global_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "serum_proteomics_global_counts",
      count_desc = "Count data from global proteomics of serum",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "serum_proteomics_global_rowfeature",
      rowfeature_desc = "RowFeature data from global proteomics of serum"
    ),
    "serum_olink" = list(
      dir = serum_olink_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "serum_olink_counts",
      count_desc = "Count data from Olink assay of serum",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "serum_olink_rowfeature",
      rowfeature_desc = "RowFeature data from Olink assay of serum"
    ),
    "nasal_viralload" = list(
      dir = nasal_viralload_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "nasal_viralload_counts",
      count_desc = "Count data for nasal viral load",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "nasal_viralload_rowfeature",
      rowfeature_desc = "RowFeature data for nasal viral load"
    ),
    "serum_rbd_abtiters" = list(
      dir = serum_rbd_abtiters_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "serum_rbd_abtiters_counts",
      count_desc = "Count data for RBD ab-titer in serum",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "serum_rbd_abtiters_rowfeature",
      rowfeature_desc = "RowFeature data for RBD ab-titer in serum"
    ),
    "serum_sarscov2_abtiters" = list(
      dir = serum_sarscov2_abtiters_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "serum_sarscov2_abtiters_counts",
      count_desc = "Count data for SARS-CoV-2 ab-titer in serum",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "serum_sarscov2_abtiters_rowfeature",
      rowfeature_desc = "RowFeature data for SARS-CoV-2 ab-titer in serum"
    ),
    "plasma_metabolomics_global" = list(
      dir = plasma_metabolomics_global_dir,
      count_pattern = "Counts.*$",
      count_file_path = "",
      count_var_name = "plasma_metabolomics_global_counts",
      count_desc = "Count data for global metabolomics of plasma",
      rowfeature_pattern = "RowFeature.*$",
      rowfeature_file_path = "",
      rowfeature_var_name = "plasma_metabolomics_global_rowfeature",
      rowfeature_desc = "RowFeature data for global metabolomics of plasma"
    )
  )
  
  data_env[["assay_files_info"]] <- assay_files_info
  return(data_env)
}


#### Load assay datasets
load_assay_data <- function(data_env, DATA_VERSION) {
  ### Fetch path of files that match the indicated pattern in the provided directory for each assay file and load the data in the variable names indicated in var_name element of the assay_files_info list. Uses the version of data specified in DATA_VERSION.
  
  assay_files_info <- data_env[[ "assay_files_info" ]]
  clinical_data <- data_env[[ "clinical_data" ]]
  
  for (f in names(assay_files_info)) {
    f_info <- assay_files_info[[f]]
    
    # Loading count data
    f_info$count_file_path <-
      fetch_file_names(DATA_VERSION, f_info$dir, f_info$count_pattern)
    
    if (is.null(f_info$count_file_path)) {
      warning(
        "\"",
        f_info$count_desc,
        "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
      )
      
    } else {
      ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
      message("Loading ",
              f_info$count_desc,
              " from ",
              f_info$count_file_path,
              "\n")
      data <-
        read_data_files(f_info$count_file_path,
                        row.names = 1,
                        header = TRUE)
      data <-
        data %>% filter(row.names(.) %in% clinical_data$sample_id)
      data_env[[ f_info$count_var_name ]] <- data
      
    }
    
    
    
    # Loading rowfeature data
    if (!is.null(f_info$rowfeature_pattern)) {
      f_info$rowfeature_file_path <-
        fetch_file_names(DATA_VERSION, f_info$dir, f_info$rowfeature_pattern)
      
      if (is.null(f_info$rowfeature_file_path)) {
        warning(
          "\"",
          f_info$rowfeature_desc,
          "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
        )
        
      } else {
        ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
        message(
          "Loading ",
          f_info$rowfeature_desc,
          " from ",
          f_info$rowfeature_file_path,
          "\n"
        )
        data <-
          read_data_files(f_info$rowfeature_file_path,
                          row.names = 1,
                          header = TRUE)
        data_env[[ f_info$rowfeature_var_name ]] <- data
        
      }
    }
    
    
    assay_files_info[[f]] <- f_info
    
  }
  
  data_env[[ "assay_files_info" ]] <- assay_files_info
  
  return(data_env)
}


#### Prepare color schemes
prepare_color_schemes <- function(data_env) {
  
  # Load clinical_data table into the local environment
  clinical_data <- data_env[[ "clinical_data" ]]
  
  ### Color for NA values
  data_env[[ "color_na" ]] <- "grey"
  
  ### Color scheme for sex
  data_env[[ "color_sex" ]] <- c("Female" = "#DC2543", "Male" = "#1E94A0")
  
  ### Color scheme for sample_type
  sample_type_levels <- unique( clinical_data$sample_type )
  data_env[[ "color_sample_type" ]] <- get_categorical_color_vector( sample_type_levels )
  
  ### Color scheme for respiratory_status
  respiratory_status_levels <- unique( sort( clinical_data$respiratory_status ) )
  data_env[[ "color_respiratory_status" ]] <- get_ordinal_color_vector( respiratory_status_levels )
  
  ### Color scheme for event_type
  event_type_levels <- unique( sort( clinical_data$event_type ) )
  data_env[[ "color_event_type" ]] <- get_ordinal_color_vector( event_type_levels )
  
  ### Color scheme for race
  race_levels <- unique( sort( clinical_data$race ) )
  data_env[[ "color_race" ]] <- get_categorical_color_vector( race_levels )
  
  
  data_env[[ "color_schemes" ]] <- data.frame(
    object_name = c("color_sex", "color_sample_type", "color_respiratory_status", "color_event_type", "color_race", "color_na"),
    feature_id = c("sex", "sample_type", "respiratory_status", "event_type", "race", "NA"),
    feature_name = c("Sex", "Sample type", "Respiratory status", "Event type", "Race", "NA values")
  )
  
  return(data_env)
}


#### Generate and print sample plots
generate_sample_plots <- function() {
  
  # Distribution of samples by sex
  bar_sex <- clinical_data %>% count(sex) %>% 
    ggplot( aes(x="Sex", y=n, fill=sex) ) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual( values = color_sex ) + theme_classic() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency")
  
  
  # Distribution of samples by respiratory_status
  bar_respiratory_status <- clinical_data %>% count(respiratory_status) %>% 
    ggplot( aes(x="Respiratory status", y=n, fill = as.factor( respiratory_status ) ) ) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual( values = color_respiratory_status, na.value = color_na ) + theme_classic() + 
    guides( fill = guide_legend(title="Respiratory status") ) +
    xlab("") + ylab("Frequency")
  
  
  # Scatter plot of ITIH4 and SERPINA1 levels
  scatter_ITIH4_vs_SERPINA1 <- plasma_proteomics_targeted_counts %>% ggplot(aes(y = ITIH4, x = SERPINA1 ) ) + 
    geom_point(size = 3, shape = 21, alpha = 0.8) + 
    geom_smooth(method = "lm", alpha = .15) + 
    theme_classic() +
    xlab("SERPINA1 levels") + ylab("ITI45 levels") +
    ggtitle("ITIH4 levels vs. SERPINA1 levels")
  
  
  
  # Jitter plot of ITIH4 levels vs. respiratory_status colored by event_type
  df <- data.frame( 
    ITIH4 = plasma_proteomics_targeted_counts[,"ITIH4"], 
    respiratory_status = clinical_data[ match( row.names(plasma_proteomics_targeted_counts), clinical_data$sample_id ), "respiratory_status" ], 
    event_type = clinical_data[ match( row.names(plasma_proteomics_targeted_counts), clinical_data$sample_id ), "event_type" ] )
  
  jitter_ITIH4_vs_resp <- df %>% ggplot(aes(y = ITIH4, x = as.factor( respiratory_status ), fill = event_type ) ) + 
    geom_point(position = position_jitter(width = 0.25), size = 3, shape = 21, alpha = 0.8) + 
    scale_fill_manual( values = color_event_type ) + 
    theme_classic() +
    labs( fill = "Event type") +
    xlab("respiratory_status") + ylab("ITIH4 levels") +
    ggtitle("ITIH4 levels vs. respiratory_status")
  
  
  
  # Example of imputing missing values
  plasma_proteomics_targeted_counts_imputed = impute::impute.knn(as.matrix(plasma_proteomics_targeted_counts))
  
  na_vals <-
    is.na(
      plasma_proteomics_targeted_counts %>% pivot_longer(cols = colnames(.),  names_to = "features") %>% pull(value)
    )
  
  imputed_vals <-
    plasma_proteomics_targeted_counts_imputed$data %>% data.frame() %>% 
    pivot_longer(cols = colnames(.),  names_to = "features") %>% 
    filter(na_vals) %>% pull(value)
  
  unimputed_vals <-
    plasma_proteomics_targeted_counts %>% 
    pivot_longer(cols = colnames(.),  names_to = "features") %>% 
    filter(!na_vals) %>% pull(value)
  
  # Comparing the Non-values with the NA values after impution.
  hist_imputed_vals <- rbind(
    data.frame(label = "Non-NA values", vals = unimputed_vals),
    data.frame(label = "NA values", vals = imputed_vals)
  ) %>%
    ggplot(aes(x = vals, fill = label)) +
    geom_density(adjust = 1.5, alpha = .4) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_classic() +
    ggtitle("Distribution of Non-NA values and\n NA values after imputation") +
    labs(fill = "") +
    xlab("Expression values")
  
  
  print(ggarrange(bar_respiratory_status, jitter_ITIH4_vs_resp, scatter_ITIH4_vs_SERPINA1, hist_imputed_vals,
                  widths = c(1,1.5), labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2))
  
  
  # Generate plot:
  print(ggplot(clinical_data) +
          geom_line(aes(x = event_date, y = respiratory_status, group = participant_id),
                    alpha = .1) +
          geom_point(
            aes(
              x = event_date,
              y = respiratory_status,
              group = participant_id,
              fill = as.character(respiratory_status)
            ),
            size = 4,
            shape = 21
          ) +
          # add color scale:
          scale_fill_manual(values = color_respiratory_status, na.value = color_na) +
          ggtitle("Respiratory Status vs. Time") +
          theme_classic() +
          labs(fill = "respiratory_status"))
  
  
  # Generate heatmap:
  melted.df <-
    tidyr::pivot_longer(clinical_data, cols = respiratory_status)
  melted.df <-
    mutate(melted.df, respiratory_status = as.character(value))
  ggplot(melted.df) +
    geom_tile(aes(x = participant_id, y = event_type, fill = respiratory_status)) +
    # add color scale:
    scale_fill_manual(values = color_respiratory_status, na.value = color_na) +
    ggtitle("respiratory_status heatmap") +
    theme_bw() +
    theme(axis.text.x = element_blank())
  
}




#### Use this function to load all available datasets and return an R environment containing all data objects.
load_IMPACC_datasets <- function(DATA_VERSION, KEEP_COVID19_POS) {
  
  
  #### Define Data directories: Use the absolute path of the directory containing "current" and "legacy" directories.
  data_dirs <- c(
    bld_cytof_dir = "/data/bld-cytof",
    bld_gwas_dir = "/data/bld-gwas",
    clinical_dir = "/data/clinical",
    ea_cytof_dir = "/data/ea-cytof",
    ea_metagenomics_dir = "/data/ea-metagenomics",
    ea_transcriptomics_dir = "/data/ea-transcriptomics",
    plasma_metabolomics_global_dir = "/data/metabolomics/plasma-metabolomics-global",
    plasma_metabolomics_targeted_dir = "/data/metabolomics/plasma-metabolomics-targeted",
    serum_metabolomics_global_dir = "/data/metabolomics/serum-metabolomics-global",
    nasal_metagenomics_dir = "/data/nasal-metagenomics",
    nasal_transcriptomics_dir = "/data/nasal-transcriptomics",
    nasal_viralload_dir = "/data/nasal-viralload",
    nasal_viralseq_dir = "/data/nasal-viralseq",
    pbmc_transcriptomics_dir = "/data/pbmc-transcriptomics",
    plasma_proteomics_targeted_dir = "/data/proteomics/plasma-proteomics-targeted",
    plasma_proteomics_global_dda_dir = "/data/proteomics/plasma-proteomics-global-DDA",
    plasma_proteomics_global_dia_dir = "/data/proteomics/plasma-proteomics-global-DIA",
    serum_proteomics_global_dir = "/data/proteomics/serum-proteomics-global",
    serum_olink_dir = "/data/serum-olink",
    serum_rbd_abtiters_dir = "/data/serum-rbd-abtiters",
    serum_sarscov2_abtiters_dir = "/data/serum-sarscov2-abtiters"
  )
  
  #### Create an empty environment for storing various data objects
  data_env <- env(data_dirs = data_dirs)
  
  #### Prepare clinical data variables
  data_env <- prepare_clinical_data(data_env, DATA_VERSION, KEEP_COVID19_POS)
  
  #### Prepare info for the assay data files
  data_env <- prepare_assay_files_info(data_env)
  
  #### Fetch appropriate data files and load assay datasets
  data_env <- load_assay_data(data_env, DATA_VERSION)
  
  #### Prepare predefined color schemes
  data_env <- prepare_color_schemes(data_env)
  
  return(data_env)
  
}

