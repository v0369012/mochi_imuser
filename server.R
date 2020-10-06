library(tidyverse)
# install.packages("remotes")
# remotes::install_github("jbisanz/qiime2R")
library(qiime2R)
library(plyr)
library(kableExtra)
library(shinyFiles)
library(myip)
library(shiny)
library(shinybusy)
library(vegan)


# library(reticulate)

server <- function(session, input, output) {
  
  
  observe({
    system("sudo service nginx start")
    system("sudo chmod -R 777 /var/www/html")
  })
  
  # Home page ----------------------------------------------------------------------------------------------------- 
  output$test <- renderUI({
    source("/srv/shiny-server/ui.R")
    # return(parseDirPath(roots = c(raw_data ="/home/imuser/raw_data"), selection = input$dirs))
    Sys.setenv(PATH='/usr/lib/rstudio-server/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/imuser/miniconda3/bin:/home/imuser/miniconda3/envs/qiime2-2019.10/bin')
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    Sys.setenv(LANG="C.UTF-8")
    HTML(
      paste0(
        "<h3 style='font-weight:bolder;'>MOCHI (MicrobiOme CHaracterization Image)</h3>",
        "<p>Thanks for using MOCHI.</p>",
        "<b>MOCHI</b><span> is a 16S rRNA analytical tool for microbiota. We simplified the operation of </span>",
        "<a href='https://qiime2.org/' target='_blank'>QIIME2</a>",
        "<span> by presenting it in a web interfaced power by the R package ",
        "<a href='https://shiny.rstudio.com/' target='_blank'>Shiny</a>",
        "<span>. All you need is a simple point and click with the mouse and filling in the required parameters for obtaining beautiful plots.</span>",
        "<br>",
        "<a href='https://i.imgur.com/nVZEArv.png' target='_blank'><img src='http://",
        my_qiime_ip_port,
        "/pipeline.png' alt='pipeline' style='width:auto;height:auto;max-width:75%;margin:10px;'></a>",
        "<hr>",
        "<h4 style='font-weight:bolder;'>The advantages of MOCHI</h4>",
        "<ul style='list-style-type:decimal;'>
                    <li>Friendly user interface (UI) design: analysis with just point and click and parameter fill-ins, no programming required</li>
                    <li>Cross platform: simple set up with <a href='https://www.docker.com/' target='_blank'>Docker</a> containers which runs on Linux, Windows, or macOS.</li>
                    <li>Local computing resource: runs on your premise where your local capability is your limit, not subject to limits in the reliability and availability of computing resources over the networks</li>
                    <li>Integration of third party software: MOCHI integrates several tools such as <a href='https://github.com/marbl/Krona/wiki' target='_blank'>Krona</a>、<a href='https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php' target='_blank'>FAPROTAX</a></li>
                 </ul>"
      )
      
    )
    
    # dir_list <- list()
    # dir_list[[1]] <- is_empty(parseDirPath(roots = c(raw_data ="/home/imuser/raw_data"), selection = input$dirs))
    # return(dir_list)
    # return(my_qiime_ip)
    # return(getwd())
    # return(system("whoami", intern = T))
    # return(system("echo ~", intern = T))
    # return(system("python /home/imuser/test.py", intern = T))
    # return(system("locale -a", intern = T))
    # use_virtualenv("qiime2-2019.10")
    # return(Sys.getenv("PATH"))
    # return(system("echo $PATH", intern = T))
    # return(system("which qiime", intern = T))
    # return(system(paste(qiime_cmd, "--version"), intern = T))
    # return(system("echo $LANG", intern = T))
    # print(as.character(file.exists("/home/imuser/www/Krona_rawdata.html")))
    
  })
  
  
  shinyDirChoose(input,
                 id = 'dirs',
                 roots = c(raw_data ="/home/imuser/raw_data"))       # Browse server side dirs
  
  
  
  
  # reactive object---------------------------------------------------------------------------------------------
  ## common objects
  
  my_IP <- reactive(my_qiime_local_ip) # get the ip of local
  
  my_qiime_port <- reactive(":8099") # give the port for this tool
  
  TaxaTable <- reactive({
    
    req(input$taxonomic_table)
    req(input$sample_data)
    infile1 <- input$taxonomic_table
    
    if (is.null(infile1))
      return(NULL)
    
    validate(
      need(infile1 != "", message = F),
      need(is.matrix(read_qza(infile1$datapath)$data), message = "")
    )
    
    read_qza(input$taxonomic_table$datapath)$data
    
  }) # read the input file (.qza)
  
  TaxaTable_merge <- reactive({
    
    as_output_taxtable <- function(df_data){
      
      df_data_rownames<-row.names(df_data)
      
      df_data <- cbind(Species=df_data_rownames, df_data) %>% as.data.frame()
      
      row.names(df_data) <- NULL
      
      # remove the non-sense string
      df_data$Species<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", "", df_data$Species)
      df_data$Species<-gsub("k__|p__|c__|o__|f__|g__|s__", "", df_data$Species)
      
      
      library(tidyr)
      df_data <- df_data %>%
        separate(
          col = Species,
          into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
          sep = ";"
        )
      
      # replace "__" to "Unassigned"
      df_data<-replace(df_data, df_data=="", "Unassigned")
      df_data<-replace(df_data, df_data=="__", "Unassigned")
      
      return(as.data.frame(df_data))
    } # clean the taxatable 
    
    taxatable <- as_output_taxtable(TaxaTable())
    
    taxatable$Level <- paste(taxatable$Kingdom,
                             taxatable$Phylum,
                             taxatable$Class,
                             taxatable$Order,
                             taxatable$Family,
                             taxatable$Genus,
                             taxatable$Species,
                             sep = ";")
    
    taxatable_merge <- taxatable[, -(1:7)]
    taxatable_merge[,-ncol(taxatable_merge)] <- data.frame(sapply(taxatable_merge[,-ncol(taxatable_merge)], function(x){
      as.character(x) %>% as.numeric()
    }))
    
    taxatable_merge_result <- taxatable_merge %>% group_by(Level) %>% summarise_all(funs(sum)) %>% as.data.frame()
    rownames(taxatable_merge_result) <- taxatable_merge_result$Level
    taxatable_merge_result <- taxatable_merge_result[, -1]
    
    return(taxatable_merge_result)
  }) # Some OTUs may be unassigned or multi-matched, in this object, both are considered as unassigned
  
  Metadata <- reactive({
    
    infile2 <- input$sample_data
    
    if (is.null(infile2))
      return(NULL)
    
    validate(
      need(infile2 != "", message = F)
    )
    
    if(is.data.frame(try(read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t"), silent = T))){
      
      metadata <- read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t")
      colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "-", ".")
      
      col_vector <- apply(metadata[,1:ncol(metadata)], MARGIN = 2, FUN = function(x){length(unique(x))})
      
      if( 1 %in% col_vector){
        
        position_non1 <- which(col_vector!=1)
        
        metadata <- metadata[,position_non1]
        
      }
      
      colnames(metadata)[1] <- "SampleID"
      
      for (i in 1:ncol(metadata)) {
        
        metadata[,i] <- as.character(metadata[,i])
        metadata[,i] <- replace_na(metadata[,i], "NA")
        
      }
      
      
      
      
      return(metadata)
    }
    
    
  }) # read the input sample data file
  
  Metadata_FA <- reactive({
    
    infile2 <- input$sample_data_FA
    
    validate(
      need(infile2 != "", message = F)
    )
    
    if(is.data.frame(try(read.table(input$sample_data_FA$datapath, header = T, na.strings = "", sep = "\t"), silent = T))){
      
      metadata <- read.table(input$sample_data_FA$datapath, header = T, na.strings = "", sep = "\t")
      colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "-", ".")
      
      col_vector <- apply(metadata[,1:ncol(metadata)], MARGIN = 2, FUN = function(x){length(unique(x))})
      
      if( 1 %in% col_vector){
        
        position_non1 <- which(col_vector!=1)
        
        metadata <- metadata[,position_non1]
        
      }
      
      colnames(metadata)[1] <- "SampleID"
      
      for (i in 1:ncol(metadata)) {
        
        metadata[,i] <- as.character(metadata[,i])
        metadata[,i] <- replace_na(metadata[,i], "NA")
        
      }
      
      return(metadata)
    }
    
    
    
  }) # read the input sample data file
  
  Metadata_stats <- reactive({
    
    infile2 <- input$sample_data
    
    if (is.null(infile2))
      return(NULL)
    
    validate(
      need(infile2 != "", message = F)
    )
    
    if(is.data.frame(try(read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t"), silent = T))){
      
      metadata <- read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t")
      colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "-", ".")
      
      col_vector <- apply(metadata[,1:ncol(metadata)], MARGIN = 2, FUN = function(x){length(unique(x))})
      
      if( 1 %in% col_vector){
        
        position_non1 <- which(col_vector!=1)
        
        metadata <- metadata[,position_non1]
        
      }
      
      colnames(metadata)[1] <- "SampleID"
      
      for (i in 1:ncol(metadata)) {
        
        metadata[,i] <- as.character(metadata[,i])
        metadata[,i] <- replace_na(metadata[,i], "NA")
        
      }
      
      OKstats_col <- lapply(colnames(metadata)[-1], function(x){
        if(sum(table(metadata[, x])<=1)==0){
          return(x)
        }
      }) %>% unlist()
      
      metadata <- metadata[, c("SampleID", OKstats_col)]
      
      return(metadata)
    }
  })
  
  output$word_metadata_NA_1 <- renderText({
    
    if("NA" %in% Metadata_stats()[,input$metadata_barplot]){
      NA_position <- which(Metadata()[,input$metadata_barplot]=="NA")
      NA_sampleid <- Metadata()[,1][NA_position]
      print(paste("Your sample id", NA_sampleid, "has missing value and therefore could be removed."))
    }
  })
  
  output$word_metadata_NA_2 <- renderText({
    
    if("NA" %in% Metadata()[,input$metadata_hm]){
      NA_position <- which(Metadata()[,input$metadata_hm]=="NA")
      NA_sampleid <- Metadata()[,1][NA_position]
      print(paste("Your sample id", NA_sampleid, "has missing value and therefore could be removed."))
    }
  })
  
  output$word_metadata_NA_3 <- renderText({
    
    if("NA" %in% Metadata()[,input$metadata_alpha]){
      NA_position <- which(Metadata()[,input$metadata_alpha]=="NA")
      NA_sampleid <- Metadata()[,1][NA_position]
      print(paste("Your sample id", NA_sampleid, "has missing value and therefore could be removed."))
    }
  })
  
  # output$word_metadata_NA_4 <- renderText({
  #   
  #   if("NA" %in% Metadata()[,input$metadata_beta]){
  #     NA_position <- which(Metadata()[,input$metadata_beta]=="NA")
  #     NA_sampleid <- Metadata()[,1][NA_position]
  #     print(paste("Your sample id", NA_sampleid, "has missing value and therefore could be removed."))
  #   }
  # })
  
  
  # Button change by the metadata -------------------------------------------------------------------------------------
  
  observe({
    
    newchoices <- colnames(Metadata())
    newchoices_stats <- colnames(Metadata_stats())
    # newchoices_FA <- colnames(Metadata_FA())
    updateRadioButtons(session, "metadata1", choices = newchoices, inline = T)
    updateRadioButtons(session, "metadata_alpha", choices = newchoices_stats[-1], inline = T)
    updateRadioButtons(session, "metadata_beta", choices = newchoices_stats[-1], inline = T)
    updateRadioButtons(session, "metadata_hm", choices = newchoices_stats, inline = T)
    updateRadioButtons(session, "metadata_barplot", choices = newchoices_stats, inline = T)
    updateRadioButtons(session, "metadata_phylo_alpha", choices = newchoices_stats[-1], inline = T)
    updateRadioButtons(session, "metadata_phylo_beta", choices = newchoices_stats[-1], inline = T)
    # updateSelectInput(session, "metadata8", choices = newchoices[-(1:2)])
    updateRadioButtons(session, "metadata_ANCOM", choices = newchoices_stats[-1], inline = T)
    # updateRadioButtons(session, "metadata10", choices = newchoices_FA[-(1:2)], inline = T)
    
  })
  
  # Check the sample data input
  observe({
    
    req(input$sample_data)
    
    library(shinyBS)
    
    if(is.data.frame(try(read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t" ), silent = T))){
      
      sample_data <- read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t" )
      # sample_data <- readr::read_tsv(input$sample_data$datapath, col_names = T, comment = "#", na = "NA")
      '%!in%' <- function(x,y)!('%in%'(x,y))
      
      if(colnames(sample_data[1]) %!in% c("sample.id", "sample-id", "SampleID")) {
        createAlert(session, 
                    anchorId = "sample_data_alert", 
                    alertId = "sampleAlert", 
                    title = "Oops!",
                    content = "Please check your sample data. Your first column names should be 'sample.id' or 'sample-id' or 'SampleID'.", 
                    append = T,
                    style = "danger")
      } else {
        closeAlert(session, "sampleAlert")
        
      }
      
    }else{
      createAlert(session, 
                  anchorId = "sample_data_alert", 
                  alertId = "sampleAlert", 
                  title = "Oops!",
                  content = "Please check your sample data format.", 
                  append = T, 
                  style = "danger")
    }
    
    
  })
  
  # give the problem message
  observe({
    
    req(input$sample_data)
    library(shinyBS)
    
    if(is.data.frame(try(read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t"), silent = T))){
      
      metadata <- read.table(input$sample_data$datapath, header = T, na.strings = "", sep = "\t" )
      colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "-", ".")
      
      col_vector <- apply(metadata[,2:ncol(metadata)], MARGIN = 2, FUN = function(x){length(unique(x))})
      
      position_1 <- which(col_vector==1)
      
      
      noOKstats_col <- lapply(colnames(metadata)[-1], function(x){
        if(sum(table(metadata[, x])<=1)!=0){
          return(x)
        }
      }) %>% unlist()
      
      noOK_col <- c(names(col_vector)[position_1], noOKstats_col) %>% paste(collapse = ", ")
      
      word_content <- paste("Your column '", noOK_col, "' would be ignored because don't have more than 2 categories or sample size of every category is less than 2.")
      
      
      if(length(position_1)!=0 | length(noOK_col)!=0) {
        createAlert(session, 
                    anchorId = "sample_data_alert", 
                    alertId = "sampleAlert", 
                    title = "Warning!",
                    content = word_content, 
                    append = FALSE,
                    style = "warning")
        
      } else {
        closeAlert(session, "sampleAlert")
        
      }
      
    }
    
    
    
  })
  
  
  # Check the taxonomic table input
  observe({
    
    req(input$taxonomic_table)
    
    
    library(shinyBS)
    
    if(!is.matrix(read_qza(input$taxonomic_table$datapath)$data)) {
      createAlert(session, 
                  anchorId = "taxatable_alert", 
                  alertId = "taxaAlert", 
                  title = "Oops!",
                  content = "Please check your input taxonomic table file.", 
                  append = T,
                  style = "danger")
    } else {
      closeAlert(session, "taxaAlert")
      
    }
  })
  
  # Check the seq input
  observe({
    
    req(input$rep_seq__dada2_upload)
    
    library(shinyBS)
    
    if(class(read_qza(input$rep_seq__dada2_upload$datapath)$data)!="DNAStringSet") {
      createAlert(session, 
                  anchorId = "seq_alert", 
                  alertId = "seqAlert", 
                  title = "Oops!",
                  content = "Please check your input sequences file.", 
                  append = T,
                  style = "danger")
    } else {
      closeAlert(session, "seqAlert")
      
    }
  })
  
  observe({
    newchoices_FA <- colnames(Metadata_FA())
    updateRadioButtons(session, "metadata10", choices = newchoices_FA[-1], inline = T)
  })
  
  observe({
    
    as_output_taxtable <- function(df_data){
      
      df_data_rownames<-row.names(df_data)
      
      df_data <- cbind(Species=df_data_rownames, df_data) %>% as.data.frame()
      
      row.names(df_data) <- NULL
      
      # remove the non-sense string
      df_data$Species<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", "", df_data$Species)
      df_data$Species<-gsub("k__|p__|c__|o__|f__|g__|s__", "", df_data$Species)
      
      
      library(tidyr)
      df_data <- df_data %>%
        separate(
          col = Species,
          into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
          sep = ";"
        )
      
      # replace "__" to "Unassigned"
      df_data<-replace(df_data, df_data=="", "Unassigned")
      df_data<-replace(df_data, df_data=="__", "Unassigned")
      
      return(as.data.frame(df_data))
    } # clean the taxatable 
    
    taxatable <- as_output_taxtable(TaxaTable())
    
    val <- length(unique(taxatable[, input$select2]))
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "integer", value = val,
                      min = 1, max = val, step = 1)
  })
  
  observe({
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    min_sampling_depth <- colSums(TaxaTable()) %>% min()
    
    # This will change the value of input$inText, based on x
    updateTextInput(session, "sampling_depth", value = min_sampling_depth)
    
  })
  
  
  observeEvent(req(input$sample_data, input$taxonomic_table), {
    
    if(is.matrix(read_qza(input$taxonomic_table$datapath)$data)) {
      
      shinyjs::toggle("taxatable_ui")
      
      shinyjs::toggle("taxabarplot_ui")
      
      shinyjs::toggle("taxaheatmap_ui") 
      
      shinyjs::toggle("krona_ui")
      
      shinyjs::toggle("alpha_ui")
      
      shinyjs::toggle("beta_ui")
      
      shinyjs::toggle("phylo_ui")
      
      shinyjs::toggle("ancom_ui")
      
    }else{
      shinyjs::hide("taxatable_ui")
      
      shinyjs::hide("taxabarplot_ui")
      
      shinyjs::hide("taxaheatmap_ui") 
      
      shinyjs::hide("krona_ui")
      
      shinyjs::hide("alpha_ui")
      
      shinyjs::hide("beta_ui")
      
      shinyjs::hide("phylo_ui")
      
      shinyjs::hide("ancom_ui")
    }
    
  })
  
  observeEvent(input$ANCOM_start, {
    shinyjs::toggle("ancom_output_ui")
    shinyjs::toggle("ancom_table_download")
  })
  
  observeEvent(req(input$phylogenetic_tree, input$table_dada2_upload, input$rep_seq_dada2_upload), {
    
    if(file.exists("/home/imuser/qiime_output/core-metrics-results/faith_pd_vector.qza")){
      shinyjs::toggle("phylo_output_ui")
    }
    
    
  })
  
  # functional analysis input file
  
  ## check sample data input
  observe({
    
    req(input$sample_data_FA)
    
    library(shinyBS)
    
    if(is.data.frame(try(read.table(input$sample_data_FA$datapath, header = T, na.strings = "", sep = "\t" ), silent = T))){
      
      sample_data <- read.table(input$sample_data_FA$datapath, header = T, na.strings = "", sep = "\t" )
      # sample_data <- readr::read_tsv(input$sample_data$datapath, col_names = T, comment = "#", na = "NA")
      '%!in%' <- function(x,y)!('%in%'(x,y))
      
      if(colnames(sample_data[1]) %!in% c("sample.id", "sample-id")) {
        createAlert(session, 
                    anchorId = "sample_data_alert_FA", 
                    alertId = "sampleAlert_FA", 
                    title = "Oops!",
                    content = "Please check your sample data. Your first column names should be 'sample.id' or 'sample-id'.", 
                    append = T,
                    style = "danger")
      } else {
        closeAlert(session, "sampleAlert_FA")
        
      }
      
    }else{
      createAlert(session, 
                  anchorId = "sample_data_alert_FA", 
                  alertId = "sampleAlert_FA", 
                  title = "Oops!",
                  content = "Please check your sample data format.", 
                  append = T, 
                  style = "danger")
    }
    
    
  })
  
  ## give the problem message
  observe({
    
    req(input$sample_data_FA)
    library(shinyBS)
    
    if(is.data.frame(try(read.table(input$sample_data_FA$datapath, header = T, na.strings = "", sep = "\t" ), silent = T))){
      
      metadata <- read.table(input$sample_data_FA$datapath, header = T, na.strings = "", sep = "\t" )
      colnames(metadata) <- stringr::str_replace_all(colnames(metadata), "-", ".")
      
      col_vector <- apply(metadata[,1:ncol(metadata)], MARGIN = 2, FUN = function(x){length(unique(x))})
      
      position_1 <- which(col_vector==1)
      
      
      noOKstats_col <- lapply(colnames(metadata)[-1], function(x){
        if(sum(table(metadata[, x])<=1)!=0){
          return(x)
        }
      }) %>% unlist()
      
      noOK_col <- c(names(col_vector)[position_1], noOKstats_col) %>% paste(collapse = ", ")
      
      word_content <- paste("Your column '", noOK_col, "' would be ignored because don't have more than 2 categories or sample size of every category is less than 2.")
      
      if(length(position_1)!=0 | length(noOK_col)!=0) {
        createAlert(session, 
                    anchorId = "sample_data_alert_FA", 
                    alertId = "sampleAlert_FA", 
                    title = "Warning!",
                    content = word_content, 
                    append = FALSE,
                    style = "warning")
        
      } else {
        closeAlert(session, "sampleAlert_FA")
        
      }
      
    }
    
    
    
  })
  
  ## check taxatable input
  observe({
    
    req(input$taxonomic_table_FA)
    
    
    library(shinyBS)
    
    if(!is.matrix(read_qza(input$taxonomic_table_FA$datapath)$data)) {
      createAlert(session, 
                  anchorId = "taxatable_alert_FA", 
                  alertId = "taxaAlert_FA", 
                  title = "Oops!",
                  content = "Please check your input taxonomic table file.", 
                  append = T,
                  style = "danger")
    } else {
      closeAlert(session, "taxaAlert_FA")
      
    }
  })
  
  
  observeEvent(req(input$sample_data_FA, input$taxonomic_table_FA, input$function_analysis), {
    
    shinyjs::toggle("func_table_ID")
  })
  
  
  # Sequences preprocessing -----
  
  # Demultiplexed (single)-------------------------------------------------------------------------------------------------------
  
  sample_file_name <- reactive({
    
    inFile <- input$sample_data
    
    if (is.null(inFile))
      return(NULL)
    
    return(inFile$name)
  })
  
  observeEvent(input$demultiplexed_single_ends, {
    
    start_time <- Sys.time()
    
    Sys.setenv(LANG="C.UTF-8")
    
    # showModal(modalDialog(title = "Running demultiplexed for single end ...", "Waiting for a moment", footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    raw_data_path_list <- list()
    raw_data_path_list[[1]] <- parseDirPath(roots = c(raw_data ="/home/imuser/raw_data"), selection = input$dirs)
    
    Sys.setenv(PATH='/usr/lib/rstudio-server/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/imuser/miniconda3/bin:/home/imuser/miniconda3/envs/qiime2-2019.10/bin')
    file.remove("/home/imuser/qiime_output/demux_single_end.qza", "/home/imuser/qiime_output/demux_single_end.qzv")
    file.remove("", "/home/imuser/qiime_output/demux_single_end.qzv")
    
    system(paste(qiime_cmd, "tools import --type", "'SampleData[SequencesWithQuality]'", "--input-path", raw_data_path_list[[1]],
                 "--input-format 'CasavaOneEightSingleLanePerSampleDirFmt'" ,'--output-path /home/imuser/qiime_output/demux_single_end.qza'))
    
    system(paste(qiime_cmd, 'demux summarize --i-data /home/imuser/qiime_output/demux_single_end.qza --o-visualization /home/imuser/qiime_output/demux_single_end.qzv'))
    # viewer_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime_2_ll_quick_viewer'
    # system('kill -9 $(lsof -t -i:8080 -sTCP:LISTEN)')
    # system(paste(viewer_cmd, '--filename /home/imuser/qiime_output/demux_single_end.qzv &'))
    
    unlink("/home/imuser/qiime_output/demux_single_unzip/new_dirname", recursive = T)
    unlink("/var/www/html/demux_single_unzip/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/demux_single_end.qzv /home/imuser/qiime_output/demux_single_end.zip")
    system("unzip -d /home/imuser/qiime_output/demux_single_unzip /home/imuser/qiime_output/demux_single_end.qzv")
    unzip_dirnames <- list.files("/home/imuser/qiime_output/demux_single_unzip", full.names = T)
    system(paste0("mv ", unzip_dirnames, " /home/imuser/qiime_output/demux_single_unzip/new_dirname"))
    system("cp -rf /home/imuser/qiime_output/demux_single_unzip /var/www/html/")
    # unlink("/home/imuser/qiime_output/demux_single_unzip/new_dirname", recursive = T)
    
    
    # if(is_empty(raw_data_path_list)){
    #   output$word_demultiplexed_single_end <- renderText({
    #     print('Please select the directory.')
    #   })
    #   
    # }else if(file.exists('/home/imuser/qiime_output/demux_single_end.qzv')){
    #   output$word_demultiplexed_single_end <- renderText({
    #     print('Successfully')
    #   })
    #   
    # }else{
    #   output$word_demultiplexed_single_end <- renderText({
    #     print('Error')
    #   })
    # }
    
    # removeModal()
    remove_modal_spinner()
    
    end_time <- system.time()
    
    spent_time <- format(round(end_time-start_time, digits = 2))
    
    if(file.exists('/home/imuser/qiime_output/demux_single_end.qzv')){
      showModal(modalDialog(title = strong("Successful!"), 
                            HTML(
                              paste0(
                                "The process took ", spent_time, ". ",
                                "You can click the button ", strong('View') ," to inspect the result.")
                            ), 
                            footer = NULL, easyClose = T, size = "l"))
    }else{
      showModal(modalDialog(title = strong("Error!", style = "color: red"), 
                            "Please check your files.", 
                            footer = NULL, easyClose = T, size = "l"))
    }
    
    
  })
  
  # Demultiplexed (paired)-------------------------------------------------------------------------------------------------------
  
  observeEvent(input$demultiplexed_paired_ends, {
    
    start_time <- Sys.time()
    
    Sys.setenv(LANG="C.UTF-8")
    
    # showModal(modalDialog(title = "Running demultiplexed for paired end ...", "Waiting for a moment",footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    raw_data_path <- parseDirPath(roots = c(raw_data ="/home/imuser/raw_data"), selection = input$dirs)
    
    Sys.setenv(PATH='/usr/lib/rstudio-server/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/imuser/miniconda3/bin:/home/imuser/miniconda3/envs/qiime2-2019.10/bin')
    file.remove("/home/imuser/qiime_output/demux_paired_end.qza", "/home/imuser/qiime_output/demux_paired_end.qzv")
    
    system(paste(qiime_cmd, "tools import --type", "'SampleData[PairedEndSequencesWithQuality]'", "--input-path", raw_data_path,
                 "--input-format 'CasavaOneEightSingleLanePerSampleDirFmt'" ,'--output-path /home/imuser/qiime_output/demux_paired_end.qza'))
    
    system(paste(qiime_cmd, 'demux summarize --i-data /home/imuser/qiime_output/demux_paired_end.qza --o-visualization /home/imuser/qiime_output/demux_paired_end.qzv'))
    
    
    unlink("/home/imuser/qiime_output/demux_paired_unzip/new_dirname", recursive = T)
    unlink("/var/www/html/demux_paired_unzip/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/demux_paired_end.qzv /home/imuser/qiime_output/demux_paired_end.zip")
    system("unzip -d /home/imuser/qiime_output/demux_paired_unzip /home/imuser/qiime_output/demux_paired_end.qzv")
    unzip_dirnames <- list.files("/home/imuser/qiime_output/demux_paired_unzip", full.names = T)
    system(paste0("mv ", unzip_dirnames, " /home/imuser/qiime_output/demux_paired_unzip/new_dirname"))
    system("cp -rf /home/imuser/qiime_output/demux_paired_unzip /var/www/html/")
    
    
    # removeModal()
    remove_modal_spinner()
    
    end_time <- Sys.time()
    spent_time <- format(round(end_time-start_time, digits = 2))
    
    if(file.exists('/home/imuser/qiime_output/demux_paired_end.qzv')){
      showModal(modalDialog(title = strong("Successful!"), 
                            HTML(
                              paste0(
                                "The process took ", spent_time, ". ",
                                "You can click the button ", strong('View') ," to inspect the result.")
                            ), 
                            footer = NULL, easyClose = T, size = "l"))
      
    }else{
      showModal(modalDialog(title = strong("Error!", style = "color: red"), 
                            "Please check your files.", 
                            footer = NULL, easyClose = T, size = "l"))
    }
    
  })
  
  
  
  # Denoising_single -----------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$denoising_single, {
    
    start_time <- Sys.time()
    # showModal(modalDialog(title = "Running denoising for single end ...", "Waiting for a moment", footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    
    file.remove("/home/imuser/qiime_output/rep-seqs-dada2_single.qza", 
                "/home/imuser/qiime_output/table-dada2_single.qza",
                "/home/imuser/qiime_output/stats-dada2_single.qza",
                "/home/imuser/qiime_output/rep-seqs-dada2_single.qzv",
                "/home/imuser/qiime_output/table-dada2_single.qzv",
                "/home/imuser/qiime_output/stats-dada2_single.qzv")
    
    if(is.null(input$sample_data_single$datapath)){
      add_metadata_table <- ""
      add_metadata_rarefaction <- ""
    }else{
      add_metadata_table <- paste("--m-sample-metadata-file", input$sample_data_single$datapath)
      add_metadata_rarefaction <- paste("--m-metadata-file", input$sample_data_single$datapath)
    }
    
    system(paste(qiime_cmd, "dada2 denoise-single --i-demultiplexed-seqs /home/imuser/qiime_output/demux_single_end.qza", 
                 '--p-trim-left', input$trim_left_single, 
                 '--p-trunc-len', input$trunc_len_single,
                 "--p-trunc-q", input$qvalue_single,
                 "--p-n-threads", input$threads_single,
                 "--p-min-fold-parent-over-abundance", input$chimera_single,
                 "--p-n-reads-learn", input$n_reads_single,
                 "--o-representative-sequences /home/imuser/qiime_output/rep-seqs-dada2_single.qza",
                 '--o-table /home/imuser/qiime_output/table-dada2_single.qza', 
                 '--o-denoising-stats /home/imuser/qiime_output/stats-dada2_single.qza'))
    
    system(paste(qiime_cmd, 'metadata tabulate --m-input-file /home/imuser/qiime_output/stats-dada2_single.qza --o-visualization /home/imuser/qiime_output/stats-dada2_single.qzv'))
    # system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table-dada2.qza --o-visualization /home/imuser/qiime_output/table-dada2.qzv --m-sample-metadata-file', paste0("/home/imuser/qiime_output/",sample_file_name())))
    system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table-dada2_single.qza',
                 add_metadata_table,
                 '--o-visualization /home/imuser/qiime_output/table-dada2_single.qzv'
                 # , '--m-sample-metadata-file', input$sample_data_single$datapath
    )
    )
    system(paste(qiime_cmd, 'feature-table tabulate-seqs --i-data /home/imuser/qiime_output/rep-seqs-dada2_single.qza --o-visualization /home/imuser/qiime_output/rep-seqs-dada2_single.qzv'))
    
    unlink("/home/imuser/qiime_output/denoise_single_stats/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_single_stats/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/stats-dada2.qzv /home/imuser/qiime_output/stats-dada2.zip")
    system("unzip -d /home/imuser/qiime_output/denoise_single_stats /home/imuser/qiime_output/stats-dada2_single.qzv")
    unzip_dirnames_stats <- list.files("/home/imuser/qiime_output/denoise_single_stats", full.names = T)
    system(paste0("mv ", unzip_dirnames_stats, " /home/imuser/qiime_output/denoise_single_stats/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_single_stats/ /var/www/html/")
    
    unlink("/home/imuser/qiime_output/denoise_single_table/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_single_table/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/table-dada2.qzv /home/imuser/qiime_output/table-dada2.zip")
    system("unzip -d /home/imuser/qiime_output/denoise_position_table /home/imuser/qiime_output/table-dada2_single.qzv")
    unzip_dirnames_table <- list.files("/home/imuser/qiime_output/denoise_single_table", full.names = T)
    system(paste0("mv ", unzip_dirnames_table, " /home/imuser/qiime_output/denoise_table/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_single_table/ /var/www/html/")
    
    unlink("/home/imuser/qiime_output/denoise_single_seqs/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_single_seqs/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/rep-seqs-dada2.qzv /home/imuser/qiime_output/rep-seqs-dada2.zip")
    system("unzip -d /home/imuser/qiime_output/denoise_single_seqs /home/imuser/qiime_output/rep-seqs-dada2_single.qzv")
    unzip_dirnames_seqs <- list.files("/home/imuser/qiime_output/denoise_single_seqs", full.names = T)
    system(paste0("mv ", unzip_dirnames_seqs, " /home/imuser/qiime_output/denoise_single_seqs/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_single_seqs/ /var/www/html/")
    
    # alpha-rarefaction
    max_depth <- read_qza("/home/imuser/qiime_output/table-dada2.qza")[["data"]] %>% colSums() %>% min()
    system(paste(qiime_cmd, "phylogeny align-to-tree-mafft-fasttree --i-sequences /home/imuser/qiime_output/rep-seqs-dada2_single.qza", 
                 "--p-n-threads", input$threads_single,
                 "--o-alignment /home/imuser/qiime_output/aligned-rep-seqs-dada2_single.qza --o-masked-alignment /home/imuser/qiime_output/masked-aligned-rep-seqs-dada2_single.qza",
                 "--o-tree /home/imuser/qiime_output/unrooted-tree_single.qza --o-rooted-tree /home/imuser/qiime_output/rooted-tree_single.qza"))
    
    system(paste(qiime_cmd, 'diversity alpha-rarefaction --i-table /home/imuser/qiime_output/table-dada2_single.qza',
                 '--p-max-depth', max_depth,
                 '--i-phylogeny /home/imuser/qiime_output/rooted-tree_single.qza',
                 add_metadata_rarefaction,
                 '--o-visualization /home/imuser/qiime_output/rarefaction-dada2_single.qzv'))
    unlink("/home/imuser/qiime_output/denoise_single_rarefaction/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_single_rarefaction/new_dirname", recursive = T)
    system("unzip -d /home/imuser/qiime_output/denoise_single_rarefaction /home/imuser/qiime_output/rarefaction-dada2_single.qzv")
    unzip_dirnames_seqs <- list.files("/home/imuser/qiime_output/denoise_single_rarefaction", full.names = T)
    system(paste0("mv ", unzip_dirnames_seqs, " /home/imuser/qiime_output/denoise_single_rarefaction/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_single_rarefaction/ /var/www/html/")
    
    # removeModal()
    remove_modal_spinner()
    end_time <- Sys.time()
    spent_time <- format(round(end_time-start_time, digits = 2))
    
    if(all.equal(file.exists(c('/home/imuser/qiime_output/stats-dada2_single.qzv', 
                               '/home/imuser/qiime_output/table-dada2_single.qzv', 
                               '/home/imuser/qiime_output/rep-seqs-dada2_single.qzv',
                               '/home/imuser/qiime_output/rarefaction-dada2_single.qzv')), 
                 c(T, T, T, T))){
      showModal(modalDialog(title = strong("Denoising succecessfully!"), 
                            HTML(
                              paste0(
                                "The process took ", spent_time, ". ",
                                "You can inspect the results!")
                            ), 
                            footer = NULL, easyClose = T, size = "l"))
    }else{
      showModal(modalDialog(title = strong("Error!", style = "color: red"),
                            "Please check your files.", 
                            footer = NULL, easyClose = T, size = "l"))
    }
    
  })
  
  
  
  observeEvent(input$my_cores_single, {
    
    # output$message_mycore_single_position <- renderUI({
    my_cores <- parallel::detectCores()
    str1 <- paste0("The number of your threads is ", my_cores, ".")
    str2 <- " If 0 is provided, all available cores will be used. [default: all threads - 2]"
    #   HTML(paste(str1, str2, sep = "</br>"))
    # })
    
    showModal(modalDialog(title = strong("Message"),
                          HTML(paste(str1, str2, sep = "</br>")), 
                          footer = NULL, easyClose = T, size = "l"))
    
  })
  
  
  
  observeEvent(input$Q_cores_single, {
    # output$message_thread_single_position <- renderText(
    #   "The number of threads to use for multithreaded processing.")
    showModal(modalDialog(title = strong("Message"),
                          "The number of threads to use for multithreaded processing.", 
                          footer = NULL, easyClose = T, size = "l"))
  })
  
  observeEvent(input$word_chimera_single, {
    # output$message_thread_single_position <- renderText(
    #   "The number of threads to use for multithreaded processing.")
    showModal(modalDialog(title = strong("Message"),
                          "The minimum abundance of potential parents of a
                         sequence being tested as chimeric, expressed as a
                         fold-change versus the abundance of the sequence
                         being tested. Values should be greater than or equal
                         to 1 (i.e. parents should be more abundant than the
                         sequence being tested).", 
                          footer = NULL, easyClose = T, size = "l"))
  })
  
  observeEvent(input$Q_learn_reads_single, {
    # output$message_learn_single_position <- renderText(
    #   " The number of reads to use when training the error model. Smaller numbers will result in a shorter run time but a less reliable error model. [default: 1000000]")
    showModal(modalDialog(title = strong("Message"),
                          "The number of reads to use when training the error model. Smaller numbers will result in a shorter run time but a less reliable error model. [default: 1000000]", 
                          footer = NULL, easyClose = T, size = "l"))
  })
  
  
  # Denoising_paired -----------------------------------------------------------------------------------------------------------------------  
  observeEvent(input$denoising_paired, {
    
    start_time <- Sys.time()
    # showModal(modalDialog(title = "Running denoising for paired end ...", "Waiting for a moment",footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    
    file.remove("/home/imuser/qiime_output/rep-seqs-dada2_paired.qza", 
                "/home/imuser/qiime_output/table-dada2_paired.qza",
                "/home/imuser/qiime_output/stats-dada2_paired.qza",
                "/home/imuser/qiime_output/rep-seqs-dada2_paired.qzv",
                "/home/imuser/qiime_output/table-dada2_paired.qzv",
                "/home/imuser/qiime_output/stats-dada2_paired.qzv")
    
    if(is.null(input$sample_data_paired$datapath)){
      add_metadata_table <- ""
      add_metadata_rarefaction <- ""
    }else{
      add_metadata_table <- paste("--m-sample-metadata-file", input$sample_data_paired$datapath)
      add_metadata_rarefaction <- paste("--m-metadata-file", input$sample_data_paired$datapath)
    }
    
    system(paste(qiime_cmd, "dada2 denoise-paired --i-demultiplexed-seqs /home/imuser/qiime_output/demux_paired_end.qza", 
                 '--p-trim-left-f', input$trim_left_f_paired,
                 '--p-trim-left-r', input$trim_left_r_paired,
                 '--p-trunc-len-f', input$trunc_len_f_paired,
                 '--p-trunc-len-r', input$trunc_len_r_paired,
                 "--p-trunc-q", input$qvalue_paired,
                 "--p-n-threads", input$threads_paired,
                 "--p-min-fold-parent-over-abundance", input$chimera_paired,
                 "--p-n-reads-learn", input$n_reads_paired,
                 "--o-representative-sequences /home/imuser/qiime_output/rep-seqs-dada2_paired.qza" ,
                 '--o-table /home/imuser/qiime_output/table-dada2_paired.qza', 
                 '--o-denoising-stats /home/imuser/qiime_output/stats-dada2_paired.qza'))
    
    system(paste(qiime_cmd, 'metadata tabulate --m-input-file /home/imuser/qiime_output/stats-dada2_paired.qza',
                 '--o-visualization /home/imuser/qiime_output/stats-dada2_paired.qzv'))
    # system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table-dada2.qza --o-visualization /home/imuser/qiime_output/table-dada2.qzv --m-sample-metadata-file', paste0("/home/imuser/qiime_output/",sample_file_name())))
    system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table-dada2_paired.qza',
                 add_metadata_table,
                 '--o-visualization /home/imuser/qiime_output/table-dada2_paired.qzv'
                 # , '--m-sample-metadata-file', input$sample_data_paired_position$datapath
    ))
    system(paste(qiime_cmd, 'feature-table tabulate-seqs --i-data /home/imuser/qiime_output/rep-seqs-dada2_paired.qza --o-visualization /home/imuser/qiime_output/rep-seqs-dada2_paired.qzv'))
    
    unlink("/home/imuser/qiime_output/denoise_paired_stats/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_paired_stats/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/stats-dada2.qzv /home/imuser/qiime_output/stats-dada2.zip")
    system("unzip -d /home/imuser/qiime_output/denoise_paired_stats /home/imuser/qiime_output/stats-dada2_paired.qzv")
    unzip_dirnames_stats <- list.files("/home/imuser/qiime_output/denoise_paired_stats", full.names = T)
    system(paste0("mv ", unzip_dirnames_stats, " /home/imuser/qiime_output/denoise_paired_stats/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_paired_stats/ /var/www/html/")
    
    unlink("/home/imuser/qiime_output/denoise_paired_table/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_paired_table/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/table-dada2.qzv /home/imuser/qiime_output/table-dada2.zip")
    system("unzip -d /home/imuser/qiime_output/denoise_paired_table /home/imuser/qiime_output/table-dada2_paired.qzv")
    unzip_dirnames_table <- list.files("/home/imuser/qiime_output/denoise_paired_table", full.names = T)
    system(paste0("mv ", unzip_dirnames_table, " /home/imuser/qiime_output/denoise_paired_table/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_paired_table/ /var/www/html/")
    
    unlink("/home/imuser/qiime_output/denoise_paired_seqs/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_paired_seqs/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/rep-seqs-dada2.qzv /home/imuser/qiime_output/rep-seqs-dada2.zip")
    system("unzip -d /home/imuser/qiime_output/denoise_paired_seqs /home/imuser/qiime_output/rep-seqs-dada2_paired.qzv")
    unzip_dirnames_seqs <- list.files("/home/imuser/qiime_output/denoise_paired_seqs", full.names = T)
    system(paste0("mv ", unzip_dirnames_seqs, " /home/imuser/qiime_output/denoise_paired_seqs/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_paired_seqs/ /var/www/html/")
    
    # alpha-rarefaction
    max_depth <- read_qza("/home/imuser/qiime_output/table-dada2_paired.qza")[["data"]] %>% colSums() %>% min()
    system(paste(qiime_cmd, "phylogeny align-to-tree-mafft-fasttree --i-sequences /home/imuser/qiime_output/rep-seqs-dada2_paired.qza", 
                 "--p-n-threads", input$threads_paired,
                 "--o-alignment /home/imuser/qiime_output/aligned-rep-seqs-dada2_paired.qza --o-masked-alignment /home/imuser/qiime_output/masked-aligned-rep-seqs-dada2_paired.qza",
                 "--o-tree /home/imuser/qiime_output/unrooted-tree_paired.qza --o-rooted-tree /home/imuser/qiime_output/rooted-tree_paired.qza"))
    
    system(paste(qiime_cmd, 'diversity alpha-rarefaction --i-table /home/imuser/qiime_output/table-dada2_paired.qza',
                 '--p-max-depth', max_depth,
                 '--i-phylogeny /home/imuser/qiime_output/rooted-tree_paired.qza',
                 add_metadata_rarefaction,
                 '--o-visualization /home/imuser/qiime_output/rarefaction-dada2_paired.qzv'))
    unlink("/home/imuser/qiime_output/denoise_paired_rarefaction/new_dirname", recursive = T)
    unlink("/var/www/html/denoise_paired_rarefaction/new_dirname", recursive = T)
    system("unzip -d /home/imuser/qiime_output/denoise_paired_rarefaction /home/imuser/qiime_output/rarefaction-dada2_paired.qzv")
    unzip_dirnames_seqs <- list.files("/home/imuser/qiime_output/denoise_paired_rarefaction", full.names = T)
    system(paste0("mv ", unzip_dirnames_seqs, " /home/imuser/qiime_output/denoise_paired_rarefaction/new_dirname"))
    system("cp -r /home/imuser/qiime_output/denoise_paired_rarefaction/ /var/www/html/")
    
    
    
    # removeModal()
    remove_modal_spinner()
    
    end_time <- Sys.time()
    spent_time <- format(round(end_time-start_time, digits = 2))
    
    if(all.equal(file.exists(c('/home/imuser/qiime_output/stats-dada2_paired.qzv', 
                               '/home/imuser/qiime_output/table-dada2_paired.qzv', 
                               '/home/imuser/qiime_output/rep-seqs-dada2_paired.qzv',
                               '/home/imuser/qiime_output/rarefaction-dada2_paired.qzv')), 
                 c(T, T, T, T))){
      # output$word_denoising_paired_position <- renderText({
      #   print('Denoising successfully')
      # })
      showModal(modalDialog(title = strong("Denoising succecessfully!"), 
                            HTML(
                              paste0(
                                "The process took ", spent_time, ". ",
                                "You can inspect the results!")
                            ), 
                            footer = NULL, easyClose = T, size = "l"))
    }else{
      # output$word_denoising_paired_position <- renderText({
      #   print('Error!')
      # })
      showModal(modalDialog(title = strong("Error!", style = "color: red"),
                            "Please check your files.", 
                            footer = NULL, easyClose = T, size = "l"))
    }
    
  })
  
  
  
  observeEvent(input$my_cores_paired, {
    
    my_cores <- parallel::detectCores()
    str1 <- paste0("The number of your threads is ", my_cores, ".")
    str2 <- " If 0 is provided, all available cores will be used. [default: all threads - 2]"
    
    showModal(modalDialog(title = strong("Message"),
                          HTML(paste(str1, str2, sep = "</br>")), 
                          footer = NULL, easyClose = T, size = "l"))
    
  })
  
  observeEvent(input$Q_cores_paired, {
    # output$message_thread_paired_position <- renderText(
    #   " The number of threads to use for multithreaded processing.")
    showModal(modalDialog(title = strong("Message"),
                          "The number of threads to use for multithreaded processing.", 
                          footer = NULL, easyClose = T, size = "l"))
  })
  
  observeEvent(input$word_chimera_paired, {
    # output$message_thread_single_position <- renderText(
    #   "The number of threads to use for multithreaded processing.")
    showModal(modalDialog(title = strong("Message"),
                          "The minimum abundance of potential parents of a
                         sequence being tested as chimeric, expressed as a
                         fold-change versus the abundance of the sequence
                         being tested. Values should be greater than or equal
                         to 1 (i.e. parents should be more abundant than the
                         sequence being tested).", 
                          footer = NULL, easyClose = T, size = "l"))
  })
  
  observeEvent(input$Q_learn_reads_paired, {
    # output$message_learn_paired_position <- renderText(
    #   " The number of reads to use when training the error model. Smaller numbers will result in a shorter run time but a less reliable error model. [default: 1000000]")
    showModal(modalDialog(title = strong("Message"),
                          " The number of reads to use when training the error model. Smaller numbers will result in a shorter run time but a less reliable error model. [default: 1000000]", 
                          footer = NULL, easyClose = T, size = "l"))
  })
  
  
  
  # Clustering -----------------------------------------------------------------------------------------------
  observeEvent(input$clustering, {
    
    output$word_clustering <- renderText({
      
      showModal(modalDialog(title = "Running denoising for single end ...", "Waiting for a moment",footer = NULL))
      
      
      qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
      
      if(input$seqs_type == "Single end"){
        
        system(paste(qiime_cmd, "vsearch dereplicate-sequences --i-sequences",'/home/imuser/qiime_output/demux_single_end.qza', 
                     "--o-dereplicate-table /home/imuser/qiime_output/table_dereplicate.qza",
                     "--o-dereplicate-sequences /home/imuser/qiime_output/rep_seqs_dereplicate.qza"))
      }else{
        
        system(paste(qiime_cmd, "vsearch dereplicate-sequences --i-sequences",'/home/imuser/qiime_output/demux_paired_end.qza', 
                     "--p-threads", input$threads_clustering,
                     "--o-dereplicate-table /home/imuser/qiime_output/table_dereplicate.qza",
                     "--o-dereplicate-sequences /home/imuser/qiime_output/rep_seqs_dereplicate.qza"))
      }
      
      
      # system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table_dereplicate.qza --o-visualization /home/imuser/qiime_output/table_dereplicate.qzv --m-sample-metadata-file', input$sample_data_paired$datapath))
      # system(paste(qiime_cmd, 'feature-table tabulate-seqs --i-data /home/imuser/qiime_output/rep_seqs_dereplicate.qza --o-visualization /home/imuser/qiime_output/rep_seqs_dereplicate.qzv'))
      
      system(paste(qiime_cmd, "vsearch cluster-features-de-novo --i-table /home/imuser/qiime_output/table_dereplicate.qza",
                   "--i-sequences /home/imuser/qiime_output/rep-seqs_dereplicate.qza",
                   "--i-perc-identity", input$clustering_percent,
                   "--o-clustered-table /home/imuser/qiime_output/table_dn_percent.qza",
                   "--o-clustered-sequences /home/imuser/qiime_output/rep_seqs_dn_percent.qza"))
      
      
      system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table_dn_percent.qza --o-visualization /home/imuser/qiime_output/table_dn_percent.qzv --m-sample-metadata-file', input$sample_data_paired_position$datapath))
      system(paste(qiime_cmd, 'feature-table tabulate-seqs --i-data /home/imuser/qiime_output/rep_seqs_dn_percent.qza --o-visualization /home/imuser/qiime_output/rep_seqs_dn_percent.qzv'))
      
      if(all.equal(file.exists(c('/home/imuser/qiime_output/table_dn_percent.qzv', '/home/imuser/qiime_output/rep_seqs_dn_percent.qzv')), c(T, T, T))){
        output$word_clustering <- renderText({
          print('Clustering successfully')
        })
      }else{
        output$word_clustering <- renderText({
          print('Error!')
        })
      }
      
      removeModal()
      
      
    })
  })
  
  
  
  
  ## OTU picking ------------------------------------------------------------------------------------------------
  
  # observeEvent(input$OTU_picking, {
  #     
  #     showModal(modalDialog(title = "Running OTU picking ...", "Waiting for a moment", footer = NULL))
  #         
  #         qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
  #         system(paste(qiime_cmd, 'vsearch cluster-features-de-novo --i-table /home/imuser/qiime_output/table-dada2.qza --i-sequences /home/imuser/qiime_output/rep-seqs-dada2.qza --p-perc-identity', 
  #                      input$OTU_identity,
  #                      '--o-clustered-table /home/imuser/qiime_output/table-dada2-dn.qza --o-clustered-sequences /home/imuser/qiime_output/rep-seqs-dada2-dn.qza'))
  #         
  #         system(paste(qiime_cmd, 'feature-table summarize --i-table /home/imuser/qiime_output/table-dada2-dn.qza --o-visualization /home/imuser/qiime_output/table-dada2-dn.qzv'))
  #         system(paste(qiime_cmd, 'feature-table tabulate-seqs --i-data /home/imuser/qiime_output/rep-seqs-dada2-dn.qza --o-visualization /home/imuser/qiime_output/rep-seqs-dada2-dn.qzv'))
  #         
  #         if(file.exists(c('/home/imuser/qiime_output/table-dada2-dn.qzv', '/home/imuser/qiime_output/rep-seqs-dada2-dn.qzv'))){
  #             output$word_OTUs <- renderText({
  #                 print('Picking successfully')
  #                 })
  #         }else{
  #             output$word_OTUs <- renderText({
  #                 print('Error!')
  #             })
  #         }
  #         
  #         
  #     removeModal()
  # })
  # 
  
  
  
  # Taxonomic Analysis ------------------------------------------------------------------------------------------
  
  
  output$out_f <- renderUI({
    
    if (input$primer_f == "other"){
      isolate({
        text_list <- list()
        text_list[[1]] <- textInput(inputId = "primer_f_manu", label = "Give the forward primer sequences")
        return(text_list)
      })
    }else{
      isolate({
        text_list <- list()
        text_list[[1]] <- p("")
        return(text_list)
      })
    }
    
    
  })
  
  output$out_r <- renderUI({
    
    if (input$primer_r == "other"){
      isolate({
        text_list <- list()
        text_list[[1]] <- textInput(inputId = "primer_r_manu", label = "Give the reverse primer sequences")
        return(text_list)
      })
    }else{
      isolate({
        text_list <- list()
        text_list[[1]] <- p("")
        return(text_list)
      })
    }
    
    
  })
  
  
  observeEvent(input$show_primer, {
    # output$mk_taxa <- renderUI({
    output$taxonomy_output <- renderDataTable({
      
      # addResourcePath(prefix = "text_files",directoryPath = "/home/imuser/")
      # # file.copy(from = "/home/imuser/text_files/Primer_table.html", to = "/var/www/html/", overwrite = T)
      # # withMathJax(includeMarkdown((path = "/home/imuser/text_files/Primer_table.Rmd")))
      # 
      # # tags$iframe(src="Primer_table.html")
      # # HTML(markdown::markdownToHTML(knit('/home/imuser/text_files/Primer_table.Rmd', quiet = F)))
      # rmarkdown::render('/home/imuser/text_files/Primer_table.Rmd')
      primer_table <- data.frame(Primer=c("8F", "27F", "CC [F]", "357F", "515F", "533F", "16S.1100.F16", "1237F", 
                                          "519R", "CD [R]", "907R", "1100R", "1391R", "1492R (l)", "1492R (s)"), 
                                 Sequences=c("AGAGTTTGATCCTGGCTCAG", "AGAGTTTGATCMTGGCTCAG", "CCAGACTCCTACGGGAGGCAGC", "CTCCTACGGGAGGCAGCAG",
                                             "GTGCCAGCMGCCGCGGTAA", "GTGCCAGCAGCCGCGGTAA", "CAACGAGCGCAACCCT", "GGGCTACACACGYGCWAC",
                                             "GWATTACCGCGGCKGCTG", "CTTGTGCGGGCCCCCGTCAATTC", "CCGTCAATTCMTTTRAGTTT", "AGGGTTGCGCTCGTTG",
                                             "GACGGGCGGTGTGTRCA", "GGTTACCTTGTTACGACTT", "ACCTTGTTACGACTT"),
                                 Target_Group=c(rep("Universal", 11), "Bacterial", rep("Universal", 3)),
                                 Reference=c("Turner et al. 1999", "Lane et al. 1991", "Rudi et al. 1997", "Turner et al. 1999", "Turner et al. 1999",
                                             "Weisburg et al. 1991", "Turner et al. 1999", "Turner et al. 1999", "Turner et al. 1999", "Rudi et al. 1997",
                                             "Lane et al. 1991", "Turner et al. 1999", "Turner et al. 1999", "Turner et al. 1999", "Lane et al. 1991"))
      colnames(primer_table)[2] <- "Sequences (5'-3')"
      colnames(primer_table)[3] <- "Target group"
      
      return(primer_table)
      
    })
  })
  
  
  
  observeEvent(input$start_training, {
    
    start_time <- Sys.time()
    
    # showModal(modalDialog(title = "Running taxonomic analysis ...", "Waiting for a moment", footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    
    
    
    database_list <- list(Silva_99=c(list.files(path = "/home/imuser/taxa_database/silva/rep_set_16S_only/99", full.names = T), "/home/imuser/taxa_database/silva/taxonomy/16S_only/99/taxonomy_7_levels.txt"),
                          Silva_97=c(list.files(path = "/home/imuser/taxa_database/silva/rep_set_16S_only/97", full.names = T), "/home/imuser/taxa_database/silva/taxonomy/16S_only/97/taxonomy_7_levels.txt"),
                          Silva_94=c(list.files(path = "/home/imuser/taxa_database/silva/rep_set_16S_only/94", full.names = T), "/home/imuser/taxa_database/silva/taxonomy/16S_only/94/taxonomy_7_levels.txt"),
                          Silva_90=c(list.files(path = "/home/imuser/taxa_database/silva/rep_set_16S_only/90", full.names = T), "/home/imuser/taxa_database/silva/taxonomy/16S_only/90/taxonomy_7_levels.txt"),
                          GreenGene_99=c("/home/imuser/taxa_database/greengenes/rep_set/99_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/99_otu_taxonomy.txt"),
                          GreenGene_97=c("/home/imuser/taxa_database/greengenes/rep_set/97_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/97_otu_taxonomy.txt"),
                          GreenGene_94=c("/home/imuser/taxa_database/greengenes/rep_set/94_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/94_otu_taxonomy.txt"),
                          GreenGene_91=c("/home/imuser/taxa_database/greengenes/rep_set/91_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/91_otu_taxonomy.txt"),
                          GreenGene_88=c("/home/imuser/taxa_database/greengenes/rep_set/88_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/88_otu_taxonomy.txt"),
                          GreenGene_85=c("/home/imuser/taxa_database/greengenes/rep_set/85_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/85_otu_taxonomy.txt"),
                          GreenGene_82=c("/home/imuser/taxa_database/greengenes/rep_set/82_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/82_otu_taxonomy.txt"),
                          GreenGene_79=c("/home/imuser/taxa_database/greengenes/rep_set/79_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/79_otu_taxonomy.txt"),
                          GreenGene_76=c("/home/imuser/taxa_database/greengenes/rep_set/76_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/76_otu_taxonomy.txt"),
                          GreenGene_73=c("/home/imuser/taxa_database/greengenes/rep_set/73_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/73_otu_taxonomy.txt"),
                          GreenGene_70=c("/home/imuser/taxa_database/greengenes/rep_set/70_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/70_otu_taxonomy.txt"),
                          GreenGene_67=c("/home/imuser/taxa_database/greengenes/rep_set/67_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/67_otu_taxonomy.txt"),
                          GreenGene_64=c("/home/imuser/taxa_database/greengenes/rep_set/64_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/64_otu_taxonomy.txt"),
                          GreenGene_61=c("/home/imuser/taxa_database/greengenes/rep_set/61_otus.fasta", "/home/imuser/taxa_database/greengenes/taxonomy/61_otu_taxonomy.txt")
    )
    
    #  transform database to .qza
    rm("/home/imuser/qiime_output/identity_otus.qza", "/home/imuser/qiime_output/ref-taxonomy-7.qza")
    system(paste(qiime_cmd, "tools import --type 'FeatureData[Sequence]' --input-path", database_list[[input$select_database]][1], 
                 "--output-path /home/imuser/qiime_output/identity_otus.qza"))
    system(paste(qiime_cmd, "tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path", database_list[[input$select_database]][2], 
                 "--output-path /home/imuser/qiime_output/ref-taxonomy-7.qza"))
    
    primer_list <- list("8F"="AGAGTTTGATCCTGGCTCAG",
                        "27F"="AGAGTTTGATCMTGGCTCAG",
                        "CC [F]"="CCAGACTCCTACGGGAGGCAGC",
                        "357F"="CTCCTACGGGAGGCAGCAG",
                        "515F"="GTGCCAGCMGCCGCGGTAA",
                        "533F"="GTGCCAGCAGCCGCGGTAA",
                        "16S.1100.F16"="CAACGAGCGCAACCCT",
                        "1237F"="GGGCTACACACGYGCWAC",
                        "519R"="GWATTACCGCGGCKGCTG",
                        "CD [R]"="CTTGTGCGGGCCCCCGTCAATTC",
                        "907R"="CCGTCAATTCMTTTRAGTTT",
                        "1100R"="AGGGTTGCGCTCGTTG",
                        "1391R"="GACGGGCGGTGTGTRCA",
                        "1492R (l)"="GGTTACCTTGTTACGACTT",
                        "1492R (s)"="ACCTTGTTACGACTT" 
    )
    
    
    rm("/home/imuser/qiime_output/ref-seqs.qza")
    if(input$primer_f != "other" & input$primer_r != "other"){
      system(paste(qiime_cmd, "feature-classifier extract-reads --i-sequences /home/imuser/qiime_output/identity_otus.qza --p-f-primer", primer_list[[input$primer_f]], 
                   "--p-r-primer", primer_list[[input$primer_r]], 
                   # "--p-trunc-len", input$trunc_length, 
                   "--p-min-length", input$min_length, 
                   "--p-max-length", input$max_length, 
                   "--p-n-jobs", input$n_jobs, 
                   "--o-reads /home/imuser/qiime_output/ref-seqs.qza"))
    }
    
    if(input$primer_f == "other" & input$primer_r != "other"){
      system(paste(qiime_cmd, "feature-classifier extract-reads --i-sequences /home/imuser/qiime_output/identity_otus.qza --p-f-primer", input$primer_f_manu, 
                   "--p-r-primer", primer_list[[input$primer_r]], 
                   # "--p-trunc-len", input$trunc_length, 
                   "--p-min-length", input$min_length, 
                   "--p-max-length", input$max_length, 
                   "--p-n-jobs", input$n_jobs, 
                   "--o-reads /home/imuser/qiime_output/ref-seqs.qza"))
    }
    
    if(input$primer_f != "other" & input$primer_r == "other"){
      system(paste(qiime_cmd, "feature-classifier extract-reads --i-sequences /home/imuser/qiime_output/identity_otus.qza --p-f-primer", primer_list[[input$primer_f]], 
                   "--p-r-primer", input$primer_r_manu, 
                   # "--p-trunc-len", input$trunc_length, 
                   "--p-min-length", input$min_length, 
                   "--p-max-length", input$max_length, 
                   "--p-n-jobs", input$n_jobs, 
                   "--o-reads /home/imuser/qiime_output/ref-seqs.qza"))
    }
    
    if(input$primer_f == "other" & input$primer_r == "other"){
      system(paste(qiime_cmd, "feature-classifier extract-reads --i-sequences /home/imuser/qiime_output/identity_otus.qza --p-f-primer", input$primer_f_manu, 
                   "--p-r-primer", input$primer_r_manu, 
                   # "--p-trunc-len", input$trunc_length, 
                   "--p-min-length", input$min_length, 
                   "--p-max-length", input$max_length, 
                   "--p-n-jobs", input$n_jobs, 
                   "--o-reads /home/imuser/qiime_output/ref-seqs.qza"))
    }
    
    rm("/home/imuser/qiime_output/classifier.qza")
    system(paste(qiime_cmd, "feature-classifier fit-classifier-naive-bayes --i-reference-reads /home/imuser/qiime_output/ref-seqs.qza --i-reference-taxonomy /home/imuser/qiime_output/ref-taxonomy-7.qza --o-classifier /home/imuser/qiime_output/classifier.qza"))
    
    rm("/home/imuser/qiime_output/taxonomy.qza")
    if(input$seqs_type == "Single end"){
      system(paste(qiime_cmd, "feature-classifier classify-sklearn --i-classifier /home/imuser/qiime_output/classifier.qza --i-reads /home/imuser/qiime_output/rep-seqs-dada2_single.qza", 
                   "--p-n-jobs", input$n_jobs,
                   "--o-classification /home/imuser/qiime_output/taxonomy.qza"))
    }else{
      system(paste(qiime_cmd, "feature-classifier classify-sklearn --i-classifier /home/imuser/qiime_output/classifier.qza --i-reads /home/imuser/qiime_output/rep-seqs-dada2_paired.qza", 
                   "--p-n-jobs", input$n_jobs,
                   "--o-classification /home/imuser/qiime_output/taxonomy.qza"))
    }
    
    # system(paste(qiime_cmd, "feature-classifier classify-sklearn --i-classifier /home/imuser/qiime_output/classifier.qza --i-reads /home/imuser/qiime_output/rep-seqs-dada2.qza", 
    #              "--p-n-jobs", input$n_jobs,
    #              "--o-classification /home/imuser/qiime_output/taxonomy.qza"))
    
    rm("/home/imuser/qiime_output/taxonomy.qzv")
    system(paste(qiime_cmd, "metadata tabulate --m-input-file /home/imuser/qiime_output/taxonomy.qza --o-visualization /home/imuser/qiime_output/taxonomy.qzv"))
    
    rm("/home/imuser/qiime_output/taxatable7.qza")
    if(input$seqs_type == "Single end"){
      system(paste(qiime_cmd, "taxa collapse --i-table /home/imuser/qiime_output/table-dada2_single.qza",
                   "--i-taxonomy /home/imuser/qiime_output/taxonomy.qza --p-level 7 --o-collapsed-table /home/imuser/qiime_output/taxatable7.qza"))
    }else{
      system(paste(qiime_cmd, "taxa collapse --i-table /home/imuser/qiime_output/table-dada2_paired.qza",
                   "--i-taxonomy /home/imuser/qiime_output/taxonomy.qza --p-level 7 --o-collapsed-table /home/imuser/qiime_output/taxatable7.qza"))
    }
    
    
    system("rm -r /home/imuser/qiime_output/taxonomy_unzip/new_dirname")
    system("rm -r /home/imuser/var/www/html/taxonomy_unzip/new_dirname")
    system("unzip -d /home/imuser/qiime_output/taxonomy_unzip /home/imuser/qiime_output/taxonomy.qzv")
    unzip_dirnames_taxa <- list.files("/home/imuser/qiime_output/taxonomy_unzip", full.names = T)
    system(paste0("mv ", unzip_dirnames_taxa, " /home/imuser/qiime_output/taxonomy_unzip/new_dirname"))
    system("cp -r /home/imuser/qiime_output/taxonomy_unzip/ /var/www/html/")
    
    
    
    # removeModal()
    remove_modal_spinner()
    
    end_time <- Sys.time()
    spent_time <- format(round(end_time-start_time, digits = 2))
    
    if(file.exists("/home/imuser/qiime_output/taxonomy.qzv")){
      # output$word_training <- renderText({
      #   print('Training successfully')
      # })
      showModal(modalDialog(title = strong("Taxonomic analysis has been finished!"), 
                            HTML(
                              paste0(
                                "The process took ", spent_time, ". ",
                                "You can inspect the results!")
                            ), 
                            footer = NULL, easyClose = T, size = "l"))
    }else{
      # output$word_training <- renderText({
      #   print('Error')
      # })
      showModal(modalDialog(title = strong("Error!", style = "color: red"),
                            "Please check your files.", 
                            footer = NULL, easyClose = T, size = "l"))
    }
    
  })
  
  
  
  output$taxatable_download <- downloadHandler(
    
    filename <-"taxonomic_table.qza",
    
    content = function(file){
      file.copy("/home/imuser/qiime_output/taxatable7.qza", file)
    }
    
  )
  
  output$table_dada2_download <- downloadHandler(
    filename = "table_for_phylo.qza",
    
    content = function(file){
      if(input$seqs_type == "Single end"){
        file.copy("/home/imuser/qiime_output/table-dada2_single.qza", file)
      }else{
        file.copy("/home/imuser/qiime_output/table-dada2_paired.qza", file)
      }
      
    }
  )
  
  output$rep_seq_dada2_download <- downloadHandler(
    filename = "rep_seqs_for_phylo.qza",
    
    content = function(file){
      if(input$seqs_type == "Single end"){
        file.copy("/home/imuser/qiime_output/rep-seqs-dada2_single.qza", file)
      }else{
        file.copy("/home/imuser/qiime_output/rep-seqs-dada2_paired.qza", file)
      }
      
    }
  )
  
  
  
  
  # Data analysis ----
  
  # common reactive objects ----
  TaxaTable_output <- reactive({
    
    # revise the species names 
    as_output_taxtable <- function(df_data){
      
      df_data_rownames<-row.names(df_data)
      
      df_data <- cbind(Species=df_data_rownames, df_data) %>% as.data.frame()
      
      row.names(df_data) <- NULL
      
      # remove the non-sense string
      df_data$Species<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", "", df_data$Species)
      df_data$Species<-gsub("k__|p__|c__|o__|f__|g__|s__", "", df_data$Species)
      
      
      library(tidyr)
      df_data <- df_data %>%
        separate(
          col = Species,
          into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
          sep = ";"
        )
      
      # replace "__" to "Unassigned"
      df_data<-replace(df_data, df_data=="", "Unassigned")
      df_data<-replace(df_data, df_data=="__", "Unassigned")
      
      return(as.data.frame(df_data))
    }
    
    
    as_taxtable_perFeature<-function(taxatable_data, metadata_data){
      
      
      metadata_feature_name <- colnames(metadata_data)
      
      
      
      metadata_feature_summary <- lapply(1:length(metadata_data), function(i){
        metadata_data[,i] %>% unique() %>% as.character()
      })
      names(metadata_feature_summary) <- metadata_feature_name
      
      
      k<-function(x, taxatable_data_k, metadata_data_k){
        
        sample_<-lapply(1:length(metadata_feature_summary[[x]]), function(i){
          
          metadata_data_k[,1][which(metadata_data_k[,x]==metadata_feature_summary[[x]][i])] %>% as.character()
          
        })
        
        names(sample_) <- metadata_feature_summary[[x]]
        
        
        taxatable_data_ <- lapply(1:length(metadata_feature_summary[[x]]), function(i){
          
          rowSums(as.data.frame(taxatable_data_k[,sample_[[i]]]))
          
        })
        
        taxatable_data_<-as.data.frame(taxatable_data_)
        colnames(taxatable_data_)<-names(sample_)
        
        return(taxatable_data_)
      }
      
      
      all_taxatable_perFeature <- sapply(1:length(metadata_feature_name), function(i){
        
        k(i, taxatable_data, metadata_data)
        
      })
      
      names(all_taxatable_perFeature) <- metadata_feature_name
      return(all_taxatable_perFeature)
      
    }
    
    
    TaxaTable_data_list <- as_taxtable_perFeature(TaxaTable(), Metadata())
    
    TaxaTable_data_list_output <- as_output_taxtable(TaxaTable_data_list[[input$metadata1]])
    
    colnames(TaxaTable_data_list_output)[1:7] <- c(
      paste0('Kingdom', " (K=", length(unique(TaxaTable_data_list_output[, 1])), ")"),
      paste0('Phylum', " (K=", nrow(unique(TaxaTable_data_list_output[, 1:2])), ")"),
      paste0('Class', " (K=", nrow(unique(TaxaTable_data_list_output[, 1:3])), ")"),
      paste0('Order', " (K=", nrow(unique(TaxaTable_data_list_output[, 1:4])), ")"),
      paste0('Family', " (K=", nrow(unique(TaxaTable_data_list_output[, 1:5])), ")"),
      paste0('Genus', " (K=", nrow(unique(TaxaTable_data_list_output[, 1:6])), ")"),
      paste0('Species', " (K=", nrow(unique(TaxaTable_data_list_output[, 1:7])), ")")
    )
    
    
    
    a <- as_output_taxtable(TaxaTable())
    ft_names <- input$metadata1 
    b <- list()
    b <- lapply(as.character(unique(Metadata()[, ft_names])), function(i){
      
      b_metadata_select <- Metadata()[Metadata()[, ft_names] %in% i,] %>% select("SampleID") 
      b_metadata_select <- b_metadata_select[, "SampleID"] %>% as.character()
      b_select <- a %>% select(b_metadata_select)
      
    })
    
    k <- a[, 1:7]
    for (i in 1:length(unique(Metadata()[, ft_names]))) {
      
      k <- cbind(k, b[[i]])
      
    }
    
    
    c <- list()
    c <- lapply(unique(Metadata()[, ft_names]), function(i){
      
      c_metadata_select <- Metadata()[Metadata()[, ft_names] %in% i,] %>% select("SampleID") 
      c_metadata_select_length <- c_metadata_select[, "SampleID"] %>% as.character() %>% length()
      
    })
    
    names(c) <- unique(Metadata()[, ft_names])
    
    
    sketch = htmltools::withTags(table(style="text-align:center;",
                                       class = 'display',
                                       col(style="width:50%"),
                                       thead(
                                         tr(
                                           th(rowspan = 2, paste0('Kingdom', " (K=", length(unique(k[, 1])), ")")),
                                           th(rowspan = 2, paste0('Phylum', " (K=", nrow(unique(k[, 1:2])), ")")),
                                           th(rowspan = 2, paste0('Class', " (K=", nrow(unique(k[, 1:3])), ")")),
                                           th(rowspan = 2, paste0('Order', " (K=", nrow(unique(k[, 1:4])), ")")),
                                           th(rowspan = 2, paste0('Family', " (K=", nrow(unique(k[, 1:5])), ")")),
                                           th(rowspan = 2, paste0('Genus', " (K=", nrow(unique(k[, 1:6])), ")")),
                                           th(rowspan = 2, paste0('Species', " (K=", nrow(unique(k[, 1:7])), ")")),
                                           # th(colspan = 8, 'gut'),
                                           # th(colspan = 8, 'left palm'),
                                           # th(colspan = 9, 'right palm'),
                                           # th(colspan = 9, 'tongue')
                                           lapply(1:length(unique(Metadata()[, ft_names])), function(i){
                                             
                                             th(colspan = c[[i]], paste0(unique(Metadata()[, ft_names])[i]), " ( N=", c[[i]], ")")
                                           })
                                         ),
                                         tr(
                                           lapply(colnames(k)[-(1:7)], th)
                                         )
                                       )
    ))
    
    
    
    
    if(input$metadata1=="SampleID"){
      
      return(DT::datatable(TaxaTable_data_list_output, rownames = F))
    }else{
      # use rownames = FALSE here because we did not generate a cell for row names in the header
      return(DT::datatable(k, container = sketch, rownames = F))
    }
    
  }) # show the taxatable in Data analysis
  
  alpha_diversity_table <- reactive({
    
    as_alpha_diversity_table <- function(taxatable_data){
      
      alpha_diversity_richness<-sapply(1:ncol(taxatable_data), function(i){
        nrow(subset(as.data.frame(taxatable_data), 
                    as.data.frame(taxatable_data)[,i]>0))
      })
      
      library(fossil)
      alpha_diversity_Choa1<-sapply(1:ncol(taxatable_data), function(i){
        chao1(taxatable_data[,i])
      })
      
      alpha_diversity_ACE<-sapply(1:ncol(taxatable_data), function(i){
        ACE(taxatable_data[,i])
      })
      
      library(vegan)
      alpha_diversity_Shannon<-sapply(1:ncol(taxatable_data), function(i){
        diversity(taxatable_data[,i], index = "shannon")
      })
      
      alpha_diversity_Simpsom<-sapply(1:ncol(taxatable_data), function(i){
        diversity(taxatable_data[,i], index = "simpson")
      })
      
      alpha_diversity_InvSimpson<-sapply(1:ncol(taxatable_data), function(i){
        diversity(taxatable_data[,i], index = "invsimpson")
      })
      
      
      ShannonEvenness<-function(data_vector){
        
        shannon_diversity<-diversity(data_vector, index = "shannon")
        species_number<-length(data_vector)
        Hmax<-log(species_number)
        J<-shannon_diversity/Hmax
        print(J)
        
      }
      
      alpha_diversity_ShannonEvenness<-sapply(1:ncol(taxatable_data), function(i){
        ShannonEvenness(taxatable_data[,i])
      })
      
      
      SimpsonEvenness<-function(data_vector){
        
        D<-diversity(data_vector, index = "simpson")
        species_number<-length(data_vector)
        E<-(1/D)/species_number
        print(E)
        
      }
      
      alpha_diversity_SimpsonEveness<-sapply(1:ncol(taxatable_data), function(i){
        SimpsonEvenness(taxatable_data[,i])
      })
      
      GoodCoverage<-function(data_vector){
        
        n<-length(which(data_vector==1))
        N<-length(which(data_vector!=0))
        C<-(1-(n/N))
        print(C)
        
      }
      
      alpha_diversity_GoodCoverage<-sapply(1:ncol(taxatable_data), function(i){
        GoodCoverage(taxatable_data[,i])
      })
      
      
      alpha_diversity_table<-data.frame(Sample=colnames(taxatable_data),
                                        Richness=alpha_diversity_richness,
                                        Chao1=alpha_diversity_Choa1,
                                        ACE=alpha_diversity_ACE,
                                        Shannon_diverstiy=alpha_diversity_Shannon,
                                        Simpon_diversity=alpha_diversity_Simpsom,
                                        InvSimpson_diversity=alpha_diversity_InvSimpson,
                                        Shannon_evenness=alpha_diversity_ShannonEvenness,
                                        Simpson_evenness=alpha_diversity_SimpsonEveness,
                                        Goods_coverage=alpha_diversity_GoodCoverage
      )
      
      
      alpha_diversity_table[,2:ncol(alpha_diversity_table)]<-
        round(alpha_diversity_table[,2:ncol(alpha_diversity_table)], 4)  
      
      return(alpha_diversity_table)
    }
    
    return(as_alpha_diversity_table(TaxaTable_merge()))
    
  })
  
  alpha_anova_boxplot <- reactive({
    
    A_diversity <- alpha_diversity_table()
    
    colnames(A_diversity)[1] <- colnames(Metadata_stats())[1]
    A_diversity_metadata <- merge(A_diversity, Metadata_stats(), by= colnames(Metadata_stats())[1])
    
    A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)], function(x){
      
      A_diversity_metadata_merge <- merge(A_diversity,
                                          Metadata_stats()[, c(colnames(A_diversity_metadata)[1], x)],
                                          by= colnames(Metadata_stats())[1])
    })
    
    names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)]
    
    library(reshape2)
    i <- which(colnames(Metadata_stats())==input$metadata_alpha)-1
    A_diversity_metadata_list_melt <- melt(A_diversity_metadata_list[[i]],
                                           id.vars=c("SampleID",
                                                     names(A_diversity_metadata_list)[i]))
    
    colnames(A_diversity_metadata_list_melt)[2] <- "feature_name"
    
    A_diversity_metadata_list_melt <- filter(A_diversity_metadata_list_melt, feature_name != "NA")
    
    index <- input$select_diversity
    A_diversity_metadata_list_melt_forplot <- subset(A_diversity_metadata_list_melt,
                                                     variable==index)
    anova_result <- aov(value ~ feature_name, A_diversity_metadata_list_melt_forplot)
    anova_summary <- summary(anova_result)
    # tukey_result <- agricolae::HSD.test(anova_result, "feature_name", group = T)
    # group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
    
    stat_box_data <- function(y, upper_limit = max(A_diversity_metadata_list_melt_forplot$value) * 1.15) {
      return( 
        data.frame(
          y = 0.95 * upper_limit,
          label = paste('n =', length(y))
        )
      )
    }
    
    ggplot(A_diversity_metadata_list_melt_forplot, aes(x = feature_name, y = value)) +
      # geom_text(data = data.frame(),
      #           aes(x = rownames(group_data), y = max(A_diversity_metadata_list_melt_forplot$value) + 1, label = group_data$groups),
      #           col = 'black',
      #           size = 5) +
      geom_boxplot() +
      ggtitle("Alpha diversity") +
      xlab(input$metadata_alpha) +
      ylab(input$select_diversity)+
      labs(caption = paste0("p value of ANOVA = ", round(anova_summary[[1]][[5]][1], 4))) + theme(text = element_text(size = 15)) + stat_summary(fun.data = stat_box_data,
                                                                                                                                                 geom = "text", 
                                                                                                                                                 hjust = 0.5,
                                                                                                                                                 vjust = 0.9)
  })
  
  alpha_anova_tukey <- reactive({
    
    A_diversity <- alpha_diversity_table()
    
    colnames(A_diversity)[1] <- colnames(Metadata_stats())[1]
    A_diversity_metadata <- merge(A_diversity, Metadata_stats(), by= colnames(Metadata_stats())[1])
    
    A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)], function(x){
      
      A_diversity_metadata_merge <- merge(A_diversity,
                                          Metadata_stats()[, c(colnames(A_diversity_metadata)[1], x)],
                                          by= colnames(Metadata_stats())[1])
    })
    
    names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)]
    
    library(reshape2)
    i <- which(colnames(Metadata_stats())==input$metadata_alpha)-1
    A_diversity_metadata_list_melt <- melt(A_diversity_metadata_list[[i]],
                                           id.vars=c("SampleID",
                                                     names(A_diversity_metadata_list)[i]))
    
    colnames(A_diversity_metadata_list_melt)[2] <- "feature_name"
    
    A_diversity_metadata_list_melt <- filter(A_diversity_metadata_list_melt, feature_name != "NA")
    
    index <- input$select_diversity
    A_diversity_metadata_list_melt_forplot <- subset(A_diversity_metadata_list_melt,
                                                     variable==index) 
    
    anova_result <- aov(value ~ feature_name, A_diversity_metadata_list_melt_forplot)
    
    tukey_result <- TukeyHSD(anova_result, "feature_name")[[1]]
    tukey_result_table <- data.frame(comparisons = rownames(tukey_result), tukey_result)
    tukey_result_table[,"comparisons"] <- stringr::str_replace_all(tukey_result_table[,"comparisons"], "-", " / ")
    tukey_result_table_separate <- separate(data = tukey_result_table, col = "comparisons", into = c("Group_A", "Group_B"), sep = " / ")
    order_by_count <- as.data.frame(table(tukey_result_table_separate[,1]))[,1]
    tukey_result_table_separate_order <- tukey_result_table_separate[rev(order(factor(tukey_result_table_separate$Group_A, levels = order_by_count))),]
    colnames(tukey_result_table_separate_order)[c(1,2,3,6)] <- c("Group A", "Group B", "Diff", "P value")
    return(tukey_result_table_separate_order[, c(1,2,3,6)])
    
  })
  
  alpha_KW_boxplot <- reactive({
    
    A_diversity <- alpha_diversity_table()
    
    colnames(A_diversity)[1] <- colnames(Metadata_stats())[1]
    A_diversity_metadata <- merge(A_diversity, Metadata_stats(), by= colnames(Metadata_stats())[1])
    
    A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)], function(x){
      
      A_diversity_metadata_merge <- merge(A_diversity,
                                          Metadata_stats()[, c(colnames(A_diversity_metadata)[1], x)],
                                          by= colnames(Metadata_stats())[1])
    })
    
    names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)]
    
    library(reshape2)
    i <- which(colnames(Metadata_stats())==input$metadata_alpha)-1
    A_diversity_metadata_list_melt <- melt(A_diversity_metadata_list[[i]],
                                           id.vars=c("SampleID",
                                                     names(A_diversity_metadata_list)[i]))
    
    colnames(A_diversity_metadata_list_melt)[2] <- "feature_name"
    
    A_diversity_metadata_list_melt <- filter(A_diversity_metadata_list_melt, feature_name != "NA")
    
    index <- input$select_diversity
    A_diversity_metadata_list_melt_forplot <- subset(A_diversity_metadata_list_melt,
                                                     variable==index)
    KW_result <- kruskal.test(value ~ feature_name, A_diversity_metadata_list_melt_forplot)
    
    stat_box_data <- function(y, upper_limit = max(A_diversity_metadata_list_melt_forplot$value) * 1.15) {
      return( 
        data.frame(
          y = 0.95 * upper_limit,
          label = paste('n =', length(y))
        )
      )
    }
    
    ggplot(A_diversity_metadata_list_melt_forplot, aes(x = feature_name, y = value)) +
      # geom_text(data = data.frame(),
      #           aes(x = rownames(group_data), y = max(A_diversity_metadata_list_melt_forplot$value) + 1, label = group_data$groups),
      #           col = 'black',
      #           size = 5) +
      geom_boxplot() +
      ggtitle("Alpha diversity") +
      xlab(input$metadata_alpha) +
      ylab(input$select_diversity)+
      labs(caption=paste0("p value of KW-test = ", round(KW_result$p.value, 4))) + theme(text = element_text(size = 15)) + stat_summary(fun.data = stat_box_data,
                                                                                                                                        geom = "text", 
                                                                                                                                        hjust = 0.5,
                                                                                                                                        vjust = 0.9)
  })
  
  alpha_KW_Dunn <- reactive({
    
    A_diversity <- alpha_diversity_table()
    
    colnames(A_diversity)[1] <- colnames(Metadata_stats())[1]
    A_diversity_metadata <- merge(A_diversity, Metadata_stats(), by= colnames(Metadata_stats())[1])
    
    A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)], function(x){
      
      A_diversity_metadata_merge <- merge(A_diversity,
                                          Metadata_stats()[, c(colnames(A_diversity_metadata)[1], x)],
                                          by= colnames(Metadata_stats())[1])
    })
    
    names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[11:ncol(A_diversity_metadata)]
    
    library(reshape2)
    i <- which(colnames(Metadata_stats())==input$metadata_alpha)-1
    A_diversity_metadata_list_melt <- melt(A_diversity_metadata_list[[i]],
                                           id.vars=c("SampleID",
                                                     names(A_diversity_metadata_list)[i]))
    
    colnames(A_diversity_metadata_list_melt)[2] <- "feature_name"
    
    A_diversity_metadata_list_melt <- filter(A_diversity_metadata_list_melt, feature_name != "NA")
    
    index <- input$select_diversity
    A_diversity_metadata_list_melt_forplot <- subset(A_diversity_metadata_list_melt,
                                                     variable==index)
    
    Dunn_result <- dunn.test::dunn.test(A_diversity_metadata_list_melt_forplot$value, A_diversity_metadata_list_melt_forplot$feature_name)
    Dunn_table <- data.frame(comparisons=Dunn_result$comparisons,
                             Z=Dunn_result$Z,
                             Pvalue= Dunn_result$P)
    Dunn_table[,"comparisons"] <- stringr::str_replace_all(Dunn_table[,"comparisons"], "-", " / ")
    Dunn_table_separate <- separate(data = Dunn_table, col = "comparisons", into = c("Group_A", "Group_B"), sep = " / ")
    order_by_count <- as.data.frame(table(Dunn_table_separate[,1]))[,1]
    Dunn_table_separate_order <- Dunn_table_separate[order(factor(Dunn_table_separate$Group_A, levels = order_by_count)),]
    colnames(Dunn_table_separate_order) <- c("Group A", "Group B", "Z", "P value")
    return(Dunn_table_separate_order)
    
  })
  
  
  BetaTable_bray <- reactive({
    
    Bray_df_data<-vegan::vegdist(t(TaxaTable_merge()), method = "bray") %>% as.matrix() %>% as.data.frame()
    return(Bray_df_data)
    
  })
  
  Beta_dsmx_heatmap <- reactive({
    
    beta_dsmx <- BetaTable_bray() %>% as.matrix()
    plot_ly(x=colnames(beta_dsmx),
            y=rownames(beta_dsmx),
            z=beta_dsmx,
            # colors = palette(50),
            type = "heatmap") %>% layout(title="Beta diversity distance matrix heatmap", xaxis=list(tickangle=45))
  })
  
  BetaPlot <- reactive({
    
    nonNA_position <- which(Metadata_stats()[, input$metadata_beta]!="NA")
    taxatable_beta <- TaxaTable_merge()[, nonNA_position]
    metadata_beta <- Metadata_stats()[nonNA_position,]
    
    library(vegan)
    Bray_df_data<-vegdist(t(taxatable_beta), method = "bray")
    
    
    # PCA
    taxatable_beta_percent <- t(t(taxatable_beta)/rowSums(t(taxatable_beta)))
    pca_Bray_df_data <- prcomp(t(taxatable_beta_percent))
    PCA_rowname <- pca_Bray_df_data$x %>% row.names()
    sample_original_names <- pca_Bray_df_data$x %>% row.names()
    
    update_rownames <- function(feature_name, metadata, i){
      
      names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
      names_order <- which(names_TF==T)
      
      sample_names <- metadata[,1][names_order] %>% as.character()
      return(sample_names)
    }
    
    names_list <- lapply(colnames(metadata_beta)[1:ncol(metadata_beta)], function(i){
      
      sapply(1:length(unique(metadata_beta[, i])), function(j){
        
        update_rownames(i, metadata_beta,j)
        
        
      })
      
    })
    
    names(names_list) <- colnames(metadata_beta)[1:ncol(metadata_beta)]
    
    for (i in 1:length(names(names_list))) {
      names(names_list[[i]]) <- unique(metadata_beta[,names(names_list)[i]]) %>% as.character()
    }
    
    metadata_beta_arrange <- arrange(metadata_beta, SampleID)
    
    PCA_rowname_arrange <- metadata_beta_arrange[, input$metadata_beta]
    
    
    library(ggplot2)
    pca_Bray_df_data_plot <- data.frame(sample=PCA_rowname_arrange, 
                                        PC1=pca_Bray_df_data$x[,1],
                                        PC2=pca_Bray_df_data$x[,2])
    pc_prop <- pca_Bray_df_data$sdev^2/sum(pca_Bray_df_data$sdev^2)
    
    library(ggrepel)
    pca_Bray_df_data_plot_gg <- ggplot(data = pca_Bray_df_data_plot, aes(x=PC1, y=PC2, label=sample_original_names, color=sample))+
      geom_point(size=1.5)+
      geom_text_repel(show.legend = FALSE)+
      xlab(paste("PC1 (", round(pc_prop[1], 2)*100, "%", ")", sep = ""))+
      ylab(paste("PC2 (", round(pc_prop[2], 2)*100, "%", ")", sep = ""))+
      geom_vline(xintercept =0, linetype="dotted")+
      geom_hline(yintercept = 0, linetype="dotted")+
      theme_bw()+
      ggtitle("PCA plot")+
      scale_colour_discrete(input$metadata_beta) + theme(text = element_text(size = 15)) 
    
    
    # PCoA
    library(ape)
    pcoa_Bray_df_data<-pcoa(Bray_df_data)
    
    # Get the name of every point in PCoA
    # Change the names by features
    PCoA_rowname <- pcoa_Bray_df_data$vectors %>% row.names()
    
    sample_original_names <- pcoa_Bray_df_data$vectors %>% row.names()
    
    update_rownames <- function(feature_name, metadata, i){
      
      names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
      names_order <- which(names_TF==T)
      
      sample_names <- metadata[,1][names_order] %>% as.character()
      return(sample_names)
    }
    
    names_list <- lapply(colnames(metadata_beta)[1:ncol(metadata_beta)], function(i){
      
      sapply(1:length(unique(metadata_beta[, i])), function(j){
        
        update_rownames(i, metadata_beta,j)
        
        
      })
      
    })
    
    names(names_list) <- colnames(metadata_beta)[1:ncol(metadata_beta)]
    
    for (i in 1:length(names(names_list))) {
      names(names_list[[i]]) <- unique(metadata_beta[,names(names_list)[i]]) %>% as.character()
    }
    
    metadata_beta_arrange <- arrange(metadata_beta, SampleID)
    
    PCoA_rowname_arrange <- metadata_beta_arrange[, input$metadata_beta]
    
    
    
    library(ggplot2)
    pcoa_Bray_df_data_plot <- data.frame(sample=PCoA_rowname_arrange, 
                                         PC1=pcoa_Bray_df_data$vectors[,1],
                                         PC2=pcoa_Bray_df_data$vectors[,2])
    
    library(ggrepel)
    pcoa_Bray_df_data_plot_gg<-ggplot(data = pcoa_Bray_df_data_plot, aes(x=PC1, y=PC2, label=sample_original_names, color=sample))+
      geom_point(size=1.5)+
      geom_text_repel(show.legend = FALSE)+
      xlab(paste("PC1 (", round(pcoa_Bray_df_data$values[1,2],2)*100, "%", ")", sep = ""))+
      ylab(paste("PC2 (", round(pcoa_Bray_df_data$values[2,2],2)*100, "%", ")", sep = ""))+
      geom_vline(xintercept =0, linetype="dotted")+
      geom_hline(yintercept = 0, linetype="dotted")+
      theme_bw()+
      ggtitle("PCoA plot")+
      scale_colour_discrete(input$metadata_beta) + theme(text = element_text(size = 15)) 
    
    # NMDS
    comm_bray <- vegdist(t(taxatable_beta))
    metaMDS_beta_df_data<-metaMDS(comm_bray, distance = "bray")
    NMDS_beta_df_data<-data.frame(NMDS1=metaMDS_beta_df_data$points[,1], NMDS2=metaMDS_beta_df_data$points[,2])
    NMDS_rowname <- PCoA_rowname_arrange
    NMDS_beta_df_data_plot<-data.frame(NMDS_beta_df_data, sample=NMDS_rowname)
    
    NMDS_beta_df_data_plot_gg<-ggplot(data = NMDS_beta_df_data_plot, aes(x=NMDS1, y=NMDS2, label=sample_original_names,color=sample))+
      geom_point(size=1.5)+
      geom_text_repel(show.legend = FALSE)+
      xlab("NMDS1")+
      ylab("NMDS2")+
      geom_vline(xintercept =0, linetype="dotted")+
      geom_hline(yintercept = 0, linetype="dotted")+
      theme_bw()+
      labs(title="NMDS plot", caption=paste("stress=", as.character(round(metaMDS_beta_df_data$stress, 4)), sep = ""))+
      #labs(caption = "A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation.")+
      scale_colour_discrete(input$metadata_beta) + theme(text = element_text(size = 15)) 
    
    
    if(input$sep == "PCoA" & input$beta_cluster == F) {
      return(pcoa_Bray_df_data_plot_gg)
    }
    else if(input$sep == "PCoA" & input$beta_cluster == T){
      return(pcoa_Bray_df_data_plot_gg + stat_ellipse(type = "t"))
    }
    else if(input$sep == "NMDS" & input$beta_cluster == F){
      return(NMDS_beta_df_data_plot_gg)
    }
    else if(input$sep == "NMDS" & input$beta_cluster == T){
      return(NMDS_beta_df_data_plot_gg + stat_ellipse(type = "t"))
    }
    else if(input$sep == "PCA" & input$beta_cluster == F){
      return(pca_Bray_df_data_plot_gg)
    }else{
      return(pca_Bray_df_data_plot_gg + stat_ellipse(type = "t"))
    }
    
    
  })
  
  output$NMDS_stress <- renderText({
    
    if(input$sep=="NMDS"){
      return("A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, > 0.2 is good/ok, and stress > 0.3 provides a poor representation.")
    }
    
  })
  # Taxatable examples for download (From QIIME2)------------------------------------------------------------------
  
  output$downloadData <- downloadHandler(
    
    filename <-"example_taxonomic_table.qza",
    
    content = function(file){
      file.copy("/home/imuser/example_files/taxtable.qza", file)
    },
    
    contentType = "application/qza"
    
  )
  
  # Metadata examples for download(From QIIME2)--------------------------------------------------------------------
  
  output$downloadMetaData <- downloadHandler(
    
    filename <-"example_metadata.tsv",
    
    content = function(file){
      file.copy("/home/imuser/example_files/sample-metadata.tsv", file)
    },
    
    contentType = "application/tsv"
    
  )
  
  
  
  # Taxatable_output---------------------------------------------------------------------------------------------------
  
  output$contents <- renderDataTable({
    
    return(TaxaTable_output())
    
  })
  
  # output$taxatable_summary <- renderText({
  #   
  #   # taxatable_summary <- data.frame(
  #   #   Kingdom = length(unique(TaxaTable_output()[, "Kingdom"])),
  #   #   Phylum = length(unique(TaxaTable_output()[, "Phylum"])),
  #   #   Class = length(unique(TaxaTable_output()[, "Class"])),
  #   #   Order = length(unique(TaxaTable_output()[, "Order"])),
  #   #   Family = length(unique(TaxaTable_output()[, "Family"])),
  #   #   Genus = length(unique(TaxaTable_output()[, "Genus"])),
  #   #   Species = length(unique(TaxaTable_output()[, "Species"]))
  #   # )
  #   # return(taxatable_summary)
  #   
  #   a <- kable(Metadata(), caption = "Hello") %>% kable_styling("striped", full_width = F) %>% add_header_above(header = c(" "=3, "A"=3, "B"=4))
  #   return(a)
  # })
  
  
  # Taxatable for download------------------------------------------------------------------------------------
  
  TaxaTable_forDL <- reactive({
    
    as_output_taxtable <- function(df_data){
      
      df_data_rownames<-row.names(df_data)
      
      df_data <- cbind(Species=df_data_rownames, df_data) %>% as.data.frame()
      
      row.names(df_data) <- NULL
      
      # remove the non-sense string
      # df_data$Species<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", "", df_data$Species)
      # df_data$Species<-gsub("k__|p__|c__|o__|f__|g__|s__", "", df_data$Species)
      
      
      library(tidyr)
      df_data <- df_data %>%
        separate(
          col = Species,
          into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
          sep = ";"
        )
      
      # replace "__" to "Unassigned"
      # df_data<-replace(df_data, df_data=="", "Unassigned")
      # df_data<-replace(df_data, df_data=="__", "Unassigned")
      
      return(as.data.frame(df_data))
    }
    a <- TaxaTable() %>% as_output_taxtable
    return(a)
  })
  
  output$downloadTaxaTable<-downloadHandler(
    
    
    filename = "TaxaTable_result.csv",
    content = function(file) {
      
      write.csv(TaxaTable_forDL(), file, row.names = FALSE)
      
    }
  )
  
  # Heatmap matrix for download---------------------------------------------------------------------------------------
  
  output$downloadHMmatrix<-downloadHandler(
    
    filename = "Heatmap_matrix.csv",
    content = function(file) {
      
      options(scipen=999)# not shown by scientific sign
      
      microbioHeatmap_matrix<-function(taxtable_data, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
        
        library(tidyr)
        species_taxtable_data<-row.names(taxtable_data) 
        taxtable_data<-as_tibble(taxtable_data) 
        taxtable_data<-add_column(taxtable_data,Species=species_taxtable_data , .before = 1) 
        
        taxtable_data<-gather(taxtable_data,
                              key = "Sample_ID",
                              value = "read_count",
                              colnames(taxtable_data[,-1]))
        
        taxtable_data<-separate(data = taxtable_data,
                                col = "Species",
                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                sep = ";")
        
        taxtable_data_Level<-taxtable_data[, c(Level,"Sample_ID","read_count")]
        
        colnames(taxtable_data_Level)[1]<-"Levels"
        
        taxtable_data_Level<-aggregate(read_count~Levels+Sample_ID,taxtable_data_Level,FUN = sum)
        
        taxtable_data_Level_spread<-spread(taxtable_data_Level,key = Sample_ID,value = read_count)
        
        row.names(taxtable_data_Level_spread)<-taxtable_data_Level_spread$Levels
        taxtable_data_Level_spread<-taxtable_data_Level_spread[,-1]
        
        taxtable_data_Level_spread<-t(taxtable_data_Level_spread) 
        taxtable_data_Level_spread_prop<-taxtable_data_Level_spread/rowSums(taxtable_data_Level_spread)
        taxtable_data_Level_spread_prop<-t(taxtable_data_Level_spread_prop)
        
        # remove non-sense string
        rownames(taxtable_data_Level_spread_prop)<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                        "", 
                                                        rownames(taxtable_data_Level_spread_prop))
        
        
        
        rownames(taxtable_data_Level_spread_prop)<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                        "", 
                                                        rownames(taxtable_data_Level_spread_prop))
        
        
        
        # replace "__" to "Unassigned"
        rownames(taxtable_data_Level_spread_prop)<-
          replace(rownames(taxtable_data_Level_spread_prop), 
                  rownames(taxtable_data_Level_spread_prop)=="__", 
                  "Unassigned")
        
        
        return(taxtable_data_Level_spread_prop)
      }
      
      if(input$select == "Kingdom") {
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Kingdom"), file, row.names = T)
      } else if (input$select == "Phylum"){
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Phylum"), file, row.names = T)
      } else if(input$select == "Class"){
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Class"), file, row.names = T)
      } else if(input$select == "Order"){
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Order"), file, row.names = T)
      } else if(input$select == "Family"){
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Family"), file, row.names = T)
      } else if(input$select == "Genus"){
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Genus"), file, row.names = T)
      }
      else {
        write.csv(microbioHeatmap_matrix(TaxaTable_merge(),"Species"), file, row.names = T)
      }
    }
  )
  
  
  # Interactive heatmap----------------------------------------------------------------------------------------
  
  output$crimeplot <- renderPlotly({
    
    options(scipen=999)# not shown by scientific sign
    
    microbioHeatmap_ly<-function(taxtable_data, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      
      library(tidyr)
      species_taxtable_data<-row.names(taxtable_data) 
      taxtable_data<-as_tibble(taxtable_data) 
      taxtable_data<-add_column(taxtable_data,Species=species_taxtable_data , .before = 1) 
      
      taxtable_data<-gather(taxtable_data,
                            key = "Sample_ID",
                            value = "read_count",
                            colnames(taxtable_data[,-1]))
      
      taxtable_data<-separate(data = taxtable_data,
                              col = "Species",
                              into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                              sep = ";")
      
      taxtable_data_Level<-taxtable_data[, c(Level,"Sample_ID","read_count")]
      
      colnames(taxtable_data_Level)[1]<-"Levels"
      
      taxtable_data_Level<-aggregate(read_count~Levels+Sample_ID,taxtable_data_Level,FUN = sum)
      
      taxtable_data_Level_spread<-spread(taxtable_data_Level,key = Sample_ID,value = read_count)
      
      row.names(taxtable_data_Level_spread)<-taxtable_data_Level_spread$Levels
      taxtable_data_Level_spread<-taxtable_data_Level_spread[,-1]
      
      taxtable_data_Level_spread<-t(taxtable_data_Level_spread) 
      taxtable_data_Level_spread_prop<-taxtable_data_Level_spread/rowSums(taxtable_data_Level_spread)
      taxtable_data_Level_spread_prop<-t(taxtable_data_Level_spread_prop)
      
      #去除物種前面的代碼
      rownames(taxtable_data_Level_spread_prop)<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                      "", 
                                                      rownames(taxtable_data_Level_spread_prop))
      
      
      
      rownames(taxtable_data_Level_spread_prop)<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                      "", 
                                                      rownames(taxtable_data_Level_spread_prop))
      
      
      
      #將"__"變為"Unassigned"
      rownames(taxtable_data_Level_spread_prop)<-
        replace(rownames(taxtable_data_Level_spread_prop), 
                rownames(taxtable_data_Level_spread_prop)=="__", 
                "Unassigned")
      
      taxtable_data_Level_spread_prop_logtrf<-log(taxtable_data_Level_spread_prop+0.01)
      
      #使sample_ID隨按鈕變動
      y_name <-taxtable_data_Level_spread_prop_logtrf %>% colnames()
      
      sample_original_names <- taxtable_data_Level_spread_prop_logtrf %>% colnames()
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        library(dplyr)
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(Metadata_stats())[1:ncol(Metadata_stats())], function(i){
        
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          
          update_rownames(i, Metadata_stats(),j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(Metadata_stats())[1:ncol(Metadata_stats())]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      
      library('plotly')
      plot_ly(x = Metadata_stats()[,input$metadata_hm],
              y = rownames(taxtable_data_Level_spread_prop_logtrf),
              z = taxtable_data_Level_spread_prop_logtrf,
              type = "heatmap"
      ) %>% layout(xaxis=list(tickangle=45))
      
      # heatmaply::heatmaply(
      #   x = taxtable_data_Level_spread_prop_logtrf
      #   
      # ) %>% layout(xaxis=list(tickangle=45))
      
    }
    
    microbioHeatmap_subly<-function(taxtable_data, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      
      library(tidyr)
      species_taxtable_data<-row.names(taxtable_data) 
      taxtable_data<-as_tibble(taxtable_data) 
      taxtable_data<-add_column(taxtable_data,Species=species_taxtable_data , .before = 1) 
      
      taxtable_data<-gather(taxtable_data,
                            key = "Sample_ID",
                            value = "read_count",
                            colnames(taxtable_data[,-1]))
      
      taxtable_data<-separate(data = taxtable_data,
                              col = "Species",
                              into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                              sep = ";")
      
      taxtable_data_Level<-taxtable_data[, c(Level,"Sample_ID","read_count")]
      
      colnames(taxtable_data_Level)[1]<-"Levels"
      
      taxtable_data_Level<-aggregate(read_count~Levels+Sample_ID,taxtable_data_Level,FUN = sum)
      
      taxtable_data_Level_spread<-spread(taxtable_data_Level,key = Sample_ID,value = read_count)
      
      row.names(taxtable_data_Level_spread)<-taxtable_data_Level_spread$Levels
      taxtable_data_Level_spread<-taxtable_data_Level_spread[,-1]
      
      taxtable_data_Level_spread<-t(taxtable_data_Level_spread) 
      taxtable_data_Level_spread_prop<-taxtable_data_Level_spread/rowSums(taxtable_data_Level_spread)
      taxtable_data_Level_spread_prop<-t(taxtable_data_Level_spread_prop)
      
      # remove non-sense string 
      rownames(taxtable_data_Level_spread_prop)<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                      "", 
                                                      rownames(taxtable_data_Level_spread_prop))
      
      
      
      rownames(taxtable_data_Level_spread_prop)<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                      "", 
                                                      rownames(taxtable_data_Level_spread_prop))
      
      
      
      # replace "__" to "Unassigned"
      rownames(taxtable_data_Level_spread_prop)<-
        replace(rownames(taxtable_data_Level_spread_prop), 
                rownames(taxtable_data_Level_spread_prop)=="__", 
                "Unassigned")
      
      taxtable_data_Level_spread_prop_logtrf<-log(taxtable_data_Level_spread_prop+0.01)
      
      #  Make sample_ID change by button
      y_name <-taxtable_data_Level_spread_prop_logtrf %>% colnames()
      
      sample_original_names <- taxtable_data_Level_spread_prop_logtrf %>% colnames()
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        library(dplyr)
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(Metadata_stats())[1:ncol(Metadata_stats())], function(i){
        
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          
          update_rownames(i, Metadata_stats(),j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(Metadata_stats())[1:ncol(Metadata_stats())]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      
      taxtable_data_Level_spread_prop_logtrf_list <- list()
      taxtable_data_Level_spread_prop_logtrf_list <- sapply(colnames(Metadata_stats()), function(i){
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          taxtable_data_Level_spread_prop_logtrf_list[[i]][[j]] <- taxtable_data_Level_spread_prop_logtrf[,names_list[[i]][[j]]]
        })
        
      })
      
      for (i in 1:length(names(names_list))) {
        names(taxtable_data_Level_spread_prop_logtrf_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      Heatmap_plotly <- function(heatmap_list, feature) {
        library('plotly')
        heatmap_plot_list <- list()
        heatmap_plot_list <- lapply(1:length(heatmap_list[[feature]]), function(i){
          
          if(i==1){
            
            plot_ly(x= colnames(heatmap_list[[feature]][[i]]),
                    y= rownames(heatmap_list[[feature]][[i]]),
                    z = heatmap_list[[feature]][[i]],
                    zmin = min(unlist(heatmap_list[[feature]])),
                    zmax = max(unlist(heatmap_list[[feature]])),
                    type = "heatmap",
                    colorbar=list(title="scale")
            ) %>% layout(title=feature,xaxis=list(title=names(heatmap_list[[feature]])[i], tickangle=45)) %>% add_trace(name=names(names_list[[feature]])[i], showscale=F)
            
          }else{
            
            plot_ly(x= colnames(heatmap_list[[feature]][[i]]),
                    y= rownames(heatmap_list[[feature]][[i]]),
                    z = heatmap_list[[feature]][[i]],
                    type = "heatmap",
                    showscale = F,
                    colorbar=list(title=names(heatmap_list[[feature]])[i])
            ) %>% layout(title=feature,xaxis=list(title=names(heatmap_list[[feature]])[i], tickangle=45)) %>% add_trace(name=names(names_list[[feature]])[i], showscale=F)
          }
          
        })
        
        return(heatmap_plot_list)
      }
      
      z_heatmap_plotly <- Heatmap_plotly(taxtable_data_Level_spread_prop_logtrf_list, input$metadata_hm)
      
      nonNA_position <- which(unique(Metadata_stats()[, input$metadata_hm])!="NA")
      z_heatmap_plotly <- z_heatmap_plotly[nonNA_position]
      
      return(subplot(z_heatmap_plotly, shareY = T, titleX = T))
      
    }
    
    if (input$metadata_hm=="SampleID"){
      return(microbioHeatmap_ly(TaxaTable_merge(), input$select))
    }
    return(microbioHeatmap_subly(TaxaTable_merge(), input$select))
    
    
  })
  
  
  # Interactive taxonomic barplot-----------------------------------------------------------------------------------------
  
  
  output$barplot <- renderPlotly({
    
    plot_LeveltoSamples<-function(taxtable, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), topN){
      
      species_taxtable<-row.names(taxtable) # get rowname
      barplot_taxa_table_data<-as_tibble(taxtable) # rawnames will be removed after transform to tibble
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) # Add new column
      # Percentage
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      # To make data tidy
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[,2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      # Seperate Species names by taxon
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      # Replace "__" to "Unassigned"
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      # Get Level
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      # it's convinient to make function by change column
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      #factor
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID,barplot_taxa_table_data_percent_Level,FUN = sum)
      
      #
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels != "Unassigned")
      
      TopN_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        topN <- topN
        species_name <- a$Levels %>% as.character()
        species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
        species_name_topN <- species_name_sort[1:topN]
        return(species_name_topN)
      })
      
      names(TopN_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      others_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        topN <- topN
        species_name <- a$Levels %>% as.character()
        species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
        species_name_topN <- species_name_sort[-c(1:topN)]
        return(species_name_topN)
      })
      names(others_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      others_list_1 <- lapply(1:length(others_list), function(i){
        c("Unassigned", others_list[[i]])
      })
      names(others_list_1) <- names(others_list)
      
      others_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% others_list_1[[i]]) %>% filter(Sample_ID==names(others_list_1)[i])
      })
      names(others_taxtable_list) <- names(others_list_1)
      
      others_taxtable_list_sum <- lapply(1:length(others_taxtable_list), function(i){
        data.frame(Levels = "Others", 
                   Sample_ID = names(others_taxtable_list)[i], 
                   read_percentage = sum(others_taxtable_list[[i]][,3]))
      })
      
      topN_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% TopN_list[[i]]) %>% filter(Sample_ID==names(TopN_list)[i])
      })
      
      append_taxtable_list <- lapply(1:length(others_taxtable_list_sum), function(i){
        bind_rows(topN_taxtable_list[[i]], others_taxtable_list_sum[[i]])
      })
      
      append_taxtable_list_all <- bind_rows(append_taxtable_list)
      #
      
      y_name <-append_taxtable_list_all$Sample_ID 
      
      sample_original_names <- append_taxtable_list_all$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(Metadata_stats())[1:ncol(Metadata_stats())], function(i){
        
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          
          update_rownames(i, Metadata_stats(),j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(Metadata_stats())[1:ncol(Metadata_stats())]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      
      for (k in 1:length(names(names_list[[input$metadata_barplot]]))) {
        
        
        y_name <- mapvalues(y_name, 
                            from = names_list[[input$metadata_barplot]][[k]], 
                            to = as.character(
                              rep(
                                unique(
                                  Metadata_stats()[, input$metadata_barplot])[k], 
                                length(names_list[[input$metadata_barplot]][[k]])
                              )
                            )
        )}
      
      if (unique(append_taxtable_list_all$Sample_ID == y_name)) {
        Samples_ID <- append_taxtable_list_all$Sample_ID
      }else{
        Samples_ID <- paste0( y_name, "__",append_taxtable_list_all$Sample_ID)
      }
      
      factor_levels <- append_taxtable_list_all$Levels %>% unique()
      others_position <- which(factor_levels=="Others")
      factor_levels_new <- c(as.character(factor_levels[others_position]), as.character(factor_levels[-others_position]))
      factor_levels_new <- factor(factor_levels_new, levels = factor_levels_new)
      append_taxtable_list_all$Levels <- factor(append_taxtable_list_all$Levels,
                                                levels = factor_levels_new
      )
      
      number_nameslist <- Metadata_stats()[, 1] %>% length()
      
      append_taxtable_list_all$read_percentage <- append_taxtable_list_all$read_percentage*100
      
      # plot
      library(viridis)
      p_Level<-ggplot()+geom_bar(data = append_taxtable_list_all,
                                 aes(x=Samples_ID,
                                     y=read_percentage,
                                     fill=Levels),
                                 colour="White",
                                 size=0.1,
                                 stat = "identity")+theme_bw()+theme(legend.title=element_blank())+ labs(x=paste0("Sample ID", " (n=", number_nameslist,")"), 
                                                                                                         y="Relative abundance (%)",
                                                                                                         fill="")
      
      library("plotly")
      ggplotly(p_Level) %>% layout(xaxis=list(tickangle=45))
    }
    
    plot_LeveltoSamples_noUnassigned<-function(taxtable, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), topN){
      
      topN_species_list <- lapply(1:ncol(taxtable), function(i){
        
        sort(taxtable[,i], decreasing = T)[1:topN] %>% names()
      })
      
      topN_species_union <- unique(unlist(topN_species_list))
      
      new_taxtable <- taxtable[topN_species_union,]
      
      species_taxtable<-row.names(new_taxtable) # get rowname
      barplot_taxa_table_data<-as_tibble(new_taxtable) # rawnames will be removed after transform to tibble
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) # Add new column
      # Percentage
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      # To make data tidy
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[, 2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      # Seperate Species names by taxon
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      # Replace "__" to "Unassigned"
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      # Get Level
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      # it's convinient to make function by change column
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      # factor
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID, barplot_taxa_table_data_percent_Level, FUN = sum)
      
      # remove unassigned
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels!="Unassigned")
      
      # compute percentage again
      h_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(i){
        
        filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == i)
        
      })
      
      for (i in 1:length(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID))) {
        
        h_list[[i]][, "read_percentage"] <- h_list[[i]][, "read_percentage"]/sum(h_list[[i]][, "read_percentage"])
        
      }
      
      
      barplot_taxa_table_data_percent_Level_noUnassigned_percentage <- bind_rows(h_list)
      
      
      y_name <-barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID 
      
      sample_original_names <- barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(Metadata_stats())[1:ncol(Metadata_stats())], function(i){
        
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          
          update_rownames(i, Metadata_stats(),j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(Metadata_stats())[1:ncol(Metadata_stats())]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      
      for (k in 1:length(names(names_list[[input$metadata_barplot]]))) {
        
        
        y_name <- mapvalues(y_name, 
                            from = names_list[[input$metadata_barplot]][[k]], 
                            to = as.character(
                              rep(
                                unique(
                                  Metadata_stats()[, input$metadata_barplot])[k], 
                                length(names_list[[input$metadata_barplot]][[k]])
                              )
                            )
        )}
      
      if (unique(barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID == y_name)) {
        Samples_ID <- barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID
      }else{
        Samples_ID <- paste0( y_name, "__", barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID)
      }
      
      # plot
      library(viridis)
      p_Level<-ggplot()+geom_bar(data = barplot_taxa_table_data_percent_Level_noUnassigned_percentage,
                                 aes(x=Samples_ID,
                                     y=read_percentage,
                                     fill=Levels),
                                 colour="White",
                                 size=0.1,
                                 stat = "identity")+theme_bw()+theme(legend.title=element_blank())
      # + scale_y_continuous(labels = scales::percent_format())
      
      library("plotly")
      ggplotly(p_Level) %>% layout(xaxis=list(tickangle=45))
    }
    
    plot_LeveltoSamples_sub<-function(taxtable, metadata,features, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), topN){
      
      # topN_species_list <- lapply(1:ncol(taxtable), function(i){
      #   
      #   sort(taxtable[,i], decreasing = T)[1:topN] %>% names()
      # })
      # 
      # topN_species_union <- unique(unlist(topN_species_list))
      # 
      # new_taxtable <- taxtable[topN_species_union,]
      
      species_taxtable<-row.names(taxtable) 
      barplot_taxa_table_data<-as_tibble(taxtable) 
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) 
      
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[,2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID,barplot_taxa_table_data_percent_Level,FUN = sum)
      
      #
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels != "Unassigned")
      
      TopN_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        topN <- topN
        species_name <- a$Levels %>% as.character()
        species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
        species_name_topN <- species_name_sort[1:topN]
        return(species_name_topN)
      })
      
      names(TopN_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      species_union <- unique(unlist(TopN_list))
      species_Levels <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Levels) %>% as.character()
      species_diff <- setdiff(species_Levels, species_union)
      
      # others_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
      #   a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
      #   topN <- topN
      #   species_name <- a$Levels %>% as.character()
      #   species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
      #   species_name_topN <- species_name_sort[-c(1:topN)]
      #   return(species_name_topN)
      # })
      '%!in%' <- function(x,y)!('%in%'(x,y))
      others_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        b <- filter(a, Levels %!in% species_union)
        b$Levels <- as.character(b$Levels)
        return(b)
      })
      names(others_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      others_list_1 <- lapply(1:length(others_list), function(i){
        c("Unassigned", others_list[[i]]$Levels)
      })
      names(others_list_1) <- names(others_list)
      
      others_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% others_list_1[[i]]) %>% filter(Sample_ID==names(others_list_1)[i])
      })
      names(others_taxtable_list) <- names(others_list_1)
      
      others_taxtable_list_sum <- lapply(1:length(others_taxtable_list), function(i){
        data.frame(Levels = "Others", 
                   Sample_ID = names(others_taxtable_list)[i], 
                   read_percentage = sum(others_taxtable_list[[i]][,"read_percentage"]))
      })
      
      # topN_taxtable_list <- lapply(1:length(others_list_1), function(i){
      #   filter(barplot_taxa_table_data_percent_Level, Levels %in% TopN_list[[i]]) %>% filter(Sample_ID==names(TopN_list)[i])
      # })
      topN_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% species_union) %>% filter(Sample_ID==names(TopN_list)[i])
      })
      
      append_taxtable_list <- lapply(1:length(others_taxtable_list_sum), function(i){
        bind_rows(topN_taxtable_list[[i]], others_taxtable_list_sum[[i]])
      })
      
      append_taxtable_list_all <- bind_rows(append_taxtable_list)
      #
      
      
      y_name <-append_taxtable_list_all$Sample_ID 
      
      sample_original_names <- append_taxtable_list_all$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(metadata)[1:ncol(metadata)], function(i){
        
        sapply(1:length(unique(metadata[, i])), function(j){
          
          update_rownames(i, metadata,j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(metadata[,names(names_list)[i]]) %>% as.character()
      }
      
      
      factor_levels <- append_taxtable_list_all$Levels %>% unique()
      others_position <- which(factor_levels=="Others")
      factor_levels_new <- c(as.character(factor_levels[others_position]), as.character(factor_levels[-others_position]))
      factor_levels_new <- factor(factor_levels_new, levels = factor_levels_new)
      append_taxtable_list_all$Levels <- factor(append_taxtable_list_all$Levels,
                                                levels = factor_levels_new
      )
      
      
      p_Level_plot_sub <- function(barplot_taxa_table_data_percent_Level, features, features_sub){
        
        
        # barplot_taxa_table_data_percent_Level_sub <- filter(barplot_taxa_table_data_percent_Level, Sample_ID %in% names_list[[features]][[features_sub]]) 
        sub_position <- which(barplot_taxa_table_data_percent_Level[, "Sample_ID"] %in% names_list[[features]][[features_sub]])
        barplot_taxa_table_data_percent_Level_sub <- barplot_taxa_table_data_percent_Level[sub_position,]
        
        barplot_taxa_table_data_percent_Level_sub <- barplot_taxa_table_data_percent_Level_sub %>% group_by(Levels)
        
        number_nameslist <- length(names_list[[features]][[features_sub]]) %>% as.character()
        
        barplot_taxa_table_data_percent_Level_sub$read_percentage <- barplot_taxa_table_data_percent_Level_sub$read_percentage*100
        
        library(viridis)
        
        p_Level_plot<-ggplot()+geom_bar(data = barplot_taxa_table_data_percent_Level_sub,
                                        aes(x=Sample_ID,
                                            y=read_percentage,
                                            fill=Levels),
                                        colour="White",
                                        size=0.1,
                                        # show.legend = T,
                                        stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x=paste0(features_sub, " (n=", number_nameslist,")"), 
                                                                                                                                          y="Relative abundance (%)",
                                                                                                                                          fill="")
        
        
      }
      
      p_Level_plot_list <- list()
      p_Level_plot_list <- lapply(names(names_list), function(i){
        lapply(names(names_list[[i]]), function(j){
          
          p_Level_plot_sub(append_taxtable_list_all, i, j) 
          
        })
        
      })
      
      names(p_Level_plot_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(p_Level_plot_list))) {
        names(p_Level_plot_list[[i]]) <- unique(metadata[,names(p_Level_plot_list)[i]]) %>% as.character()
      }
      
      for (i in 1:length(names(names_list))) {
        for (j in 2:length(names(names_list[[i]]))) {
          p_Level_plot_list[[i]][[j]] <- p_Level_plot_list[[i]][[j]]
        }
        
      }
      
      for (i in 2:length(p_Level_plot_list[[features]])) {
        
        p_Level_plot_list[[features]][[i]] <- style(p_Level_plot_list[[features]][[i]], showlegend = F)
        
      }
      
      nonNA_position <- which(names(p_Level_plot_list[[features]])!="NA")
      p_Level_plot_list[[features]] <- p_Level_plot_list[[features]][nonNA_position]
      
      subplot(p_Level_plot_list[[features]],
              shareY = T,
              titleX = T)
      
    }
    
    plot_LeveltoSamples_sub_noUnassigned<-function(taxtable, metadata,features, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), topN){
      
      topN_species_list <- lapply(1:ncol(taxtable), function(i){
        
        sort(taxtable[,i], decreasing = T)[1:topN] %>% names()
      })
      
      topN_species_union <- unique(unlist(topN_species_list))
      
      new_taxtable <- taxtable[topN_species_union,]
      
      species_taxtable<-row.names(new_taxtable) 
      barplot_taxa_table_data<-as_tibble(new_taxtable) 
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) 
      
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[,2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID,barplot_taxa_table_data_percent_Level,FUN = sum)
      
      # remove unassigned
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels!="Unassigned")
      
      # compute percentage again
      h_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(i){
        
        filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == i)
        
      })
      
      for (i in 1:length(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID))) {
        
        h_list[[i]][, "read_percentage"] <- h_list[[i]][, "read_percentage"]/sum(h_list[[i]][, "read_percentage"])
        
      }
      
      
      barplot_taxa_table_data_percent_Level_noUnassigned_percentage <- bind_rows(h_list)
      
      
      y_name <-barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID 
      
      sample_original_names <- barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(metadata)[1:ncol(metadata)], function(i){
        
        sapply(1:length(unique(metadata[, i])), function(j){
          
          update_rownames(i, metadata,j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(metadata[,names(names_list)[i]]) %>% as.character()
      }
      
      
      
      
      
      p_Level_plot_sub <- function(barplot_taxa_table_data_percent_Level_noUnassigned_percentage, features, features_sub){
        
        barplot_taxa_table_data_percent_Level_sub <- subset(barplot_taxa_table_data_percent_Level_noUnassigned_percentage, Sample_ID %in% names_list[[features]][[features_sub]])  
        
        barplot_taxa_table_data_percent_Level_sub <- barplot_taxa_table_data_percent_Level_sub %>% group_by(Levels)
        
        
        library(viridis)
        
        p_Level_plot<-ggplot()+geom_bar(data = barplot_taxa_table_data_percent_Level_sub,
                                        aes(x=Sample_ID,
                                            y=read_percentage,
                                            fill=Levels),
                                        colour="White",
                                        size=0.1,
                                        # show.legend = T,
                                        stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x=features_sub, fill="")
        
        
      }
      
      p_Level_plot_list <- list()
      p_Level_plot_list <- lapply(names(names_list), function(i){
        lapply(names(names_list[[i]]), function(j){
          
          p_Level_plot_sub(barplot_taxa_table_data_percent_Level_noUnassigned_percentage, i, j) 
          
        })
        
      })
      
      names(p_Level_plot_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(p_Level_plot_list))) {
        names(p_Level_plot_list[[i]]) <- unique(metadata[,names(p_Level_plot_list)[i]]) %>% as.character()
      }
      
      for (i in 1:length(names(names_list))) {
        for (j in 2:length(names(names_list[[i]]))) {
          p_Level_plot_list[[i]][[j]] <- p_Level_plot_list[[i]][[j]]
        }
        
      }
      
      for (i in 2:length(p_Level_plot_list[[features]])) {
        
        p_Level_plot_list[[features]][[i]] <- style(p_Level_plot_list[[features]][[i]], showlegend = F)
        
      }
      
      nonNA_position <- which(names(p_Level_plot_list[[features]])!="NA")
      p_Level_plot_list[[features]] <- p_Level_plot_list[[features]][nonNA_position]
      
      subplot(p_Level_plot_list[[features]],
              shareY = T,
              titleX = T)
      
    }
    
    # if (input$metadata_barplot=="SampleID" | input$metadata_barplot=="barcode.sequence"){
    #   if (input$without_unassigned == F){
    #     taxtable <- read_qza(input$taxonomic_table$datapath)[["data"]]
    #     return(plot_LeveltoSamples(taxtable, input$select2, input$integer))
    #   }else{
    #     taxtable <- read_qza(input$taxonomic_table$datapath)[["data"]]
    #     return(plot_LeveltoSamples_noUnassigned(taxtable, input$select2, input$integer))
    #   }
    #     
    # }else{
    #   
    #   if(input$without_unassigned == F){
    #     taxtable <- read_qza(input$taxonomic_table$datapath)[["data"]]
    #     return(plot_LeveltoSamples_sub(taxtable, Metadata_stats(),input$metadata_barplot, input$select2, input$integer))
    #   }else{
    #     taxtable <- read_qza(input$taxonomic_table$datapath)[["data"]]
    #     return(plot_LeveltoSamples_sub_noUnassigned(taxtable, Metadata_stats(),input$metadata_barplot, input$select2, input$integer))
    #   }
    #     
    # }
    
    if (input$metadata_barplot=="SampleID"){
      
      # taxtable <- read_qza(input$taxonomic_table$datapath)[["data"]]
      return(plot_LeveltoSamples(TaxaTable(), input$select2, input$integer))
      
    }else{
      
      # taxtable <- read_qza(input$taxonomic_table$datapath)[["data"]]
      return(plot_LeveltoSamples_sub(TaxaTable(), Metadata_stats(),input$metadata_barplot, input$select2, input$integer))
      
    }
    
  })
  
  
  # Interactive taxonomic barplot for download-----------------------------------------------------------------------------
  
  barplot_download <- reactive({
    
    plot_LeveltoSamples_save<-function(taxtable, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), topN){
      
      species_taxtable<-row.names(taxtable) # get rowname
      barplot_taxa_table_data<-as_tibble(taxtable) # rawnames will be removed after transform to tibble
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) # Add new column
      # Percentage
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      # To make data tidy
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[,2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      # Seperate Species names by taxon
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      # Replace "__" to "Unassigned"
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      # Get Level
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      # it's convinient to make function by change column
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      #factor
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID,barplot_taxa_table_data_percent_Level,FUN = sum)
      
      #
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels != "Unassigned")
      
      TopN_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        topN <- topN
        species_name <- a$Levels %>% as.character()
        species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
        species_name_topN <- species_name_sort[1:topN]
        return(species_name_topN)
      })
      
      names(TopN_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      others_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        topN <- topN
        species_name <- a$Levels %>% as.character()
        species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
        species_name_topN <- species_name_sort[-c(1:topN)]
        return(species_name_topN)
      })
      names(others_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      others_list_1 <- lapply(1:length(others_list), function(i){
        c("Unassigned", others_list[[i]])
      })
      names(others_list_1) <- names(others_list)
      
      others_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% others_list_1[[i]]) %>% filter(Sample_ID==names(others_list_1)[i])
      })
      names(others_taxtable_list) <- names(others_list_1)
      
      others_taxtable_list_sum <- lapply(1:length(others_taxtable_list), function(i){
        data.frame(Levels = "Others", 
                   Sample_ID = names(others_taxtable_list)[i], 
                   read_percentage = sum(others_taxtable_list[[i]][,3]))
      })
      
      topN_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% TopN_list[[i]]) %>% filter(Sample_ID==names(TopN_list)[i])
      })
      
      append_taxtable_list <- lapply(1:length(others_taxtable_list_sum), function(i){
        bind_rows(topN_taxtable_list[[i]], others_taxtable_list_sum[[i]])
      })
      
      append_taxtable_list_all <- bind_rows(append_taxtable_list)
      #
      
      y_name <-append_taxtable_list_all$Sample_ID 
      
      sample_original_names <- append_taxtable_list_all$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(Metadata_stats())[1:ncol(Metadata_stats())], function(i){
        
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          
          update_rownames(i, Metadata_stats(),j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(Metadata_stats())[1:ncol(Metadata_stats())]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      
      for (k in 1:length(names(names_list[[input$metadata_barplot]]))) {
        
        
        y_name <- mapvalues(y_name, 
                            from = names_list[[input$metadata_barplot]][[k]], 
                            to = as.character(
                              rep(
                                unique(
                                  Metadata_stats()[, input$metadata_barplot])[k], 
                                length(names_list[[input$metadata_barplot]][[k]])
                              )
                            )
        )}
      
      if (unique(append_taxtable_list_all$Sample_ID == y_name)) {
        Samples_ID <- append_taxtable_list_all$Sample_ID
      }else{
        Samples_ID <- paste0( y_name, "__",append_taxtable_list_all$Sample_ID)
      }
      
      factor_levels <- append_taxtable_list_all$Levels %>% unique()
      others_position <- which(factor_levels=="Others")
      factor_levels_new <- c(as.character(factor_levels[others_position]), as.character(factor_levels[-others_position]))
      factor_levels_new <- factor(factor_levels_new, levels = factor_levels_new)
      append_taxtable_list_all$Levels <- factor(append_taxtable_list_all$Levels,
                                                levels = factor_levels_new
      )
      
      number_nameslist <- Metadata_stats()[, 1] %>% length()
      
      append_taxtable_list_all$read_percentage <- append_taxtable_list_all$read_percentage*100
      
      # plot
      library(viridis)
      p_Level<-ggplot()+geom_bar(data = append_taxtable_list_all,
                                 aes(x=Samples_ID,
                                     y=read_percentage,
                                     fill=Levels),
                                 colour="White",
                                 size=0.1,
                                 stat = "identity")+theme_bw()+theme(legend.title=element_blank())+ labs(x=paste0("Sample ID", " (n=", number_nameslist,")"), 
                                                                                                         y="Relative abundance (%)",
                                                                                                         fill="")
      
      p_Level + guides(fill = guide_legend(nrow = 40, byrow = TRUE))
    }
    
    return(plot_LeveltoSamples_save(TaxaTable_merge(), input$select2, input$integer))
    
  })
  
  barplot_noUnassigned_download <- reactive({
    
    plot_LeveltoSamples_noUnassigned_save<-function(taxtable, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      
      species_taxtable<-row.names(taxtable) # get rowname
      barplot_taxa_table_data<-as_tibble(taxtable) # rawnames will be removed after transform to tibble
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) # Add new column
      # Percentage
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      # To make data tidy
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[, 2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      # Seperate Species names by taxon
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      # Replace "__" to "Unassigned"
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      # Get Level
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      # it's convinient to make function by change column
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      # factor
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID, barplot_taxa_table_data_percent_Level, FUN = sum)
      
      # remove unassigned
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels!="Unassigned")
      
      # compute percentage again
      h_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(i){
        
        filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == i)
        
      })
      
      for (i in 1:length(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID))) {
        
        h_list[[i]][, "read_percentage"] <- h_list[[i]][, "read_percentage"]/sum(h_list[[i]][, "read_percentage"])
        
      }
      
      
      barplot_taxa_table_data_percent_Level_noUnassigned_percentage <- bind_rows(h_list)
      
      
      y_name <-barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID 
      
      sample_original_names <- barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(Metadata_stats())[1:ncol(Metadata_stats())], function(i){
        
        sapply(1:length(unique(Metadata_stats()[, i])), function(j){
          
          update_rownames(i, Metadata_stats(),j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(Metadata_stats())[1:ncol(Metadata_stats())]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(Metadata_stats()[,names(names_list)[i]]) %>% as.character()
      }
      
      
      for (k in 1:length(names(names_list[[input$metadata_barplot]]))) {
        
        
        y_name <- mapvalues(y_name, 
                            from = names_list[[input$metadata_barplot]][[k]], 
                            to = as.character(
                              rep(
                                unique(
                                  Metadata_stats()[, input$metadata_barplot])[k], 
                                length(names_list[[input$metadata_barplot]][[k]])
                              )
                            )
        )}
      
      if (unique(barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID == y_name)) {
        Samples_ID <- barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID
      }else{
        Samples_ID <- paste0( y_name, "__", barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID)
      }
      
      # plot
      library(viridis)
      p_Level<-ggplot()+geom_bar(data = barplot_taxa_table_data_percent_Level_noUnassigned_percentage,
                                 aes(x=Samples_ID,
                                     y=read_percentage,
                                     fill=Levels),
                                 colour="White",
                                 size=0.1,
                                 stat = "identity")+theme_bw()+theme(legend.title=element_blank())+theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 
      # + scale_y_continuous(labels = scales::percent_format())
      
      p_Level + guides(fill = guide_legend(nrow = 40, byrow = TRUE))
    }
    
    return(plot_LeveltoSamples_noUnassigned_save(TaxaTable_merge(), input$select2))
  })
  
  barplot_sub_download <- reactive({
    
    plot_LeveltoSamples_sub<-function(taxtable, metadata, features, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), topN){
      
      # topN_species_list <- lapply(1:ncol(taxtable), function(i){
      #   
      #   sort(taxtable[,i], decreasing = T)[1:topN] %>% names()
      # })
      # 
      # topN_species_union <- unique(unlist(topN_species_list))
      # 
      # new_taxtable <- taxtable[topN_species_union,]
      
      species_taxtable<-row.names(taxtable) 
      barplot_taxa_table_data<-as_tibble(taxtable) 
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) 
      
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[,2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID,barplot_taxa_table_data_percent_Level,FUN = sum)
      
      #
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels != "Unassigned")
      
      TopN_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        topN <- topN
        species_name <- a$Levels %>% as.character()
        species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
        species_name_topN <- species_name_sort[1:topN]
        return(species_name_topN)
      })
      
      names(TopN_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      species_union <- unique(unlist(TopN_list))
      species_Levels <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Levels) %>% as.character()
      species_diff <- setdiff(species_Levels, species_union)
      
      # others_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
      #   a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
      #   topN <- topN
      #   species_name <- a$Levels %>% as.character()
      #   species_name_sort <- species_name[order(a$read_percentage, decreasing = T)]
      #   species_name_topN <- species_name_sort[-c(1:topN)]
      #   return(species_name_topN)
      # })
      '%!in%' <- function(x,y)!('%in%'(x,y))
      others_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(x){
        a <- filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == x)
        b <- filter(a, Levels %!in% species_union)
        b$Levels <- as.character(b$Levels)
        return(b)
      })
      names(others_list) <- unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID)
      
      others_list_1 <- lapply(1:length(others_list), function(i){
        c("Unassigned", others_list[[i]]$Levels)
      })
      names(others_list_1) <- names(others_list)
      
      others_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% others_list_1[[i]]) %>% filter(Sample_ID==names(others_list_1)[i])
      })
      names(others_taxtable_list) <- names(others_list_1)
      
      others_taxtable_list_sum <- lapply(1:length(others_taxtable_list), function(i){
        data.frame(Levels = "Others", 
                   Sample_ID = names(others_taxtable_list)[i], 
                   read_percentage = sum(others_taxtable_list[[i]][,"read_percentage"]))
      })
      
      # topN_taxtable_list <- lapply(1:length(others_list_1), function(i){
      #   filter(barplot_taxa_table_data_percent_Level, Levels %in% TopN_list[[i]]) %>% filter(Sample_ID==names(TopN_list)[i])
      # })
      topN_taxtable_list <- lapply(1:length(others_list_1), function(i){
        filter(barplot_taxa_table_data_percent_Level, Levels %in% species_union) %>% filter(Sample_ID==names(TopN_list)[i])
      })
      
      append_taxtable_list <- lapply(1:length(others_taxtable_list_sum), function(i){
        bind_rows(topN_taxtable_list[[i]], others_taxtable_list_sum[[i]])
      })
      
      append_taxtable_list_all <- bind_rows(append_taxtable_list)
      #
      
      
      y_name <-append_taxtable_list_all$Sample_ID 
      
      sample_original_names <- append_taxtable_list_all$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(metadata)[1:ncol(metadata)], function(i){
        
        sapply(1:length(unique(metadata[, i])), function(j){
          
          update_rownames(i, metadata,j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(metadata[,names(names_list)[i]]) %>% as.character()
      }
      
      
      factor_levels <- append_taxtable_list_all$Levels %>% unique()
      others_position <- which(factor_levels=="Others")
      factor_levels_new <- c(as.character(factor_levels[others_position]), as.character(factor_levels[-others_position]))
      factor_levels_new <- factor(factor_levels_new, levels = factor_levels_new)
      append_taxtable_list_all$Levels <- factor(append_taxtable_list_all$Levels,
                                                levels = factor_levels_new
      )
      
      
      p_Level_plot_sub <- function(barplot_taxa_table_data_percent_Level, features, features_sub){
        
        
        # barplot_taxa_table_data_percent_Level_sub <- filter(barplot_taxa_table_data_percent_Level, Sample_ID %in% names_list[[features]][[features_sub]]) 
        sub_position <- which(barplot_taxa_table_data_percent_Level[, "Sample_ID"] %in% names_list[[features]][[features_sub]])
        barplot_taxa_table_data_percent_Level_sub <- barplot_taxa_table_data_percent_Level[sub_position,]
        
        barplot_taxa_table_data_percent_Level_sub <- barplot_taxa_table_data_percent_Level_sub %>% group_by(Levels)
        
        number_nameslist <- length(names_list[[features]][[features_sub]]) %>% as.character()
        
        barplot_taxa_table_data_percent_Level_sub$read_percentage <- barplot_taxa_table_data_percent_Level_sub$read_percentage*100
        
        library(viridis)
        
        p_Level_plot<-ggplot()+geom_bar(data = barplot_taxa_table_data_percent_Level_sub,
                                        aes(x=Sample_ID,
                                            y=read_percentage,
                                            fill=Levels),
                                        colour="White",
                                        size=0.1,
                                        # show.legend = T,
                                        stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x=paste0(features_sub, " (n=", number_nameslist,")"), 
                                                                                                                                          y="Relative abundance (%)",
                                                                                                                                          fill="")
        
        
      }
      
      p_Level_plot_list <- list()
      p_Level_plot_list <- lapply(names(names_list), function(i){
        lapply(names(names_list[[i]]), function(j){
          
          p_Level_plot_sub(append_taxtable_list_all, i, j) 
          
        })
        
      })
      
      names(p_Level_plot_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(p_Level_plot_list))) {
        names(p_Level_plot_list[[i]]) <- unique(metadata[,names(p_Level_plot_list)[i]]) %>% as.character()
      }
      
      for (i in 1:length(names(names_list))) {
        for (j in 2:length(names(names_list[[i]]))) {
          p_Level_plot_list[[i]][[j]] <- p_Level_plot_list[[i]][[j]]
        }
        
      }
      
      for (i in 1:(length(p_Level_plot_list[[features]])-1)) {
        
        p_Level_plot_list[[features]][[i]] <- p_Level_plot_list[[features]][[i]] + theme(legend.position="none")
        
      }
      
      nonNA_position <- which(names(p_Level_plot_list[[features]])!="NA")
      p_Level_plot_list[[features]] <- p_Level_plot_list[[features]][nonNA_position]
      
      egg::ggarrange(plots = p_Level_plot_list[[features]], nrow=1)
      
    }
    
    return(plot_LeveltoSamples_sub(taxtable = TaxaTable_merge(), metadata = Metadata_stats(), features = input$metadata_barplot, Level = input$select2, topN = input$integer))
  })
  
  barplot_noUnassigned_sub_download <- reactive({
    
    plot_LeveltoSamples_sub_noUnassigned<-function(taxtable, metadata,features, Level=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      
      species_taxtable<-row.names(taxtable) 
      barplot_taxa_table_data<-as_tibble(taxtable) 
      barplot_taxa_table_data<-add_column(barplot_taxa_table_data,Species=species_taxtable) 
      
      barplot_taxa_table_data_percent<-t(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)])/
                                           rowSums(t(barplot_taxa_table_data[,-ncol(barplot_taxa_table_data)]))) %>% as_tibble()
      barplot_taxa_table_data_percent<-add_column(barplot_taxa_table_data_percent, Species_name=barplot_taxa_table_data$Species, .before = 1)
      
      
      library(tidyr)
      barplot_taxa_table_data_percent<-gather(barplot_taxa_table_data_percent,
                                              key = "Sample_ID",
                                              value = "read_percentage",
                                              colnames(barplot_taxa_table_data_percent[,2:ncol(barplot_taxa_table_data_percent)]))
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      barplot_taxa_table_data_percent$Species_name<-gsub("k__|p__|c__|o__|f__|g__|s__", 
                                                         "", 
                                                         barplot_taxa_table_data_percent$Species_name)
      
      
      
      barplot_taxa_table_data_percent<-separate(data = barplot_taxa_table_data_percent,
                                                col = "Species_name",
                                                into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                                sep = ";")
      
      
      
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="__",
                                               "Unassigned")
      
      barplot_taxa_table_data_percent<-replace(barplot_taxa_table_data_percent,
                                               barplot_taxa_table_data_percent=="",
                                               "Unassigned")
      
      
      barplot_taxa_table_data_percent_Level<-barplot_taxa_table_data_percent[,c(Level,"Sample_ID","read_percentage")]
      
      
      colnames(barplot_taxa_table_data_percent_Level)[1]<-"Levels"
      
      
      barplot_taxa_table_data_percent_Level$Levels<-factor(barplot_taxa_table_data_percent_Level$Levels,
                                                           levels = unique(c("Unassigned",barplot_taxa_table_data_percent_Level$Levels))
      )
      
      barplot_taxa_table_data_percent_Level<-aggregate(read_percentage~Levels+Sample_ID,barplot_taxa_table_data_percent_Level,FUN = sum)
      
      # remove unassigned
      barplot_taxa_table_data_percent_Level_noUnassigned <- filter(barplot_taxa_table_data_percent_Level, Levels!="Unassigned")
      
      # compute percentage again
      h_list <- lapply(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID), function(i){
        
        filter(barplot_taxa_table_data_percent_Level_noUnassigned, Sample_ID == i)
        
      })
      
      for (i in 1:length(unique(barplot_taxa_table_data_percent_Level_noUnassigned$Sample_ID))) {
        
        h_list[[i]][, "read_percentage"] <- h_list[[i]][, "read_percentage"]/sum(h_list[[i]][, "read_percentage"])
        
      }
      
      
      barplot_taxa_table_data_percent_Level_noUnassigned_percentage <- bind_rows(h_list)
      
      
      y_name <-barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID 
      
      sample_original_names <- barplot_taxa_table_data_percent_Level_noUnassigned_percentage$Sample_ID
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- sapply(colnames(metadata)[1:ncol(metadata)], function(i){
        
        sapply(1:length(unique(metadata[, i])), function(j){
          
          update_rownames(i, metadata,j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(metadata[,names(names_list)[i]]) %>% as.character()
      }
      
      
      
      
      
      p_Level_plot_sub <- function(barplot_taxa_table_data_percent_Level_noUnassigned_percentage, features, features_sub){
        
        barplot_taxa_table_data_percent_Level_sub <- subset(barplot_taxa_table_data_percent_Level_noUnassigned_percentage, Sample_ID %in% names_list[[features]][[features_sub]])  
        
        barplot_taxa_table_data_percent_Level_sub <- barplot_taxa_table_data_percent_Level_sub %>% group_by(Levels)
        
        
        library(viridis)
        
        p_Level_plot<-ggplot()+geom_bar(data = barplot_taxa_table_data_percent_Level_sub,
                                        aes(x=Sample_ID,
                                            y=read_percentage,
                                            fill=Levels),
                                        colour="White",
                                        size=0.1,
                                        show.legend = T,
                                        stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x=features_sub, fill="")
        
        p_Level_plot <- p_Level_plot + guides(fill = guide_legend(nrow = 40, byrow = TRUE))
        
      }
      
      p_Level_plot_list <- list()
      p_Level_plot_list <- lapply(names(names_list), function(i){
        lapply(names(names_list[[i]]), function(j){
          
          p_Level_plot_sub(barplot_taxa_table_data_percent_Level_noUnassigned_percentage, i, j) 
          
        })
        
      })
      
      names(p_Level_plot_list) <- colnames(metadata)[1:ncol(metadata)]
      
      for (i in 1:length(names(p_Level_plot_list))) {
        names(p_Level_plot_list[[i]]) <- unique(metadata[,names(p_Level_plot_list)[i]]) %>% as.character()
      }
      
      for (i in 1:length(names(names_list))) {
        for (j in 2:length(names(names_list[[i]]))) {
          p_Level_plot_list[[i]][[j]] <- p_Level_plot_list[[i]][[j]]
        }
        
      }
      
      for (i in 1:length(names(names_list))) {
        
        for (j in 1:(length(names(names_list[[i]]))-1)){
          
          p_Level_plot_list[[i]][[j]] <- p_Level_plot_list[[i]][[j]] + theme(legend.position="none")
        }
        
      }
      
      for (i in 1:length(names(names_list))) {
        
        for (j in 2:length(names(names_list[[i]]))){
          
          p_Level_plot_list[[i]][[j]] <- p_Level_plot_list[[i]][[j]] + theme(axis.title.y=element_blank())
        }
        
      }
      
      nonNA_position <- which(names(p_Level_plot_list[[features]])!="NA")
      p_Level_plot_list[[features]] <- p_Level_plot_list[[features]][nonNA_position]
      
      egg::ggarrange(plots = p_Level_plot_list[[features]], nrow=1)
      
    }
    
    return(plot_LeveltoSamples_sub_noUnassigned(taxtable = TaxaTable_merge(), metadata = Metadata_stats(), features = input$metadata_barplot, Level = input$select2))
  })
  
  
  output$download_barplot <- downloadHandler(
    
    filename = function(){
      paste0("barplot_", input$metadata_barplot, "_", input$select2, ".png")
    },
    content = function(file){
      
      
      if (input$metadata_barplot=="SampleID"){
        # if (input$without_unassigned == F){
        ggsave(file, plot = barplot_download(), width = 80, height = 40, units = "cm")
        # }else{
        #   ggsave(file, plot = barplot_noUnassigned_download(), width = 80, height = 40, units = "cm")
        # }
        
      }else{
        
        # if(input$without_unassigned == F){
        
        
        ggsave(file, plot = barplot_sub_download(), width = 80, height = 40, units = "cm")
        
        
        # }else{
        #     ggsave(file, plot = barplot_noUnassigned_sub_download(), width = 80, height = 40, units = "cm")
        # }
        
        
      }
    }
  )
  
  # Alpha diversity Table---------------------------------------------------------------------------------------
  
  output$contents2<-renderDataTable({
    
    return(alpha_diversity_table())
    
  })
  
  
  # Alpha diversity table for download--------------------------------------------------------------------------
  
  output$downloadAlpha<-downloadHandler(
    
    filename = "AlphaDiversityTable_result.csv",
    content = function(file) {
      
      write.csv(alpha_diversity_table(), file, row.names = FALSE)
      
    }
  )
  
  # Alpha diversity boxplot--------------------------------------------------------------------------------------
  
  output$alpha_boxplot <- renderPlot({
    
    if(input$select_stat=="ANOVA"){
      return(alpha_anova_boxplot())
    }else{
      return(alpha_KW_boxplot())
    }
    
    
  })
  
  # Alpha diversity statistics----------------------------------------------------------------------------------
  
  output$post_test_type <- renderText({
    
    if(input$select_stat=="ANOVA"){
      return("Tukey test")
    }else{
      return("Dunn test")
    } 
    
  })
  
  output$post_test <- renderTable({
    
    if(input$select_stat=="ANOVA"){
      return(alpha_anova_tukey())
    }else{
      return(alpha_KW_Dunn())
    }
    
  })
  
  # Alpha diversity statistics for download----------------------------------------------------------------------------------
  
  output$Alpha_posttest_DL <- downloadHandler(
    
    filename = function(){
      
      if(input$select_stat=="ANOVA"){
        paste0("Alpha_diversity_Tukey_", input$metadata_alpha, "_",input$select_diversity,".csv")
      }else{
        paste0("Alpha_diversity_Dunn_", input$metadata_alpha, "_",input$select_diversity,".csv")
      }
      
    },
    
    content = function(file){
      
      if(input$select_stat=="ANOVA"){
        write.csv(alpha_anova_tukey(), file, row.names = F)
      }else{
        write.csv(alpha_KW_Dunn(), file, row.names = F)
      }
      
    }
  )
  # Alpha diversity boxplot for download---------------------------------------------------------------------------------
  
  output$downloadAlphaBoxPlot <- downloadHandler(
    
    filename = function(){
      paste0("Alpha_diversity_boxplot_", input$metadata_alpha, "_",input$select_diversity,".png")
    },
    content = function(file){
      
      ggsave(file, plot = alpha_anova_boxplot())
      
    }
  )
  # Beta diversity distance matrix heatmap-------------------------------------------------------------------------------------------------
  
  output$beta_dsmx_hm <- renderPlotly({
    
    return(Beta_dsmx_heatmap())
  })
  
  # Beta diversity distance matrix heatmap for download------------------------------------------------------------------------------------
  
  output$download_beta_dsmx_hm <- downloadHandler(
    
    filename = "Beta_dsmx_data.csv",
    content = function(file){
      write.csv(BetaTable_bray(), file, row.names = T)
    }
  )
  
  # Produce PCoA and NMDS
  output$betaplot<-renderPlot({
    
    return(BetaPlot())
    
  })
  
  
  # PCoA and NMDS plot for download---------------------------------------------------------------------------------
  
  output$downloadBetaPlot<-downloadHandler(
    
    filename = function(){
      if(input$sep == "PCA"){
        paste0("BetaPlot_PCA_", input$metadata_beta,".png")
      }
      else if(input$sep == "PCoA") {
        paste0("BetaPlot_PCoA_", input$metadata_beta,".png")
      }
      else {
        paste0("BetaPlot_NMDS_", input$metadata_beta, ".png")
      }
    },
    content = function(file){
      
      return(ggsave(file, plot = BetaPlot()))
      
    }
  )
  
  # Permanova table------------------------------------------------------------------------------------------------
  
  Permanova_table <- reactive({
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
    
    adonis_result_table_list <- lapply(colnames(nonNA_metadata)[-1], function(x){
      
      adonis(bray_df~nonNA_metadata[,x], nonNA_metadata, permutations = 999)[["aov.tab"]]
      
    })
    
    names(adonis_result_table_list) <- colnames(nonNA_metadata)[-1]
    
    for (i in colnames(nonNA_metadata)[-1]) {
      
      rownames(adonis_result_table_list[[i]])[1] <- i
    }
    
    adonis_result_table_list_show <- adonis_result_table_list[[input$metadata_beta]][1, c(5,6)]
    colnames(adonis_result_table_list_show) <- c("R^2", "P value")
    
    
    
    adonis_result_table_list_show[, "R^2"] <- adonis_result_table_list_show[, "R^2"] %>% round(4) %>% as.character()
    adonis_result_table_list_show[, "P value"] <- adonis_result_table_list_show[, "P value"] %>% round(4) %>% as.character()
    
    return(adonis_result_table_list_show)
  })
  
  output$Permanova_title <- renderText({
    
    if(is.data.frame(Permanova_table())){
      
      return("PERMANOVA")
    }else{
      return(NULL)
    }
    
  })
  
  output$permanova_table <- renderTable({
    
    return(Permanova_table())
    
  }
  #, caption = "PERMANOVA", caption.placement = getOption("xtable.caption.placement", "top")
  )
  
  # download permanova table
  output$download_permanova <- downloadHandler(
    
    filename = "PerMANOVA_table.csv",
    content = function(file) {
      
      write.csv(Permanova_table(), file, row.names = FALSE)
      
    }
  )
  
  # When length(group_names)<=2, hide the download button of pair table
  observe({
    
    req(input$sample_data, input$taxonomic_table)
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    group_names <- unique(nonNA_metadata[, input$metadata_beta])
    
    if(length(group_names)<=2){
      shinyjs::hide("download_permanova_pair")
      shinyjs::hide("download_ANOSIM_pair")
      shinyjs::hide("download_MRPP_pair")
    }
    
  })
  
  
  Permanova_pair_table <- reactive({
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
    
    group_names <- unique(nonNA_metadata[,input$metadata_beta])
    
    if(length(group_names)>2){
      
      sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
        
        sample_cb_meta <- nonNA_metadata[,c(1,x)]
        
      })
      names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
      
      
      taxatable_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
        
        feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
        
        sample_f1_f2 <- sample_list[[input$metadata_beta]] %>% filter(sample_list[[input$metadata_beta]][,input$metadata_beta] %in% feature1_feature2)
        return(nonNA_taxtable[, sample_f1_f2[,"SampleID"]])
        
      })
      
      adonis_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
        
        bray_df_pair <- vegdist(t(taxatable_f1_f2_list[[k]]), method = "bray") %>% as.matrix() %>% as.data.frame()
        metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(bray_df_pair))
        adonis_result_pair <- adonis(bray_df_pair~metadata_pair[,input$metadata_beta], 
                                     metadata_pair, permutations = 999) 
        
      })
      
      pair_names <- c()
      df_pair<- t(combn(group_names, 2))
      for (i in 1:ncol(combn(group_names, 2))) {
        
        pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
        
      }
      names(adonis_result_pair_list) <- pair_names
      
      adonis_result_pair_list_table <- as.character(adonis_result_pair_list[[1]][["aov.tab"]][1,])
      for (i in 1:(ncol(combn(group_names, 2))-1)) {
        adonis_result_pair_list_table <- rbind(adonis_result_pair_list_table, 
                                               as.character(adonis_result_pair_list[[i+1]][["aov.tab"]][1,]))
      }
      
      colnames(adonis_result_pair_list_table) <- c("Df", "Sums Of Sqs", "Mean Sqs", "F.Model", "R^2", "P value")
      adonis_result_pair_list_table <- cbind(Comparisons=pair_names, adonis_result_pair_list_table)
      
      for (i in 3:7) {
        
        adonis_result_pair_list_table[, i] <- adonis_result_pair_list_table[, i] %>% as.numeric() %>% round(4) %>% as.character()
        
      }
      
      return(adonis_result_pair_list_table[, c(1, 6, 7)])
      
    }else{
      return(NULL)
    }
    
  })
  
  
  
  output$Permanova_pair_title <- renderText({
    
    if(is.null(Permanova_pair_table())){
      
      return(NULL)
    }else{
      return("PERMANOVA pair")
    }
    
    
  })
  
  output$permanova_pair_table <- renderTable({
    
    return(Permanova_pair_table())
  })
  
  
  # download permanova pair table
  output$download_permanova_pair <- downloadHandler(
    
    filename = "PerMANOVA_pair_table.csv",
    content = function(file) {
      
      write.csv(Permanova_pair_table(), file, row.names = FALSE)
      
    }
  )
  
  # ANOSIM table------------------------------------------------------------------------------------------------
  
  ANOSIM_table <- reactive({
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
    anosim_result <- anosim(bray_df, nonNA_metadata[, input$metadata_beta], permutations = 999)
    anosim_result_table <- data.frame(Significance=anosim_result$signif, statistic.R=anosim_result$statistic)
    anosim_result_table <- anosim_result_table[,c(2,1)]
    colnames(anosim_result_table) <- c("R", "P value")
    
    anosim_result_table$R <- anosim_result_table$R %>% as.numeric() %>% round(4) %>% as.character()
    anosim_result_table$"P value" <- anosim_result_table$"P value" %>% as.numeric() %>% round(4) %>% as.character()
    
    return(anosim_result_table)
    
  })
  
  output$ANOSIM_title <- renderText({
    
    if(is.data.frame(ANOSIM_table())){
      
      return("ANOSIM")
    }else{
      return(NULL)
    }
    
    
  })
  
  output$ANOSIM_table <- renderTable({
    
    return(ANOSIM_table())
  })
  
  # download ANOSIM table
  output$download_ANOSIM <- downloadHandler(
    
    filename = "ANOSIM_table.csv",
    content = function(file) {
      
      write.csv(ANOSIM_table(), file, row.names = FALSE)
      
    }
  )
  
  ANOSIM_pair_table <- reactive({
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
    
    group_names <- unique(nonNA_metadata[,input$metadata_beta])
    
    if(length(group_names)>2){
      
      sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
        
        sample_cb_meta <- nonNA_metadata[,c(1,x)]
        
      })
      names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
      
      
      taxatable_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
        
        feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
        
        sample_f1_f2 <- sample_list[[input$metadata_beta]] %>% filter(sample_list[[input$metadata_beta]][,input$metadata_beta] %in% feature1_feature2)
        return(TaxaTable_merge()[, sample_f1_f2[,"SampleID"]])
        
      })
      
      anosim_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
        
        bray_df_pair <- vegdist(t(taxatable_f1_f2_list[[k]]), method = "bray") %>% as.matrix() %>% as.data.frame()
        metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(bray_df_pair))
        anosim_result_pair <- anosim(bray_df_pair, metadata_pair[,input$metadata_beta],
                                     metadata_pair, permutations = 999)
        
      })
      
      pair_names <- c()
      df_pair <- t(combn(group_names, 2))
      for (i in 1:ncol(combn(group_names, 2))) {
        
        pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
        
      }
      names(anosim_result_pair_list) <- pair_names
      
      anosim_result_pair_list_signif <- c()
      for (i in 1:ncol(combn(group_names, 2))) {
        
        anosim_result_pair_list_signif[i] <- anosim_result_pair_list[[i]][["signif"]]
      }
      
      anosim_result_pair_list_R <- c()
      for (i in 1:ncol(combn(group_names, 2))) {
        
        anosim_result_pair_list_R[i] <- anosim_result_pair_list[[i]][["statistic"]]
      }
      
      anosim_result_pair_list_table <- data.frame(Comparisons=pair_names,
                                                  R=anosim_result_pair_list_R,
                                                  "P_value"=anosim_result_pair_list_signif)
      
      anosim_result_pair_list_table$R <- anosim_result_pair_list_table$R %>% as.numeric() %>% round(4) %>% as.character()
      anosim_result_pair_list_table$P_value <- anosim_result_pair_list_table$P_value %>% as.numeric() %>% round(4) %>% as.character()
      colnames(anosim_result_pair_list_table)[3] <- "P value"
      
      return(anosim_result_pair_list_table)
      
    }else{
      return(NULL)
    }
    
  })
  
  output$ANOSIM_pair_title <- renderText({
    
    if(is.null(ANOSIM_pair_table())){
      
      return(NULL)
    }else{
      return("ANOSIM pair")
    }
  })
  
  output$ANOSIM_pair_table <- renderTable({
    
    return(ANOSIM_pair_table())
  })
  
  # download ANOSIM pair table
  output$download_ANOSIM_pair <- downloadHandler(
    
    filename = "ANOSIM_pair_table.csv",
    content = function(file) {
      
      write.csv(ANOSIM_pair_table(), file, row.names = FALSE)
      
    }
  )
  
  # MRPP table------------------------------------------------------------------------------------------------
  
  MRPP_table <- reactive({
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
    mrpp_result <- mrpp(bray_df, nonNA_metadata[, input$metadata_beta], permutations = 999)
    mrpp_result_table <- data.frame(Group = 'all', 
                                    Distance = 'Bray-Curtis', 
                                    A = mrpp_result$A,            
                                    Observe.delta = mrpp_result$delta,            
                                    Expect.delta = mrpp_result$E.delta,            
                                    P.value = mrpp_result$Pvalue)
    mrpp_result_table_show <- mrpp_result_table[,3:6]
    colnames(mrpp_result_table_show)[4] <- "P value"
    
    for (i in 1:4) {
      
      mrpp_result_table_show[, i] <- mrpp_result_table_show[, i] %>% round(4) %>% as.character()
      
    }
    
    return(mrpp_result_table_show)
  })
  
  output$MRPP_title <- renderText({
    
    if(is.data.frame(MRPP_table())){
      
      return("MRPP")
    }else{
      return(NULL)
    }
  })
  
  output$MRPP_table <- renderTable({
    
    return(MRPP_table())
    
  })
  
  # download MRPP table
  output$download_MRPP <- downloadHandler(
    
    filename = "MRPP_pair_table.csv",
    content = function(file) {
      
      write.csv(MRPP_table(), file, row.names = FALSE)
      
    }
  )
  
  MRPP_pair_table <- reactive({
    
    nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    
    nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
    nonNA_metadata <- Metadata_stats()[nonNA_position, ]
    
    bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
    group_names <- unique(nonNA_metadata[,input$metadata_beta])
    
    if(length(group_names)>2){
      
      sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
        
        sample_cb_meta <- nonNA_metadata[,c(1,x)]
        
      })
      names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
      
      
      taxatable_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
        
        feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
        
        sample_f1_f2 <- sample_list[[input$metadata_beta]] %>% filter(sample_list[[input$metadata_beta]][,input$metadata_beta] %in% feature1_feature2)
        return(TaxaTable_merge()[, sample_f1_f2[,"SampleID"]])
        
      })
      
      MRPP_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
        
        bray_df_pair <- vegdist(t(taxatable_f1_f2_list[[k]]), method = "bray") %>% as.matrix() %>% as.data.frame()
        metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(bray_df_pair))
        MRPP_result_pair <- mrpp(bray_df_pair, metadata_pair[,input$metadata_beta], permutations = 999)
        
      })
      
      pair_names <- c()
      df_pair <- t(combn(group_names, 2))
      for (i in 1:ncol(combn(group_names, 2))) {
        
        pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
        
      }
      names(MRPP_result_pair_list) <- pair_names
      
      MRPP_result_pair_list_Pvalue <- c()
      for (i in 1:ncol(combn(group_names, 2))) {
        MRPP_result_pair_list_Pvalue[i] <- MRPP_result_pair_list[[i]][["Pvalue"]]
      }
      
      MRPP_result_pair_list_A <- c()
      for (i in 1:ncol(combn(group_names, 2))) {
        MRPP_result_pair_list_A[i] <- MRPP_result_pair_list[[i]][["A"]]
      }
      
      MRPP_result_pair_list_delta <- c()
      for (i in 1:ncol(combn(group_names, 2))) {
        MRPP_result_pair_list_delta[i] <- MRPP_result_pair_list[[i]][["delta"]]
      }
      
      MRPP_result_pair_list_Edelta <- c()
      for (i in 1:ncol(combn(group_names, 2))) {
        MRPP_result_pair_list_Edelta[i] <- MRPP_result_pair_list[[i]][["E.delta"]]
      }
      
      MRPP_result_pair_list_table <- data.frame(Comparisons=pair_names,
                                                A=MRPP_result_pair_list_A,
                                                delta=MRPP_result_pair_list_delta,
                                                E.delta=MRPP_result_pair_list_Edelta,
                                                P_value=MRPP_result_pair_list_Pvalue)
      
      for (i in 2:5) {
        
        MRPP_result_pair_list_table[, i] <- MRPP_result_pair_list_table[, i] %>% round(4) %>% as.character()
        
      }
      
      colnames(MRPP_result_pair_list_table)[5] <- "P value"
      
      return(MRPP_result_pair_list_table)
    }else{
      return(NULL)
    }
  })
  
  output$MRPP_pair_title <- renderText({
    
    if(is.null(MRPP_pair_table())){
      
      return(NULL)
    }else{
      return("MRPP pair")
    }
  })
  
  output$MRPP_pair_table <- renderTable({
    
    return(MRPP_pair_table())
    
    output$download_MRPP_pair <- downloadHandler(
      
      filename = "MRPP_pair_table.csv",
      content = function(file) {
        
        write.csv(MRPP_pair_table(), file, row.names = FALSE)
        
      }
    )
    
  })
  
  # Krona ------------------------------------------------------------------------------------------------------
  
  # observeEvent(input$run_krona, {
  # 
  #   Sys.setenv(TERM="xterm")
  # 
  #   # showModal(modalDialog(title = "Running Krona ...", "Waiting for a moment", footer = NULL))
  #   show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
  # 
  #   Krona_data<-function(taxtable, sampleID){
  #     taxtable_0<-taxtable[, sampleID]
  #     taxtable_1<-data.frame(taxtable_0, Species=row.names(taxtable))
  #     colnames(taxtable_1)[1]<-sampleID
  #     taxtable_2<-as_tibble(taxtable_1)
  #     taxtable_3<-separate(data = taxtable_2, col = "Species", into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = ";")
  #     write.table(taxtable_3,file = paste("taxtable_forKrona", sampleID, sep = "_"), quote = F, sep = "\t", row.names = F, col.names = F)
  #   }
  # 
  #   for (i in colnames(TaxaTable())) {
  #     Krona_data(TaxaTable(), i)
  #   }
  # 
  #   filename_forKrona <- list.files(pattern = "taxtable_forKrona") %>% paste(collapse = " ")
  #   system(paste('ktImportText -o /home/imuser/Krona_rawdata.html', filename_forKrona))
  #   file.remove(list.files(pattern = "taxtable_forKrona"))
  #   # filesstrings::file.move("/home/imuser/Krona_rawdata.html", "/home/imuser/www", overwrite = T)
  #   # filesstrings::file.move(files = c("/home/imuser/www/Krona_rawdata.html", "home/imuser/www/iframe.html"),
  #   #                         destinations = "/srv/shiny-server/www", overwrite = T)
  #   # file.remove(list.files(path = "/srv/shiny-server/www/", full.names = T))
  #   # file.copy(from = "/home/imuser/www/iframe.html", to = "/srv/shiny-server/www/", overwrite = T, copy.mode = T)
  #   # file.copy(from = "/home/imuser/www/Krona_rawdata.html", to = "/srv/shiny-server/www/", overwrite = T, copy.mode = T)
  # 
  #   # system("rm /srv/shiny-server/www/*")
  #   # system("cp /home/imuser/www/iframe.html /srv/shiny-server/www/")
  #   system("cp /home/imuser/www/Krona_rawdata.html /var/www/html/")
  #   # system("cp /home/imuser/www/Krona_rawdata.html /srv/shiny-server/www/")
  # 
  #   output$krona_output <- renderUI(includeHTML(path = "/var/www/html/iframe_krona.html"))
  # 
  # 
  #   # removeModal()
  #   remove_modal_spinner()
  # })
  
  output$krona_output <- renderUI({
    
    Krona_data<-function(taxtable, sampleID){
      taxtable_0<-taxtable[, sampleID]
      taxtable_1<-data.frame(taxtable_0, Species=row.names(taxtable))
      colnames(taxtable_1)[1]<-sampleID
      taxtable_2<-as_tibble(taxtable_1)
      taxtable_3<-separate(data = taxtable_2, col = "Species", into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = ";")
      setwd("/home/imuser/Krona_files/")
      write.table(taxtable_3,file = paste("taxtable_forKrona", sampleID, sep = "_"), quote = F, sep = "\t", row.names = F, col.names = F)
    }
    
    for (i in colnames(TaxaTable())) {
      Krona_data(TaxaTable(), i)
    }
    
    filename_forKrona <- list.files(pattern = "taxtable_forKrona") %>% paste(collapse = " ")
    system(paste('ktImportText -o /home/imuser/Krona_files/Krona_rawdata.html', filename_forKrona))
    file.remove(list.files(pattern = "taxtable_forKrona"))
    
    file.copy(from = "/home/imuser/Krona_files/Krona_rawdata.html",
              to = "/var/www/html/",
              overwrite = T)
    
    write_file(x = paste0('<iframe src="http://', my_qiime_ip, my_qiime_port(),'/Krona_rawdata.html" height=800px width=100% align = middle frameborder = 0></iframe>'),
               path = "/home/imuser/Krona_files/iframe_krona.html",
               append = F)
    file.copy(from = "/home/imuser/Krona_files/iframe_krona.html",
              to = "/var/www/html/",
              overwrite = T)
    
    
    includeHTML(path = paste0("http://",
                              my_qiime_ip, my_qiime_port(),
                              "/iframe_krona.html"))
    
  })
  
  options(shiny.sanitize.errors = T)
  options(shiny.trace = F)
  
  # Krona download ---------------------------------------------------------------------------------------------
  output$download_krona <- downloadHandler(
    
    
    filename = "Krona_results.zip",
    # filename = "selected_sampleID.tsv",
    
    content = function(file){
      system("zip -r /home/imuser/Krona.zip /home/imuser/Krona_files/*")
      file.copy("/home/imuser/Krona.zip", file)
      # file.copy("/srv/shiny-server/www/Krona_rawdata.html", file, overwrite = T)
      # file.copy("/home/imuser/selected_sampleID.tsv", file)
    }
    
  )
  
  
  # Phylogenetic diversity analysis ---------------------------------------------------------------------------
  
  output$example_rep_seqs <- downloadHandler(
    
    filename = "example_rep_seqs.qza",
    
    content = function(file){
      file.copy(from = "/home/imuser/example_files/rep-seqs.qza", to = file)
    },
    contentType = "application/qza"
  )
  
  output$example_feature_table <- downloadHandler(
    
    filename = "feature_table_forPhylo_example.qza",
    
    content = function(file){
      file.copy(from = "/home/imuser/example_files/table.qza", to = file)
    },
    contentType = "application/qza"
  )
  
  observeEvent(req(input$phylogenetic_tree, input$table_dada2_upload, input$rep_seq_dada2_upload), {
    
    # showModal(modalDialog(title = "Running making tree ...", "Waiting for a moment",footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    
    system(paste(qiime_cmd, "phylogeny align-to-tree-mafft-fasttree --i-sequences", input$rep_seq_dada2_upload$datapath, 
                 "--p-n-threads", input$threads_phylogenetic,
                 "--o-alignment /home/imuser/qiime_output/aligned-rep-seqs-dada2.qza --o-masked-alignment /home/imuser/qiime_output/masked-aligned-rep-seqs-dada2.qza",
                 "--o-tree /home/imuser/qiime_output/unrooted-tree.qza --o-rooted-tree /home/imuser/qiime_output/rooted-tree.qza"))
    
    write.table(Metadata_stats(), file='/home/imuser/metadata.tsv', quote=FALSE, sep='\t', row.names = F)
    system("rm -r /home/imuser/qiime_output/core-metrics-results/")
    system(paste(qiime_cmd, "diversity core-metrics-phylogenetic --i-phylogeny /home/imuser/qiime_output/rooted-tree.qza",
                 "--i-table", input$table_dada2_upload$datapath,
                 "--p-sampling-depth ", input$sampling_depth,
                 "--m-metadata-file /home/imuser/metadata.tsv", 
                 "--output-dir /home/imuser/qiime_output/core-metrics-results"))
    
    faith_PD <- reactive({
      read_qza("/home/imuser/qiime_output/core-metrics-results/faith_pd_vector.qza")[["data"]]
    })
    
    output$contents4 <- renderDataTable({
      a <- faith_PD()
      b <- data.frame(SampleID=rownames(a),
                      FaithPD=a[,1])
      return(b)
    })
    
    faithPD_boxplot_anova <- reactive({
      
      as_faithPD_boxplot <- function(metadata, feature)
        
      {
        
        A_diversity <- faith_PD()
        
        A_diversity <- cbind(sample.id=rownames(A_diversity), A_diversity)
        colnames(A_diversity)[1] <- colnames(metadata)[1]
        A_diversity_metadata <- merge(A_diversity, metadata, by= colnames(metadata)[1])
        
        A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)], function(x){
          
          A_diversity_metadata_merge <- merge(A_diversity, 
                                              metadata[, c(colnames(A_diversity_metadata)[1], x)], 
                                              by= colnames(metadata)[1])
        })
        
        names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)]
        
        for (i in 1:length(A_diversity_metadata_list)) {
          
          colnames(A_diversity_metadata_list[[i]])[3] <- "feature_name"
        }
        
        A_diversity_metadata_list[[feature]] <- filter(A_diversity_metadata_list[[feature]], feature_name != "NA")
        
        anova_result <- aov(faith_pd ~ feature_name, A_diversity_metadata_list[[feature]])
        anova_summary <- summary(anova_result)
        # tukey_result <- agricolae::HSD.test(anova_result, "feature_name", group = T)
        # group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
        options(scipen=999)
        
        stat_box_data <- function(y, upper_limit = max(A_diversity_metadata_list[[feature]]$faith_pd) * 1.15) {
          return( 
            data.frame(
              y = 0.95 * upper_limit,
              label = paste('n =', length(y))
            )
          )
        }
        
        ggplot(A_diversity_metadata_list[[feature]], aes(x = feature_name, y = faith_pd)) +
          # geom_text(data = data.frame(),
          #           aes(x = rownames(group_data), y = max(A_diversity_metadata_list[[feature]]$faith_pd) + 1, label = group_data$groups),
          #           col = 'black',
          #           size = 5) +
          geom_boxplot() +
          ggtitle("Phylogenetic Alpha diversity") +
          xlab(input$metadata_phylo_beta) +
          ylab("Faith_PD")+
          labs(caption=paste0("p value of ANOVA = ", round(anova_summary[[1]][[5]][1], 4))) + theme(text = element_text(size = 15)) + stat_summary(fun.data = stat_box_data,
                                                                                                                                                   geom = "text", 
                                                                                                                                                   hjust = 0.5,
                                                                                                                                                   vjust = 0.9)
        
      }
      
      
      
      return(as_faithPD_boxplot(Metadata_stats(), input$metadata_phylo_alpha
      )
      )
      
      
    })
    
    faithPD_boxplot_KWtest <- reactive({
      
      as_faithPD_boxplot <- function(metadata, feature)
        
      {
        
        A_diversity <- faith_PD()
        
        A_diversity <- cbind(sample.id=rownames(A_diversity), A_diversity)
        colnames(A_diversity)[1] <- colnames(metadata)[1]
        A_diversity_metadata <- merge(A_diversity, metadata, by= colnames(metadata)[1])
        
        A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)], function(x){
          
          A_diversity_metadata_merge <- merge(A_diversity, 
                                              metadata[, c(colnames(A_diversity_metadata)[1], x)], 
                                              by= colnames(metadata)[1])
        })
        
        names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)]
        
        for (i in 1:length(A_diversity_metadata_list)) {
          
          colnames(A_diversity_metadata_list[[i]])[3] <- "feature_name"
        }
        
        A_diversity_metadata_list[[feature]] <- filter(A_diversity_metadata_list[[feature]], feature_name!="NA")
        
        KW_result <- kruskal.test(faith_pd ~ feature_name, A_diversity_metadata_list[[feature]])
        # tukey_result <- agricolae::HSD.test(anova_result, "feature_name", group = T)
        # group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
        options(scipen=999)
        
        stat_box_data <- function(y, upper_limit = max(A_diversity_metadata_list[[feature]]$faith_pd) * 1.15) {
          return( 
            data.frame(
              y = 0.95 * upper_limit,
              label = paste('n =', length(y))
            )
          )
        }
        
        ggplot(A_diversity_metadata_list[[feature]], aes(x = feature_name, y = faith_pd)) +
          # geom_text(data = data.frame(),
          #           aes(x = rownames(group_data), y = max(A_diversity_metadata_list[[feature]]$faith_pd) + 1, label = group_data$groups),
          #           col = 'black',
          #           size = 5) +
          geom_boxplot() +
          ggtitle("Phylogenetic Alpha diversity") +
          xlab(input$metadata_phylo_beta) +
          ylab("Faith_PD")+
          labs(caption=paste0("p value of KW-test = ", round(KW_result$p.value, 4))) + theme(text = element_text(size = 15)) + stat_summary(fun.data = stat_box_data,
                                                                                                                                            geom = "text", 
                                                                                                                                            hjust = 0.5,
                                                                                                                                            vjust = 0.9) 
        
      }
      
      
      
      return(as_faithPD_boxplot(Metadata_stats(), input$metadata_phylo_alpha
      )
      )
      
      
    })
    
    
    
    output$faith_PD_boxplot <- renderPlot({
      
      if(input$select_stat_phylo=="ANOVA"){
        return(faithPD_boxplot_anova())
      }else{
        return(faithPD_boxplot_KWtest())
      }
      
    })
    
    output$download_faithPD_boxplot <- downloadHandler(
      
      filename = function(){
        paste0("faithPD_boxplot_", input$metadata_phylo_alpha, ".png")
      },
      content = function(file){
        
        ggsave(file, plot = faithPD_boxplot())
        
      })
    
    
    faith_PD_post_test_tukey <- reactive({
      
      A_diversity <- faith_PD()
      
      A_diversity <- cbind(sample.id=rownames(A_diversity), A_diversity)
      colnames(A_diversity)[1] <- colnames(Metadata_stats())[1]
      A_diversity_metadata <- merge(A_diversity, Metadata_stats(), by= colnames(Metadata_stats())[1])
      
      A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)], function(x){
        
        A_diversity_metadata_merge <- merge(A_diversity, 
                                            Metadata_stats()[, c(colnames(A_diversity_metadata)[1], x)], 
                                            by= colnames(Metadata_stats())[1])
      })
      
      names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)]
      
      for (i in 1:length(A_diversity_metadata_list)) {
        
        colnames(A_diversity_metadata_list[[i]])[3] <- "feature_name"
      }
      
      A_diversity_metadata_list[[input$metadata_phylo_alpha]] <- filter(A_diversity_metadata_list[[input$metadata_phylo_alpha]], feature_name != "NA")
      
      anova_result <- aov(faith_pd ~ feature_name, A_diversity_metadata_list[[input$metadata_phylo_alpha]])
      tukey_result <- TukeyHSD(anova_result, "feature_name")[[1]]
      tukey_result_table <- data.frame(comparisons = rownames(tukey_result), tukey_result)
      tukey_result_table[,"comparisons"] <- stringr::str_replace_all(tukey_result_table[,"comparisons"], "-", " / ")
      tukey_result_table_separate <- separate(data = tukey_result_table, col = "comparisons", into = c("Group_A", "Group_B"), sep = " / ")
      order_by_count <- as.data.frame(table(tukey_result_table_separate[,1]))[,1]
      tukey_result_table_separate_order <- tukey_result_table_separate[rev(order(factor(tukey_result_table_separate$Group_A, levels = order_by_count))),]
      colnames(tukey_result_table_separate_order)[c(1,2,3,6)] <- c("Group A", "Group B", "Diff", "P value")
      return(tukey_result_table_separate_order[,c(1,2,3,6)])
    })
    
    faith_PD_post_test_Dunn <- reactive({
      
      A_diversity <- faith_PD()
      
      A_diversity <- cbind(sample.id=rownames(A_diversity), A_diversity)
      colnames(A_diversity)[1] <- colnames(Metadata_stats())[1]
      A_diversity_metadata <- merge(A_diversity, Metadata_stats(), by= colnames(Metadata_stats())[1])
      
      A_diversity_metadata_list <- lapply(colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)], function(x){
        
        A_diversity_metadata_merge <- merge(A_diversity, 
                                            Metadata_stats()[, c(colnames(A_diversity_metadata)[1], x)], 
                                            by= colnames(Metadata_stats())[1])
      })
      
      names(A_diversity_metadata_list) <- colnames(A_diversity_metadata)[3:ncol(A_diversity_metadata)]
      
      for (i in 1:length(A_diversity_metadata_list)) {
        
        colnames(A_diversity_metadata_list[[i]])[3] <- "feature_name"
      }
      
      A_diversity_metadata_list[[input$metadata_phylo_alpha]] <- filter(A_diversity_metadata_list[[input$metadata_phylo_alpha]], feature_name!="NA")
      
      Dunn_result <- dunn.test::dunn.test(A_diversity_metadata_list[[input$metadata_phylo_alpha]]$faith_pd, A_diversity_metadata_list[[input$metadata_phylo_alpha]]$feature_name)
      Dunn_table <- data.frame(comparisons=Dunn_result$comparisons,
                               Z=Dunn_result$Z,
                               Pvalue= Dunn_result$P)
      Dunn_table[,"comparisons"] <- stringr::str_replace_all(Dunn_table[,"comparisons"], "-", " / ")
      Dunn_table_separate <- separate(data = Dunn_table, col = "comparisons", into = c("Group_A", "Group_B"), sep = " / ")
      order_by_count <- as.data.frame(table(Dunn_table_separate[,1]))[,1]
      Dunn_table_separate_order <- Dunn_table_separate[order(factor(Dunn_table_separate$Group_A, levels = order_by_count)),]
      colnames(Dunn_table_separate_order) <- c("Group A", "Group B", "Z", "P value")
      return(Dunn_table_separate_order)
    })
    
    output$post_test_type_phylo <- renderText({
      
      if(input$select_stat_phylo=="ANOVA"){
        return("Tukey test")
      }else{
        return("Dunn test")
      }
    })
    
    output$post_test_phylo <- renderTable({
      
      if(input$select_stat_phylo=="ANOVA"){
        return(faith_PD_post_test_tukey())
      }else{
        return(faith_PD_post_test_Dunn())
      }
    })
    
    
    output$unif_dm_hm <- renderPlotly({
      
      if(input$UnW_or_W=="Unweighted"){
        unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix()
        plot_ly(x=colnames(unW_unifrac_dm),
                y=rownames(unW_unifrac_dm),
                z=unW_unifrac_dm,
                # colors = palette(50),
                type = "heatmap") %>% layout(title="Unweighted unifrac distance matrix heatmap", xaxis=list(tickangle=45))
      }else{
        W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix()
        plot_ly(x=colnames(W_unifrac_dm),
                y=rownames(W_unifrac_dm),
                z=W_unifrac_dm,
                # colors = palette(50),
                type = "heatmap") %>% layout(title="Weighted unifrac distance matrix heatmap", xaxis=list(tickangle=45)) 
      }
    })
    
    
    output$download_unif_dm<-downloadHandler(
      
      filename = "unifrac_distance_matrix.csv",
      content = function(file) {
        
        write.csv(as.matrix(read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]]), file, row.names = T)
        
      }
    )
    
    unW_unif_pcoa_plot <- reactive({
      
      #PCoA
      unW_unifrac_dm_pcoa_qiime <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_pcoa_results.qza")[["data"]]
      unW_unifrac_dm_pcoa_qiime_forplot <- unW_unifrac_dm_pcoa_qiime$Vectors[,1:3]
      unW_unifrac_dm_pcoa_qiime_forplot <-merge(Metadata_stats(), unW_unifrac_dm_pcoa_qiime_forplot, by="SampleID") %>% as_tibble()
      
      unW_unifrac_dm_pcoa_qiime_forplot_table_list <- lapply(colnames(Metadata_stats()), function(i){
        
        unW_unifrac_dm_pcoa_qiime_forplot_table <- unW_unifrac_dm_pcoa_qiime_forplot[, c("SampleID", i, "PC1", "PC2")]
      })
      
      names(unW_unifrac_dm_pcoa_qiime_forplot_table_list) <- colnames(Metadata_stats())
      
      # Make all feature name to "feature"
      for (i in 1:length(unW_unifrac_dm_pcoa_qiime_forplot_table_list)) {
        names(unW_unifrac_dm_pcoa_qiime_forplot_table_list[[i]])[2] <- "feature"
      }
      
      
      unW_unifrac_dm_pcoa_qiime_plot_list <- lapply(1:length(unW_unifrac_dm_pcoa_qiime_forplot_table_list), function(i){
        
        ggplot(data = unW_unifrac_dm_pcoa_qiime_forplot_table_list[[i]], 
               aes(x=PC1, y=PC2, label=SampleID, color=feature))+
          geom_point(size=1.5)+
          ggrepel::geom_text_repel(show.legend = F)+
          xlab(paste0("PC1 (", round(unW_unifrac_dm_pcoa_qiime$ProportionExplained[1],2)*100, "%)"))+
          ylab(paste0("PC2 (", round(unW_unifrac_dm_pcoa_qiime$ProportionExplained[2],2)*100, "%)"))+
          geom_vline(xintercept =0, linetype="dotted")+
          geom_hline(yintercept = 0, linetype="dotted")+
          theme_bw()+
          ggtitle("Unweighted unifrac PCoA plot")+
          scale_colour_discrete(names(unW_unifrac_dm_pcoa_qiime_forplot_table_list)[[i]])
      })
      
      names(unW_unifrac_dm_pcoa_qiime_plot_list) <- colnames(Metadata_stats())
      return(unW_unifrac_dm_pcoa_qiime_plot_list[[input$metadata_phylo_beta]])
      
    })
    
    # unW_unif_pca_plot <- reactive({
    # 
    #   #PCA
    #   nonNA_position <- which(Metadata_stats()[, input$metadata_phylo_beta]!="NA")
    #   taxatable_beta <- TaxaTable_merge()[, nonNA_position]
    #   metadata_beta <- Metadata_stats()[nonNA_position,]
    #   colnames(metadata_beta)[1] <- "SampleID"
    #   
    #   # unW_unifrac_dm_qiime <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]]
    #   pca_Bray_df_data <- prcomp(unW_unifrac_dm_qiime)
    #   PCA_rowname <- as.matrix(unW_unifrac_dm_qiime) %>% rownames()
    #   
    #   metadata_beta_unW_unifrac <- filter(metadata_beta, SampleID %in% PCA_rowname)
    #   
    #   sample_original_names <- pca_Bray_df_data$scores %>% row.names()
    #   
    #   update_rownames <- function(feature_name, metadata, i){
    #     
    #     names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
    #     names_order <- which(names_TF==T)
    #     
    #     sample_names <- metadata[,1][names_order] %>% as.character()
    #     return(sample_names)
    #   }
    #   
    #   names_list <- lapply(colnames(metadata_beta_unW_unifrac)[1:ncol(metadata_beta_unW_unifrac)], function(i){
    #     
    #     sapply(1:length(unique(metadata_beta[, i])), function(j){
    #       
    #       update_rownames(i, metadata_beta,j)
    #       
    #       
    #     })
    #     
    #   })
    #   
    #   names(names_list) <- colnames(metadata_beta_unW_unifrac)[1:ncol(metadata_beta_unW_unifrac)]
    #   
    #   for (i in 1:length(names(names_list))) {
    #     names(names_list[[i]]) <- unique(metadata_beta_unW_unifrac[,names(names_list)[i]]) %>% as.character()
    #   }
    #   
    #   samples_unW <- as.matrix(unW_unifrac_dm_qiime) %>% rownames()
    #   colnames(metadata_beta_unW_unifrac)[1] <- "SampleID"
    #   # metadata_beta_unW <- filter(metadata_beta, SampleID %in% samples_unW)
    #   
    #   # metadata_beta_unW_arrange <- arrange(metadata_beta_unW_unifrac, "SampleID")
    #   # PCA_rowname_arrange <- metadata_beta_unW_arrange[, input$metadata_phylo_beta]
    #   
    #   
    #   
    #   
    #   library(ggplot2)
    #   pca_Bray_df_data_plot <- data.frame(sample=PCA_rowname, 
    #                                       PC1=pca_Bray_df_data$scores[,1],
    #                                       PC2=pca_Bray_df_data$scores[,2])
    #   pc_prop <- pca_Bray_df_data$sdev^2/sum(pca_Bray_df_data$sdev^2)
    #   
    #   library(ggrepel)
    #   pca_Bray_df_data_plot_gg <- ggplot(data = pca_Bray_df_data_plot, aes(x=PC1, y=PC2, label=sample_unW, color=sample))+
    #     geom_point(size = 1.5)+
    #     geom_text_repel(show.legend = FALSE)+
    #     xlab(paste("PC1 (", round(pc_prop[1], 2)*100, "%", ")", sep = ""))+
    #     ylab(paste("PC2 (", round(pc_prop[2], 2)*100, "%", ")", sep = ""))+
    #     geom_vline(xintercept = 0, linetype="dotted")+
    #     geom_hline(yintercept = 0, linetype="dotted")+
    #     theme_bw()+
    #     ggtitle("PCA plot")+
    #     scale_colour_discrete(input$metadata_phylo_beta) + theme(text = element_text(size = 15)) 
    #   
    #   
    # })
    
    unW_unif_nmds_plot <- reactive({
      
      #NMDS
      nonNA_position <- which(Metadata_stats()[, input$metadata_phylo_beta]!="NA")
      taxatable_beta <- TaxaTable_merge()[, nonNA_position]
      metadata_beta <- Metadata_stats()[nonNA_position,]
      colnames(metadata_beta)[1] <- "SampleID"
      
      unW_unifrac_dm_qiime <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]]
      sample_original_names <- as.matrix(unW_unifrac_dm_qiime) %>% rownames()
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- lapply(colnames(metadata_beta)[1:ncol(metadata_beta)], function(i){
        
        sapply(1:length(unique(metadata_beta[, i])), function(j){
          
          update_rownames(i, metadata_beta,j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(metadata_beta)[1:ncol(metadata_beta)]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(metadata_beta[,names(names_list)[i]]) %>% as.character()
      }
      
      metadata_beta_filter <- filter(metadata_beta, SampleID %in% sample_original_names)
      
      metadata_beta_filter_arrange <- arrange(metadata_beta_filter, SampleID)
      
      NMDS_rowname <- metadata_beta_filter_arrange[, input$metadata_phylo_beta]
      
      metaMDS_beta_df_data <- metaMDS(unW_unifrac_dm_qiime, distance = "bray")
      NMDS_beta_df_data <- data.frame(NMDS1=metaMDS_beta_df_data$points[,1], NMDS2=metaMDS_beta_df_data$points[,2])
      # NMDS_rowname <- NMDS_rowname_arrange
      NMDS_beta_df_data_plot <- data.frame(NMDS_beta_df_data, sample=NMDS_rowname)
      
      NMDS_beta_df_data_plot_gg <- ggplot(data = NMDS_beta_df_data_plot, aes(x=NMDS1, y=NMDS2, label=sample_original_names, color=sample))+
        geom_point(size=1.5)+
        geom_text_repel(show.legend = FALSE)+
        xlab("NMDS1")+
        ylab("NMDS2")+
        geom_vline(xintercept = 0, linetype = "dotted")+
        geom_hline(yintercept = 0, linetype = "dotted")+
        theme_bw()+
        labs(title="Unweighted unifrac NMDS plot", caption=paste("stress=", as.character(round(metaMDS_beta_df_data$stress, 4)), sep = ""))+
        #labs(caption = "A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation.")+
        scale_colour_discrete(input$metadata_phylo_beta) + theme(text = element_text(size = 15)) 
      
      return(NMDS_beta_df_data_plot_gg)
    })
    
    W_unif_pcoa_plot <- reactive({
      
      #PCoA
      W_unifrac_dm_pcoa_qiime <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_pcoa_results.qza")[["data"]]
      W_unifrac_dm_pcoa_qiime_forplot <- W_unifrac_dm_pcoa_qiime$Vectors[,1:3]
      W_unifrac_dm_pcoa_qiime_forplot <-merge(Metadata_stats(), W_unifrac_dm_pcoa_qiime_forplot, by="SampleID") %>% as_tibble()
      
      W_unifrac_dm_pcoa_qiime_forplot_table_list <- lapply(colnames(Metadata_stats()), function(i){
        
        W_unifrac_dm_pcoa_qiime_forplot_table <- W_unifrac_dm_pcoa_qiime_forplot[, c("SampleID", i, "PC1", "PC2")]
      })
      
      names(W_unifrac_dm_pcoa_qiime_forplot_table_list) <- colnames(Metadata_stats())
      
      # Make all feature name to "feature"
      for (i in 1:length(W_unifrac_dm_pcoa_qiime_forplot_table_list)) {
        names(W_unifrac_dm_pcoa_qiime_forplot_table_list[[i]])[2] <- "feature"
      }
      
      
      W_unifrac_dm_pcoa_qiime_plot_list <- lapply(1:length(W_unifrac_dm_pcoa_qiime_forplot_table_list), function(i){
        
        ggplot(data = W_unifrac_dm_pcoa_qiime_forplot_table_list[[i]], 
               aes(x=PC1, y=PC2, label=SampleID, color=feature))+
          geom_point(size=1.5)+
          ggrepel::geom_text_repel(show.legend = F)+
          xlab(paste0("PC1 (", round(W_unifrac_dm_pcoa_qiime$ProportionExplained[1],2)*100, "%)"))+
          ylab(paste0("PC2 (", round(W_unifrac_dm_pcoa_qiime$ProportionExplained[2],2)*100, "%)"))+
          geom_vline(xintercept =0, linetype="dotted")+
          geom_hline(yintercept = 0, linetype="dotted")+
          theme_bw()+
          ggtitle("Weighted unifrac PCoA plot")+
          scale_colour_discrete(names(W_unifrac_dm_pcoa_qiime_forplot_table_list)[[i]])
      })
      
      names(W_unifrac_dm_pcoa_qiime_plot_list) <- colnames(Metadata_stats())
      
      return(W_unifrac_dm_pcoa_qiime_plot_list[[input$metadata_phylo_beta]])
      
      
    })
    
    W_unif_nmds_plot <- reactive({
      
      nonNA_position <- which(Metadata_stats()[, input$metadata_phylo_beta]!="NA")
      taxatable_beta <- TaxaTable_merge()[, nonNA_position]
      metadata_beta <- Metadata_stats()[nonNA_position,]
      colnames(metadata_beta)[1] <- "SampleID"
      
      W_unifrac_dm_qiime <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]]
      sample_original_names <- as.matrix(W_unifrac_dm_qiime) %>% rownames()
      
      update_rownames <- function(feature_name, metadata, i){
        
        names_TF <- metadata[, feature_name]==unique(metadata[, feature_name])[i]
        names_order <- which(names_TF==T)
        
        sample_names <- metadata[,1][names_order] %>% as.character()
        return(sample_names)
      }
      
      names_list <- lapply(colnames(metadata_beta)[1:ncol(metadata_beta)], function(i){
        
        sapply(1:length(unique(metadata_beta[, i])), function(j){
          
          update_rownames(i, metadata_beta,j)
          
          
        })
        
      })
      
      names(names_list) <- colnames(metadata_beta)[1:ncol(metadata_beta)]
      
      for (i in 1:length(names(names_list))) {
        names(names_list[[i]]) <- unique(metadata_beta[,names(names_list)[i]]) %>% as.character()
      }
      
      metadata_beta_filter <- filter(metadata_beta, SampleID %in% sample_original_names)
      
      metadata_beta_filter_arrange <- arrange(metadata_beta_filter, SampleID)
      
      NMDS_rowname <- metadata_beta_filter_arrange[, input$metadata_phylo_beta]
      
      
      metaMDS_beta_df_data <- metaMDS(W_unifrac_dm_qiime, distance = "bray")
      NMDS_beta_df_data <- data.frame(NMDS1=metaMDS_beta_df_data$points[,1], NMDS2=metaMDS_beta_df_data$points[,2])
      # NMDS_rowname <- as.matrix(W_unifrac_dm_qiime) %>% rownames()
      NMDS_beta_df_data_plot <- data.frame(NMDS_beta_df_data, sample=NMDS_rowname)
      
      NMDS_beta_df_data_plot_gg <- ggplot(data = NMDS_beta_df_data_plot, aes(x=NMDS1, y=NMDS2, label=sample_original_names,color=sample))+
        geom_point(size=1.5)+
        geom_text_repel(show.legend = FALSE)+
        xlab("NMDS1")+
        ylab("NMDS2")+
        geom_vline(xintercept = 0, linetype = "dotted")+
        geom_hline(yintercept = 0, linetype = "dotted")+
        theme_bw()+
        labs(title="Weighted unifrac NMDS plot", caption=paste("stress=", as.character(round(metaMDS_beta_df_data$stress, 4)), sep = ""))+
        #labs(caption = "A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation.")+
        scale_colour_discrete(input$metadata_phylo_beta) + theme(text = element_text(size = 15))
      
      return(NMDS_beta_df_data_plot_gg)
    })
    
    
    output$unW_unif_ordination <- renderPlot({
      
      if(input$UnW_or_W_phylo=="Unweighted" & input$phylo_cluster == F & input$ordination_phylo == "PCoA"){
        
        return(unW_unif_pcoa_plot())
        
      }else if (input$UnW_or_W_phylo=="Unweighted" & input$phylo_cluster == F & input$ordination_phylo == "NMDS"){
        
        return(unW_unif_nmds_plot())
        
      }else if (input$UnW_or_W_phylo=="Unweighted" & input$phylo_cluster == T & input$ordination_phylo == "PCoA"){
        
        return(unW_unif_pcoa_plot() + stat_ellipse(type = "t"))
        
      }else if (input$UnW_or_W_phylo=="Unweighted" & input$phylo_cluster == T & input$ordination_phylo == "NMDS"){
        
        return(unW_unif_nmds_plot() + stat_ellipse(type = "t"))
        
      }else if(input$UnW_or_W_phylo=="Weighted" & input$phylo_cluster ==F & input$ordination_phylo == "PCoA"){
        
        return(W_unif_pcoa_plot())
        
      }else if(input$UnW_or_W_phylo=="Weighted" & input$phylo_cluster ==F & input$ordination_phylo == "NMDS"){
        
        return(W_unif_nmds_plot())
        
      }else if(input$UnW_or_W_phylo=="Weighted" & input$phylo_cluster ==T & input$ordination_phylo == "PCoA"){
        
        return(W_unif_pcoa_plot() + stat_ellipse(type = "t"))
        
      }else{
        
        return(W_unif_nmds_plot() + stat_ellipse(type = "t"))
        
      }
      
      
    })
    
    output$download_unif_plot <- downloadHandler(
      
      filename = function(){
        paste0("Unifrac_plot_", input$UnW_or_W_phylo, "_", input$metadata_phylo_beta, ".png")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          ggsave(file, plot = unW_unif_pcoa_plot())
        }else{
          ggsave(file, plot = W_unif_pcoa_plot())
        }
        
      }
      
    )
    
    
    Permanova_table_phylo_unW <- reactive({
      
      nonNA_position <- which(Metadata_stats()[, input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(unW_unifrac_dm))
      
      adonis_result_table_list <- lapply(colnames(nonNA_metadata)[-1], function(x){
        
        adonis(unW_unifrac_dm ~ nonNA_metadata[,x], nonNA_metadata, permutations = 999)[["aov.tab"]]
        
      })
      
      names(adonis_result_table_list) <- colnames(nonNA_metadata)[-1]
      
      for (i in colnames(nonNA_metadata)[-1]) {
        
        rownames(adonis_result_table_list[[i]])[1] <- i
      }
      
      adonis_result_table_list_show <- adonis_result_table_list[[input$metadata_phylo_beta]][1, c(5,6)]
      colnames(adonis_result_table_list_show) <- c("R^2", "P value")
      
      
      
      adonis_result_table_list_show[, "R^2"] <- adonis_result_table_list_show[, "R^2"] %>% round(4) %>% as.character()
      adonis_result_table_list_show[, "P value"] <- adonis_result_table_list_show[, "P value"] %>% round(4) %>% as.character()
      
      return(adonis_result_table_list_show)
    })
    
    Permanova_table_phylo_W <- reactive({
      
      nonNA_position <- which(Metadata_stats()[, input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(W_unifrac_dm))
      
      adonis_result_table_list <- lapply(colnames(nonNA_metadata)[-1], function(x){
        
        adonis(W_unifrac_dm ~ nonNA_metadata[,x], nonNA_metadata, permutations = 999)[["aov.tab"]]
        
      })
      
      names(adonis_result_table_list) <- colnames(nonNA_metadata)[-1]
      
      for (i in colnames(nonNA_metadata)[-1]) {
        
        rownames(adonis_result_table_list[[i]])[1] <- i
      }
      
      adonis_result_table_list_show <- adonis_result_table_list[[input$metadata_phylo_beta]][1, c(5,6)]
      colnames(adonis_result_table_list_show) <- c("R^2", "P value")
      
      
      
      adonis_result_table_list_show[, "R^2"] <- adonis_result_table_list_show[, "R^2"] %>% round(4) %>% as.character()
      adonis_result_table_list_show[, "P value"] <- adonis_result_table_list_show[, "P value"] %>% round(4) %>% as.character()
      
      return(adonis_result_table_list_show)
    })
    
    # Show permanova title
    output$Permanova_title_phylo <- renderText({
      
      if(is.null(Permanova_table_phylo_unW()) && is.null(Permanova_table_phylo_W())){
        
        return(NULL)
      }else{
        return("PerMANOVA")
      }
      
      
    })
    
    # Show permanova table
    output$permanova_table_phylo <- renderTable({
      
      permanova_table_phylo_list <- list(Permanova_table_phylo_unW(), Permanova_table_phylo_W())
      names(permanova_table_phylo_list) <- c("Unweighted", "Weighted")
      return(permanova_table_phylo_list[[input$UnW_or_W_phylo]])
      
    })
    
    # Download permanova table
    output$download_permanova_phylo <- downloadHandler(
      
      filename = function(){
        paste0("UniFrac_PerMANOVA_table_", input$UnW_or_W_phylo, ".csv")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          write_csv(Permanova_table_phylo_unW(), file)
        }else{
          write_csv(Permanova_table_phylo_W(), file)
        }
      }
    )
    
    
    # When length(group_names)<=2, hide the download button of pair table
    observe({
      
      req(input$sample_data, input$taxonomic_table)
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      group_names <- unique(nonNA_metadata[, input$metadata_phylo_beta])
      
      if(length(group_names)<=2){
        shinyjs::hide("download_permanova_pair_phylo")
        shinyjs::hide("download_ANOSIM_pair_phylo")
        shinyjs::hide("download_MRPP_pair_phylo")
      }
      
    })
    
    Permanova_pair_table_phylo_unW <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(unW_unifrac_dm))
      
      group_names <- unique(nonNA_metadata[,input$metadata_phylo_beta])
      
      if(length(group_names)>2){
        
        sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
          
          sample_cb_meta <- nonNA_metadata[,c(1,x)]
          
        })
        names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
        
        
        unW_unifrac_dm_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
          
          feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
          sample_f1_f2 <- sample_list[[input$metadata_phylo_beta]] %>% filter(sample_list[[input$metadata_phylo_beta]][,input$metadata_phylo_beta] %in% feature1_feature2)
          unW_unifrac_ds <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]]
          return(usedist::dist_subset(unW_unifrac_ds, sample_f1_f2[,"SampleID"]))
        })
        
        
        adonis_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
          
          unW_unifrac_dm_pair <- unW_unifrac_dm_f1_f2_list[[k]] %>% as.matrix() %>% as.data.frame()
          metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(unW_unifrac_dm_pair))
          adonis_result_pair <- adonis(unW_unifrac_dm_pair ~ metadata_pair[, input$metadata_phylo_beta], metadata_pair, permutations = 999)
          
        })
        
        pair_names <- c()
        df_pair<- t(combn(group_names, 2))
        for (i in 1:ncol(combn(group_names, 2))) {
          
          pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
          
        }
        names(adonis_result_pair_list) <- pair_names
        
        adonis_result_pair_list_table <- as.character(adonis_result_pair_list[[1]][["aov.tab"]][1,])
        for (i in 1:(ncol(combn(group_names, 2))-1)) {
          adonis_result_pair_list_table <- rbind(adonis_result_pair_list_table, 
                                                 as.character(adonis_result_pair_list[[i+1]][["aov.tab"]][1,]))
        }
        
        colnames(adonis_result_pair_list_table) <- c("Df", "Sums Of Sqs", "Mean Sqs", "F.Model", "R^2", "P value")
        adonis_result_pair_list_table <- cbind(Comparisons=pair_names, adonis_result_pair_list_table)
        
        for (i in 3:7) {
          
          adonis_result_pair_list_table[, i] <- adonis_result_pair_list_table[, i] %>% as.numeric() %>% round(4) %>% as.character()
          
        }
        
        return(as.data.frame(adonis_result_pair_list_table[, c(1, 6, 7)]))
        
      }else{
        return(NULL)
      }
      
    })
    
    Permanova_pair_table_phylo_W <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(W_unifrac_dm))
      
      group_names <- unique(nonNA_metadata[,input$metadata_phylo_beta])
      
      if(length(group_names)>2){
        
        sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
          
          sample_cb_meta <- nonNA_metadata[,c(1,x)]
          
        })
        names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
        
        
        W_unifrac_dm_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
          
          feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
          sample_f1_f2 <- sample_list[[input$metadata_phylo_beta]] %>% filter(sample_list[[input$metadata_phylo_beta]][,input$metadata_phylo_beta] %in% feature1_feature2)
          W_unifrac_ds <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]]
          return(usedist::dist_subset(W_unifrac_ds, sample_f1_f2[,"SampleID"]))
        })
        
        
        adonis_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
          
          W_unifrac_dm_pair <- W_unifrac_dm_f1_f2_list[[k]] %>% as.matrix() %>% as.data.frame()
          metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(W_unifrac_dm_pair))
          adonis_result_pair <- adonis(W_unifrac_dm_pair ~ metadata_pair[, input$metadata_phylo_beta], metadata_pair, permutations = 999)
          
        })
        
        pair_names <- c()
        df_pair<- t(combn(group_names, 2))
        for (i in 1:ncol(combn(group_names, 2))) {
          
          pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
          
        }
        names(adonis_result_pair_list) <- pair_names
        
        adonis_result_pair_list_table <- as.character(adonis_result_pair_list[[1]][["aov.tab"]][1,])
        for (i in 1:(ncol(combn(group_names, 2))-1)) {
          adonis_result_pair_list_table <- rbind(adonis_result_pair_list_table, 
                                                 as.character(adonis_result_pair_list[[i+1]][["aov.tab"]][1,]))
        }
        
        colnames(adonis_result_pair_list_table) <- c("Df", "Sums Of Sqs", "Mean Sqs", "F.Model", "R^2", "P value")
        adonis_result_pair_list_table <- cbind(Comparisons=pair_names, adonis_result_pair_list_table)
        
        for (i in 3:7) {
          
          adonis_result_pair_list_table[, i] <- adonis_result_pair_list_table[, i] %>% as.numeric() %>% round(4) %>% as.character()
          
        }
        
        return(as.data.frame(adonis_result_pair_list_table[, c(1, 6, 7)]))
        
      }else{
        return(NULL)
      }
      
    })
    
    # Show permanova pair title
    output$Permanova_pair_title_phylo <- renderText({
      
      if(is.null(Permanova_pair_table_phylo_unW()) && is.null(Permanova_pair_table_phylo_W())){
        
        return(NULL)
      }else{
        return("PerMANOVA pair")
      }
      
      
    })
    
    # Show permanova pair table
    output$permanova_pair_table_phylo <- renderTable({
      
      permanova_pair_table_phylo_list <- list(Permanova_pair_table_phylo_unW(), Permanova_pair_table_phylo_W())
      names(permanova_pair_table_phylo_list) <- c("Unweighted", "Weighted")
      return(permanova_pair_table_phylo_list[[input$UnW_or_W_phylo]])
      
    })
    
    # Download permanova pair table
    output$download_permanova_pair_phylo <- downloadHandler(
      
      filename = function(){
        paste0("Unifrac_PerMANOVA_pair_table_", input$UnW_or_W_phylo, ".csv")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          write_csv(Permanova_pair_table_phylo_unW(), file)
        }else{
          write_csv(Permanova_pair_table_phylo_W(), file)
        }
      }
    )
    
    
    ANOSIM_table_phylo_unW <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(unW_unifrac_dm))
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      anosim_result <- anosim(unW_unifrac_dm, nonNA_metadata[, input$metadata_phylo_beta], permutations = 999)
      anosim_result_table <- data.frame(Significance=anosim_result$signif, statistic.R=anosim_result$statistic)
      anosim_result_table <- anosim_result_table[,c(2,1)]
      colnames(anosim_result_table) <- c("R", "P value")
      
      anosim_result_table$R <- anosim_result_table$R %>% as.numeric() %>% round(4) %>% as.character()
      anosim_result_table$"P value" <- anosim_result_table$"P value" %>% as.numeric() %>% round(4) %>% as.character()
      
      return(anosim_result_table)
      
    })
    ANOSIM_table_phylo_W <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(W_unifrac_dm))
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      anosim_result <- anosim(W_unifrac_dm, nonNA_metadata[, input$metadata_phylo_beta], permutations = 999)
      anosim_result_table <- data.frame(Significance=anosim_result$signif, statistic.R=anosim_result$statistic)
      anosim_result_table <- anosim_result_table[,c(2,1)]
      colnames(anosim_result_table) <- c("R", "P value")
      
      anosim_result_table$R <- anosim_result_table$R %>% as.numeric() %>% round(4) %>% as.character()
      anosim_result_table$"P value" <- anosim_result_table$"P value" %>% as.numeric() %>% round(4) %>% as.character()
      
      return(anosim_result_table)
      
    })
    ANOSIM_pair_table_phylo_unW <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(unW_unifrac_dm))
      
      group_names <- unique(nonNA_metadata[,input$metadata_phylo_beta])
      
      if(length(group_names)>2){
        
        sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
          
          sample_cb_meta <- nonNA_metadata[,c(1,x)]
          
        })
        
        names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
        
        
        unW_unifrac_dm_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
          
          feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
          sample_f1_f2 <- sample_list[[input$metadata_phylo_beta]] %>% filter(sample_list[[input$metadata_phylo_beta]][,input$metadata_phylo_beta] %in% feature1_feature2)
          unW_unifrac_ds <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]]
          return(usedist::dist_subset(unW_unifrac_ds, sample_f1_f2[,"SampleID"]))
        })
        
        anosim_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
          
          unW_unifrac_dm_pair <- unW_unifrac_dm_f1_f2_list[[k]] %>% as.matrix() %>% as.data.frame()
          metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(unW_unifrac_dm_pair))
          adonis_result_pair <- anosim(unW_unifrac_dm_pair, metadata_pair[, input$metadata_phylo_beta], metadata_pair, permutations = 999)
          
        })
        
        pair_names <- c()
        df_pair <- t(combn(group_names, 2))
        for (i in 1:ncol(combn(group_names, 2))) {
          
          pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
          
        }
        names(anosim_result_pair_list) <- pair_names
        
        anosim_result_pair_list_signif <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          
          anosim_result_pair_list_signif[i] <- anosim_result_pair_list[[i]][["signif"]]
        }
        
        anosim_result_pair_list_R <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          
          anosim_result_pair_list_R[i] <- anosim_result_pair_list[[i]][["statistic"]]
        }
        
        anosim_result_pair_list_table <- data.frame(Comparisons=pair_names,
                                                    R=anosim_result_pair_list_R,
                                                    "P_value"=anosim_result_pair_list_signif)
        
        anosim_result_pair_list_table$R <- anosim_result_pair_list_table$R %>% as.numeric() %>% round(4) %>% as.character()
        anosim_result_pair_list_table$P_value <- anosim_result_pair_list_table$P_value %>% as.numeric() %>% round(4) %>% as.character()
        colnames(anosim_result_pair_list_table)[3] <- "P value"
        
        return(anosim_result_pair_list_table)
        
      }else{
        return(NULL)
      }
      
    })
    ANOSIM_pair_table_phylo_W <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(W_unifrac_dm))
      
      group_names <- unique(nonNA_metadata[,input$metadata_phylo_beta])
      
      if(length(group_names)>2){
        
        sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
          
          sample_cb_meta <- nonNA_metadata[,c(1,x)]
          
        })
        
        names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
        
        
        W_unifrac_dm_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
          
          feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
          sample_f1_f2 <- sample_list[[input$metadata_phylo_beta]] %>% filter(sample_list[[input$metadata_phylo_beta]][,input$metadata_phylo_beta] %in% feature1_feature2)
          W_unifrac_ds <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]]
          return(usedist::dist_subset(W_unifrac_ds, sample_f1_f2[,"SampleID"]))
        })
        
        anosim_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
          
          W_unifrac_dm_pair <- W_unifrac_dm_f1_f2_list[[k]] %>% as.matrix() %>% as.data.frame()
          metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(W_unifrac_dm_pair))
          anosim_result_pair <- anosim(W_unifrac_dm_pair, metadata_pair[,input$metadata_phylo_beta], metadata_pair, permutations = 999)
          
        })
        
        pair_names <- c()
        df_pair <- t(combn(group_names, 2))
        for (i in 1:ncol(combn(group_names, 2))) {
          
          pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
          
        }
        names(anosim_result_pair_list) <- pair_names
        
        anosim_result_pair_list_signif <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          
          anosim_result_pair_list_signif[i] <- anosim_result_pair_list[[i]][["signif"]]
        }
        
        anosim_result_pair_list_R <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          
          anosim_result_pair_list_R[i] <- anosim_result_pair_list[[i]][["statistic"]]
        }
        
        anosim_result_pair_list_table <- data.frame(Comparisons=pair_names,
                                                    R=anosim_result_pair_list_R,
                                                    "P_value"=anosim_result_pair_list_signif)
        
        anosim_result_pair_list_table$R <- anosim_result_pair_list_table$R %>% as.numeric() %>% round(4) %>% as.character()
        anosim_result_pair_list_table$P_value <- anosim_result_pair_list_table$P_value %>% as.numeric() %>% round(4) %>% as.character()
        colnames(anosim_result_pair_list_table)[3] <- "P value"
        
        return(anosim_result_pair_list_table)
        
      }else{
        return(NULL)
      }
      
    })
    
    # Show ANOSIM title
    output$ANOSIM_title_phylo <- renderText({
      
      if(is.null(ANOSIM_table_phylo_unW()) && is.null(ANOSIM_table_phylo_W())){
        
        return(NULL)
      }else{
        return("ANOSIM")
      }
      
      
    })
    
    # Show ANOSIM table
    output$ANOSIM_table_phylo <- renderTable({
      
      ANOSIM_table_phylo_list <- list(ANOSIM_table_phylo_unW(), ANOSIM_table_phylo_W())
      names(ANOSIM_table_phylo_list) <- c("Unweighted", "Weighted")
      return(ANOSIM_table_phylo_list[[input$UnW_or_W_phylo]])
      
    })
    
    # Download ANOSIM table
    output$download_ANOSIM_phylo <- downloadHandler(
      
      filename = function(){
        paste0("Unifrac_ANOSIM_table_", input$UnW_or_W_phylo, ".csv")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          write_csv(ANOSIM_table_phylo_unW(), file)
        }else{
          write_csv(ANOSIM_table_phylo_W(), file)
        }
      }
    )
    
    # Show ANOSIM pair title
    output$ANOSIM_pair_title_phylo <- renderText({
      
      if(is.null(ANOSIM_pair_table_phylo_unW()) && is.null(ANOSIM_pair_table_phylo_W())){
        
        return(NULL)
      }else{
        return("ANOSIM pair")
      }
      
    })
    
    # Show ANOSIM pair table
    output$ANOSIM_pair_table_phylo <- renderTable({
      
      ANOSIM_pair_table_phylo_list <- list(ANOSIM_pair_table_phylo_unW(), ANOSIM_pair_table_phylo_W())
      names(ANOSIM_pair_table_phylo_list) <- c("Unweighted", "Weighted")
      return(ANOSIM_pair_table_phylo_list[[input$UnW_or_W_phylo]])
      
    })
    
    
    
    
    # Download ANOSIM pair table
    output$download_ANOSIM_pair_phylo <- downloadHandler(
      
      filename = function(){
        paste0("Unifrac_ANOSIM_pair_table_", input$UnW_or_W_phylo, ".csv")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          write_csv(ANOSIM_pair_table_phylo_unW(), file)
        }else{
          write_csv(ANOSIM_pair_table_phylo_W(), file)
        }
      }
    )
    
    
    MRPP_table_phylo_unW <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(unW_unifrac_dm))
      
      mrpp_result <- mrpp(unW_unifrac_dm, nonNA_metadata[, input$metadata_phylo_beta], permutations = 999)
      mrpp_result_table <- data.frame(Group = 'all', 
                                      Distance = 'Bray-Curtis', 
                                      A = mrpp_result$A,            
                                      Observe.delta = mrpp_result$delta,            
                                      Expect.delta = mrpp_result$E.delta,            
                                      P.value = mrpp_result$Pvalue)
      mrpp_result_table_show <- mrpp_result_table[,3:6]
      colnames(mrpp_result_table_show)[4] <- "P value"
      
      for (i in 1:4) {
        
        mrpp_result_table_show[, i] <- mrpp_result_table_show[, i] %>% round(4) %>% as.character()
        
      }
      
      return(mrpp_result_table_show)
    })
    MRPP_table_phylo_W <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      # bray_df <- vegdist(t(nonNA_taxtable), method = "bray") %>% as.matrix() %>% as.data.frame()
      W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(W_unifrac_dm))
      
      mrpp_result <- mrpp(W_unifrac_dm, nonNA_metadata[, input$metadata_phylo_beta], permutations = 999)
      mrpp_result_table <- data.frame(Group = 'all', 
                                      Distance = 'Bray-Curtis', 
                                      A = mrpp_result$A,            
                                      Observe.delta = mrpp_result$delta,            
                                      Expect.delta = mrpp_result$E.delta,            
                                      P.value = mrpp_result$Pvalue)
      mrpp_result_table_show <- mrpp_result_table[,3:6]
      colnames(mrpp_result_table_show)[4] <- "P value"
      
      for (i in 1:4) {
        
        mrpp_result_table_show[, i] <- mrpp_result_table_show[, i] %>% round(4) %>% as.character()
        
      }
      
      return(mrpp_result_table_show)
    })
    MRPP_pair_table_phylo_unW <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      unW_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(unW_unifrac_dm))
      
      group_names <- unique(nonNA_metadata[,input$metadata_phylo_beta])
      
      if(length(group_names)>2){
        
        sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
          
          sample_cb_meta <- nonNA_metadata[,c(1,x)]
          
        })
        names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
        
        unW_unifrac_dm_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
          
          feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
          sample_f1_f2 <- sample_list[[input$metadata_phylo_beta]] %>% filter(sample_list[[input$metadata_phylo_beta]][,input$metadata_phylo_beta] %in% feature1_feature2)
          unW_unifrac_ds <- read_qza("/home/imuser/qiime_output/core-metrics-results/unweighted_unifrac_distance_matrix.qza")[["data"]]
          return(usedist::dist_subset(unW_unifrac_ds, sample_f1_f2[,"SampleID"]))
          
        })
        
        MRPP_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
          
          unW_unifrac_dm_pair <- unW_unifrac_dm_f1_f2_list[[k]] %>% as.matrix() %>% as.data.frame()
          metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(unW_unifrac_dm_pair))
          MRPP_result_pair <- mrpp(unW_unifrac_dm_pair, metadata_pair[, input$metadata_phylo_beta], permutations = 999)
          
        })
        
        pair_names <- c()
        df_pair <- t(combn(group_names, 2))
        for (i in 1:ncol(combn(group_names, 2))) {
          
          pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
          
        }
        names(MRPP_result_pair_list) <- pair_names
        
        MRPP_result_pair_list_Pvalue <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_Pvalue[i] <- MRPP_result_pair_list[[i]][["Pvalue"]]
        }
        
        MRPP_result_pair_list_A <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_A[i] <- MRPP_result_pair_list[[i]][["A"]]
        }
        
        MRPP_result_pair_list_delta <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_delta[i] <- MRPP_result_pair_list[[i]][["delta"]]
        }
        
        MRPP_result_pair_list_Edelta <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_Edelta[i] <- MRPP_result_pair_list[[i]][["E.delta"]]
        }
        
        MRPP_result_pair_list_table <- data.frame(Comparisons=pair_names,
                                                  A=MRPP_result_pair_list_A,
                                                  delta=MRPP_result_pair_list_delta,
                                                  E.delta=MRPP_result_pair_list_Edelta,
                                                  P_value=MRPP_result_pair_list_Pvalue)
        
        for (i in 2:5) {
          
          MRPP_result_pair_list_table[, i] <- MRPP_result_pair_list_table[, i] %>% round(4) %>% as.character()
          
        }
        
        colnames(MRPP_result_pair_list_table)[5] <- "P value"
        
        return(MRPP_result_pair_list_table)
      }else{
        return(NULL)
      }
    })
    MRPP_pair_table_phylo_W <- reactive({
      
      nonNA_position <- which(Metadata_stats()[,input$metadata_phylo_beta] != "NA")
      nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
      
      nonNA_taxtable <- TaxaTable_merge()[,nonNA_sampleid]
      nonNA_metadata <- Metadata_stats()[nonNA_position, ]
      
      W_unifrac_dm <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]] %>% as.matrix() %>% as.data.frame()
      
      nonNA_metadata <- filter(nonNA_metadata, SampleID %in% rownames(W_unifrac_dm))
      
      group_names <- unique(nonNA_metadata[,input$metadata_phylo_beta])
      
      if(length(group_names)>2){
        
        sample_list <- lapply(2:ncol(nonNA_metadata), function(x){
          
          sample_cb_meta <- nonNA_metadata[,c(1,x)]
          
        })
        names(sample_list) <- colnames(nonNA_metadata)[2:ncol(nonNA_metadata)]
        
        W_unifrac_dm_f1_f2_list <- lapply(1:ncol(combn(group_names, 2)), function(i){
          
          feature1_feature2 <- t(combn(group_names, 2))[i,] %>% as.character()
          sample_f1_f2 <- sample_list[[input$metadata_phylo_beta]] %>% filter(sample_list[[input$metadata_phylo_beta]][,input$metadata_phylo_beta] %in% feature1_feature2)
          W_unifrac_ds <- read_qza("/home/imuser/qiime_output/core-metrics-results/weighted_unifrac_distance_matrix.qza")[["data"]]
          return(usedist::dist_subset(W_unifrac_ds, sample_f1_f2[,"SampleID"]))
          
        })
        
        MRPP_result_pair_list <- lapply(1:ncol(combn(group_names, 2)), function(k){
          
          W_unifrac_dm_pair <- W_unifrac_dm_f1_f2_list[[k]] %>% as.matrix() %>% as.data.frame()
          metadata_pair <- nonNA_metadata %>% filter(SampleID %in% colnames(W_unifrac_dm_pair))
          MRPP_result_pair <- mrpp(W_unifrac_dm_pair, metadata_pair[, input$metadata_phylo_beta], permutations = 999)
          
        })
        
        pair_names <- c()
        df_pair <- t(combn(group_names, 2))
        for (i in 1:ncol(combn(group_names, 2))) {
          
          pair_names[i] <- paste(df_pair[i,1], df_pair[i,2], sep = " - ")
          
        }
        names(MRPP_result_pair_list) <- pair_names
        
        MRPP_result_pair_list_Pvalue <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_Pvalue[i] <- MRPP_result_pair_list[[i]][["Pvalue"]]
        }
        
        MRPP_result_pair_list_A <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_A[i] <- MRPP_result_pair_list[[i]][["A"]]
        }
        
        MRPP_result_pair_list_delta <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_delta[i] <- MRPP_result_pair_list[[i]][["delta"]]
        }
        
        MRPP_result_pair_list_Edelta <- c()
        for (i in 1:ncol(combn(group_names, 2))) {
          MRPP_result_pair_list_Edelta[i] <- MRPP_result_pair_list[[i]][["E.delta"]]
        }
        
        MRPP_result_pair_list_table <- data.frame(Comparisons=pair_names,
                                                  A=MRPP_result_pair_list_A,
                                                  delta=MRPP_result_pair_list_delta,
                                                  E.delta=MRPP_result_pair_list_Edelta,
                                                  P_value=MRPP_result_pair_list_Pvalue)
        
        for (i in 2:5) {
          
          MRPP_result_pair_list_table[, i] <- MRPP_result_pair_list_table[, i] %>% round(4) %>% as.character()
          
        }
        
        colnames(MRPP_result_pair_list_table)[5] <- "P value"
        
        return(MRPP_result_pair_list_table)
      }else{
        return(NULL)
      }
    })
    
    # Show MRPP title
    output$MRPP_title_phylo <- renderText({
      
      if(is.null(MRPP_table_phylo_unW()) && is.null(MRPP_table_phylo_W())){
        
        return(NULL)
      }else{
        return("ANOSIM")
      }
      
      
    })
    
    # Show MRPP table
    output$MRPP_table_phylo <- renderTable({
      
      MRPP_table_phylo_list <- list(MRPP_table_phylo_unW(), MRPP_table_phylo_W())
      names(MRPP_table_phylo_list) <- c("Unweighted", "Weighted")
      return(MRPP_table_phylo_list[[input$UnW_or_W_phylo]])
      
    })
    
    # Download MRPP table
    output$download_MRPP_phylo <- downloadHandler(
      
      filename = function(){
        paste0("Unifrac_MRPP_table_", input$UnW_or_W_phylo, ".csv")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          write_csv(MRPP_table_phylo_unW(), file)
        }else{
          write_csv(MRPP_table_phylo_W(), file)
        }
      }
    )
    
    # Show MRPP pair title
    output$MRPP_pair_title_phylo <- renderText({
      
      if(is.null(MRPP_pair_table_phylo_unW()) && is.null(MRPP_pair_table_phylo_W())){
        
        return(NULL)
      }else{
        return("MRPP pair")
      }
      
      
    })
    
    # Show MRPP pair table
    output$MRPP_pair_table_phylo <- renderTable({
      
      MRPP_pair_table_phylo_list <- list(MRPP_pair_table_phylo_unW(), MRPP_pair_table_phylo_W())
      names(MRPP_pair_table_phylo_list) <- c("Unweighted", "Weighted")
      return(MRPP_pair_table_phylo_list[[input$UnW_or_W_phylo]])
      
    })
    
    # Download MRPP pair table
    output$download_MRPP_pair_phylo <- downloadHandler(
      
      filename = function(){
        paste0("Unifrac_MRPP_pair_table_", input$UnW_or_W_phylo, ".csv")
      },
      content = function(file){
        if(input$UnW_or_W_phylo=="Unweighted"){
          write_csv(MRPP_pair_table_phylo_unW(), file)
        }else{
          write_csv(MRPP_pair_table_phylo_W(), file)
        }
      }
    )
    
    # removeModal()
    remove_modal_spinner()
    
  })
  
  
  
  
  
  
  
  
  # ANCOM ------------------------------------------------------------------------------------------------------
  observeEvent(input$ANCOM_start, {
    
    # showModal(modalDialog(title = "Running ANCOM ...", "Waiting for a moment", footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    Sys.setenv(LANG="C.UTF-8")
    file.remove("/home/imuser/qiime_output/ancom_comparison.qzv")
    # selected_metadata <- Metadata_stats() %>% filter(Metadata_stats()[,input$metadata8] == input$metadata8_factor) %>% as_tibble()
    # colnames(selected_metadata)[1] <- "SampleID"
    # write.table(selected_metadata[,1], file='/home/imuser/selected_sampleID.tsv', quote=FALSE, sep='\t', row.names = F)
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    
    nonNA_position <- which(Metadata_stats()[, input$metadata_ANCOM]!="NA")
    nonNA_sampleid <- Metadata_stats()[,1][nonNA_position]
    nonNA_metadata <- Metadata_stats()[nonNA_position,]
    
    nonNA_metadata_categorical <- apply(as.matrix(nonNA_metadata[,2:ncol(nonNA_metadata)]), MARGIN = 2, FUN = function(x){
      if(sum(grepl("^[A-Za-z]+", x, perl = T)) > 0){
        return(x)
      }else{
        x <- paste0("[ ", x, " ]")
      }
      
    })
    
    # for (r in 1:nrow(nonNA_metadata)) {
    #   for (c in 2:ncol(nonNA_metadata)) {
    #     
    #     x <- nonNA_metadata[r,c]
    #     if(sum(grepl("^[A-Za-z]+", x, perl = T)) > 0){
    #       return(x)
    #     }else{
    #       x <- paste0("[ ", x, " ]")
    #     }
    #   }
    #   
    # }
    
    
    nonNA_metadata_1stcol <- nonNA_metadata[,1] %>% as.data.frame()
    colnames(nonNA_metadata_1stcol) <- "SampleID"
    nonNA_metadata_categorical <- as.data.frame(nonNA_metadata_categorical)
    colnames(nonNA_metadata[,2:ncol(nonNA_metadata)])
    nonNA_metadata_categorical <- cbind(nonNA_metadata_1stcol, nonNA_metadata_categorical)
    colnames(nonNA_metadata_categorical) <- c("SampleID", colnames(nonNA_metadata)[-1])
    
    write.table(nonNA_metadata_categorical, file='/home/imuser/nonNA_metadata_categorical.tsv', quote = F, sep='\t', row.names = F)
    
    # system("sed -i '2i\ #q2:types  categorical	categorical	categorical	categorical	categorical	categorical	categorical	categorical' /home/imuser/nonNA_metadata.tsv")
    
    # # add categorical line
    # categorical_word <- paste(rep("categorical", ncol(Metadata_stats())), collapse = " ")
    # type_word <- paste("#q2:types", categorical_word)
    # a <- read.csv('/home/imuser/nonNA_metadata.tsv', stringsAsFactors = F, header = F)
    # b <- rbind(type_word, a) %>% as_tibble()
    # colnames(b) <- b[2,]
    # c <- b[-2,]
    # write.table(c, file = "/home/imuser/nonNA_metadata_addline.tsv", sep = "/t", row.names = F, quote = F, col.names = T)
    
    file.copy(from = input$taxonomic_table$datapath, to = "/home/imuser/upload_taxtable.qza", overwrite = T)
    nonNA_sampleid_1 <- c("SampleID", nonNA_sampleid) # for qiime2 reading format
    write.table(x = nonNA_sampleid_1, file = "/home/imuser/nonNA_sampleid.tsv",quote = F, row.names = F, col.names = F)
    
    system(paste(qiime_cmd, "feature-table filter-samples", "--i-table /home/imuser/upload_taxtable.qza", "--m-metadata-file /home/imuser/nonNA_sampleid.tsv", "--o-filtered-table /home/imuser/qiime_output/nonNA_table.qza"))
    system(paste(qiime_cmd, "composition add-pseudocount --i-table", "/home/imuser/qiime_output/nonNA_table.qza", "--o-composition-table /home/imuser/qiime_output/comp_table.qza"))
    system(paste(qiime_cmd, "composition ancom --i-table /home/imuser/qiime_output/comp_table.qza --m-metadata-file", '/home/imuser/nonNA_metadata_categorical.tsv',
                 "--m-metadata-column", input$metadata_ANCOM,"--o-visualization /home/imuser/qiime_output/ancom_comparison.qzv"))
    
    unlink("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname", recursive = T)
    # system("cp /home/imuser/qiime_output/ancom_comparison.qzv /home/imuser/qiime_output/ancom_comparison.zip")
    system("unzip -d /home/imuser/qiime_output/ancom_comparison_unzip /home/imuser/qiime_output/ancom_comparison.qzv")
    ancom_unzip_dirname <- list.files("/home/imuser/qiime_output/ancom_comparison_unzip", full.names = T)
    system(paste0("mv ", ancom_unzip_dirname, " /home/imuser/qiime_output/ancom_comparison_unzip/new_dirname"))
    
    output$word_ancom_plotly <- renderUI({
      
      h3("ANCOM Volcano Plot", 
         style = "color: black;top: 10px;")
    })
    
    output$ancom_plotly <- renderPlotly({
      
      ancom_data <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/data.tsv", sep = "\t", header = T)
      ancom_sig <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/ancom.tsv", sep = "\t", header = T)
      names(ancom_sig)[1] <- "id"
      
      ancom_merge <- merge(x = ancom_data, y = ancom_sig[, c(1,3)], by = "id")
      
      plot_ly(data = ancom_merge,
              x = ~ clr,
              y = ~ W,
              color = ~ Reject.null.hypothesis,
              colors = c("grey","red")[1:length(unique(ancom_sig$Reject.null.hypothesis))],
              type = "scatter",
              hoverinfo = "text",
              hovertext = paste("Species:", ancom_merge$id,
                                "<br> clr:", ancom_merge$clr,
                                "<br> W:", ancom_merge$W) 
      ) %>% layout(showlegend = T)
      
    })
    
    
    output$annotation_ancom <- renderUI({
      HTML("<p>The W value is essentially a count of the number of sub-hypotheses that have passed for a given species.<br>The clr represent log-fold change relative to the average microbe.</p>")
    })
    
    output$word_ancom_table <- renderUI({
      
      h3("ANCOM statistical results (Speceis with significant w value)", 
         style = "color: black;top: 10px;")
    })
    
    output$ancom_sig <- renderDataTable({
      
      ancom_sig <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/ancom.tsv", sep = "\t", header = T)
      names(ancom_sig)[1] <- "Speceis"
      
      ancom_sig_true <- filter(ancom_sig, Reject.null.hypothesis == "True")
      return(ancom_sig_true[,1:2])
      
    })
    
    
    
    # removeModal()
    remove_modal_spinner()
    
    if (file.exists("/home/imuser/qiime_output/ancom_comparison.qzv")){
      
      # output$word_ANCOM <- renderText(print("ANCOM successfully!"))
      # showModal(modalDialog(title = strong("ANCOM has been finished."), 
      #                       "You can inspect the results!", 
      #                       footer = NULL, easyClose = T, size = "l"))
      
    }else if (is.character(Metadata_stats()[, input$metadata_ANCOM])==F){
      
      # output$word_ANCOM <- renderText(print("Error!The factor should be categorical data."))
      showModal(modalDialog(title = strong("Error!", style = "color: red"),
                            "The factor should be categorical data.", 
                            footer = NULL, easyClose = T, size = "l"))
    }else{
      showModal(modalDialog(title = strong("Error!", style = "color: red"),
                            "Please check your files.", 
                            footer = NULL, easyClose = T, size = "l"))
    }
    
    
  })
  
  
  output$ancom_plot_download <- downloadHandler(
    filename = "ANCOM_plot.png",
    content = function(file){
      ancom_data <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/data.tsv", sep = "\t", header = T)
      ancom_sig <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/ancom.tsv", sep = "\t", header = T)
      names(ancom_sig)[1] <- "id"
      
      ancom_merge <- merge(x = ancom_data, y = ancom_sig[, c(1,3)], by = "id")
      
      k <- ggplot(data = ancom_merge, 
                  aes(clr, W, col = Reject.null.hypothesis))+geom_point()+scale_colour_discrete("Reject null hypothesis")+theme_bw()
      
      ggsave(file, plot = k, width = 80, height = 40, units = "cm")
    }
  )
  
  output$ancom_table_download <- downloadHandler(
    filename = "ANCOM_table.csv",
    content = function(file){
      
      ancom_data <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/data.tsv", sep = "\t", header = T)
      ancom_sig <- read.table("/home/imuser/qiime_output/ancom_comparison_unzip/new_dirname/data/ancom.tsv", sep = "\t", header = T)
      names(ancom_sig)[1] <- "id"
      
      ancom_merge <- merge(x = ancom_data, y = ancom_sig[, c(1,3)], by = "id")
      write_csv(ancom_merge, file, col_names = T)
    }
  )
  
  
  
  
  # Function Analysis-------------------------------------------------------------------------------------------
  observeEvent(input$function_analysis, {
    
    # showModal(modalDialog(title = "Running FAPROTAX ...", "Waiting for a moment", footer = NULL))
    show_modal_spinner(spin = "circle", color = "#317EAC", text = "Please wait...")
    
    qiime_cmd <- '/home/imuser/miniconda3/envs/qiime2-2019.10/bin/qiime'
    
    system(paste(qiime_cmd, "tools export --input-path", input$taxonomic_table_FA$datapath, "--output-path /home/imuser/qiime_output/exported-feature-table7"))
    system("/home/imuser/miniconda3/envs/python-2.7/bin/python2.7 /home/imuser/FAPROTAX_1.2.1/collapse_table.py --force -i /home/imuser/qiime_output/exported-feature-table7/feature-table.biom -o /home/imuser/FAPROTAX_output/func-table7.biom -g /home/imuser/FAPROTAX_1.2.1/FAPROTAX.txt -r /home/imuser/FAPROTAX_output/report7-record.txt --out_groups2records_table /home/imuser/FAPROTAX_output/groups2record.biom")
    system(paste(qiime_cmd, "tools import --input-path /home/imuser/FAPROTAX_output/func-table7.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path /home/imuser/qiime_output/func-table7.qza"))
    system(paste(qiime_cmd, "tools import --input-path /home/imuser/FAPROTAX_output/groups2record.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path /home/imuser/qiime_output/groups2record.qza"))
    
    # removeModal() 
    remove_modal_spinner()
    
    
    output$func_table_BY_sampleid <- renderDataTable({
      
      func_table_BY_sampleid <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
      
      TF_all0 <- apply(func_table_BY_sampleid, 1, function(x) !all(x==0))
      
      func_table_BY_sampleid_filtered <- func_table_BY_sampleid[TF_all0,]
      
      func_table_BY_sampleid_filtered_tibble <- cbind("Type"=rownames(func_table_BY_sampleid_filtered), func_table_BY_sampleid_filtered) %>% as_tibble()
      
      return(func_table_BY_sampleid_filtered_tibble)
      
    })
    
    
    # output$func_table_BY_speciesname <- renderDataTable({
    #   
    #   func_table_BY_sampleid <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
    #   
    #   TF_all0 <- apply(func_table_BY_sampleid, 1, function(x) !all(x==0))
    #   
    #   func_table_BY_sampleid_filtered <- func_table_BY_sampleid[TF_all0,]
    #   
    #   func_name_filtered <- rownames(func_table_BY_sampleid_filtered)
    #   
    #   
    #   func_table_BY_speciesname <- read_qza("/home/imuser/qiime_output/groups2record.qza")[["data"]]
    #   
    #   func_table_BY_speciesname_filtered <- func_table_BY_speciesname[,func_name_filtered]
    #   
    #   func_table_BY_speciesname_filtered_tibble <- cbind("Species names"=rownames(func_table_BY_speciesname_filtered), func_table_BY_speciesname_filtered) %>% as_tibble()
    #   
    #   return(func_table_BY_speciesname_filtered_tibble)
    # })
    
    
    # output$Function_Heatmap <- renderPlotly({
    #   
    #   func_table_forHM <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
    #   TF_all0 <- apply(func_table_forHM, 1, function(x) !all(x==0))
    #   func_table_forHM_filtered <- func_table_forHM[TF_all0,]
    #   palette <- colorRampPalette(c("darkblue", "blue", "lightblue1",
    #                                 "green","yellow", "red", "darkred"))
    #   set.seed(9876)
    #   
    #   a <- list(tickangle = 45)
    #   plot_ly(x=colnames(func_table_forHM_filtered),
    #           y=rownames(func_table_forHM_filtered),
    #           z=func_table_forHM_filtered, 
    #           colors = palette(50),
    #           type="heatmap") %>% layout(xaxis = a, showlegend = FALSE)
    #   
    # })
    output$Function_barplot <- renderPlot({
      
      func_table_BY_sampleid <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
      
      TF_all0 <- apply(func_table_BY_sampleid, 1, function(x) !all(x==0))
      
      func_table_BY_sampleid_filtered <- func_table_BY_sampleid[TF_all0,]
      
      func_table_BY_sampleid_filtered_tibble <- cbind("Type"=rownames(func_table_BY_sampleid_filtered), func_table_BY_sampleid_filtered) %>% as_tibble()
      
      func_table_BY_sampleid_filtered_tibble[,-1] <- apply(func_table_BY_sampleid_filtered_tibble[,-1], MARGIN = 1, FUN = as.numeric)
      
      # df_barplot <- data.frame(
      #   Type = func_table_BY_sampleid_filtered_tibble[,1],
      #   reads = rowSums(func_table_BY_sampleid_filtered_tibble[,-1]),
      #   treatment = rep(c("low", "hight"), length(func_table_BY_sampleid_filtered_tibble[,1])),
      #   stringsAsFactors = F
      # )
      
      df_barplot <- gather(func_table_BY_sampleid_filtered_tibble, Sample_ID, reads, -Type)
      
      df_barplot$feature <- mapvalues(df_barplot$Sample_ID,
                                      from = Metadata_FA()[,1],
                                      to = as.character(Metadata_FA()[, input$metadata10]))
      
      ggplot(df_barplot,
             aes(y = reads, fill = feature, x = Type)) + geom_bar(stat = "identity", position = "dodge", alpha = 1, width = .8) + coord_flip() + theme_bw() + guides(fill=guide_legend(title=input$metadata10)) + ggtitle("Function prediction")
      
    })
    
    output$function_report <- renderUI({
      
      a <- read_table("/home/imuser/FAPROTAX_output/report7-record.txt") %>% as.data.frame()
      a_report <- a[104:106,2]
      a_report[1] <- str_replace_all(a_report[1], pattern = "record", replacement = "OTU")
      a_report[1] <- str_replace_all(a_report[1], pattern = "group", replacement = "functional group")
      a_report[2] <- str_replace_all(a_report[2], pattern = "record", replacement = "OTU")
      a_report[2] <- str_replace_all(a_report[2], pattern = "group", replacement = "functional group")
      a_report[2] <- str_remove(a_report[2], pattern = "\\(leftovers\\)")
      a_report[3] <- str_replace_all(a_report[3], pattern = "record", replacement = "OTU")
      a_report[3] <- str_replace_all(a_report[3], pattern = "group", replacement = "functional group")
      
      HTML(paste("<h4 style='margin-bottom: -10px'>Summary</h4>",a_report[1], a_report[2] ,a_report[3], sep = "<br/>"))
      
      
    })
    
    
  })
  
  
  
  # output$word_FA1 <- renderText({
  #   
  #   if(is.null(input$taxonomic_table_FA$datapath) == F){
  #     
  #     return("")
  #     
  #   }else{
  #     return("Please upload the file from the last steps of the sequnecing preprocessing.")
  #   }
  # })
  
  
  
  
  # output$word_FA2 <- renderText({
  #   
  #   if(file.exists("/home/imuser/qiime_output/func-table7.qza")){
  #     
  #     return("")
  #     
  #   }else{
  #     
  #     return("Please finish the steps of the sequnecing preprocessing.")
  #   }
  # })
  
  
  
  
  # output$word_FA3 <- renderText({
  #   
  #   if(file.exists("/home/imuser/qiime_output/func-table7.qza")){
  #     
  #     return("")
  #     
  #   }else{
  #     
  #     return("Please finish the steps of the sequnecing preprocessing.")
  #   }
  # })
  
  
  
  
  
  output$func_table_ID<-downloadHandler(
    
    filename = "function_table_bySampleID.csv",
    content = function(file) {
      
      func_table_BY_sampleid <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
      
      TF_all0 <- apply(func_table_BY_sampleid, 1, function(x) !all(x==0))
      
      func_table_BY_sampleid_filtered <- func_table_BY_sampleid[TF_all0,]
      
      func_table_BY_sampleid_filtered_tibble <- cbind("Type"=rownames(func_table_BY_sampleid_filtered), func_table_BY_sampleid_filtered) %>% as_tibble()
      
      write.csv(func_table_BY_sampleid_filtered_tibble, file, row.names = F)
      
    }
  )
  
  output$func_table_Sp<-downloadHandler(
    
    filename = "function_table_bySpeciesName.csv",
    content = function(file) {
      
      func_table_BY_sampleid <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
      
      TF_all0 <- apply(func_table_BY_sampleid, 1, function(x) !all(x==0))
      
      func_table_BY_sampleid_filtered <- func_table_BY_sampleid[TF_all0,]
      
      func_name_filtered <- rownames(func_table_BY_sampleid_filtered)
      
      
      func_table_BY_speciesname <- read_qza("/home/imuser/qiime_output/groups2record.qza")[["data"]]
      
      func_table_BY_speciesname_filtered <- func_table_BY_speciesname[,func_name_filtered]
      
      func_table_BY_speciesname_filtered_tibble <- cbind("Species names"=rownames(func_table_BY_speciesname_filtered), func_table_BY_speciesname_filtered) %>% as_tibble()
      
      write.csv(func_table_BY_speciesname_filtered_tibble, file, row.names = F)
      
    }
  )
  
  
  output$FA_plot_download <- downloadHandler(
    
    filename = "Functional_analysis_plot.png",
    
    content = function(file){
      
      func_table_BY_sampleid <- read_qza("/home/imuser/qiime_output/func-table7.qza")[["data"]]
      
      TF_all0 <- apply(func_table_BY_sampleid, 1, function(x) !all(x==0))
      
      func_table_BY_sampleid_filtered <- func_table_BY_sampleid[TF_all0,]
      
      func_table_BY_sampleid_filtered_tibble <- cbind("Type"=rownames(func_table_BY_sampleid_filtered), func_table_BY_sampleid_filtered) %>% as_tibble()
      
      func_table_BY_sampleid_filtered_tibble[,-1] <- apply(func_table_BY_sampleid_filtered_tibble[,-1], MARGIN = 1, FUN = as.numeric)
      
      df_barplot <- gather(func_table_BY_sampleid_filtered_tibble, Sample_ID, reads, -Type)
      
      df_barplot$feature <- mapvalues(df_barplot$Sample_ID,
                                      from = Metadata_FA()[,1],
                                      to = as.character(Metadata_FA()[, input$metadata10]))
      
      FA_plot <- ggplot(df_barplot,
                        aes(y = reads, fill = feature, x = Type)) + geom_bar(stat = "identity", position = "dodge", alpha = 1, width = .8) + coord_flip() + theme_bw() + guides(fill=guide_legend(title=input$metadata10))
      
      ggsave(file, plot = FA_plot, width = 80, height = 40, units = "cm")
    }
  )
  
  
  # Tutorial ---------------------------------------------------------------------------------------------------
  # output$tutorial <- renderUI({
  #     # p <- includeHTML("/home/imuser/text_files/tutorial.html")
  #     # return(HTML(p))
  #   HTML(markdown::markdownToHTML(knit("/home/imuser/text_files/tutorial.Rmd", quiet = T)))
  # })
  
}

