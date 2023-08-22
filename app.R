if (!require(shiny)) install.packages('shiny')
library(shiny)
if (!require(ids)) install.packages('ids')
library(ids)
if (!require(shinyFiles)) install.packages('shinyFiles')
library(shinyFiles)
if (!require(readr)) install.packages('readr')
library(readr)
if (!require(DT)) install.packages('DT')
library(DT)
if (!require("TmCalculator")) install.packages("TmCalculator")
library("TmCalculator")
 if (!require("Peptides")) install.packages("Peptides", dependencies=TRUE)
library("Peptides")
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("DT")) install.packages("DT")
library(DT)
if (!require("bslib")) install.packages("bslib", dependencies=TRUE)
library(bslib)
if (!require("stringr")) install.packages("stringr")
library(stringr)
if (!require("openxlsx")) install.packages("openxlsx")


# Create the Construct Design data frame
Construct_Design_df <- data.frame(
  Position = character(0),
  Target_Name = character(0),
  Vector_Name = character(0),
  Template_DNA_ID = character(0),
  N_Terminal_Residue = character(0),
  C_Terminal_Residue = character(0),
  N_Terminal_Extension_DNA = character(0),
  C_Terminal_Extension_DNA = character(0),
  Forward_Primer_Name = character(0),
  Reverse_Primer_Name = character(0),
  PCR_Product_Name = character(0),
  PCR_Product_Sequence = character(0),
  PCR_Product_Length = integer(0),
  N_Terminal_Residue_ = character(0),
  C_Terminal_Residue_ = character(0),
  Construct_Name = character(0),
  Construct_Coding_DNA_Sequence = character(0),
  Construct_Protein_Sequence = character(0),
  Construct_Protein_Sequence_no_tags = character(0),
  Expected_Mass = integer(0),
  Predicted_PI = integer(0),
  Expected_Mass_no_tags = integer(0),
  Predicted_PI_no_tags = integer(0)
)

# Create the Primer data frame
Primer_df <- data.frame(
  Primer_name = character(0),
  Primer_sequence = character(0),
  Primer_length = integer(0),
  Template_DNA_ID = character(0),
  Amplification_position = character(0),
  Tm_whole_primer = integer(0),
  Tm_matching_sequence = integer(0),
  CG_ratio = integer(0))

#positions in the plate
plate_position <- c(
  "A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12",
  "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B09", "B10", "B11", "B12",
  "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12",
  "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12",
  "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09", "E10", "E11", "E12",
  "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12",
  "G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", "G11", "G12",
  "H01", "H02", "H03", "H04", "H05", "H06", "H07", "H08", "H09", "H10", "H11", "H12"
)


# Define a function to translate a DNA codon sequence into a protein sequence
translate_codons <- function(codon_sequence) {
  codon_table <- list(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  
  
  codon_triplets <- strsplit(codon_sequence, "(?<=\\G...)", perl=TRUE)[[1]]
  
  protein_sequence <- sapply(codon_triplets, function(codon) {
    translated <- codon_table[[codon]]
    if (!is.null(translated)) {
      return(translated)
    } else {
      return("?")  # Placeholder for unknown codons
    }
  })
  
  return(paste(protein_sequence, collapse = ""))
}


# Define a function to generate a list of amino acid positions and symbols
generate_amino_acid_list <- function(protein_sequence) {
  codon_table <- list(
    "F" = "Phe", "L" = "Leu", "I" = "Ile", "M" = "Met",
    "V" = "Val", "S" = "Ser", "P" = "Pro", "T" = "Thr",
    "A" = "Ala", "Y" = "Tyr", "*" = "Ter", "H" = "His",
    "Q" = "Gln", "N" = "Asn", "K" = "Lys", "D" = "Asp",
    "E" = "Glu", "C" = "Cys", "W" = "Trp", "R" = "Arg",
    "G" = "Gly"
  )
  
  amino_acid_list <- list()
  for (i in 1:nchar(protein_sequence)) {
    amino_acid <- substr(protein_sequence, i, i)
    if (amino_acid %in% names(codon_table)) {
      amino_acid_name <- codon_table[[amino_acid]]
      amino_acid_list[[i]] <- paste(i, "-", amino_acid_name)
    }
  }
  
  return(amino_acid_list)
}

    #Reverse complementary sequence
    reverse_complement <- function(seq) {
      seq <- chartr("ATCG", "TAGC", seq)
      return(rev(strsplit(seq, "")[[1]]))
    }

  # Function to calculate the Tm of a DNA sequence
  calculate_tm <- function(sequence) {
    tm <- Tm_NN(sequence, Na = 50, dnac1 = 500, nn_table = "DNA_NN4", ambiguous = TRUE, dNTPs = 25)
    return(tm$Tm)
  }
  
  cg_percentage <- function(sequence) {
    sequence <- toupper(sequence)  # Convert to uppercase for consistency
    c_count <- sum(strsplit(sequence, "")[[1]] == "C")
    g_count <- sum(strsplit(sequence, "")[[1]] == "G")
    total_bases <- nchar(sequence)
    cg_percentage <- ((c_count + g_count) / total_bases) * 100
    return(cg_percentage)
  }
  
  # Function to extract forward primer sequence with x nucleotides
  extract_forward_primer <- function(sequence, start_pos, num_nucleotides) {
    start_pos <- max(1, start_pos) # Ensure start_pos is not negative
    end_pos <- min(nchar(sequence), start_pos + num_nucleotides - 1)
    
    forward_primer <- substr(sequence, start_pos, end_pos)
    
    return(forward_primer)
  }
  
  # Function to extract reverse primer sequence with x nucleotides and calculate reverse complement
  extract_reverse_primer <- function(sequence, end_pos, num_nucleotides) {
    end_pos <- max(1, end_pos) # Ensure end_pos is not negative
    start_pos <- max(1, end_pos - num_nucleotides + 1)
    
    reverse_primer <- substr(sequence, start_pos, end_pos)
    
    reverse_complementary_sequence <- paste(reverse_complement(reverse_primer), collapse = "")
    
    return(reverse_complementary_sequence)
  }
  
  #Function for merging overhanging sequences 
  merge_overlaping_sequences <- function(five, three, stream) {
    overlap <- ""  # Initialize the overlap variable
    
    if (grepl(three, five)) {
      seqq_five <- substr(five, 1, ((gregexec(three, five)[[1]][1])+(nchar(three))-1))
      seqq_three <- substr(three, ((gregexec(three, five)[[1]][1])+(nchar(three))-1),nchar(five))
    } else {
      for (s in 1:(min(nchar(five), nchar(three)))) {
        test_overlap <- substr(three, 1, s)  # Chopped sequence preparation
        start <- nchar(five) - (s-1) 
        end <- nchar(five)
        five_test <- substr(five, start, end)
        
        if (test_overlap == five_test) {
          overlap <- test_overlap
        } 
      }
      
      if (nchar(overlap) > 0) {
        start_overlap <- nchar(overlap) + 1
        three_nooverlap <- substr(three, start_overlap, nchar(three))
        seqq_five <- paste(five, three_nooverlap, sep = "")
        seqq_three <- paste(five, three_nooverlap, sep = "")
      } else {
        seqq_five <- paste(five, three, sep = "")  # No overlap found, concatenate as-is
        seqq_three <- paste(five, three, sep = "")  # No overlap found, concatenate as-is
      }
    }
    if(stream == "up"){
    return(seqq_five)
    }
    if(stream == "down"){
      return(seqq_three)
    }
  }
  
  
  
  # Function to prepare extensions
  prep_extension <- function(vector, amplicon, forward_primer,reverse_primer, df) {
    f_extension <- toupper(df[df$Vector.Name == vector, ]$Forward.Primer.5..Extension)
    r_extension <- toupper(df[df$Vector.Name == vector, ]$Reverse.Primer.5..Extension)
    vec5extension<- toupper(df[df$Vector.Name == vector, ]$Vector.5..Sequence)
    vec3extension<- toupper(df[df$Vector.Name == vector, ]$Vector.3..Sequence)
    primerf <- forward_primer
    first_three_f <- substr(primerf, 1, 3)
    if (!(first_three_f == "ATG")) {
      atg_extension <- df[df$Vector.Name == vector, ]$Forward.Primer.5..Conditional.Extension
      f_extension <- paste(f_extension, atg_extension, sep = "")
    }
    primerr <- reverse_primer
    first_three_r <- substr(primerr, 1, 3)
    if (!(first_three_r %in% c("TTA", "TCA", "CTA"))) {
      tca_extension <- df[df$Vector.Name == vector, ]$Reverse.Primer.5..Conditional.Extension
      r_extension <- paste(r_extension, tca_extension, sep = "")
    }
    
    reverse_complementary_reverse_expansion_ <- paste(reverse_complement(r_extension), collapse = "")
    full_amplicon_<- paste(f_extension,amplicon,reverse_complementary_reverse_expansion_, sep = "")
    
    fiveprimesequence <- merge_overlaping_sequences(five = vec5extension, three = f_extension, stream = "up")
    aaa5<- paste("this is five:",fiveprimesequence)

    threeprimesequence <- merge_overlaping_sequences(five = reverse_complementary_reverse_expansion_ ,three = vec3extension, stream = "down")
    aaa3<- paste("this is three:",threeprimesequence)

   
    wholeCDS <- paste(fiveprimesequence,amplicon,threeprimesequence,sep = "")
    
    startCDS_translation <- (df[df$Vector.Name == vector, ]$Vector.CDS.Start.Index)
    
    translated_amplicon <- translate_codons(substr(wholeCDS,startCDS_translation,nchar(wholeCDS)))
    
    final_position <- nchar(wholeCDS)
    
    if(grepl("\\*", translated_amplicon)){
    final_position <- (((gregexpr("\\*", translated_amplicon)[[1]][1])*3)-1)
    }
    
    startCDS_translation_termination <- substr(wholeCDS,startCDS_translation,(startCDS_translation+ final_position))
    
    final_translated_protein <- translate_codons(startCDS_translation_termination)
    
    final_translated_protein_mw <- mw(final_translated_protein)
    
    final_translated_protein_PI <- pI(final_translated_protein)
    
    n_therminaltag_position <-(df[df$Vector.Name == vector, ]$Vector.N.Terminal.Tag.Length)

    c_therminaltag_position <-(df[df$Vector.Name == vector, ]$Vector.C.Terminal.Tag.Length)

    final_translated_protein_notag <- substr(final_translated_protein,n_therminaltag_position,(nchar(final_translated_protein)-c_therminaltag_position))

    final_translated_protein_notag_mw <- mw(final_translated_protein_notag)
    
    final_translated_protein_notag_pi <- pI(final_translated_protein_notag)
    
    
    
    extended_packge <- list(full_f_primer = paste(f_extension,forward_primer,sep = ""), extension_f=f_extension,
                             full_r_primer = paste(r_extension,reverse_primer,sep = ""), extension_r=r_extension,
                             full_amplicon = full_amplicon_, product_length = nchar(full_amplicon_),
                            constructCDS = startCDS_translation_termination, protein_with_tag = final_translated_protein,
                            protein_without_tag = final_translated_protein_notag, mw_with_tag = final_translated_protein_mw, 
                            mw_without_tag = final_translated_protein_notag_mw, pi_with_tag = final_translated_protein_PI,
                            pi_without_tag = final_translated_protein_notag_pi, f_tm = calculate_tm(forward_primer),
                            r_tm = calculate_tm(reverse_primer)
                            )
    
    return(extended_packge) 
  }
  

  
  
  ui <- fluidPage(
    tags$head(
      tags$style(
        HTML(".custom-title-panel { height: 100px; background-color: #ffffff; }")
      )
    ),
    theme = bs_theme(bootswatch = "journal"),
    titlePanel(
      div(class = "custom-title-panel",
          h1("LIC construct and primer design tool"),
          div(
            style = "position: absolute; top: 20px; right: 0;",
            img(src = "CQMEDlogo.png", width = 300, height = 70)
          )
      )
    ),
      mainPanel(width = 12,
        tabsetPanel(
          tabPanel("Step-1: Template DNA info",
                   br(),
                   fluidRow(
                     column(3, textInput("author", "Project's author")),
                     column(3, textInput("plate_id", "Plate")),
                     column(3, textInput("targetName", "Target Name:")),
                     column(3, textInput("templateId", "Template DNA ID:"))
                   ),
                   textAreaInput("templateSequence", "Template DNA Sequence: (Paste DNA sequence here)",
                                 value = "ATGACGCGTATATCGGCTAGCGCGTATATATTAGCGCGAGTAGCGGCGATTCGCGGATTCGGATGCGATGCTAGCGCGATCGATCG", 
                                 width = "100%"),
                   fluidRow(
                     column(4, textInput("proteinSequence", "Reference Protein Sequence:")),
                     column(4, textInput("notes", "Notes (optional):")),
                     column(4, numericInput("num", label = "Number of protein variables", value = 1))
                   ),
                   br(),
                   h2("Step-2: Select Vectors for this project"),
                   uiOutput("vector_choices")
          ),
          tabPanel("Step-3: Setting up the constructs",
                   verbatimTextOutput(outputId = "protein"),
                   uiOutput("interaction_slider"),
                   uiOutput("selectedValues"),
                   uiOutput("vecchoices"),
                   br(),
                   br(),
                   fluidRow(
                      column (4,actionButton("makeprimer", "Step-4: Make Primer")),
                      column (4,downloadButton("downloadExcel", "Step-5: Download Construct File")),
                      column (4, downloadButton("downloadScarab", "Step-6: Download Scarab File"))
                   ),
                   br(),
                   br(),
                   br(),
                   div(
                     style = "text-align: center;",  # Center the content
                     h2("Download oligo forms of the desired company")
                   ),
                   br(),
                   fluidRow(
                   column (3, offset = 1, downloadButton("downloadThermo",
                                  label = div(
                                    img(src = "thermo_logo.png", width = 120, height = 30),  # Replace with your logo's path
                                    ""
                                  )
                                )),
                   column (3, offset = 1, downloadButton("downloadExxtend",
                                  label = div(
                                    img(src = "exxtend_logo.png", width = 120, height = 30),  # Replace with your logo's path
                                    ""
                                  )
                   )),
                   column (3, offset = 1, downloadButton("downloadCSVidt",
                                  label = div(
                                    img(src = "idt_logo.png", width = 120, height = 30),  # Replace with your logo's path
                                    ""
                                  )
                                  )
                   ))
          ),
          tabPanel("Construct Design",
                   dataTableOutput("data_table_construct")
          ),
          tabPanel("Primer",
                   dataTableOutput("data_table_primer")
          ),
          tabPanel("Vector Spreadsheet",
                   tags$iframe(style="height:800px; width:100%", src="Vector_CQMED.pdf")
          )
        )
      )
    )
  
  
  vector_PCR_product_name <- function(target_name, n) {
    prefix <- target_name
    vector <- paste0(prefix, "-", sprintf("ab%03d", 1:n))
    return(vector)
  }
  
  vector_construct_name <- function(target_name, n) {
    prefix <- target_name
    vector <- paste0(prefix, "-", sprintf("cb%03d", 1:n))
    return(vector)
  }
  
  vector_primerf_name <- function(target_name, n) {
    prefix <- target_name
    vector <- paste0(prefix, "-", sprintf("fb%03d", 1:n))
    return(vector)
  }
  
  vector_primerr_name <- function(target_name, n) {
    prefix <- target_name
    vector <- paste0(prefix, "-", sprintf("rb%03d", 1:n))
    return(vector)
  }
  
  get_primer_name <- function(primer_sequence) {
    primer_sequence_clean <- str_trim(primer_sequence)
    primer_name_return <- Primer_df[Primer_df$Primer_sequence == primer_sequence_clean,]$Primer_name
    return(primer_name_return)
  }
  
#--------------------------------------------------------------------------------------------------------------------------------
# Define server
server <- function(input, output) {
  
  shinyDirChoose(input, "outputFolder", roots = c(wd = getwd()))
  
  output$CQMEDlogo <- renderImage({
    image_path <- "CQMEDlogo.png"  # Relative path within the www directory
    list(src = image_path, alt = "CQMEDlogo")
    }, deleteFile = FALSE)
  
  # Create a reactive expression to collect the data
  collected_data <- reactive({
    data <- data.frame(
      Author = input$author,
      Target_Name = input$targetName,
      Template_DNA_ID = input$templateId,
      Template_DNA_Sequence = input$templateSequence,
      Reference_Protein_Sequence = input$proteinSequence,
      Notes = input$notes
    )
    return(data)
  })
  
  #Create a reactive value to store the current state of your data
  rv_construct_design <- reactiveVal(Construct_Design_df)
  
  rv_primer <- reactiveVal(Primer_df)
  
  
  # Process the template sequence
  processed_template_sequence <- reactive({
    toupper(input$templateSequence)  # Convert to upper case
  })
  
  # Generate a list of amino acid positions and three-letter codes
  protein_aminoacids <- reactive({
    translate_codons(processed_template_sequence())
  })
  
  amino_acid_list_selection <- reactive({
    generate_amino_acid_list(protein_aminoacids())
  })
  
  tamanho <- reactive({
    length(amino_acid_list_selection())
  })
  
  output$protein <- renderText(paste("Protein sequence", protein_aminoacids(), sep="\n"))
  
  output$interaction_slider <- renderUI({
    num_sliders <- input$num
    slider_inputs <- lapply(1:num_sliders, function(i) {
      sliderInput(
        inputId = paste0("slider", i),  # Unique input ID for each slider
        label = paste("Select Range", i, "for PCR cloning:"),
        step = 1,
        min = 1,
        max = tamanho(),
        value = c(1, 1),
        dragRange = TRUE,
        width = "95%"
      )
    })
    do.call(tagList, slider_inputs)
  })
  
  # Display selected values of all sliders
  output$selectedValues <- renderUI({
    num_sliders <- input$num
    slider_value_text <- lapply(1:num_sliders, function(i) {
      selected_range <- input[[paste0("slider", i)]]
      paste("For construct", i, "you've selected the aa range from: ",
            amino_acid_list_selection()[[selected_range[1]]],
            " to ",
            amino_acid_list_selection()[[selected_range[2]]], sep = " ")
    })
    do.call(tagList, lapply(slider_value_text, function(sentence) {
      HTML(paste(sentence, "<br>"))  # Wrap each sentence in HTML <br> tag
    }))
  })

  # Load the dataframe from Template_vectors_SGC.csv
  Template_vectors_SGC <- reactive({
    read.csv("Template_vectors_SGC.csv", stringsAsFactors = FALSE)
  })  
    
  output$vector_choices <- renderUI({
    vector_options <- unique(Template_vectors_SGC()$Vector.Name)
    checkboxGroupInput("selected_vectors", "",
                       choices = vector_options,
                       inline = T,
                       width = "100%")
  })
  
  output$vecchoices <- renderUI({
    num_sliders <- input$num
    vector_choosen <- input$selected_vectors
    vector_construct <- lapply(1:num_sliders, function(i) {
      checkboxGroupInput(paste0("selected_vectors_construct",i), paste0("Select Vectors for construct",i),
                         choices = vector_choosen,
                         inline = T,
                         width = "100%"
      )
    })
    do.call(tagList, vector_construct)
  })
  
  # Create a reactive expression for the sequences
  extracted_sequences <- eventReactive(input$makeprimer, {
    num_sliders <- input$num
    matchprimer <- lapply(1:num_sliders, function(i) {
      selected_range <- input[[paste0("slider", i)]]
      seq <- processed_template_sequence()
      start_pos <- ((selected_range[1] - 1) * 3) + 1
      num_forward <- 18
      final_pos <- (selected_range[2] * 3)
      num_backward <- 18
      tm_threshold <- 45  # Adjust Tm threshold as needed
      amplicon_minor <- substr(seq,start_pos,final_pos) #amplicon without the primer extensions
      
      # In your server function, retrieve the selected vectors
      selected_checkboxes <- input[[paste0("selected_vectors_construct", i)]]
      num_of_contructs_for_this_vector <- length(selected_checkboxes)
      
      #Making matching forward primer
      
      forward_primer <- extract_forward_primer(seq, start_pos, num_forward)
      
      while (calculate_tm(forward_primer) <= tm_threshold && nchar(forward_primer) < nchar(seq)) {
        num_forward <- num_forward +1
        forward_primer <- extract_forward_primer(seq, start_pos,num_forward)
      }
      
      #Making matching reverse primer
      
      reverse_primer<- extract_reverse_primer(seq,final_pos,num_backward)
      
      while (calculate_tm(reverse_primer) <= tm_threshold && nchar(reverse_primer) < nchar(seq)) {
        num_backward <- num_backward - 1
        reverse_primer<- extract_reverse_primer(seq,final_pos,num_backward)
     }
      
      #prep_extension <- function(vector, amplicon, forward_primer,reverse_primer, df)
      # Apply prep_extension function to each vector in the list ok
      extended_list <- lapply(selected_checkboxes, prep_extension,amplicon = amplicon_minor, forward_primer=forward_primer, reverse_primer = reverse_primer, df = Template_vectors_SGC())
      
      extended <- bind_rows(extended_list)
      
      
      # Create the Construct Design data frame
        Construct_Design_feeder_df <- data.frame(
        Target_Name = rep(input$targetName, num_of_contructs_for_this_vector),
        Vector_Name = unlist(selected_checkboxes),
        Template_DNA_ID = rep(input$templateId, num_of_contructs_for_this_vector),
        N_Terminal_Residue = rep(paste(amino_acid_list_selection()[[selected_range[1]]]), num_of_contructs_for_this_vector),
        C_Terminal_Residue = rep(paste(amino_acid_list_selection()[[selected_range[2]]]), num_of_contructs_for_this_vector),
        N_Terminal_Extension_DNA =extended$extension_f,
        C_Terminal_Extension_DNA = extended$extension_r,
        Forward_Primer_Name = extended$full_f_primer,
        Reverse_Primer_Name = extended$full_r_primer,
        #PCR_Product_Name 
        PCR_Product_Sequence = extended$full_amplicon,
        PCR_Product_Length = extended$product_length,
        N_Terminal_Residue_ = rep(paste(amino_acid_list_selection()[[selected_range[1]]]), num_of_contructs_for_this_vector),
        C_Terminal_Residue_ = rep(paste(amino_acid_list_selection()[[selected_range[2]]]), num_of_contructs_for_this_vector),
        #Construct_Name 
        Construct_Coding_DNA_Sequence = extended$constructCDS,
        Construct_Protein_Sequence = extended$protein_with_tag,
        Construct_Protein_Sequence_no_tags = extended$protein_without_tag,
        Expected_Mass = extended$mw_with_tag,
        Predicted_PI = extended$pi_with_tag,
        Expected_Mass_no_tags = extended$mw_without_tag,
        Predicted_PI_no_tags = extended$pi_without_tag,
        Matching_tm_f = extended$f_tm,
        Matching_tm_r = extended$r_tm
      )
      
      Construct_Design_feeder_df
    })
    
    
    Construct_Design_final_df <- bind_rows(matchprimer)  # Combine the list of data frames
    
    # Order Construct_Design_final_df by the "Vector_Name" column
    Construct_Design_final_df <- Construct_Design_final_df[order(Construct_Design_final_df$Vector_Name), ]
    
    Construct_Design_final_df[1:nrow(Construct_Design_final_df),"Position"] <- plate_position[1:nrow(Construct_Design_final_df)]
    
        Construct_Design_final_df[1:nrow(Construct_Design_final_df),"PCR_Product_Name"] <- vector_PCR_product_name(input$targetName, nrow(Construct_Design_final_df))
      
    Construct_Design_final_df[1:nrow(Construct_Design_final_df),"Construct_Name"] <- vector_construct_name(input$targetName, nrow(Construct_Design_final_df))
                                              
    
    return(Construct_Design_final_df)
    
  })
  
  #
  
  # Update the DataTable when extracted_sequences() changes
  observeEvent(extracted_sequences(), {
    rv_construct_design(Construct_Design_df)
    
    column_names_extracted_data <- c(
      "Position",
      "Target_Name",
      "Vector_Name",
      "Template_DNA_ID",
      "N_Terminal_Residue",
      "C_Terminal_Residue",
      "N_Terminal_Extension_DNA",
      "C_Terminal_Extension_DNA",
      "Forward_Primer_Name",
      "Reverse_Primer_Name",
      "PCR_Product_Name",
      "PCR_Product_Sequence",
      "PCR_Product_Length",
      "N_Terminal_Residue_",
      "C_Terminal_Residue_",
      "Construct_Name",
      "Construct_Coding_DNA_Sequence",
      "Construct_Protein_Sequence",
      "Construct_Protein_Sequence_no_tags",
      "Expected_Mass",
      "Predicted_PI",
      "Expected_Mass_no_tags",
      "Predicted_PI_no_tags"
    )
    
    extracted_data <- extracted_sequences()  # Get the extracted sequences data
    
    updated_Construct_Design_df <- bind_rows(rv_construct_design(),extracted_data[,column_names_extracted_data])
    
    
    Primer_df <- data.frame(
      Primer_name = c(rep("F",nrow(updated_Construct_Design_df)),rep("R",nrow(updated_Construct_Design_df))),
      Primer_sequence = c(updated_Construct_Design_df$Forward_Primer_Name,updated_Construct_Design_df$Reverse_Primer_Name),
      Primer_length = nchar(c(updated_Construct_Design_df$Forward_Primer_Name,
                        updated_Construct_Design_df$Reverse_Primer_Name)),
      Template_DNA_ID = rep(updated_Construct_Design_df$Template_DNA_ID,2),
      Amplification_position = c(updated_Construct_Design_df$N_Terminal_Residue,updated_Construct_Design_df$C_Terminal_Residue),
      Tm_whole_primer = calculate_tm(c(updated_Construct_Design_df$Forward_Primer_Name,
                                        updated_Construct_Design_df$Reverse_Primer_Name)),
      Tm_matching_sequence = c(extracted_data$Matching_tm_f,extracted_data$Matching_tm_r),
      CG_ratio = cg_percentage(c(updated_Construct_Design_df$Forward_Primer_Name,
                                 updated_Construct_Design_df$Reverse_Primer_Name))
      )
    
    Primer_df <- Primer_df[!duplicated(Primer_df$Primer_sequence), ]
    
    Primer_df[,"Primer_name"] <- c(vector_primerf_name(input$targetName,sum(str_count(Primer_df$Primer_name, "F"))),
                                   vector_primerr_name(input$targetName,sum(str_count(Primer_df$Primer_name, "R"))))
    
    # Update the Forward and Reverse primer names in updated_Construct_Design_df
    
    for (g in 1:nrow(updated_Construct_Design_df)){
    updated_Construct_Design_df[g,"Forward_Primer_Name"] <- Primer_df[Primer_df$Primer_sequence == updated_Construct_Design_df[g,"Forward_Primer_Name"],]$Primer_name
    
    updated_Construct_Design_df[g,"Reverse_Primer_Name"] <- Primer_df[Primer_df$Primer_sequence == updated_Construct_Design_df[g,"Reverse_Primer_Name"],]$Primer_name
    }
    
    rv_construct_design(updated_Construct_Design_df)
    
    rv_primer(Primer_df)

  })
  
  # Render the DataTable using the reactive value
  output$data_table_construct <- renderDataTable({
    datatable(rv_construct_design()
    )
  })
  
  output$data_table_primer <- renderDataTable({
      datatable(rv_primer()
      )
  })
  
  output$downloadExcel <- downloadHandler(
    filename = function() {
     paste(input$author, input$targetName, input$templateId, Sys.Date(),
                                          "Template_Contruct_Design.xlsx", sep = "_")
    },
    content = function(file) {
      wb <- createWorkbook()
      addWorksheet(wb, "Template DNA")
      writeData(wb, sheet = 1, x = collected_data())
      
      addWorksheet(wb, "Construc Design")
      writeData(wb, sheet = 2, x = rv_construct_design())
      
      addWorksheet(wb, "Primer")
      writeData(wb, sheet = 3, x = rv_primer())
      
      saveWorkbook(wb, file)
    }
  )
  
  output$downloadScarab <- downloadHandler(
    filename = function() {
      paste(input$plate_id, input$author, Sys.Date(),
            "Export to Scarab.xlsx", sep = "_")
    },
    content = function(file) {
      wb <- createWorkbook()
      
      # Create a data frame
      alleles_df <- data.frame(
        Allele_ID = rv_construct_design()$PCR_Product_Name,
        Plate = rep(input$plate_id, nrow(rv_construct_design())),
        Entry_Clone_ID = rv_construct_design()$Template_DNA_ID,
        Forward_Primer_ID = rv_construct_design()$Forward_Primer_Name,
        Reverse_Primer_ID = rv_construct_design()$Reverse_Primer_Name,
        DNA_Sequence = rep("", nrow(rv_construct_design())),
        Protein_Sequence = rep("", nrow(rv_construct_design())),
        Allele_Status = rep("In Progress", nrow(rv_construct_design())),
        Location = rep("", nrow(rv_construct_design())),
        Comments = rep("", nrow(rv_construct_design())),
        ELN_Experiment = rep("", nrow(rv_construct_design())),
        Date_Record = rep("", nrow(rv_construct_design())),
        Person = rep("", nrow(rv_construct_design())),
        Plate_Well = rv_construct_design()$Position,
        DNA_Sequence_Length = rep("", nrow(rv_construct_design()))
      )
      
      for (e in 1:nrow(rv_construct_design())) {
        alleles_df[e, "Comments"] <- paste("Region:", rv_construct_design()[e, "N_Terminal_Residue"], "to", 
                                           rv_construct_design()[e, "C_Terminal_Residue"], sep = " ")
      }
      
      
      
      addWorksheet(wb, "Alleles")
      writeData(wb, sheet = "Alleles", x = alleles_df)
      
      # Create a data frame
      construct_df <- data.frame(
        Construct_ID = rv_construct_design()$Construct_Name,
        Construct_Plate_ID = rep(input$plate_id, nrow(rv_construct_design())),
        Construct_Plate_Well = rv_construct_design()$Position,
        Vector_ID = rv_construct_design()$Vector_Name,
        Allele_ID = rv_construct_design()$PCR_Product_Name,
        Construct_Status = rep("", nrow(rv_construct_design())),
        Construct_Protein_Sequence = rep("", nrow(rv_construct_design())),
        Expected_Mass = rep("", nrow(rv_construct_design())),
        Restriction_Enzyme_1_ID = rep("", nrow(rv_construct_design())),
        Restriction_Enzyme_2_ID = rep("", nrow(rv_construct_design())),
        Construct_Protein_Sequence_no_tag = rep("", nrow(rv_construct_design())),
        Expected_Mass_no_tag = rep("", nrow(rv_construct_design())),
        Construct_DNA_Sequence = rep("", nrow(rv_construct_design())),
        Construct_Location = rep("", nrow(rv_construct_design())),
        Construct_ELN_Reference = rep("", nrow(rv_construct_design())),
        Construct_Comments = rep("", nrow(rv_construct_design()))
      )
      
      for (e in 1:nrow(rv_construct_design())) {
        construct_df[e, "Construct_Comments"] <- paste("Region:", rv_construct_design()[e, "N_Terminal_Residue"], "to", 
                                           rv_construct_design()[e, "C_Terminal_Residue"], sep = " ")
      }
      
      
      addWorksheet(wb, "Constructs")
      writeData(wb, sheet = "Constructs", x = construct_df)
      
      
      saveWorkbook(wb, file)
    }
  )

  output$downloadThermo <- downloadHandler(
    filename = function() {
      paste(input$author, input$targetName, input$templateId, Sys.Date(),
            "Formulário de Oligos de DNA Thermo.xlsx", sep = "_")
    },
    content = function(file) {
      # Load the existing Excel file
      existing_wb <- loadWorkbook("Formulário de Oligos de DNA Thermo.xlsx")
      
      thermo_order_primers_df <- data.frame(
        "Nome do pesquisador" = rep(input$author, nrow(rv_primer())),
        "Nome do oligonucleotídeo" = rv_primer()$Primer_name,
        "Sequência" = rv_primer()$Primer_sequence
      )
      
      # Determine the row where you want to start adding data (e.g., row 21)
      start_row <- 29
      
      # Write the data frame starting from the specified row
      writeData(wb = existing_wb, sheet = "Formulário", x = thermo_order_primers_df, startCol = 1, startRow = start_row, colNames = FALSE)
      
      # Save the modified Excel file back to the same location
      saveWorkbook(existing_wb, file)
    }
  )
  
  output$downloadExxtend <- downloadHandler(
    filename = function() {
      paste(input$author, input$targetName, input$templateId, Sys.Date(),
            "PLANILHA EXXTEND.xlsx", sep = "_")
    },
    content = function(file) {
      # Load the existing Excel file
      existing_wb <- loadWorkbook("PLANILHA EXXTEND.xlsx")
      
      thermo_order_primers_df <- data.frame(
        "Nome do pesquisador" = rep(input$author, nrow(rv_primer())),
        "Nome do oligonucleotídeo" = rv_primer()$Primer_name,
        "Sequência" = rv_primer()$Primer_sequence
      )
      
      # Determine the row where you want to start adding data (e.g., row 21)
      start_row <- 13
      
      # Write the data frame starting from the specified row
      writeData(wb = existing_wb, sheet = "Pedido e Cotação", x = thermo_order_primers_df, startCol = 3, startRow = start_row, colNames = FALSE)
      
      # Save the modified Excel file back to the same location
      saveWorkbook(existing_wb, file)
    }
  )
  
  # Define a download handler
  output$downloadCSVidt <- downloadHandler(
    filename = function() {
      paste(input$author, input$targetName, input$templateId, Sys.Date(),
            "_IDT_Primers.csv" , sep = "") # Set the file name
    },
    content = function(file) {
      idt_order_primers_df <- data.frame(
        "C1" = rv_primer()$Primer_name,
        "C2" = rv_primer()$Primer_sequence,
        "C3" = rep("25nm", nrow(rv_primer())),
        "C4" = rep("STD", nrow(rv_primer()))
      )
      write.csv(idt_order_primers_df, file, row.names = FALSE, col.names = FALSE, quote = FALSE)  # Exclude column names
    }
  )
  
  
}

  

# Run the app
shinyApp(ui, server)
