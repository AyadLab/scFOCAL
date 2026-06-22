library(shiny)
library(shinyFiles) # <--- ADD THIS LINE HERE
library(Seurat)
library(ggplot2)
library(tibble)
library(cowplot)
library(viridis)
library(dplyr)
library(ggsci)
library(ggrepel)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(reshape2)
library(Hmisc)
library(corrplot)
library(pheatmap)
library(grid)
library(MAST)
library(shinydashboard)
library(shinythemes)
library(scales)
library(ggforce)
library(EnhancedVolcano)
library(DT)

options(shiny.maxRequestSize = 30000*1024^2) # increase limit to 15gb?

################################################################################
#
################################################################################

# # LINCS Response Signature Data
# LINCS.ResponseSigs <- read.delim(
#   file = "matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt", header = T)
# row.names(LINCS.ResponseSigs) <- LINCS.ResponseSigs$Genes
# LINCS.ResponseSigs <- na.omit(LINCS.ResponseSigs)
# newNames <- gsub("-", ".", rownames(LINCS.ResponseSigs))
# rownames(LINCS.ResponseSigs) <- newNames
# L1000_genes <- colnames(LINCS.ResponseSigs)
# L1000_compounds <- rownames(LINCS.ResponseSigs)

ui <- fluidPage(
  tags$script(src = "https://kit.fontawesome.com/070e476711.js"),
  img(src = "scFOCAL.png", width = "30%"),
  hr(),
  tags$head(tags$style(HTML(
    ".nav.nav-pills.nav-stacked > .active > a, .nav.nav-pills.nav-stacked > .active > a:hover {
      background-color: #000000;
    }
    .well {
      min-height: 20px;
      max-width: 5000px;
      padding: 19px;
      margin-bottom: 20px;
      background-color: #ffffff;
      border: 1px solid #ffffff;
      border-radius: 4px;
      -webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
      box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
      font-family: 'sans-serif', Arial Rounded MT Bold;
    }
    "
  ))),
  
  navlistPanel(well = T, fluid = T, widths = c(2, 8),
               
               # ------------------------------------------------------------------
               # OVERVIEW TAB
               # ------------------------------------------------------------------
               tabPanel(tags$div(
                 tags$i(class = "fa-sharp fa-solid fa-desktop"),
                 tags$span("- Overview"),
                 tags$style(type = "text/css", "li a{color:#000000;}")
               ), img(src = "scSynergySeq Diagram (2).png", width = "100%")
               ),
               
               # ------------------------------------------------------------------
               # RUN scFOCAL DEV TAB
               # ------------------------------------------------------------------
               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-magnifying-glass-chart"),
                 tags$span("- Run scFOCAL - dev")
               ),
               
               tabsetPanel(
                 
                 # --- 1. DATA UPLOAD ---
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-upload"),
                   tags$span("1. Data Upload"),
                   tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 br(),
                 h4("Step 1: Upload Seurat object RDS file"),
                 hr(),
                 br(),
                 br(),
                 # radioButtons("upload_source", "Select Data Source:",
                 #              choices = c("Upload from Computer" = "local",
                 #                          "Browse Server Files" = "server")),
                 # # splitLayout(
                 #   # Local Upload
                 #   conditionalPanel(
                 #     condition = "input.upload_source == 'local'",
                 #     fileInput(inputId = "seurobjRDS",
                 #               label = NULL,
                 #               buttonLabel = "Browse...",
                 #               placeholder = "No file selected",
                 #               width = NULL,
                 #               multiple = F,
                 #               accept = c(".RDS", ".Rds", ".rds"))
                 #   ),
                 #   # Server Upload
                 #   conditionalPanel(
                 #     condition = "input.upload_source == 'server'",
                 #     shinyFilesButton("file_server", "Browse Server Files", "Please select an RDS file", multiple = FALSE),
                 #     tags$div(tags$b("Selected file:"), textOutput("selected_server_file"), style = "margin-bottom: 15px; margin-top: 5px;")
                 #   )
                 # ) 
                 # ), # <-- FIXED: Closed 1. Data Upload Tab
                 splitLayout(
                   cellWidths = c("40%", "60%"),
                   
                   # LEFT SIDE: Source Selector
                   radioButtons(inputId = "upload_source", 
                                label = "Select Data Source:",
                                choices = c("Upload from Computer" = "local",
                                            "Browse Server Files" = "server")),
                   
                   # RIGHT SIDE: Dynamic Browse Buttons & Success Messages
                   tags$div(
                     
                     # --- LOCAL UPLOAD UI ---
                     conditionalPanel(
                       condition = "input.upload_source == 'local'",
                       fileInput(inputId = "seurobjRDS",
                                 label = "Choose Local File:",
                                 buttonLabel = "Browse...",
                                 placeholder = "No file selected",
                                 width = "100%",
                                 multiple = FALSE,
                                 accept = c(".RDS", ".Rds", ".rds")),
                       
                       # Local Success Message (Shows only when file input is not null)
                       conditionalPanel(
                         condition = "input.seurobjRDS != null",
                         tags$h5(icon("check-circle"), " File successfully uploaded!", style = "color: #28a745; font-weight: bold;")
                       )
                     ),
                     
                     # --- SERVER UPLOAD UI ---
                     conditionalPanel(
                       condition = "input.upload_source == 'server'",
                       tags$label("Choose Server File:"),
                       tags$br(),
                       
                       # We style this button to visually match the local "Browse..." button
                       shinyFilesButton(id = "file_server", 
                                        label = "Browse...", 
                                        title = "Please select an RDS file", 
                                        multiple = FALSE, 
                                        icon = icon("folder-open")),
                       tags$br(), tags$br(),
                       
                       # Server Success Message (shinyFiles button changes from an integer to an object when a file is picked)
                       conditionalPanel(
                         condition = "typeof input.file_server === 'object' && input.file_server !== null",
                         tags$h5(icon("check-circle"), " File successfully selected:", style = "color: #28a745; font-weight: bold;"),
                         tags$div(textOutput("selected_server_file"), style = "color: #28a745; margin-top: 5px;")
                       )
                     )
                   )
                 )
                 ),
                 # --- 2. PRE-PROCESSING ---
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-gears"),
                   tags$span("2. Pre-processing"),
                   tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 
                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    br(),
                                    p("Please first upload a seurat object as an .RDS file"),
                                    hr(),
                                    br(),
                                    p(em("Due to the size of some datasets, please allow a few minutes after upload completes for the file to be loaded into the working environment"))
                                  )
                 ),
                 
                 conditionalPanel(condition = 'output.seuratLoaded',
                                  h4("Success! Seurat object loaded.", style = "color: green;"),
                                  conditionalPanel(condition = 'output.seurat_obj_v5',
                                                   h4("Version 5 Seurat Object Detected"))
                 ),
                 
                 conditionalPanel(condition = "output.seuratLoaded",
                                  splitLayout(
                                    wellPanel(
                                      p(strong("Data Pre-processing")),
                                      hr(),
                                      br(),
                                      uiOutput(outputId = "reductionUseRDS"),
                                      uiOutput(outputId = "groupByRDS"),
                                      checkboxInput(inputId = "noNormalCellIdentsCheckbox",
                                                    label = "My data has no normal cell reference", value = FALSE),
                                      conditionalPanel(condition = "!input.noNormalCellIdentsCheckbox",
                                                       uiOutput(outputId = "normalCellIdentsRDS")
                                      ),
                                      uiOutput('cancerCellIdentsRDS'),
                                      uiOutput('subject_replicate_ident')
                                    ), 
                                    tabsetPanel(type = "tabs",
                                                tabPanel("Identities",
                                                         plotOutput(outputId = "dimplotRDS") 
                                                ),
                                                tabPanel("Subject/Replicate",
                                                         plotOutput(outputId = "dimplotRDS_subject")
                                                )
                                    )
                                  ),
                                  splitLayout(
                                    wellPanel(
                                      plotOutput(outputId = "controlHighlight")
                                    ),
                                    wellPanel(
                                      plotOutput(outputId = "diseasedHighlight")
                                    )
                                  )
                 )
                 ), # <-- Closes 2. Pre-processing Tab
                 
                 # --- 3. DISEASE SIGNATURES ---
                 tabPanel(tags$div(
                   tags$i(class = "fa-solid fa-signature"),
                   tags$span("3. Disease Signatures"),
                   tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    br(),
                                    p("Please first upload a seurat object as an .RDS file"),
                                    hr(),
                                    br(),
                                    p(em("Due to the size of some datasets, please allow a few minutes after upload completes for the file to be loaded into the working environment"))
                                  )
                 ),
                 conditionalPanel(condition = "output.seuratLoaded",
                                  wellPanel(
                                    img(src = "SlicedDiseaseSignatureCalc.png", width = "100%")
                                  ),
                                  conditionalPanel(condition = "output.controlsNotSelected",
                                                   h4("No test cell populations selected"),
                                                   br(),
                                                   p(strong("Please return to the pre-processing tab to select test cell populations")),
                                                   hr()
                                  ),
                                  conditionalPanel(condition = "output.controlsSelected",
                                                   wellPanel(
                                                     sliderInput(inputId = "cellSubsetMax_1",
                                                                 label = "Select maximum number of cells per condition to test",
                                                                 min = 100, max = 40000, step = 100, value = 500
                                                     ),
                                                     radioButtons(inputId = "diseaseSignaturesByReplicate",
                                                                  label = "Calculate overall or subject/replicate disease signatures?",
                                                                  choices = list("Overall", "By Subject/Replicate"),
                                                                  selected = c("Overall")),
                                                     actionButton(inputId = "CalcDiseaseSig",
                                                                  label = "Calculate Disease Signatures",
                                                                  icon = icon("sync"))
                                                   )
                                  )
                 ),
                 conditionalPanel(condition = "input.CalcDiseaseSig >= 1", 
                                  tabsetPanel(
                                    tabPanel(title = "Heatmap",
                                             wellPanel(
                                               plotOutput(outputId = "diseaseSigRDS"),
                                               downloadButton("diseaseSigDownload", label = "Download signature .csv"),
                                               downloadButton("diseaseSigHeatmapDownload", label = "Download Plot .pdf")
                                             )
                                    ),
                                    tabPanel(title = "Table",
                                             wellPanel(DT::dataTableOutput(outputId = "diseaseSigTable"))
                                    ),
                                    tabPanel(title = "Compare",
                                             "Under construction!"),
                                    tabPanel(title = "Reversal",
                                             wellPanel(
                                               radioButtons(inputId = "whichDiseaseSignatures",
                                                            choices = c("Whole dataset disease signatures",
                                                                        "By subject / replicate disease signatures"),
                                                            inline = F, selected = "Whole dataset disease signatures",
                                                            label = "Select which disease signature set to use."),
                                               conditionalPanel(condition = 'input.whichDiseaseSignatures == "By subject / replicate disease signatures"',
                                                                conditionalPanel(condition = 'input.calcSlicedDiseaseSigs >= 1',
                                                                                 actionButton(inputId = "reverseSlicedDiseaseSigs",
                                                                                              label = "Score compound reversal"),
                                                                                 hr()
                                                                ),
                                                                conditionalPanel(condition = 'input.calcSlicedDiseaseSigs < 1',
                                                                                 "Please calculate sliced disease signatures first."
                                                                )
                                               ),
                                               conditionalPanel(condition = 'input.whichDiseaseSignatures == "Whole dataset disease signatures"',
                                                                conditionalPanel(condition = 'input.CalcDiseaseSig >= 1',
                                                                                 actionButton(inputId = "reverseDiseaseSigs",
                                                                                              label = "Score compound reversal")
                                                                ),
                                                                conditionalPanel(condition = 'input.CalcDiseaseSig < 1',
                                                                                 "Please calculate disease signatures first."
                                                                )
                                               ),
                                               conditionalPanel(condition = 'input.reverseDiseaseSigs >= 1',
                                                                tabsetPanel(selected = "QC",
                                                                            tabPanel("QC",
                                                                                     plotOutput(outputId = "DiseaseSigReversalQC")
                                                                            ),
                                                                            tabPanel("Heatmap",
                                                                                     sliderInput(inputId = "varNum",
                                                                                                 label = "Select number of top most variable compounds to plot",
                                                                                                 min = 1, max = 2000, step = 1,
                                                                                                 value = 100),
                                                                                     plotOutput(outputId = "DiseaseSigReversalHeatmap")
                                                                            ),
                                                                            tabPanel("Barplots",
                                                                                     selectizeInput(inputId = "wcompoundToBarplot",
                                                                                                    choices = L1000_compounds,
                                                                                                    selected = "alisertib",
                                                                                                    label = "Select small molecule to plot"
                                                                                     ),
                                                                                     plotOutput(outputId = "DiseaseSigReversalBarplot")
                                                                            )
                                                                ),
                                                                downloadButton(outputId = "reversedDiseaseSigDL", label = "Download reversal score matrix")
                                               ),
                                               conditionalPanel(condition = 'input.reverseSlicedDiseaseSigs >= 1',
                                                                tabsetPanel(selected = "QC",
                                                                            tabPanel("QC",
                                                                                     plotOutput(outputId = "slicedDiseaseSigReversalQC")
                                                                            ),
                                                                            tabPanel("Heatmap",
                                                                                     sliderInput(inputId = "varNum",
                                                                                                 label = "Select number of top most variable compounds to plot",
                                                                                                 min = 1, max = 2000, step = 1,
                                                                                                 value = 100),
                                                                                     plotOutput(outputId = "slicedDiseaseSigReversalHeatmap")
                                                                            ),
                                                                            tabPanel("Barplots",
                                                                                     selectizeInput(inputId = "compoundToBarplot",
                                                                                                    choices = L1000_compounds,
                                                                                                    selected = "alisertib",
                                                                                                    label = "Select small molecule to plot"
                                                                                     ),
                                                                                     plotOutput(outputId = "slicedDiseaseSigReversalBarplot")
                                                                            )
                                                                ),
                                                                downloadButton(outputId = "reversedSlicedDiseaseSigDL", label = "Download reversal score matrix")
                                               )
                                             )
                                    )
                                  )
                 )
                 ), # <-- Closes 3. Disease Signatures Tab
                 
                 # --- 4. CELL-DRUG CONNECTIVITY ---
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-pills"),
                   tags$span("4. Cell-Drug Connectivity"),
                   tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    br(),
                                    p("Please first upload a seurat object as an .RDS file"),
                                    hr(),
                                    br(),
                                    p(em("Due to the size of some datasets, please allow a few minutes after upload completes for the file to be loaded into the working environment"))
                                  )
                 ),
                 conditionalPanel(condition = "output.seuratLoaded",
                                  conditionalPanel(condition = "output.controlsNotSelected",
                                                   conditionalPanel(condition = "!input.noNormalCellIdentsCheckbox",
                                                                    h4("No test cell populations selected"),
                                                                    br(),
                                                                    p(strong("Please return to the pre-processing tab to select test cell populations")),
                                                                    hr()
                                                   ),
                                                   conditionalPanel(condition = "input.noNormalCellIdentsCheckbox",
                                                                    "in!",
                                                                    conditionalPanel(condition="output.controlsNotSelected_noControl",
                                                                                     h4("No test cell populations selected"),
                                                                                     br(),
                                                                                     p(strong("Please return to the pre-processing tab to select test cell populations")),
                                                                                     hr()
                                                                    ),
                                                                    conditionalPanel(condition="output.controlsSelected_noControl",
                                                                                     wellPanel(
                                                                                       h4("Step 1: Calculate L1000 small molecule vs single-cell connectivity"),
                                                                                       br(),
                                                                                       splitLayout(
                                                                                         wellPanel(
                                                                                           conditionalPanel(condition = "output.uploadCorrelationMatrix",
                                                                                                            hr(),
                                                                                                            uiOutput("L1000_release_InSilico"),
                                                                                                            actionButton(inputId = "CalculateRDS_L1000_Spearman_Mat", label = "Calculate single-cell compound discordance")
                                                                                           ),
                                                                                           checkboxInput(inputId = "uploadCorrelationMatrix", label = "Upload a previously calculated Cell-Drug Connectivity file."),
                                                                                           conditionalPanel(condition = "input.uploadCorrelationMatrix",
                                                                                                            fileInput(inputId = "corrMatUpload",
                                                                                                                      label = NULL,
                                                                                                                      buttonLabel = "Browse...",
                                                                                                                      placeholder = "No file selected",
                                                                                                                      width = NULL,
                                                                                                                      multiple = F,
                                                                                                                      accept = c(".csv", ".CSV"))
                                                                                           )
                                                                                         ),
                                                                                         wellPanel(
                                                                                           conditionalPanel(condition = "output.corrMatrixCalculated",
                                                                                                            hr(),
                                                                                                            h4("Discordance calculations complete. Download to avoid recalculation"),
                                                                                                            downloadButton(outputId = "RDScorrMatDownload", label = "Download single-cell vs small molecule correlations")
                                                                                           ),
                                                                                           conditionalPanel(condition = "output.corrMatUploaded",
                                                                                                            hr(),
                                                                                                            h5("Discordance matrix upload detected!"),
                                                                                                            h5("Integrating into Seurat object..."),
                                                                                                            hr()
                                                                                           )
                                                                                         )
                                                                                       ),
                                                                                       h5("Step 1.2 (Optional): Calculate cell-drug connectivities for custom TCS input."),
                                                                                       checkboxInput(inputId = "customTCSasRef", label = "Utilize custom TCS signature?"),
                                                                                       wellPanel(
                                                                                         conditionalPanel(condition = "input.customTCSasRef",
                                                                                                          splitLayout(
                                                                                                            wellPanel(
                                                                                                              textInput(inputId = "Custom_TCS_nameInput",
                                                                                                                        label = "1.2.1: Enter a name for your uploaded TCS (avoid spaces)",
                                                                                                                        placeholder = "Custom_Drug_Signature", width = "50%"),
                                                                                                              fileInput(inputId = "customTCScsvPath",
                                                                                                                        label = "1.2.2: Upload custom response signature. (.csv)",
                                                                                                                        buttonLabel = "Browse...",
                                                                                                                        placeholder = "No file selected",
                                                                                                                        width = "50%",
                                                                                                                        multiple = F,
                                                                                                                        accept = c(".csv")),
                                                                                                              conditionalPanel(condition = 'output.TCSuploaded',
                                                                                                                               actionButton(inputId = "customTCSbrowser_button",
                                                                                                                                            label = "Browse TCS upload",
                                                                                                                                            icon = icon("table")),
                                                                                                                               hr(),
                                                                                                                               radioButtons(inputId = "TCS_scoring_method",
                                                                                                                                            label = "1.2.3: Select connectivity scoring approach:",
                                                                                                                                            choices = list("Spearman", "Singscore", "ssGSEA"),
                                                                                                                                            selected = "Spearman", inline = T),
                                                                                                                               actionButton(inputId = "CalculateCustomTCSConnectivities", label = "Calculate custom TCS single-cell connectivity"),
                                                                                                                               hr(),
                                                                                                                               conditionalPanel(condition = 'input.CalculateCustomTCSConnectivities',
                                                                                                                                                downloadButton(outputId = 'downloadCustomTCScellConnectivities', label = "Download custom TCS connectivities"))
                                                                                                              )
                                                                                                            ),
                                                                                                            wellPanel(
                                                                                                              conditionalPanel(condition = 'output.TCSuploaded',
                                                                                                                               hr(),
                                                                                                                               h5("Success! Custom TCS loaded."),
                                                                                                                               hr())
                                                                                                            )
                                                                                                          )
                                                                                         )
                                                                                       ),
                                                                                       tags$hr(style = "border-top: 3px solid black;"),
                                                                                       conditionalPanel(condition = "output.corrMatrixCalculated || output.corrMatUploaded",
                                                                                                        h4("Step 1: Select reference compound and set parameters."),
                                                                                                        br(),
                                                                                                        radioButtons(inputId = "Reference_L1000_or_Custom",
                                                                                                                     label = "Step 1.1: L1000 or Custom TCS?",
                                                                                                                     choices = list("L1000 Derived", "Custom TCS Upload"),
                                                                                                                     selected = "L1000 Derived", inline = TRUE),
                                                                                                        conditionalPanel(condition = "input.Reference_L1000_or_Custom == 'L1000 Derived'",
                                                                                                                         wellPanel(
                                                                                                                           splitLayout(
                                                                                                                             wellPanel(
                                                                                                                               selectizeInput(inputId = "referenceCompound",
                                                                                                                                              label = "Select reference compound",
                                                                                                                                              choices = L1000_compounds,
                                                                                                                                              selected = "alisertib",
                                                                                                                                              multiple = F,
                                                                                                                                              width = NULL,
                                                                                                                                              size = NULL),
                                                                                                                               sliderInput(inputId = "sigMinCutoffRDS",
                                                                                                                                           label = "Set visualization cutoff minimum:",
                                                                                                                                           value = -1, min = -1, max = 1, step = 0.1),
                                                                                                                               sliderInput(inputId = "sigMaxCutoffRDS",
                                                                                                                                           label = "Set visualization cutoff maximum:",
                                                                                                                                           value = 0, min = -1, max = 1, step = 0.1),
                                                                                                                               sliderInput(inputId = "correlationCutoff_RDS", "correlationCutoff_RDS",
                                                                                                                                           label = "Set sensitivity and resistance cut-offs [adjustable pseudodose]",
                                                                                                                                           value = c(0, 0), min = -1, max = 1, step = 0.05)
                                                                                                                             ),
                                                                                                                             plotOutput(outputId = "refSignatureUMAPRDS")
                                                                                                                           )
                                                                                                                         ),
                                                                                                                         plotOutput(outputId = "scSynergySeq1RDS"),
                                                                                                                         hr(),
                                                                                                                         br(),
                                                                                                                         wellPanel(
                                                                                                                           splitLayout(
                                                                                                                             textOutput("referenceCompound"), # FIXED TRAILING COMMA
                                                                                                                             actionButton("perturbationButton", label = "Run Perturbation Analysis", icon = icon("redo"))
                                                                                                                           ),
                                                                                                                           conditionalPanel(condition = "input.perturbationButton",
                                                                                                                                            hr(),
                                                                                                                                            h4("Success! Please navigate to the results tab."),
                                                                                                                                            hr())
                                                                                                                         )
                                                                                                        )
                                                                                       )
                                                                                     )
                                                                    )
                                                   )
                                  ),
                                  conditionalPanel(condition = "output.controlsSelected", 
                                                   wellPanel(
                                                     h4("Step 1: Calculate L1000 small molecule vs single-cell connectivity"),
                                                     br(),
                                                     splitLayout(
                                                       wellPanel(
                                                         conditionalPanel(condition = "output.uploadCorrelationMatrix",
                                                                          hr(),
                                                                          uiOutput("L1000_release_InSilico"),
                                                                          actionButton(inputId = "CalculateRDS_L1000_Spearman_Mat", label = "Calculate single-cell compound discordance")
                                                         ),
                                                         checkboxInput(inputId = "uploadCorrelationMatrix", label = "Upload a previously calculated Cell-Drug Connectivity file"),
                                                         conditionalPanel(condition = "input.uploadCorrelationMatrix",
                                                                          fileInput(inputId = "corrMatUpload",
                                                                                    label = NULL,
                                                                                    buttonLabel = "Browse...",
                                                                                    placeholder = "No file selected",
                                                                                    width = NULL,
                                                                                    multiple = F,
                                                                                    accept = c(".csv", ".CSV"))
                                                         )
                                                       ),
                                                       wellPanel(
                                                         conditionalPanel(condition = "output.corrMatrixCalculated",
                                                                          hr(), 
                                                                          h4("Discordance calculations complete. Download to avoid recalculation"),
                                                                          downloadButton(outputId = "RDScorrMatDownload", label = "Download single-cell vs small molecule correlations")
                                                         ),
                                                         conditionalPanel(condition = "output.corrMatUploaded",
                                                                          hr(),
                                                                          h5("Discordance matrix upload detected!"),
                                                                          h5("Integrating into Seurat object..."),
                                                                          hr()
                                                         )
                                                       )
                                                     ),
                                                     h5("Step 1.2 (Optional): Calculate cell-drug connectivities for custom TCS input."),
                                                     checkboxInput(inputId = "customTCSasRef", label = "Utilize custom TCS signature?"),
                                                     wellPanel(
                                                       conditionalPanel(condition = "input.customTCSasRef",
                                                                        splitLayout(
                                                                          wellPanel(
                                                                            textInput(inputId = "Custom_TCS_nameInput",
                                                                                      label = "1.2.1: Enter a name for your uploaded TCS (avoid spaces)",
                                                                                      placeholder = "Custom_Drug_Signature", width = "50%"),
                                                                            fileInput(inputId = "customTCScsvPath",
                                                                                      label = "1.2.2: Upload custom response signature. (.csv)",
                                                                                      buttonLabel = "Browse...",
                                                                                      placeholder = "No file selected",
                                                                                      width = "50%",
                                                                                      multiple = F,
                                                                                      accept = c(".csv")),
                                                                            conditionalPanel(condition = 'output.TCSuploaded',
                                                                                             actionButton(inputId = "customTCSbrowser_button",
                                                                                                          label = "Browse TCS upload",
                                                                                                          icon = icon("table")),
                                                                                             hr(),
                                                                                             radioButtons(inputId = "TCS_scoring_method",
                                                                                                          label = "1.2.3: Select connectivity scoring approach:",
                                                                                                          choices = list("Spearman", "Singscore", "ssGSEA"),
                                                                                                          selected = "Spearman", inline = T),
                                                                                             actionButton(inputId = "CalculateCustomTCSConnectivities", label = "Calculate custom TCS single-cell connectivity"),
                                                                                             hr(),
                                                                                             conditionalPanel(condition = 'input.CalculateCustomTCSConnectivities',
                                                                                                              downloadButton(outputId = 'downloadCustomTCScellConnectivities', label = "Download custom TCS connectivities"))
                                                                            )
                                                                          ),
                                                                          wellPanel(
                                                                            conditionalPanel(condition = 'output.TCSuploaded',
                                                                                             hr(),
                                                                                             h5("Success! Custom TCS loaded."),
                                                                                             hr())
                                                                          )
                                                                        )
                                                       )
                                                     )
                                                   )
                                  )
                 )
                 ), # <-- Closes 4. Cell-Drug Connectivity Tab
                 
                 # --- 5. IN SILICO PERTURBATION ---
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-magnifying-glass-chart"),
                   tags$span("5. In Silico Perturbation"),
                   tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    br(),
                                    p("Please first upload a seurat object as an .RDS file"),
                                    hr(),
                                    br(),
                                    p(em("Due to the size of some datasets, please allow a few minutes after upload completes for the file to be loaded into the working environment"))
                                  )
                 ),
                 conditionalPanel(condition = "output.seuratLoaded",
                                  conditionalPanel(condition = "output.controlsNotSelected",
                                                   h4("No test cell populations selected"),
                                                   br(),
                                                   p(strong("Please return to the pre-processing tab to select test cell populations")),
                                                   hr()
                                  ),
                                  conditionalPanel(condition = "output.controlsSelected",
                                                   hr(),
                                                   h4("'Treat' your dataset to characterize sensitive and resistant cell populations"),
                                                   br(),
                                                   tabsetPanel(
                                                     tabPanel(tags$div(
                                                       tags$i(class = "fa-sharp fa-solid fa-magnifying-glass-chart"),
                                                       tags$span("Initialize Perturbation Experiment"), 
                                                       tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                                                     ),
                                                     tags$hr(style = "border-top: 3px solid black;"),
                                                     conditionalPanel(condition = "output.corrMatrixCalculated || output.corrMatUploaded",
                                                                      h4("Step 1: Select reference compound and set parameters."),
                                                                      br(),
                                                                      radioButtons(inputId = "Reference_L1000_or_Custom",
                                                                                   label = "Step 1.1: L1000 or Custom TCS?",
                                                                                   choices = list("L1000 Derived", "Custom TCS Upload"),
                                                                                   selected = "L1000 Derived", inline = TRUE),
                                                                      conditionalPanel(condition = "input.Reference_L1000_or_Custom == 'L1000 Derived'",
                                                                                       wellPanel(
                                                                                         splitLayout(
                                                                                           wellPanel(
                                                                                             selectizeInput(inputId = "referenceCompound",
                                                                                                            label = "Select reference compound",
                                                                                                            choices = L1000_compounds,
                                                                                                            selected = "alisertib",
                                                                                                            multiple = F,
                                                                                                            width = NULL,
                                                                                                            size = NULL),
                                                                                             sliderInput(inputId = "sigMinCutoffRDS",
                                                                                                         label = "Set visualization cutoff minimum:",
                                                                                                         value = -1, min = -1, max = 1, step = 0.1),
                                                                                             sliderInput(inputId = "sigMaxCutoffRDS",
                                                                                                         label = "Set visualization cutoff maximum:",
                                                                                                         value = 0, min = -1, max = 1, step = 0.1),
                                                                                             sliderInput(inputId = "correlationCutoff_RDS", "correlationCutoff_RDS",
                                                                                                         label = "Set sensitivity and resistance cut-offs [adjustable pseudodose]",
                                                                                                         value = c(0, 0), min = -1, max = 1, step = 0.05)
                                                                                           ),
                                                                                           plotOutput(outputId = "refSignatureUMAPRDS")
                                                                                         )
                                                                                       ),
                                                                                       plotOutput(outputId = "scSynergySeq1RDS"),
                                                                                       wellPanel(
                                                                                         splitLayout(
                                                                                           wellPanel(
                                                                                             h4("Step 2: In Silico Perturbation"),
                                                                                             textOutput("referenceCompound") # FIXED TRAILING COMMA
                                                                                           ),
                                                                                           wellPanel(
                                                                                             h4(),
                                                                                             actionButton("perturbationButton", label = "Run Perturbation Analysis", icon = icon("redo"))
                                                                                           )
                                                                                         ),
                                                                                         conditionalPanel(condition = "input.perturbationButton",
                                                                                                          hr(),
                                                                                                          h4("Success! Please navigate to the results tab."),
                                                                                                          hr())
                                                                                       )
                                                                      ),
                                                                      conditionalPanel(condition = "input.Reference_L1000_or_Custom =='Custom TCS Upload'",
                                                                                       conditionalPanel(condition = "input.CalculateCustomTCSConnectivities",
                                                                                                        splitLayout(
                                                                                                          wellPanel(
                                                                                                            sliderInput(inputId = "sigMinCutoff_customTCS",
                                                                                                                        label = "Set visualization cutoff minimum:",
                                                                                                                        value = -1, min = -1, max = 1, step = 0.1),
                                                                                                            sliderInput(inputId = "sigMaxCutoff_customTCS",
                                                                                                                        label = "Set visualization cutoff maximum:",
                                                                                                                        value = 0, min = -1, max = 1, step = 0.1),
                                                                                                            sliderInput(inputId = "correlationCutoff_customTCS", "correlationCutoff_customTCS",
                                                                                                                        label = "Set sensitivity and resistance cut-offs [adjustable pseudodose]",
                                                                                                                        value = c(0, 0), min = -1, max = 1, step = 0.05)
                                                                                                          ),
                                                                                                          wellPanel(plotOutput(outputId = "CustomRefSignatureUMAPRDS"))
                                                                                                        ),
                                                                                                        wellPanel(
                                                                                                          plotOutput(outputId = "scSynergySeq_custTCS")
                                                                                                        ),
                                                                                                        wellPanel(
                                                                                                          splitLayout(
                                                                                                            textOutput("referenceCompound_customTCS"), # FIXED TRAILING COMMA
                                                                                                            actionButton("perturbationButton_customTCS", label = "Run Perturbation Analysis", icon = icon("redo"))
                                                                                                          ),
                                                                                                          conditionalPanel(condition = "input.perturbationButton_customTCS",
                                                                                                                           hr(),
                                                                                                                           h4("Success! Please navigate to the results tab."),
                                                                                                                           hr()
                                                                                                          )
                                                                                                        )
                                                                                       )
                                                                      )
                                                     ),
                                                     conditionalPanel(condition = "input.CalculateCustomTCSConnectivities == 0",
                                                                      "custom TCS connectivities not yet calculated, please return to the 'Cell-Drug Connectivity' Tab.")
                                                     ),
                                                     
                                                     tabPanel(tags$div(
                                                       tags$i(class = "fa-sharp fa-solid fa-magnifying-glass-chart"),
                                                       tags$span("Results and Analysis"), 
                                                       tags$style(type = "text/css", "li a{color:#000000; font-family: 'sans-serif', Arial Rounded MT Bold;}")
                                                     ),
                                                     conditionalPanel(condition = "input.perturbationButton == 0",
                                                                      hr()
                                                     ),
                                                     conditionalPanel(condition = "input.perturbationButton != 0",
                                                                      br(),
                                                                      h4("L1000 TCS Result: "),
                                                                      hr(),
                                                                      tabsetPanel(
                                                                        tabPanel("Resistance Signature",
                                                                                 wellPanel(
                                                                                   splitLayout(
                                                                                     wellPanel(
                                                                                       sliderInput(inputId = "cellSubsetMax_pert",
                                                                                                   label = "Select number of cells per condition to test",
                                                                                                   min = 100, max = 40000, step = 100, value = 500
                                                                                       ),
                                                                                       selectizeInput(inputId = "pertDEx_testUse",
                                                                                                      label = "Select differential expression test to use",
                                                                                                      choices = c("MAST", "roc", "wilcoxon", "deseq"),
                                                                                                      selected = "MAST"),
                                                                                       actionButton(inputId = "calcSensVsResistantDEGS",
                                                                                                    label = "Run differential expression analysis"),
                                                                                       conditionalPanel(condition = "output.resVsSensCalced",
                                                                                                        hr(),
                                                                                                        downloadButton(outputId = "resVsSensDEGSdownload",
                                                                                                                       label = "Download differential expression results")
                                                                                       )
                                                                                     ),
                                                                                     wellPanel(
                                                                                       plotOutput(outputId = "resVsSensVolcano")
                                                                                     )
                                                                                   ),
                                                                                   DT::dataTableOutput(outputId = "resVsSensTable")
                                                                                 )
                                                                        ),
                                                                        tabPanel("Proportion shift",
                                                                                 wellPanel(
                                                                                   splitLayout(
                                                                                     plotOutput(outputId = "scFOCAL_StateShiftAlluvial")
                                                                                   ),
                                                                                   wellPanel(
                                                                                     sliderInput(
                                                                                       inputId = "quantile_slider",
                                                                                       label = "Select Quantile Cutoff (%):",
                                                                                       min = 5,
                                                                                       max = 45,
                                                                                       value = 20, 
                                                                                       step = 5
                                                                                     ),
                                                                                     wellPanel(
                                                                                       style = "height: 1000px;",
                                                                                       plotOutput("main_plot", height = "100%")
                                                                                     )
                                                                                   )
                                                                                 )
                                                                        ),
                                                                        tabPanel("Combination Scoring",
                                                                                 splitLayout(
                                                                                   wellPanel(
                                                                                     plotOutput(outputId = "l1000in_limmaVolc"),
                                                                                     downloadButton(outputId = "download_limmaVolc", label = "Download volcano .pdf")
                                                                                   ),
                                                                                   wellPanel(
                                                                                     plotOutput(outputId = "rcm_overall")
                                                                                   )
                                                                                 ),
                                                                                 wellPanel(
                                                                                   plotOutput(outputId = 'combinationScoreBar'),
                                                                                   actionButton(inputId = 'combinationTableButton',
                                                                                                label = "View data in table", icon = icon("table")),
                                                                                   downloadButton(outputId = "combinationScoreDownload",
                                                                                                  label = "Download Combination Scoring Results Table")
                                                                                 )
                                                                        )
                                                                      )
                                                     ),
                                                     conditionalPanel(condition = "input.perturbationButton_customTCS",
                                                                      hr(),
                                                                      h4("Custom TCS Result"),
                                                                      hr(),
                                                                      downloadButton(outputId = "download_customTCS_sensitivityAssigner", label = "Download custom TCS sensitivity dataframe"),
                                                                      hr(),
                                                                      tabsetPanel(
                                                                        tabPanel("Resistance Signature",
                                                                                 wellPanel(
                                                                                   splitLayout(
                                                                                     wellPanel(
                                                                                       sliderInput(inputId = "cellSubsetMax_pert_customTCS",
                                                                                                   label = "Select number of cells per condition to test",
                                                                                                   min = 100, max = 40000, step = 100, value = 500
                                                                                       ),
                                                                                       selectizeInput(inputId = "pertDEx_testUse_customTCS",
                                                                                                      label = "Select differential expression test to use",
                                                                                                      choices = c("MAST", "roc", "wilcoxon", "deseq"),
                                                                                                      selected = "MAST"),
                                                                                       actionButton(inputId = "calcSensVsResistantDEGS_customTCS",
                                                                                                    label = "Run differential expression analysis"),
                                                                                       conditionalPanel(condition = "output.resVsSensCalced_customTCS",
                                                                                                        hr(),
                                                                                                        downloadButton(outputId = "resVsSensDEGSdownload_customTCS",
                                                                                                                       label = "Download differential expression results")
                                                                                       )
                                                                                     ),
                                                                                     wellPanel(
                                                                                       plotOutput(outputId = "resVsSensVolcano_customTCS")
                                                                                     )
                                                                                   ),
                                                                                   DT::dataTableOutput(outputId = "resVsSensTable_customTCS")
                                                                                 )
                                                                        ),
                                                                        tabPanel("Proportion Shift",
                                                                                 wellPanel(
                                                                                   splitLayout(
                                                                                     wellPanel(
                                                                                       plotOutput(outputId = "scFOCAL_StateShiftAlluvial_customTCS")
                                                                                     ),
                                                                                     wellPanel()
                                                                                   ),
                                                                                   wellPanel(
                                                                                     sliderInput(
                                                                                       inputId = "quantile_slider_customTCS",
                                                                                       label = "Select Quantile Cutoff (%):",
                                                                                       min = 5,
                                                                                       max = 45,
                                                                                       value = 20, 
                                                                                       step = 5
                                                                                     ),
                                                                                     wellPanel(
                                                                                       style = "height: 1000px;",
                                                                                       plotOutput("main_plot_customTCS", height = "100%")
                                                                                     )
                                                                                   )
                                                                                 )
                                                                        ),
                                                                        tabPanel("Combination Scoring",
                                                                                 hr("Rank L1000 Small Molecules for combination with your reference drug using the scFOCAL combination index."),
                                                                                 splitLayout(
                                                                                   wellPanel(
                                                                                     plotOutput(outputId = "customTCS_limmaVolc"),
                                                                                     downloadButton(outputId = "download_customTCS_limmaVolc", label = "Download volcano .pdf")
                                                                                   ),
                                                                                   wellPanel(
                                                                                     plotOutput(outputId = "customTCS_rcm_overall")
                                                                                   )
                                                                                 ),
                                                                                 wellPanel(
                                                                                   plotOutput(outputId = 'customTCS_combinationScoreBar'),
                                                                                   actionButton(inputId = 'customTCS_combinationTableButton',
                                                                                                label = "View data in table", icon = icon("table")),
                                                                                   downloadButton(outputId = "combinationScoreDownload_customTCS",
                                                                                                  label = "Download Combination Scoring Results Table")
                                                                                 )
                                                                        )
                                                                      )
                                                     )
                                                     )
                                                   )
                                  )
                 )
                 )
               )
               ), # <-- Closes 5. In Silico Perturbation Tab
               
               # ------------------------------------------------------------------
               # TUTORIALS TAB
               # ------------------------------------------------------------------
               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-question"),
                 tags$span("- Tutorials")
               ), 
               h1(strong("Tutorials")),
               hr(),
               br(),
               p(strong("Welcome to scFOCAL! This tutorial will walk you through the scFOCAL workflow.")),
               h2(strong("Data Upload")),
               p("Upload data as a Seurat object that contains expression data from RNAseq. File format should be .RDS or .rds. Please wait for a success message before moving to Step 2. This may take a few minutes."),
               h2(strong("Data Pre-Processing")),
               p("Select cell population paramaters for further study. Identify which cell population(s) you would like to set as controls and which cell population(s) you would like to set as transformed. These selections will affect how disease signatures are calculated (refer to next step)."),
               h2(strong("Disease and Reversal Signatures")),
               p("Now, you can calculate transcriptomic disease signatures compared to healthy controls. Use the slider bar to select how many cells you would like to test per condition. Click 'Calculate Disease Signatures'. Selecting more cells to test will give a more accurate signature, but it may take an extended period of time so please be patient. Once the heatmap appears, you will be able to download the disease signature matrix using the 'Download signature .csv' button. You can also download the heatmap as a pdf."),
               p("To study how compounds reverse the disease signature of each cell subtype, navigate to the 'Reversal' tab. Select how many compounds you would like to study. Once the heatmap appears, you will be able to download the corresponding matrix that provides reversal values for every cell subtype with each compound tested. To do this, click the 'Download reversal score matrix' button."),
               p("Within the reversal tab, you can use the barplot tab to graphically compare reversal scores of cell populations treated with a single compound. Select a compound of interest using the dropdown menu located under 'Select a small molecule to plot'."),
               h2(strong("In Silico Perturbation")),
               p("To calculate a correlation matrix between L1000 compounds and RNAseq data, first select which release of the L1000 compound set you would like to use. This can be done with the dropdown menu found under 'Select L1000 release'. Then to calculate the correlation matrix, click on the 'Calculate single-cell compound discordance' button. When calculation is complete, you can download the matrix. Since this is a computationally-intense calculation, it may take extended time. Therefore, we highly recommend that you download the correlation matrix to avoid recalculation. Additionally, you have the option of uploading your own correlation matrix if you have already calculated it. If you would like to upload your own correlation matrix, please ensure that it is in .CSV or .csv format."),
               p("You can now select a compound of interest for synergy analysis. Under 'Step 2: Select reference compound and set parameters', select the compound to which you would like to identify synergistic compounds. Additionally, select computational parameters."),
               p("Click on 'Run Perturbation Analysis'. You will receive a Success message when analysis is complete, at which point you can navigate to the Results tab to study the analysis.")
               ),
               
               # ------------------------------------------------------------------
               # INSTALL R PACKAGE TAB
               # ------------------------------------------------------------------
               tabPanel(tags$div(
                 tags$i(class = "fa-brands fa-github"),
                 tags$span("- Install R Package")
               ), 
               h1("scFOCAL R Package"),
               hr(),
               br(),
               p(strong("Due to the likely large size of many single-cell RNA seq datasets, the scFOCAL shiny app is locally deployable via the scFOCAL R Package. Please read the license agreement below and agree to access install instructions.")),
               h1("License Agreement"),
               br(), 
               p(strong("1. The Board of Trustees of the Georgetown University (?Georgetown?) provides scFOCAL software and code (?Service?) free of charge for non-commercial use only. Use of the Service by any commercial entity for any purpose, including research, is prohibited.")),
               p(strong("2. By using the Service, you agree to be bound by the terms of this Agreement. Please read it carefully.")),
               p(strong("3. You agree not to use the Service for commercial advantage, or in the course of for-profit activities. You agree not to use the Service on behalf of any organization that is not a non-profit organization. Commercial entities wishing to use this Service should contact Georgetown University?s Office of Technology Licensing.")),
               p(strong("4. THE SERVICE IS OFFERED ?AS IS?, AND, TO THE EXTENT PERMITTED BY LAW, GEORGETOWN MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED. GEORGETOWN SHALL NOT BE LIABLE FOR ANY CLAIMS OR DAMAGES WITH RESPECT TO ANY LOSS OR OTHER CLAIM BY YOU OR ANY THIRD PARTY ON ACCOUNT OF, OR ARISING FROM THE USE OF THE SERVICE. YOU HEREBY AGREE TO DEFEND AND INDEMNIFY GEORGETOWN, ITS TRUSTEES, EMPLOYEES, OFFICERS, STUDENTS, AGENTS, FACULTY, REPRESENTATIVES, AND VOLUNTEERS (?GEORGETOWN INDEMNITEES?) FROM ANY LOSS OR CLAIM ASSERTED AGAINST GEORGETOWN INDEMNITEES ARISING FROM YOUR USE OF THE SERVICE.")),
               p(strong("5. All rights not expressly granted to you in this Agreement are reserved and retained by GEORGETOWN or its licensors or content providers. This Agreement provides no license under any patent.")),
               p(strong("6. You agree that this Agreement and any dispute arising under it is governed by the laws of the District of Columbia, United States of America, applicable to agreements negotiated, executed, and performed within the DISTRICT OF COLUMBIA")),
               p(strong("7. Subject to your compliance with the terms and conditions set forth in this Agreement, GEORGETOWN grants you a revocable, non-exclusive, non-transferable right to access and make use of the Service.")),
               br(),
               p(em("Do you accept the terms and conditions in this agreement?")),
               fluidRow(column(width = 4, actionButton(inputId = "acceptAgreement", label = "Accept")), column(width = 4, actionButton(inputId = "rejectAgreement", label = "Reject and close"))),
               conditionalPanel(condition = "input.acceptAgreement",
                                "Please scroll down for instructions."),
               br(),
               conditionalPanel(condition = "input.acceptAgreement",
                                p(strong("Thank You!")),
                                h2("How to install the scFOCAL R package"),
                                p(strong("The scFOCAL R Package can be installed from github using devtools.")),
                                p(code("devtools::install_github('AyadLab/scFOCAL')")),
                                p(strong("Once installed, load the scFOCAL library.")),
                                p(code("library(scFOCAL)")),
                                p(strong("Launch the scFOCAL shiny UI using the RunscFOCAL() function.")),
                                p(code("scFOCAL::runscFOCAL()"))
               )
               ),
               
               # ------------------------------------------------------------------
               # CONTACT TAB
               # ------------------------------------------------------------------
               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-envelope"),
                 tags$span("- Contact")
               ), 
               
               tags$h4(
                 a(href = "https://github.com/AyadLab/scFOCAL",
                   "For issues: Please report them on the scFOCAL git."),
                 align = "center"
               ),
               hr(),
               splitLayout(
                 wellPanel(
                   style = "text-align: center;",
                   br(),
                   p(strong("Robert K. Suter, PhD")),
                   p(em("Assistant Professor")),
                   p(em("Department of Oncology")),
                   p(em("Lombardi Comprehensive Cancer Center")),
                   p(em("Georgetown University Medical Center")),
                   hr(),
                   img(src = "rks_headshot_2022.png", height = "160px"),
                   br(), br(),
                   a(href="mailto:rks82@georgetown.edu", "Contact"),
                   hr(), 
                   tags$h4(
                     a(href="https://suterlab.com", "SuterLab.com"),
                     align = "center"
                   ),
                   hr()
                 ),
                 
                 wellPanel(
                   style = "text-align: center;",
                   br(),
                   p(strong("Nagi G. Ayad, PhD")),
                   p(em("Professor")),
                   p(em("Department of Oncology")),
                   p(em("Lombardi Comprehensive Cancer Center")),
                   p(em("Georgetown University Medical Center")),
                   hr(),
                   img(src = "nagi_headshot.png", height = "160px"),
                   br(), br(),
                   a(href="mailto:na853@georgetown.edu", "Contact"),
                   hr(),
                   tags$h4(
                     a(href="https://sites.google.com/georgetown.edu/ayadlab/home", "AyadLab.com"),
                     align = "center"
                   ),
                   hr()
                 )
               )
               )
  )
)
