library(shiny)
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

# LINCS Response Signature Data
LINCS.ResponseSigs <- read.delim(
  file = "matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt", header = T)
row.names(LINCS.ResponseSigs) <- LINCS.ResponseSigs$Genes
LINCS.ResponseSigs <- na.omit(LINCS.ResponseSigs)
newNames <- gsub("-", ".", rownames(LINCS.ResponseSigs))
rownames(LINCS.ResponseSigs) <- newNames
L1000_genes <- colnames(LINCS.ResponseSigs)
L1000_compounds <- rownames(LINCS.ResponseSigs)

ui <- fluidPage(
  tags$script(src = "https://kit.fontawesome.com/070e476711.js"),
  img(src = "ISOSCELES.png", width = 700, height = 100),
  hr(),
  # br(),
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

                            "))),
  navlistPanel(well = T, fluid = T, widths = c(3, 7),
               # tags$head(tags$style(HTML(".tab-content { height: 83vh; overflow-y: auto !important; }" ))),
               tabPanel(tags$div(
                 tags$i(class = "fa-sharp fa-solid fa-desktop"),
                 tags$span("- Overview"),
                 tags$style(type = "text/css", "li a{color:#000000;}")
               ), img(src = "scSynergySeq Diagram (2).png", width = 500, height = 500)
               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-magnifying-glass-chart"),
                 tags$span("- Run ISOSCELES")
               ), #put contents for actual application here

               tabsetPanel(
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-upload"),
                   tags$span("1. Data Upload"),
                   tags$style(type = "text/css", "li a{color:#000000; font-samily: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 br(),
                 h4("Upload Seurat object RDS file"),
                 hr(),
                 br(),
                 br(),
                 splitLayout(
                   fileInput(inputId = "seurobjRDS",
                             label = NULL,
                             buttonLabel = "Browse...",
                             placeholder = "No file selected",
                             width = NULL,
                             multiple = F,
                             accept = c(".RDS", ".Rds", ".rds")),
                   conditionalPanel(condition = 'output.seuratLoaded',
                                    h4("Success! Seurat object loaded."))
                 ),
                 br(),
                 hr(),
                 p(em("Depending on file size, after the file upload above is complete, it will still take some additional time for the data to be loaded into the processing environment. Please be patient."))
                 ),
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-gears"),
                   tags$span("2. Pre-processing"),
                   tags$style(type = "text/css", "li a{color:#000000; font-samily: 'sans-serif', Arial Rounded MT Bold;}")
                 ),

                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    # hr(),
                                    br(),
                                    p("Please first upload a seurat object as an .RDS file"),
                                    hr(),
                                    br(),
                                    p(em("Due to the size of some datasets, please allow a few minutes after upload completes for the file to be loaded into the working environment"))
                                  )
                 ),

                 conditionalPanel(condition = "output.seuratLoaded",
                                  splitLayout(
                                    wellPanel(
                                      p(strong("Data Pre-processing")),
                                      hr(),
                                      br(),
                                      uiOutput(outputId = "reductionUseRDS"),
                                      uiOutput(outputId = "groupByRDS"),
                                      uiOutput(outputId = "normalCellIdentsRDS"),
                                      uiOutput('cancerCellIdentsRDS')),
                                    tabsetPanel(type = "tabs",
                                                tabPanel(
                                                  "2D",
                                                  plotOutput(outputId = "dimplotRDS") # I'm guessing this is also a problem using the conditional panels now..
                                                ),
                                                tabPanel(
                                                  "3D",
                                                  "Under Construction!"
                                                ))
                                  ),
                                  splitLayout(
                                    wellPanel(
                                      # p(strong("Control Cells")),
                                      plotOutput(outputId = "controlHighlight")
                                    ),
                                    wellPanel(
                                      # p(strong("Diseased Cells")),
                                      plotOutput(outputId = "diseasedHighlight")
                                    )

                                  )
                 )

                 ),
                 tabPanel(tags$div(
                   tags$i(class = "fa-solid fa-signature"),
                   tags$span("3. Disease Signatures"),
                   tags$style(type = "text/css", "li a{color:#000000; font-samily: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    br(),
                                    p("Please first upload a seurat object as an .RDS file"),
                                    hr(),
                                    br(),
                                    p(em("Due to the size of some datasets, please allow a few minutes after upload completes for the file to be loaded into the working environment"))
                                  ),
                 ),
                 conditionalPanel(condition = "output.seuratLoaded",

                                  wellPanel(
                                    img(src = "SlicedDiseaseSignatureCalc.png", width = 650, height = 220),
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
                                                                 label = "Select number of cells per condition to test",
                                                                 min = 100, max = 40000, step = 100, value = 500
                                                     ),
                                                     actionButton(inputId = "CalcMarkers",
                                                                  label = "Identify Cluster Markers",
                                                                  icon = icon("sync")),
                                                     actionButton(inputId = "CalcDiseaseSig",
                                                                  label = "Calculate Disease Signatures",
                                                                  icon = icon("sync"))
                                                   )

                                  ),


                                  conditionalPanel(condition = "input.CalcDiseaseSig >= 1", # Change this condition to detect completion of diffexp testing
                                                   tabsetPanel(

                                                     tabPanel(title = "Heatmap",

                                                              wellPanel(
                                                                plotOutput(outputId = "diseaseSigRDS"),
                                                                downloadButton("diseaseSigDownload", label = "Download signature .csv"),
                                                                downloadButton("diseaseSigHeatmapDownload", label = "Download Plot .pdf")
                                                              )#,
                                                              # wellPanel(tableOutput(outputId = "diseaseSigTable"))

                                                     ),

                                                     tabPanel(title = "Table",
                                                              wellPanel(DT::dataTableOutput(outputId = "diseaseSigTable")
                                                              )
                                                     ),

                                                     tabPanel(title = "Compare",
                                                              "Under construction!"),

                                                     tabPanel(title = "Reversal",
                                                              wellPanel(
                                                                radioButtons(inputId = "whichDiseaseSignatures",
                                                                             choices = c("Whole dataset disease signatures"),
                                                                             inline = F, selected = "Whole dataset disease signatures",
                                                                             label = "Select which disease signature set to use."),
                                                                conditionalPanel(condition = 'input.whichDiseaseSignatures == "Sliced dataset disease signatures"',
                                                                                 conditionalPanel(condition = 'input.calcSlicedDiseaseSigs >= 1', # Need a message to show that you need to go back...
                                                                                                  actionButton(inputId = "reverseSlicedDiseaseSigs",
                                                                                                               label = "Score compound reversal")
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

                                                                # If disease signatures haven't been calculated yet, need alternative.
                                                              )

                                                     )
                                                   ))

                 )

                 ),
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-pills"),
                   tags$span("4. In Silico Perturbation"),
                   tags$style(type = "text/css", "li a{color:#000000; font-samily: 'sans-serif', Arial Rounded MT Bold;}")
                 ),

                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    # hr(),
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
                                                   wellPanel(
                                                     h4("Step 1: Calculate L1000 small molecule vs single-cell correlations"),
                                                     br(),
                                                     splitLayout(
                                                       wellPanel(
                                                         conditionalPanel(condition = "output.uploadCorrelationMatrix",
                                                                          hr(),
                                                                          uiOutput("L1000_release_InSilico"),
                                                                          actionButton(inputId = "CalculateRDS_L1000_Spearman_Mat", label = "Calculate single-cell compound discordance")
                                                         ),
                                                         checkboxInput(inputId = "uploadCorrelationMatrix", label = "Alternatively, upload precalculated correlation file", ),
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
                                                                          hr(), # this part needs fixing... should align with left side
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
                                                     # br(),
                                                     hr(),
                                                     conditionalPanel(condition = "output.corrMatrixCalculated || output.corrMatUploaded",
                                                                      h4("Step 2: Select reference compound and set parameters."),
                                                                      br(),

                                                                      wellPanel(
                                                                        splitLayout(
                                                                          wellPanel(
                                                                            selectizeInput(inputId = "referenceCompound",
                                                                                           label = "Select reference compound (searchable...)",
                                                                                           choices = L1000_compounds,
                                                                                           selected = "alisertib",
                                                                                           multiple = F,
                                                                                           width = NULL,
                                                                                           size = NULL),
                                                                            # wellPanel(
                                                                            sliderInput(inputId = "sigMinCutoffRDS",
                                                                                        label = "Set visualization cutoff minimum:",
                                                                                        value = -1, min = -1, max = 1, step = 0.1),
                                                                            sliderInput(inputId = "sigMaxCutoffRDS",
                                                                                        label = "Set visualization cutoff maximum:",
                                                                                        value = 0, min = -1, max = 1, step = 0.1),
                                                                            sliderInput(inputId = "correlationCutoff_RDS", "correlationCutoff_RDS",
                                                                                        label = "Set sensitivity and resistance cut-offs [adjustable pseudodose]",
                                                                                        value = c(0, 0), min = -1, max = 1, step = 0.05)
                                                                            # )
                                                                          ),
                                                                          plotOutput(outputId = "refSignatureUMAPRDS")
                                                                        )
                                                                      ),
                                                                      wellPanel(
                                                                        plotOutput(outputId = "scSynergySeq1RDS")
                                                                      ),
                                                                      hr(),
                                                                      br(),
                                                                      wellPanel(
                                                                        splitLayout(
                                                                          textOutput("referenceCompound", ),
                                                                          actionButton("perturbationButton", label = "Run Perturbation Analysis", icon = icon("redo"))
                                                                        ),
                                                                        conditionalPanel(condition = "input.perturbationButton",
                                                                                         hr(),
                                                                                         h4("Success! Please navigate to the results tab."),
                                                                                         hr())
                                                                      )

                                                     )
                                                     # textOutput('finalMatrixText')
                                                     # sliderInput(inputId = "cellSubsetMax_1",
                                                     #             label = "Select number of cells per condition to test",
                                                     #             min = 100, max = 40000, step = 100, value = 500
                                                     # )
                                                   )

                                  )

                 )
                 ),
                 tabPanel(tags$div(
                   tags$i(class = "fa-sharp fa-solid fa-magnifying-glass-chart"),
                   tags$span("5. Results"),
                   tags$style(type = "text/css", "li a{color:#000000; font-samily: 'sans-serif', Arial Rounded MT Bold;}")
                 ),
                 # results section here
                 conditionalPanel(condition = "output.seuratNotLoaded",
                                  wellPanel(
                                    h4("No data detected."),
                                    # hr(),
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
                                                   # what goes next here?
                                                   conditionalPanel(condition = "input.perturbationButton == 0",
                                                                    "not pressed"
                                                   ),
                                                   conditionalPanel(condition = "input.perturbationButton != 0",
                                                                    h4("Explore predicted expression and pharmacotranscriptomic effects"),
                                                                    hr(),
                                                                    tabsetPanel(
                                                                      tabPanel("Gene expression",
                                                                               wellPanel(
                                                                                 splitLayout(
                                                                                   wellPanel(
                                                                                     sliderInput(inputId = "cellSubsetMax_pert",
                                                                                                 label = "Select number of cells per condition to test",
                                                                                                 min = 100, max = 40000, step = 100, value = 500
                                                                                     ),
                                                                                     selectizeInput(inputId = "pertDEx_testUse", label = "Select differential expression test to use", choices = c("MAST", "roc", "wilcoxon", "deseq"), selected = "MAST"),
                                                                                     actionButton(inputId = "calcSensVsResistantDEGS", label = "Run differential expression analysis"),
                                                                                     conditionalPanel(condition = "output.resVsSensCalced",
                                                                                                      hr(),
                                                                                                      downloadButton(outputId = "resVsSensDEGSdownload", label = "Download differential expression results")
                                                                                     )

                                                                                   ),
                                                                                   wellPanel(
                                                                                     plotOutput(outputId = "resVsSensVolcano")
                                                                                   )
                                                                                 ),
                                                                                 DT::dataTableOutput(outputId = "resVsSensTable")
                                                                               )
                                                                      ),
                                                                      tabPanel("Pharmacotranscriptomics",
                                                                               wellPanel(
                                                                                 splitLayout(
                                                                                   wellPanel(
                                                                                     sliderInput(inputId = "synergyPredictorSlider",
                                                                                                 label = "select number of compounds to label",
                                                                                                 value = 30, min = 1, max = 1468, step = 1),
                                                                                   ),
                                                                                   wellPanel(
                                                                                     h4("How to interpret this plot:"),
                                                                                     img(src = "perturbationCompoundFigure.png", width = 350, height = 130),
                                                                                     hr()
                                                                                   )
                                                                                 )
                                                                               ),
                                                                               plotOutput(outputId = "scSynergySeq3_RDS")
                                                                      ),
                                                                      tabPanel("Proportion shift",
                                                                               # sliderInput(inputId = "nCompoundSlider_RDS",
                                                                               #             label = "Select number of compounds to visualize",
                                                                               #             value = 30, min = 1, max = 1468, step = 1),
                                                                               # plotOutput(outputId = "scSynergySeq2_RDS")
                                                                               wellPanel(
                                                                                 splitLayout(
                                                                                   plotOutput(outputId = "ISOSCELES_StateShiftAlluvial"),
                                                                                   "Plot of normalized shift here"
                                                                                 )
                                                                               )
                                                                      )
                                                                    )
                                                   )
                                  )
                 )
                 )
               )
               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-question"),
                 tags$span("- Tutorials")
               ), #put contents for tutorial here
               h1(strong("Tutorials")),
               hr(),
               br(),
               p(strong("Welcome to ISOSCELES! This tutorial will walk you through the ISOSCELES workflow.")),
               h2(strong("Data Upload")),
               p("Upload data as a Seurat object that contains expression data from RNAseq. File format should be .RDS or .rds. Please wait for a success message before moving to Step 2.
                 This may take a few minutes."
               ),
               h2(strong("Data Pre-Processing")),
               p("Select cell population paramaters for further study. Identify which cell population(s) you would like to set as controls and which cell population(s) you would like
                to set as transformed. These selections will affect how disease signatures are calculated (refer to next step)."),
               h2(strong("Disease and Reversal Signatures")),
               p("Now, you can calculate transcriptomic disease signatures compared to healthy controls. Use the slider bar to select how many cells you would like to test per condition. Click 'Calculate Disease Signatures'.
                Selecting more cells to test will give a more accurate signature, but it may take an extended period of time so please be patient. Once the heatmap appears, you will be able
                to download the disease signature matrix using the 'Download signature .csv' button. You can also download the heatmap as a pdf."),
               p("To study how compounds reverse the disease signature of each cell subtype, navigate to the 'Reversal' tab. Select how many compounds you would like to study. Once the heatmap appears, you will be able to
                download the corresponding matrix that provides reversal values for every cell subtype with each compound tested. To do this, click the 'Download reversal score matrix' button."),
               p("Within the reversal tab, you can use the barplot tab to graphically compare reversal scores of cell populations treated with a single compound. Select a compound of interest using the dropdown menu located under
                'Select a small molecule to plot'."),
               h2(strong("In Silico Perturbation")),
               p("To calculate a correlation matrix between L1000 compounds and RNAseq data, first select which release of the L1000 compound set you would like to use. This can be done with the dropdown menu found under 'Select L1000'
                release'. Then to calculate the correlation matrix, click on the 'Calculate single-cell compound discordance' button. When calculation is complete, you can download the matrix. Since this is a computationally-intense calculation,
                it may take extended time. Therefore, we highly recommend that you download the correlation matrix to avoid recalculation. Additionally, you have the option of uploading your own correlation matrix if you have already calculated it.
                If you would like to upload your own correlation matrix, please ensure that it is in .CSV or .csv format."),
               p("You can now select a compound of interest for synergy analysis. Under 'Step 2: Select reference compound and set parameters', select the compound to which you would like to identify synergistic compounds. Additionally,
                select computational parameters."),
               p("Click on 'Run Perturbation Analysis'. You will receive a Success message when analysis is complete, at which point you can navigate to the Results tab to study the analysis."),
               ),

               # tabPanel(tags$div(
               #   tags$i(class="fa-sharp fa-solid fa-question"),
               #   tags$span("- FAQ")
               # ), #put contents for FAQ here
               # ),

               tabPanel(tags$div(
                 tags$i(class = "fa-brands fa-github"),
                 tags$span("- Install R Package")
               ), # put contents for R package installation here
               h1("ISOSCELES R Package"),
               hr(),
               br(),
               p(strong("Due to the likely large size of many single-cell RNA seq datasets, the ISOSCELES shiny app is locally deployable via the ISOSCELES R Package. Please read the license agreement below and agree to access install instructions.")),
               h1("License Agreement"),
               # hr(),
               br(), # Note this is a draft... I dont know if this is right. Not a lawyer!
               p(strong("1. The Board of Trustees of the Georgetown University (“Georgetown”) provides ISOSCELES software and code (“Service”) free of charge for non-commercial use only. Use of the Service by any commercial entity for any purpose, including research, is prohibited.")),
               p(strong("2. By using the Service, you agree to be bound by the terms of this Agreement. Please read it carefully.")),
               p(strong("3. You agree not to use the Service for commercial advantage, or in the course of for-profit activities. You agree not to use the Service on behalf of any organization that is not a non-profit organization. Commercial entities wishing to use this Service should contact Georgetown University’s Office of Technology Licensing.")),
               p(strong("4. THE SERVICE IS OFFERED “AS IS”, AND, TO THE EXTENT PERMITTED BY LAW, GEORGETOWN MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED. GEORGETOWN SHALL NOT BE LIABLE FOR ANY CLAIMS OR DAMAGES WITH RESPECT TO ANY LOSS OR OTHER CLAIM BY YOU OR ANY THIRD PARTY ON ACCOUNT OF, OR ARISING FROM THE USE OF THE SERVICE. YOU HEREBY AGREE TO DEFEND AND INDEMNIFY GEORGETOWN, ITS TRUSTEES, EMPLOYEES, OFFICERS, STUDENTS, AGENTS, FACULTY, REPRESENTATIVES, AND VOLUNTEERS (“GEORGETOWN INDEMNITEES”) FROM ANY LOSS OR CLAIM ASSERTED AGAINST GEORGETOWN INDEMNITEES ARISING FROM YOUR USE OF THE SERVICE.")),
               p(strong("5. All rights not expressly granted to you in this Agreement are reserved and retained by GEORGETOWN or its licensors or content providers. This Agreement provides no license under any patent.")),
               p(strong("6. You agree that this Agreement and any dispute arising under it is governed by the laws of the District of Columbia, United States of America, applicable to agreements negotiated, executed, and performed within the DISTRICT OF COLUMBIA")),
               p(strong("7. Subject to your compliance with the terms and conditions set forth in this Agreement, GEORGETOWN grants you a revocable, non-exclusive, non-transferable right to access and make use of the Service.")),
               # hr(),
               br(),
               p(em("Do you accept the terms and conditions in this agreement?")),
               fluidRow(column(width = 4, actionButton(inputId = "acceptAgreement", label = "Accept")), column(width = 4, actionButton(inputId = "rejectAgreement", label = "Reject and close"))),
               conditionalPanel(condition = "input.acceptAgreement",
                                "Please scroll down for instructions."),
               # hr(),
               br(),
               conditionalPanel(condition = "input.acceptAgreement",
                                p(strong("Thank You!")),
                                h2("How to install the ISOSCELES R package"),
                                p(strong("The ISOSCELES R Package can be installed from github using devtools.")),
                                p(code("devtools::install_github('AyadLab/ISOSCELES')")),
                                p(strong("Once installed, load the ISOSCELES library.")),
                                p(code("library(ISOSCELES)")),
                                p(strong("Launch the ISOSCELES shiny UI using the RunISOSCELES() function.")),
                                p(code("ISOSCELES::runISOSCELES()"))
               )
               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-info"),
                 tags$span("- Resources")
               ), #put contents for FAQ here

               wellPanel(
                 h3("Other Resources from our laboratories:"),
                 hr(),
                 # a(href="synergyseq.com", "SynergySeq"),
                 br(),
                 # splitLayout(
                 a(href="http://www.SynergySeq.com", "SynergySeq"),
                 # h5("SynergySeq: Identify drug combinations that transcriptionally reverse a disease gene expression signature."),
                 # ),
                 br(),
                 br(),
                 # splitLayout(h4("Read the publication here: "),
                 a(href= "https://www.nature.com/articles/s41467-018-07659-z", "Nature Communications 9, 5315")
               )
               #
               # )

               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-download"),
                 tags$span("- Downloads")
               ), #put contents for FAQ here
               ),

               tabPanel(tags$div(
                 tags$i(class="fa-sharp fa-solid fa-envelope"),
                 tags$span("- Contact")
               ), #put contents for Contact here #isosceles.app@gmail.com
               splitLayout(
                 wellPanel(
                   # hr(),
                   br(),
                   p(strong("Robert K. Suter, PhD")),
                   p(em("Bioinformatician")),
                   hr(),
                   img(src = "rks_headshot_2022.png", width = 185, height = 160),
                   br(), br(),
                   a(href="mailto:rks82@georgetown.edu", "Contact..."),
                   hr()
                 ),
                 wellPanel(
                   # hr(),
                   br(),
                   p(strong("Nagi G. Ayad, PhD")),
                   p(em("Principal Investigator")),
                   hr(),
                   img(src = "nagi_headshot.png", width = 92.5, height = 160),
                   br(), br(),
                   a(href="mailto:na853@georgetown.edu", "Contact..."),
                   hr()
                 )
               )
               )
  )
)
