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
library(lme4)
library(edgeR)
library(ggpubr)
library(stringr)
library(viridis)
library(ComplexHeatmap)

server <- function(input, output) {

  # L1000 data as reactive
  L1000_genes <- reactive({
    colnames(get0("LINCS.ResponseSigs", envir = asNamespace("scFOCAL")))
  })

  LINCS.ResponseSigs <- reactive({
    get0("LINCS.ResponseSigs", envir = asNamespace("scFOCAL"))
  })

  ################################################################################
  #
  ################################################################################
  # handle upload of Seurat object
  #rdsSeurat <- reactive({
   # req(input$seurobjRDS)
   # readRDS(input$seurobjRDS$datapath)
 # })
  # handle upload of Seurat object (Dynamically handles Local vs Server)
  rdsSeurat <- reactive({
    if (input$upload_source == "local") {
      req(input$seurobjRDS)
      readRDS(input$seurobjRDS$datapath)
    } else {
      req(!is.integer(input$file_server))
      file_info <- parseFilePaths(volumes, input$file_server)
      req(nrow(file_info) > 0)
      
      file_path <- as.character(file_info$datapath[1])
      req(file.exists(file_path)) 
      
      readRDS(file_path)
    }
  })

  ######
  # SeuratObject Version Detection
  seurat_obj_version <- reactive({
    if (unlist(SeuratObject::Version(rdsSeurat()))[1] < 5){
      print("Older version Seurat object detected: (server.R #81)")
    } else {
      print("V5 Seurat object detected: (server.R #83)")
      return(NULL)
    }
  })

  output$seurat_obj_v5 <- reactive({
    return(is.null(seurat_obj_version()))
  })
  outputOptions(output, 'seurat_obj_v5', suspendWhenHidden = FALSE)

  #####
  # OrthogonAL Detection - test if Seurat object is output from OrthogonAL
  isOrthogonAL <- reactive({
    RDSseurat <- rdsSeurat()
    assays <- names(RDSseurat@assays)
    print(assays)
    if ("RNA_ortho" %in% assays){
      print("Object detected is from OrthogonAL")
      "Object detected is from OrthogonAL"
    } else {
      print("not OrthogonAL object")
      return(NULL)
    }
  })

  output$OrthogonAL <- reactive({
    print(c(isOrthogonAL()))
    return(is.null(isOrthogonAL()))
  })
  outputOptions(output, 'OrthogonAL', suspendWhenHidden = FALSE)

  #####

  sl <- reactive({
    if (class(rdsSeurat()) == "Seurat"){
      return(NULL)
    }
  })

  output$seuratLoaded <- reactive({
    return(is.null(sl()))
  })
  outputOptions(output, 'seuratLoaded', suspendWhenHidden = FALSE)

  output$seuratNotLoaded <- reactive({
    return(is.null(input$seurobjRDS$datapath))
  })

  outputOptions(output, 'seuratNotLoaded', suspendWhenHidden = FALSE)
  # --- UPDATED: Dynamic check for whether a file has been loaded yet ---
  output$seuratNotLoaded <- reactive({
    if (input$upload_source == "local") {
      return(is.null(input$seurobjRDS$datapath))
    } else {
      return(is.integer(input$file_server))
    }
  })
  outputOptions(output, 'seuratNotLoaded', suspendWhenHidden = FALSE)
  # Custom TCS Upload
  ##############################################################################
  # handle upload of Seurat object
  customTCScsv <- reactive({
    req(input$customTCScsvPath)
    read.csv(file = input$customTCScsvPath$datapath, header = T)
  })

  ######
  # SeuratObject Version Detection
  customTCSvalidation <- reactive({ # need to add in functionality to test for correct csv format...
    # if (class(customTCScsv()) == "data.frame"){
    #   print("Custom TCS successfully loaded")
    # } else {
    #   print("TEST")
    #   return(NULL)
    # }
    req(input$customTCScsvPath$datapath)
    return(NULL)
  })

  output$TCSuploaded <- reactive({
    return(is.null(customTCSvalidation()))
  })
  outputOptions(output, 'TCSuploaded', suspendWhenHidden = FALSE)

  # modal - browse TCS upload

  output$customTCS_table <- DT::renderDataTable({
    DT::datatable(customTCScsv())
  })

  observeEvent(input$customTCSbrowser_button, {
    showModal(modalDialog(
      title = "Custom TCS Browser",
      tags$head(tags$style(".modal-dialog{ width:85% !important;}")),
      DT::dataTableOutput("customTCS_table"),
      br(),
      "plotSomethingHere",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  ##############################################################################

  # Handle upload of pre-calculated correlation matrix
  corrMatUpload <- reactive({
    req(input$corrMatUpload)
    read.csv(input$corrMatUpload$datapath, row.names = 1, header = T)
  })

  cm <- reactive({
    if (class(corrMatUpload()) == "data.frame"){
      return(NULL)
    }
  })

  output$corrMatUploaded <- reactive({
    return(is.null(cm()))
  })
  outputOptions(output, 'corrMatUploaded', suspendWhenHidden = FALSE)

  output$corrMatNotUploaded <- reactive({
    return(is.null(input$corrMatUpload$datapath))
  })
  outputOptions(output, 'corrMatNotUploaded', suspendWhenHidden = FALSE)

  output$L1000_release_InSilico <- renderUI({ # select from dropdown which L1000 release to use for in silico perturbation:
    selectInput('L1000_Release', #input$groupByRDS
                choices = c("2015", "2017", "2021"),
                label = "Select L1000 Release",
                selected = "2017",
                multiple = F)
  })

  output$reductionUseRDS <- renderUI({ # select from dropdown which reduction loadings to plot across app...
    RDSseurat <- rdsSeurat()
    reductionsRDS <- Reductions(RDSseurat)
    selectInput('reductionUseRDS',
                choices = unique(reductionsRDS),
                label = "Select reduction to plot data throughout app...",
                selected = "umap",
                multiple = F)
  })

  output$groupByRDS <- renderUI({ # select from dropdown which metadata slot to group cells by:
    RDSseurat <- rdsSeurat() # this can't run until upload is a valid seurat object in the environment
    metadata_slots_RDS <- colnames(RDSseurat@meta.data)
    selectInput('groupByRDS', #input$groupByRDS
                choices = unique(metadata_slots_RDS),
                label = "Select attribute to group by.",
                selected = NULL,
                multiple = F)
  })

  output$subject_replicate_ident <- renderUI({
    RDSseurat <- rdsSeurat()
    metadata_slots_RDS <- colnames(RDSseurat@meta.data)
    selectInput('subject_replicate_ident',
                choices = unique(metadata_slots_RDS),
                label = "Select subject/replicate metadata slot",
                selected = NULL,
                multiple = F)
  })

  output$normalCellIdentsRDS <- renderUI({
    RDSseurat <- rdsSeurat()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
    print(unique(Idents(RDSseurat)))
    selectInput('normalCellIdentsRDS',
                choices = unique(Idents(RDSseurat)),
                label = "Select control cell population(s).",
                selected = NULL,
                multiple = T
    )
  })

  output$dimplotRDS <- renderPlot({
    p1 <- Seurat::DimPlot(rdsSeurat(), reduction = input$reductionUseRDS, group.by = input$groupByRDS) +
      theme_void() +
      theme(plot.title = element_text(size = 18, face = "bold"))
    plot_grid(plotlist = list(p1))
  })

  output$dimplotRDS_subject <- renderPlot({
    p1 <- Seurat::DimPlot(rdsSeurat(), reduction = input$reductionUseRDS, group.by = input$subject_replicate_ident) + ###
      theme_void() +
      theme(plot.title = element_text(size = 18, face = "bold"))
    plot_grid(plotlist = list(p1))
  })

  output$cancerCellIdentsRDS <- renderUI({
    RDSseurat <- rdsSeurat()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
    selectInput('cancerCellIdentsRDS',
                choices = unique(Idents(RDSseurat)),
                label = "Which cells are transformed?",
                selected = NULL,
                multiple = T
    )
  })

  output$controlHighlight <- renderPlot({
    if (!is.null(input$normalCellIdentsRDS)){
      RDSseurat <- rdsSeurat()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      DimPlot(object = RDSseurat,
              reduction = input$reductionUseRDS,
              cells.highlight = WhichCells(RDSseurat,
                                           idents = input$normalCellIdentsRDS),
              cols.highlight = c("red")) +
        theme_void() +
        ggtitle("Control Population") +
        theme(plot.title = element_text(size = 18, face = "bold"))
    }
  })

  output$diseasedHighlight <- renderPlot({
    if (!is.null(input$cancerCellIdentsRDS)){
      RDSseurat <- rdsSeurat()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      DimPlot(object = RDSseurat,
              reduction = input$reductionUseRDS,
              cells.highlight = WhichCells(RDSseurat,
                                           idents = input$cancerCellIdentsRDS),
              cols.highlight = c("navy")) +
        theme_void() +
        ggtitle("Diseased Population") +
        theme(plot.title = element_text(size = 18, face = "bold"))
    }
  })

  # Control display of disease signature calculation controls for population id in preprocessing

  ps <- reactive({
    if (!is.null(input$normalCellIdentsRDS) == TRUE){
      if (!is.null(input$cancerCellIdentsRDS) == TRUE){
        return(TRUE)
      }
    }
    # return(!is.null(input$normalCellIdentsRDS))
  })

  ps_noControl <- reactive({
    if (input$noNormalCellIdentsCheckbox == TRUE){
      if (!is.null(input$cancerCellIdentsRDS) == TRUE){
        return(TRUE)
      }
    }
  })

  output$controlsSelected <- reactive({
    return(!is.null(ps()))
  })
  output$controlsNotSelected <- reactive({
    return(is.null(ps()))
  })

  output$controlsSelected_noControl <- reactive({
    return(!is.null(ps_noControl()))
  })
  output$controlsNotSelected_noControl <- reactive({
    return(is.null(ps_noControl()))
  })

  outputOptions(output, 'controlsSelected', suspendWhenHidden = FALSE)
  outputOptions(output, 'controlsNotSelected', suspendWhenHidden = FALSE)
  outputOptions(output, 'controlsSelected_noControl', suspendWhenHidden = FALSE)
  outputOptions(output, 'controlsNotSelected_noControl', suspendWhenHidden = FALSE)

  ##############################################################################
  #                                                                            #
  ##############################################################################

  output$splitBySelecter <- renderUI({
    selectInput('sliceBy',
                choices = unique(colnames(testdat@meta.data)),
                label = "Select how to slice your data into multiple datasets to test, i.e. by patient",
                selected = "orig.ident",
                multiple = F
    )
  })

  output$splitBySelecterRDS <- renderUI({
    RDSseurat <- rdsSeurat()
    selectInput('sliceByRDS',
                choices = unique(colnames(RDSseurat@meta.data)),
                label = "Select how to slice your data into multiple datasets to test, i.e. by patient (RDS)",
                selected = NULL,
                multiple = F
    )
  })

  # Marker Calculation
  output$markers <- renderPlot({
    if (input$CalcMarkers == T && input$dataChoice == "10X_GBM"){
      Idents(testdat) <- input$groupBy
      allMarkers <- FindAllMarkers(testdat, test.use = "MAST", only.pos = T)
      top5markers <- allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
      p1 <- Seurat::DoHeatmap(testdat, features = top5markers$gene)
      plot_grid(plotlist = list(p1))
    }
  })

  output$markersRDS <- renderPlot({
    if (input$CalcMarkers == T && input$dataChoice == "RDS_Upload"){
      RDSseurat <- rdsSeurat()
      Idents(RDSseurat) <- input$groupByRDS
      allMarkers <- FindAllMarkers(RDSseurat, test.use = "MAST", only.pos = T)
      top5markers <- allMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
      p2 <- Seurat::DoHeatmap(RDSseurat, features = top5markers$gene, slot = "data")
      plot_grid(plotlist = list(p2))
    }
  })

  # Disease Signature Calculation

  # Unsliced Disease Signatures

  data <- reactiveVal()
  data(data.frame())

  observeEvent(input$CalcDiseaseSig, {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Running MAST", value = 0)
    diseaseSigs <- data()
    tumor_markers <- data.frame()
    RDSseurat <- rdsSeurat()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)

    if (is.null(isOrthogonAL()) == TRUE){ # is this backwards? ############################################################################################# <
      print("Is null = F, non-orthogonal")
    } else {
      print("Is.null = T, OrthogonalObject")
      DefaultAssay(RDSseurat) <- c("RNA_ortho")
    }

    if (input$diseaseSignaturesByReplicate == 'By Subject/Replicate') {
      print("would calc by replicate here!!")
      # subject_replicate_ident

      for (replicate in unique(RDSseurat@meta.data[,c(input$subject_replicate_ident)])){
        print(replicate)
        RDSseurat <- SetIdent(RDSseurat, value = input$subject_replicate_ident)
        sub_obj <- subset(RDSseurat, idents = c(replicate))
        sub_obj <- SetIdent(sub_obj, value = input$groupByRDS)
        print(table(Idents(sub_obj)))
        for (i in 1:length(input$cancerCellIdentsRDS)){
          if (input$cancerCellIdentsRDS[i] %in% unique(sub_obj@meta.data[,input$groupByRDS])){

            ident_counts <- table(Idents(sub_obj))
            print(names(ident_counts))
            if (ident_counts[input$cancerCellIdentsRDS[i]] > 3) {
              # Proceed with the analysis for this identity
              cat("Processing identity:", input$cancerCellIdentsRDS[i], "with", ident_counts[input$cancerCellIdentsRDS[i]], "cells.\n")
              progress$inc(1/(length(input$cancerCellIdentsRDS)*length(unique(RDSseurat@meta.data[,c(input$subject_replicate_ident)]))),
                           detail = paste("Running Cluster ", input$cancerCellIdentsRDS[i], " of subject ", replicate))
              clusterDiseaseSig <- FindMarkers(sub_obj,
                                               ident.1 = input$cancerCellIdentsRDS[i],
                                               ident.2 = input$normalCellIdentsRDS,
                                               test.use = "MAST",
                                               only.pos = F,
                                               max.cells.per.ident = input$cellSubsetMax_1) # Can use max.cells.per.ident to set a slider for option to speed up analysis.
              clusterDiseaseSig$cluster <- input$cancerCellIdentsRDS[i]
              clusterDiseaseSig$gene <- rownames(clusterDiseaseSig)
              clusterDiseaseSig$replicate <- replicate
              tumor_markers <- rbind(tumor_markers, clusterDiseaseSig)
            } else {
              # Skip this identity due to insufficient cell count
              cat("Skipping identity:", input$cancerCellIdentsRDS[i], "with only", ident_counts[input$cancerCellIdentsRDS[i]], "cells.\n")
            }


          }

        }

      }

    }
    if (input$diseaseSignaturesByReplicate == 'Overall') {
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      for (i in 1:length(input$cancerCellIdentsRDS)){
        progress$inc(1/length(input$cancerCellIdentsRDS), detail = paste("Running Cluster", input$cancerCellIdentsRDS[i]))
        clusterDiseaseSig <- FindMarkers(RDSseurat,
                                         ident.1 = input$cancerCellIdentsRDS[i],
                                         ident.2 = input$normalCellIdentsRDS,
                                         test.use = "MAST",
                                         only.pos = F,
                                         max.cells.per.ident = input$cellSubsetMax_1
        ) # Can use max.cells.per.ident to set a slider for option to speed up analysis.
        clusterDiseaseSig$cluster <- input$cancerCellIdentsRDS[i]
        clusterDiseaseSig$gene <- rownames(clusterDiseaseSig)
        tumor_markers <- rbind(tumor_markers, clusterDiseaseSig)
      }
    }

    print(head(tumor_markers))
    data(tumor_markers)
    # }
  })

  diseaseSigs <- reactive({
    data()
  })

  observeEvent(input$CalcDiseaseSig, {
    print(head(rownames(diseaseSigs())))
  })


  # output$diseaseSig <- renderPlot({
  #   if (input$CalcDiseaseSig == T && input$dataChoice == "10X_GBM")
  #     top_dis_markers <- diseaseSigs() %>% group_by(cluster) %>% top_n(wt = avg_log2FC, n = 10)
  #   sigHeatMap <- DoHeatmap(testdat, features = top_dis_markers$gene)
  #   sigHeatMap
  # })

  output$diseaseSigRDS <- renderPlot({
    if (input$CalcDiseaseSig == T) {
      RDSseurat <- rdsSeurat()
      top_dis_markers <- diseaseSigs() %>% group_by(cluster) %>% top_n(wt = avg_log2FC, n = 10)
      sigHeatMap <- DoHeatmap(RDSseurat, features = top_dis_markers$gene, group.by = input$groupByRDS) + scale_fill_gradient2(low = muted("navy"), mid = "white", high = muted("red"))
      sigHeatMap
    }
  })

  # output$diseaseSigDownload <- downloadHandler(
  #   filename = function() {
  #     paste("DiseaseSignature", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(testdat, file)
  #   })

  output$diseaseSigTable <- DT::renderDataTable({
    diseaseSigs()
  })

  output$diseaseSigDownload <- downloadHandler(
    filename = function() {
      paste(input$cancerCellIdents, "_DiseaseSignatures.csv", sep = "")
    },
    content = function(file) {
      write.csv(diseaseSigs(), file, row.names = FALSE)
    }
  )

  ### Whole dataset signature reversal
  # Here - add in scripts to rank disordance with calculated sliced disease signatures.

  walincs <- reactiveVal()
  wblincs <- reactiveVal()
  wreversalMat <- reactiveVal()
  observeEvent(input$reverseDiseaseSigs, {

    if (input$diseaseSignaturesByReplicate == 'By Subject/Replicate') {
      diseaseSigs <- diseaseSigs()
      replicates <- unique(diseaseSigs$replicate)
      rep_final_matrix <- data.frame()
      for (rep in replicates){
        print(paste("Reversing signatures from", rep, sep = " "))
        rep_dS <- subset(diseaseSigs, diseaseSigs$replicate == rep)
        u <- reshape2::dcast(rep_dS, formula = cluster ~ gene, value.var = "avg_log2FC", fill = 0, fun.aggregate = mean)
        rownames(u) <- u$cluster
        #   u$group <- NULL
        u <- as.data.frame(t(u))

        # QC - that as total DEGs increases so does the number within L1000 geneset.
        LINCS_DEGs <- rep_dS[which(rep_dS$gene %in% L1000_genes()),]
        a <- as.vector(table(rep_dS$cluster))
        b <- as.vector(table(LINCS_DEGs$cluster))
        alincs(a)
        blincs(b)

        Final_Matrix <- data.frame()
        print(unique(LINCS_DEGs$cluster))
        for (i in 1:length(unique(LINCS_DEGs$cluster))) {
          vect1 <- LINCS_DEGs[which(LINCS_DEGs$cluster == unique(LINCS_DEGs$cluster)[[i]]), "avg_log2FC"]
          vect1_Genes <- as.character(LINCS_DEGs[which(LINCS_DEGs$cluster == unique(LINCS_DEGs$cluster)[[i]]), "gene"])
          names(vect1) <- vect1_Genes
          setdiff(vect1_Genes, colnames(LINCS.ResponseSigs()))
          LINCS.ResponseSigs.Filtered <- LINCS.ResponseSigs()[,vect1_Genes]
          Ranked_List <- apply(LINCS.ResponseSigs.Filtered, 1, function(x){
            mult <- vect1 * x
            a <- length(mult[mult>0])
            b <- length(mult[mult<0])
            if (a==0){a <- 1}
            c <- b / a
          })

          Final_Matrix <- rbind(Final_Matrix, Ranked_List)


          clusters <- unique(LINCS_DEGs$cluster)

          colnames(Final_Matrix) <- names(Ranked_List)
          rownames(Final_Matrix) <- clusters[1:i]
          rownames(Final_Matrix) <- paste0(rownames(Final_Matrix), "_", rep)
        }

        rep_final_matrix <- rbind(rep_final_matrix, Final_Matrix)
        print("final replicate reversal mat")
        print(head(rep_final_matrix))
        wreversalMat(rep_final_matrix)

      }
    } else {

      # Here - run processing needed to reverse signatures. Can we add in a progress meter as well?
      sliced_signatures2 <- diseaseSigs()
      # sliced_signatures2$group <- paste0(sliced_signatures2$slice, "_", sliced_signatures2$cluster)
      #
      u <- reshape2::dcast(sliced_signatures2, formula = cluster ~ gene, value.var = "avg_log2FC", fill = 0, fun.aggregate = mean)
      rownames(u) <- u$cluster
      #   u$group <- NULL
      u <- as.data.frame(t(u))

      # QC - that as total DEGs increases so does the number within L1000 geneset.
      LINCS_DEGs <- sliced_signatures2[which(sliced_signatures2$gene %in% L1000_genes()),]
      a <- as.vector(table(sliced_signatures2$cluster))
      b <- as.vector(table(LINCS_DEGs$cluster))
      alincs(a)
      blincs(b)

      Final_Matrix <- data.frame()
      print(unique(LINCS_DEGs$cluster))
      for (i in 1:length(unique(LINCS_DEGs$cluster))) {
        vect1 <- LINCS_DEGs[which(LINCS_DEGs$cluster == unique(LINCS_DEGs$cluster)[[i]]), "avg_log2FC"]
        vect1_Genes <- as.character(LINCS_DEGs[which(LINCS_DEGs$cluster == unique(LINCS_DEGs$cluster)[[i]]), "gene"])
        names(vect1) <- vect1_Genes
        setdiff(vect1_Genes, colnames(LINCS.ResponseSigs()))
        LINCS.ResponseSigs.Filtered <- LINCS.ResponseSigs()[,vect1_Genes]
        Ranked_List <- apply(LINCS.ResponseSigs.Filtered, 1, function(x){
          mult <- vect1 * x
          a <- length(mult[mult>0])
          b <- length(mult[mult<0])
          if (a==0){a <- 1}
          c <- b / a
        })

        Final_Matrix <- rbind(Final_Matrix, Ranked_List)


        clusters <- unique(LINCS_DEGs$cluster)

        colnames(Final_Matrix) <- names(Ranked_List)
        rownames(Final_Matrix) <- clusters[1:i]
      }

      wreversalMat(Final_Matrix)
    }



  })

  # Set up download button for reversalMat()
  output$reversedDiseaseSigDL <- downloadHandler(
    filename = function() {
      paste(unique(input$cancerCellIdents), "_DiseaseSignature", "ReversalMatrix.csv", sep = "_")
    },
    content = function(file) {
      write.csv(wreversalMat(), file, row.names = T)
    }
  )


  output$DiseaseSigReversalQC <- renderPlot({
    if (input$reverseDiseaseSigs == T){
      a <- alincs()
      b <- blincs()
      plot(a, b)
    }
  })

  # the matrix needs rownames I think.
  output$DiseaseSigReversalHeatmap <- renderPlot({
    if (input$reverseDiseaseSigs == T) {
      Final_Matrix <- wreversalMat()
      Final_Matrix <- na.omit(t(Final_Matrix))
      var <- apply(Final_Matrix, 1, var) # Calculate variance
      variableCompounds <- names(var[order(var, decreasing = T)][1:input$varNum])
      pheatmap(Final_Matrix[variableCompounds,],
               scale = "row", show_rownames = F, color = viridis(256, option = "D"), fontsize = 24)#, cluster_rows = F, cluster_cols = F)
    }
  })

  output$DiseaseSigReversalBarplot <- renderPlot({
    Final_Matrix <- wreversalMat()
    Final_Matrix$group <- rownames(Final_Matrix)

    p <- ggplot(
      Final_Matrix,
      aes_string(x = "group", y = input$wcompoundToBarplot, fill = "group")) +
      geom_col(position = "identity", colour = "black", size = 0.25) +
      scale_fill_viridis(discrete = T, option = "C") +
      ggtitle(paste("L1000 Derived Consensus Signature",
                    "_",
                    input$wcompoundToBarplot)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    p
  })


  ########################################################################################
  # Sliced disease signatures
  ########################################################################################

  data2 <- reactiveVal()
  data2(data.frame())
  data3 <- reactiveVal()
  data3(data.frame())

  observeEvent(input$calcSlicedDiseaseSigs, {


    slicedDiseaseSigs <- data2()

    tumor_markers2 <- data.frame()
    tumor_markers3 <- data.frame()

    progress4 <- shiny::Progress$new()
    progress4$set(message = "Splitting your data")
    splitData <- Seurat::SplitObject(testdat, split.by = input$sliceBy)
    progress4$close()

    progress2 <- shiny::Progress$new()
    on.exit(progress2$close())
    progress2$set(message = "Running MAST on sliced cell populations", value = 0)

    for (i in 1:length(splitData)){
      progress2$inc(1/length(splitData), detail = paste("Running dataset", FetchData(splitData[[i]], vars = input$sliceBy)[1,1]))
      currObj <- splitData[[i]]
      groupByPresent <- intersect(unique(currObj@meta.data[, input$groupBy]),
                                  input$cancerCellIdents)


      # And what about the splitting takes so long????
      progress3 <- shiny::Progress$new()
      progress3$set(message = "Running MAST")

      for (j in 1:length(groupByPresent)){
        progress3$inc(1/length(groupByPresent), detail = paste("Running Cluster", groupByPresent[j]))
        clusterDiseaseSig <- FindMarkers(object = currObj,
                                         ident.1 = groupByPresent[j],
                                         ident.2 = input$normalCellIdents,
                                         test.use = "MAST",
                                         only.pos = T,
                                         max.cells.per.ident = input$cellSubsetMax#,
                                         # logfc.threshold = -1
        ) # Can use max.cells.per.ident to set a slider for option to speed up analysis.
        clusterDiseaseSig$cluster <- groupByPresent[j]
        clusterDiseaseSig$slice <- currObj@meta.data[1, input$sliceBy]
        clusterDiseaseSig$gene <- rownames(clusterDiseaseSig)
        tumor_markers2 <- rbind(tumor_markers2, clusterDiseaseSig)
      }
      progress3$close()
      data2(tumor_markers2)
      tumor_markers3 <- rbind(tumor_markers3, tumor_markers2)
      rm(currObj)
    }
    data3(tumor_markers3)

  })

  data4 <- reactiveVal()
  data4(data.frame())
  data5 <- reactiveVal()
  data5(data.frame())

  observeEvent(input$calcSlicedDiseaseSigsRDS, {

    RDSseurat <- rdsSeurat()
    slicedDiseaseSigsRDS <- data4()

    tumor_markers2 <- data.frame()
    tumor_markers3 <- data.frame()

    progress8 <- shiny::Progress$new()
    progress8$set(message = "Splitting your data")
    Idents(RDSseurat) <- input$sliceByRDS
    splitData <- Seurat::SplitObject(RDSseurat, split.by = input$sliceByRDS)
    progress8$close()

    progress9 <- shiny::Progress$new()
    on.exit(progress9$close())
    progress9$set(message = "Running MAST on sliced cell populations", value = 0)

    for (i in 1:length(splitData)){
      progress2$inc(1/length(splitData), detail = paste("Running dataset", FetchData(splitData[[i]], vars = input$sliceByRDS)[1,1]))
      currObj <- splitData[[i]]
      print(currObj)
      groupByPresent <- intersect(unique(currObj@meta.data[, input$groupByRDS]),
                                  input$cancerCellIdentsRDS)


      # And what about the splitting takes so long????
      progress10 <- shiny::Progress$new()
      progress10$set(message = "Running MAST")
      print(groupByPresent) #
      for (j in 1:length(groupByPresent)){
        print(j)
        progress10$inc(1/length(groupByPresent), detail = paste("Running Cluster", groupByPresent[j]))
        print("Just before FindMarkers...")
        clusterDiseaseSig <- FindMarkers(object = currObj,
                                         ident.1 = groupByPresent[j],
                                         ident.2 = input$normalCellIdentsRDS,
                                         test.use = "MAST",
                                         only.pos = T,
                                         max.cells.per.ident = input$cellSubsetMaxRDS)# ,
        # logfc.threshold = -1
        # ) # Can use max.cells.per.ident to set a slider for option to speed up analysis.
        print(groupByPresent[j])
        clusterDiseaseSigRDS$cluster <- groupByPresent[j]
        print(currObj@meta.data[1, input$sliceByRDS])
        clusterDiseaseSigRDS$slice <- currObj@meta.data[1, input$sliceByRDS]
        clusterDiseaseSigRDS$gene <- rownames(clusterDiseaseSigRDS)
        tumor_markers2 <- rbind(tumor_markers2, clusterDiseaseSigRDS)
      }
      progress10$close()
      data4(tumor_markers2)
      tumor_markers3 <- rbind(tumor_markers3, tumor_markers2)
      rm(currObj)
    }
    data5(tumor_markers3)

  })

  slicedDiseaseSigs <- reactive({
    data3()
  })

  slicedDiseaseSigsRDS <- reactive({
    data5()
  })

  output$slicedDiseaseSig <- renderPlot({
    if (input$calcSlicedDiseaseSigs == T){
      top_dis_markers_sliced <- slicedDiseaseSigs() %>% group_by(slice) %>% group_by(cluster) %>% top_n(wt = avg_log2FC, n = 10)
      sigHeatMap <- DoHeatmap(testdat, features = top_dis_markers_sliced$gene)
      sigHeatMap
    }
  })

  output$slicedDiseaseSigRDS <- renderPlot({
    RDSseurat <- rdsSeurat()
    if (input$calcSlicedDiseaseSigsRDS == T){
      top_dis_markers_sliced <- slicedDiseaseSigsRDS() %>% group_by(slice) %>% group_by(cluster) %>% top_n(wt = avg_log2FC, n = 10)
      sigHeatMap <- DoHeatmap(RDSseurat, features = top_dis_markers_sliced$gene)
      sigHeatMap
    }
  })

  output$slicedDiseaseSigTable <- renderTable({
    slicedDiseaseSigs()
  })

  output$slicedDiseaseSigTableRDS <- renderTable({
    slicedDiseaseSigsRDS()
  })

  output$slicedDiseaseSigDownload <- downloadHandler(
    filename = function() {
      paste(input$sliceBy, input$cancerCellIdents, "_DiseaseSignatures.csv", sep = "")
    },
    content = function(file) {
      write.csv(slicedDiseaseSigs(), file, row.names = FALSE)
    }
  )

  output$slicedDiseaseSigDownloadRDS <- downloadHandler(
    filename = function() {
      paste(input$sliceByRDS, input$cancerCellIdentsRDS, "_DiseaseSignatures.csv", sep = "")
    },
    content = function(file) {
      write.csv(slicedDiseaseSigsRDS(), file, row.names = FALSE)
    }
  )

  # Here - add in plot output for correlation matrices.

  output$slicedDiseaseSigCorrMatPlot <- renderPlot({
    sliced_signatures <- slicedDiseaseSigs()
    sliced_signatures$group <- paste0(sliced_signatures$slice, "_", sliced_signatures$cluster)
    t <- reshape2::dcast(sliced_signatures, formula = group ~ gene, value.var = "avg_log2FC", fill = 0, fun.aggregate = mean)
    rownames(t) <- t$group
    t$group <- NULL
    t <- as.data.frame(t(t))
    res <- cor(t, method = "pearson")

    corrMat <- pheatmap(res,
                        color = viridis(256, option = "D"),
                        # height = 4,
                        # width = 4,
                        border_color = NA,
                        # cellwidth = 12,
                        # cellheight = 12,
                        silent = T)
    corrMat
  })

  # Here - add in scripts to rank disordance with calculated sliced disease signatures.

  alincs <- reactiveVal()
  blincs <- reactiveVal()
  reversalMat <- reactiveVal()
  observeEvent(input$reverseSlicedDiseaseSigs, {
    # Here - run processing needed to reverse signatures. Can we add in a progress meter as well?
    sliced_signatures2 <- slicedDiseaseSigs()
    sliced_signatures2$group <- paste0(sliced_signatures2$slice, "_", sliced_signatures2$cluster)

    u <- reshape2::dcast(sliced_signatures2, formula = group ~ gene, value.var = "avg_log2FC", fill = 0, fun.aggregate = mean)
    rownames(u) <- u$group
    u$group <- NULL
    u <- as.data.frame(t(u))

    # QC - that as total DEGs increases so does the number within L1000 geneset.
    LINCS_DEGs <- sliced_signatures2[which(sliced_signatures2$gene %in% L1000_genes()),]
    a <- as.vector(table(sliced_signatures2$group))
    b <- as.vector(table(LINCS_DEGs$group))
    alincs(a)
    blincs(b)

    Final_Matrix <- data.frame()
    for (i in 1:length(unique(LINCS_DEGs$group))) {
      vect1 <- LINCS_DEGs[which(LINCS_DEGs$group == unique(LINCS_DEGs$group)[[i]]), "avg_log2FC"]
      vect1_Genes <- as.character(LINCS_DEGs[which(LINCS_DEGs$group == unique(LINCS_DEGs$group)[[i]]), "gene"])
      names(vect1) <- vect1_Genes
      setdiff(vect1_Genes, colnames(LINCS.ResponseSigs))
      LINCS.ResponseSigs.Filtered <- LINCS.ResponseSigs[,vect1_Genes]
      Ranked_List <- apply(LINCS.ResponseSigs.Filtered, 1, function(x){
        mult <- vect1 * x
        a <- length(mult[mult>0])
        b <- length(mult[mult<0])
        if (a==0){a <- 1}
        c <- b / a
      })

      Final_Matrix <- rbind(Final_Matrix, Ranked_List)


      groups <- unique(LINCS_DEGs$group)

      colnames(Final_Matrix) <- names(Ranked_List)
      rownames(Final_Matrix) <- groups[1:i]
    }

    reversalMat(Final_Matrix)

  })

  # Set up download button for reversalMat()
  output$reversedSlicedDiseaseSigDL <- downloadHandler(
    filename = function() {
      paste(unique(input$cancerCellIdents), "_DiseaseSignature", "ReversalMatrix.csv", sep = "_")
    },
    content = function(file) {
      write.csv(reversalMat(), file, row.names = T)
    }
  )


  output$slicedDiseaseSigReversalQC <- renderPlot({
    if (input$reverseSlicedDiseaseSigs == T){
      a <- alincs()
      b <- blincs()
      plot(a, b)
    }
  })

  # the matrix needs rownames I think.
  output$slicedDiseaseSigReversalHeatmap <- renderPlot({
    if (input$reverseSlicedDiseaseSigs == T) {
      Final_Matrix <- reversalMat()
      Final_Matrix <- na.omit(t(Final_Matrix))
      var <- apply(Final_Matrix, 1, var) # Calculate variance
      variableCompounds <- names(var[order(var, decreasing = T)][1:input$varNum])
      pheatmap(Final_Matrix[variableCompounds,],
               scale = "row", show_rownames = F, color = viridis(256, option = "D"))#, cluster_rows = F, cluster_cols = F)
    }
  })

  output$slicedDiseaseSigReversalBarplot <- renderPlot({
    Final_Matrix <- reversalMat()
    Final_Matrix$group <- rownames(Final_Matrix)

    p <- ggplot(
      Final_Matrix,
      aes_string(x = "group", y = input$compoundToBarplot, fill = "group")) +
      geom_col(position = "identity", colour = "black", size = 0.25) +
      scale_fill_viridis(discrete = T, option = "C") +
      ggtitle(paste("L1000 Derived Consensus Signature",
                    "_",
                    input$compoundToBarplot)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    p
  })

  ##############################################################################
  #                                                                            #
  ##############################################################################

  # Calculate correlation matrix between L1000 compounds and single-cells in RDS upload.
  RDS_Final_CorrMat <- reactiveVal()

  observeEvent(eventExpr = input$CalculateRDS_L1000_Spearman_Mat, { # Here, we should add dropdown menu to select data slot to pull scale.data from...
    RDSseurat <- rdsSeurat()
    if (unlist(SeuratObject::Version(RDSseurat))[1] < 5) {
      if (dim(RDSseurat@assays$RNA@scale.data)[1] == 0){
        print("scale.data slot is empty - scaling data")
        RDSseurat <- ScaleData(RDSseurat, do.center = T)
      }
    } else {
      print("V5 Seurat Object Detected: server.R #853")
      if (is.null(isOrthogonAL()) == TRUE){ # is this backwards? ############################################################################################# <
        print("Is null = F, non-orthogonal")
        # if (dim(RDSseurat@assays$RNA@layers$scale.data)[1] == 0){
        #   print("scale.data slot is empty - scaling data")
        #   RDSseurat <- ScaleData(RDSseurat, do.center = T)
        # }
      } else {
        print("Is.null = T, OrthogonalObject")
        DefaultAssay(RDSseurat) <- c("RNA_ortho")
        # if (dim(RDSseurat@assays$RNA_ortho@layers$scale.data)[1] == 0){
        #   print("scale.data slot is empty - scaling data")
        #   RDSseurat <- ScaleData(RDSseurat, do.center = T)
        # }
      }
      # DefaultAssay(RDSseurat) <- "RNA_ortho"
      # if (dim(RDSseurat@assays$RNA_ortho@layers$scale.data)[1] == 0){
      #   print("scale.data slot is empty - scaling data")
      #   RDSseurat <- ScaleData(RDSseurat, do.center = T)
      # }
    }

    a <- LINCS.ResponseSigs()
    a$compound <- rownames(a)
    progress5 <- shiny::Progress$new()
    progress5$set(message = "Transposing scale.data slot")

    # if (is.null(isOrthogonAL()) == FALSE){ # is this backwards?
    #   print("Is null = F, orthogonal")
    # } else {
    #   print("Is.null = T, non-OrthogonalObject")
    #   DefaultAssay(RDSseurat) <- c("RNA_ortho")
    # }

    if (is.null(isOrthogonAL()) == TRUE){ # is this backwards? ############################################################################################# <
      print("Is null = F, non-orthogonal")
    } else {
      print("Is.null = T, OrthogonalObject")
      DefaultAssay(RDSseurat) <- c("RNA_ortho")
    }
    ##########################################################################################################################################
    if (unlist(SeuratObject::Version(RDSseurat))[1] < 5){
      print("line 909: detected as less than v5")
      # if (is.null(isOrthogonAL()) == FALSE){ # is this backwards?
      #   print("Is null = F, non-orthogonal")
      # } else {
      #   print("Is.null = T, OrthogonalObject")
      #   DefaultAssay(RDSseurat) <- c("RNA_ortho")
      # }
      total.transpose <- t(RDSseurat@assays$RNA@scale.data) # which assay to integrate with?
    } else {
      print("We're here 927")
      if (is.null(isOrthogonAL()) == TRUE){ # is this backwards?
        print("Is.null = T, non-OrthogonalObject")
        # DefaultAssay(RDSseurat) <- c("RNA_ortho")
        # total.transpose <- t(RDSseurat@assays$RNA_ortho@layers$scale.data) # which assay to integrate with?
        total.transpose <- t(GetAssayData(RDSseurat, assay = "RNA", layer = "scale.data"))
      } else {
        print("Is null = F, orthogonal")
        # total.transpose <- t(RDSseurat@assays$RNA@layers$scale.data) # which assay to integrate with?
        total.transpose <- t(GetAssayData(RDSseurat, assay = "RNA_ortho", layer = "scale.data"))
      }
    }
    ##################################################################################################################################################
    progress5$close()

    progress6 <- shiny::Progress$new()
    on.exit(progress6$close())
    progress6$set(message = "Calculating Drug-Cell Connectivity", value = 0)

    # loop through compounds
    counter <- 1
    Final_Matrix <- data.frame()
    for (i in 1:length(rownames(LINCS.ResponseSigs()))){
      cmpdToOverlay <- rownames(LINCS.ResponseSigs())[i]
      progress6$inc(1/length(rownames(LINCS.ResponseSigs())), detail = paste("Calculating Correlations Against ", cmpdToOverlay))
      cmpd <- subset(a, a$compound == cmpdToOverlay)
      cmpd$Genes <- NULL
      cmpdGenes <- colnames(cmpd)
      cmpd_overlap <- colnames(total.transpose)[which(colnames(total.transpose) %in% cmpdGenes)]
      total.transpose.cmpd <- total.transpose[,cmpd_overlap] # FileA
      cmpd_ordered <- as.numeric(as.vector(t(cmpd)))
      names(cmpd_ordered) <- colnames(cmpd)
      cmpd_ordered2 <- cmpd_ordered[cmpd_overlap] # Character list in brackets orders to fit that character list...
      # head(colnames(total.transpose.cmpd) == names(cmpd_ordered2))
      SC <- cor(cmpd_ordered2, t(total.transpose.cmpd), method = "spearman")
      Final_Matrix <- rbind(Final_Matrix, SC)
      rownames(Final_Matrix) <- rownames(LINCS.ResponseSigs())[1:counter]
      counter <- counter + 1
    }
    Final_Matrix <- na.omit(Final_Matrix)
    RDS_Final_CorrMat(Final_Matrix)
  })

  ## Custom TCS Scoring
  ##############################################################################
  customTCScorrMat <- reactiveVal()
  observeEvent(eventExpr = input$CalculateCustomTCSConnectivities, { # Here, we should add dropdown menu to select data slot to pull scale.data from...
    #### tranpose scale.data
    RDSseurat <- rdsSeurat()
    if (unlist(SeuratObject::Version(RDSseurat))[1] < 5) {
      if (dim(RDSseurat@assays$RNA@scale.data)[1] == 0){
        print("scale.data slot is empty - scaling data")
        RDSseurat <- ScaleData(RDSseurat, do.center = T)
      }
    } else {
      print("V5 Seurat Object Detected: server.R #853")
      if (is.null(isOrthogonAL()) == TRUE){ # is this backwards? ############################################################################################# <
        print("Is null = F, non-orthogonal")
        # if (dim(RDSseurat@assays$RNA@layers$scale.data)[1] == 0){
        #   print("scale.data slot is empty - scaling data")
        #   RDSseurat <- ScaleData(RDSseurat, do.center = T)
        # }
      } else {
        print("Is.null = T, OrthogonalObject")
        DefaultAssay(RDSseurat) <- c("RNA_ortho")
        # if (dim(RDSseurat@assays$RNA_ortho@layers$scale.data)[1] == 0){
        #   print("scale.data slot is empty - scaling data")
        #   RDSseurat <- ScaleData(RDSseurat, do.center = T)
        # }
      }
      # DefaultAssay(RDSseurat) <- "RNA_ortho"
      # if (dim(RDSseurat@assays$RNA_ortho@layers$scale.data)[1] == 0){
      #   print("scale.data slot is empty - scaling data")
      #   RDSseurat <- ScaleData(RDSseurat, do.center = T)
      # }
    }

    a <- LINCS.ResponseSigs()
    a$compound <- rownames(a)
    progress5 <- shiny::Progress$new()
    progress5$set(message = "Transposing scale.data slot")

    if (is.null(isOrthogonAL()) == TRUE){ # is this backwards? ############################################################################################# <
      print("Is null = F, non-orthogonal")
    } else {
      print("Is.null = T, OrthogonalObject")
      DefaultAssay(RDSseurat) <- c("RNA_ortho")
    }
    ##########################################################################################################################################
    if (unlist(SeuratObject::Version(RDSseurat))[1] < 5){
      print("line 909: detected as less than v5")
      total.transpose <- t(RDSseurat@assays$RNA@scale.data) # which assay to integrate with?
    } else {
      print("We're here 927")
      if (is.null(isOrthogonAL()) == TRUE){ # is this backwards?
        print("Is.null = T, non-OrthogonalObject")
        total.transpose <- t(GetAssayData(RDSseurat, assay = "RNA", layer = "scale.data"))
      } else {
        print("Is null = F, orthogonal")
        total.transpose <- t(GetAssayData(RDSseurat, assay = "RNA_ortho", layer = "scale.data"))
      }
    }

    progress5$close()


    ####

    custTCS <- customTCScsv()
    print(colnames(custTCS))
    custTCS <- data.frame(gene = custTCS$gene,
                          log2FoldChange = custTCS$log2FoldChange,
                          padj = custTCS$padj)
    print(head(custTCS))
    # scoring approach fork
    if (input$TCS_scoring_method == "Spearman"){
      print("Spearman startin'!")
      Final_Matrix <- data.frame()
      cmpd_genes <- custTCS$gene
      print(head(colnames(total.transpose)))
      print(head(cmpd_genes))
      print(intersect(cmpd_genes,colnames(total.transpose)))
      cmpd_overlap <- colnames(total.transpose)[which(colnames(total.transpose) %in% cmpd_genes)]
      print(length(cmpd_overlap))
      total.transpose.cmpd <- total.transpose[,cmpd_overlap]
      cmpd_ordered <- as.numeric(as.vector(custTCS$log2FoldChange))
      names(cmpd_ordered) <- cmpd_genes
      cmpd_ordered2 <- cmpd_ordered[cmpd_overlap] # Character list in brackets orders to fit that character list...
      # head(colnames(total.transpose.cmpd) == names(cmpd_ordered2))
      SC <- cor(cmpd_ordered2, t(total.transpose.cmpd), method = "spearman")
      Final_Matrix <- rbind(Final_Matrix, SC)
      # rownames(Final_Matrix) <- input$Custom_TCS_nameInput
      # counter <- counter + 1
    }
    # print(head(Final_Matrix))
    rownames(Final_Matrix) <- c(input$Custom_TCS_nameInput)
    customTCScorrMat(Final_Matrix)
  })

  # need to insert conditional statement testing null state of RDS_Final_CorrMat
  output$CustomTCScorrMatrixCalculated <- reactive({
    return(!is.null(customTCScorrMat()))
  })
  outputOptions(output, 'CustomTCScorrMatrixCalculated', suspendWhenHidden = FALSE)

  output$downloadCustomTCScellConnectivities <- downloadHandler(
    filename = function() {
      paste(input$Custom_TCS_nameInput, "_connectivities_",Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame(t(customTCScorrMat())), file, row.names = T)
    }
  )


  seurat_custom_corradded <- reactive({ ############## testing...
    req(input$CalculateCustomTCSConnectivities)
    obj <- rdsSeurat()
    Final_Matrix <- as.data.frame(t(customTCScorrMat()))
    colnames(Final_Matrix) <- c(input$Custom_TCS_nameInput)
    print("Printing custom results colnames...")
    print(colnames(Final_Matrix))
    obj <- AddMetaData(obj, metadata = Final_Matrix)
    obj
  })

  ##############################################################################

  output$finalMatrixText <- renderText({
    Final_Matrix <- as.data.frame(t(RDS_Final_CorrMat()))
    print(head(colnames(Final_Matrix)))
    # print("testing final matrix rendertext")
  })

  # need to insert conditional statement testing null state of RDS_Final_CorrMat
  output$corrMatrixCalculated <- reactive({
    return(!is.null(RDS_Final_CorrMat()))
  })
  outputOptions(output, 'corrMatrixCalculated', suspendWhenHidden = FALSE)

  # show correlation matrix calc button
  output$uploadCorrelationMatrix <- reactive({
    if (input$uploadCorrelationMatrix == FALSE)
      return(TRUE)
  })
  outputOptions(output, 'uploadCorrelationMatrix', suspendWhenHidden = FALSE)


  # connec this to the download button
  # How do we fix this so that the download is always .csv and not lag making .html...?
  output$RDScorrMatDownload <- downloadHandler(
    filename = function() {
      paste("RDS_upload_L1000_consensus_corrmat.csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame(t(RDS_Final_CorrMat())), file, row.names = T) #customTCScorrMat
    }
  )

  ##############################################################################
  #
  ##############################################################################
  # # Logic to handle upload of previously calculated correlation matrix here...


  # add metadata to seurat object...
  # seurat_corradded <- reactive({
  #     req(input$CalculateRDS_L1000_Spearman_Mat)
  #     obj <- rdsSeurat()
  #     Final_Matrix <- as.data.frame(t(RDS_Final_CorrMat()))
  #     obj <- AddMetaData(obj, metadata = Final_Matrix)
  #     obj
  # })

  seurat_corradded <- reactive({
    if (input$uploadCorrelationMatrix == TRUE){
      req(input$uploadCorrelationMatrix == TRUE)
      obj <- rdsSeurat()
      Final_Matrix <- as.data.frame(t(corrMatUpload()))
      obj <- AddMetaData(obj, metadata = Final_Matrix)
      obj
    } else {
      req(input$CalculateRDS_L1000_Spearman_Mat)
      obj <- rdsSeurat()
      Final_Matrix <- as.data.frame(t(RDS_Final_CorrMat()))
      obj <- AddMetaData(obj, metadata = Final_Matrix)
      obj
    }
  })

  # seurat_corradded <- eventReactive(input$CalculateRDS_L1000_Spearman_Mat, {
  #   req(input$CalculateRDS_L1000_Spearman_Mat)
  #   obj <- rdsSeurat()
  #   Final_Matrix <- as.data.frame(t(RDS_Final_CorrMat()))
  #   obj <- AddMetaData(obj, metadata = Final_Matrix)
  #   obj
  # })
  # #
  # seurat_corradded <- eventReactive(input$uploadCorrelationMatrix, {
  #   req(input$uploadCorrelationMatrix == TRUE)
  #   obj <- rdsSeurat()
  #   Final_Matrix <- as.data.frame(t(corrMatUpload()))
  #   obj <- AddMetaData(obj, metadata = Final_Matrix)
  #   obj
  # })

  ##############################################################################
  #
  ##############################################################################

  # testing
  output$finalMatrixAddedText <- renderText({
    obj <- seurat_corradded()
    Final_Matrix <- as.data.frame(t(RDS_Final_CorrMat()))
    print(head(rownames(Final_Matrix)))
    print(class(obj))
    print(head(obj$alisertib))
    # print("testing final matrix rendertext")
  })

  # Visualize compound consensus signature
  output$sigBar <- renderPlot({ # Here, need to filter for non-zero genes!
    compoundSignature <- LINCS.ResponseSigs()[input$referenceCompound,] # subset single-compound data
    compoundSignature$Genes <- NULL # set-up
    compoundSignature <- as.data.frame(t(compoundSignature)) # set-up
    compoundSignature$Genes <- rownames(compoundSignature) # set-up
    compoundSignature <- subset(compoundSignature, compoundSignature[input$referenceCompound] != 0) # remove zero values
    for (i in 1:length(rownames(compoundSignature))){ # Assign positive or negative values so we can color by that.
      if (compoundSignature[i, input$referenceCompound] > 0){
        compoundSignature$PosNeg[i] <- "Pos"
      }
      if (compoundSignature[i, input$referenceCompound] < 0){
        compoundSignature$PosNeg[i] <- "Neg"
      }
    }

    compoundSignature <- compoundSignature %>% arrange(desc(noquote(input$referenceCompound)))
    compoundSignature$Genes <- factor(compoundSignature$Genes, levels = compoundSignature$Genes[order(compoundSignature[input$referenceCompound])])


    p <- ggplot(
      compoundSignature,
      aes_string(x = "Genes", y = input$referenceCompound, fill = "PosNeg")) +
      geom_col(position = "identity", colour = "black", size = 0.25) +
      scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
      ggtitle(paste("L1000 Derived Consensus Signature",
                    "_",
                    input$referenceCompound)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    p
  })

  ###################################################################################################
  #            Visualize compound consensus signature overlaid on single-cell Data (DimPlot)        #
  ###################################################################################################

  # output$refSignatureUMAP <- renderPlot({ # FeaturePlot of selected compound spearman correlations
  #   if (input$dataChoice == "10X_GBM") {
  #     FeaturePlot(testdat,
  #                 features = c(input$referenceCompound),
  #                 min.cutoff = input$sigMinCutoff,
  #                 max.cutoff = input$sigMaxCutoff,
  #                 reduction = input$reductionUse,
  #                 cols = c("red", "black"),
  #                 pt.size = 1, order = F
  #     )
  #   }
  # })

  output$refSignatureUMAPRDS <- renderPlot({
    FeaturePlot(seurat_corradded(),
                features = c(input$referenceCompound),
                min.cutoff = input$sigMinCutoffRDS,
                max.cutoff = input$sigMaxCutoffRDS,
                reduction = input$reductionUseRDS,
                cols = c("red", "black"),
                pt.size = 1, order = F) +
      theme_void() +
      theme(plot.title = element_text(size = 18, face = "bold"))
  })

  output$CustomRefSignatureUMAPRDS <- renderPlot({
    FeaturePlot(seurat_custom_corradded(),
                features = c(input$Custom_TCS_nameInput),
                min.cutoff = input$sigMinCutoff_customTCS,
                max.cutoff = input$sigMaxCutoff_customTCS,
                reduction = input$reductionUseRDS,
                cols = c("red", "black"),
                pt.size = 1, order = F) +
      theme_void() +
      theme(plot.title = element_text(size = 18, face = "bold"))
  })

  # Visualization of compound spearmans across clusters

  output$refSignatureVln <- renderPlot({
    if (input$dataChoice == "10X_GBM") {
      VlnPlot(testdat,
              features = gsub("-", ".", c(input$referenceCompound)),group.by = input$groupBy, fill.by = "ident")
    }
  })

  output$refSignatureVlnRDS <- renderPlot({
    if (input$dataChoice == "RDS_Upload") {
      VlnPlot(rdsSeurat(), features = c(input$referenceCompound), group.by = input$groupByRDS, fill.by = "ident")
    }
  })

  # Visualization of cell-by-cell expression correlations for paired compounds
  # need to set up a new selection button for 2nd compound...

  output$combinationCorrelation <- renderPlot({
    if (input$dataChoice == "10X_GBM") {
      FeatureScatter(testdat, feature1 = input$referenceCompound, feature2 = input$correlationCompound, shuffle = T, group.by = input$groupBy, pt.size = 1, jitter = T)
    }
    if (input$dataChoice == "RDS_Upload") {
      FeatureScatter(rdsSeurat(), feature1 = input$referenceCompound, feature2 = input$correlationCompoundRDS, shuffle = T, group.by = input$groupByRDS, pt.size = 1, jitter = T)
    }
  })

  ####################################################################################################
  # Compound-Cell Interaction Correlation Matrix

  corrMat <- reactiveVal()
  corrMatMOA <- reactiveVal()
  # observeEvent(eventExpr = input$calculateCorrMatrix, {
  #   print("is true!")
  #   # print("spearman Plotting...")
  #   comp_cell_corres <- cor(testSpearman, method = "spearman")
  #   print("Corr mat is done...")
  #   ### Add in annotations...
  #
  #   MOA <- read.csv(file = "~/Documents/Projects/scGBM/Patient/Patient_Processing_Pipeline/20191115_relationships_parsed.csv")
  #   MOA <- MOA[,c("cmap_name", "target")]
  #   spearmanMOA <- merge(comp_cell_corres, MOA, by.x = "row.names", by.y = "cmap_name")# Add in compound targets to spearman correlation dataframe
  #
  #   anno <- data.frame(compounds = spearmanMOA$Row.names, target = spearmanMOA$target)
  #
  #   anno2 <- subset(anno, duplicated(anno$compounds) == F)
  #   duplicated(anno2)
  #   duplicated(anno2$compounds)
  #   rownames(anno2) <- anno2$compounds
  #   anno2$compounds <- NULL
  #
  #   corrMat(comp_cell_corres)
  #   corrMatMOA(anno2)
  #
  #   # targets <- c("Bromodomain")
  #   # keep <- subset(anno2, anno2$target %in% targets | rownames(anno2) %in% c("alisertib", "temozolomide"))
  #   #
  #   # anno3 <- subset(anno2, rownames(anno2) %in% rownames(keep))
  #   #
  #   # comp_cell_corrMat <- pheatmap(comp_cell_corres,
  #   #                               color = viridis(256, option = "C"), show_rownames = T, show_colnames = F,
  #   #                               height = 4,
  #   #                               width = 4,
  #   #                               border_color = NA,
  #   #                               # cellwidth = 12,
  #   #                               # cellheight = 12,
  #   #                               silent = T)
  #
  # })

  # output$combinationCorrelation <- renderPlot({
  #   if (input$calculateCorrMatrix >= 1){
  #     # # print("spearman Plotting...")
  #     # comp_cell_corres <- cor(testSpearman, method = "spearman")
  #     #
  #     # ### Add in annotations...
  #     #
  #     # MOA <- read.csv(file = "~/Documents/Projects/scGBM/Patient/Patient_Processing_Pipeline/20191115_relationships_parsed.csv")
  #     # MOA <- MOA[,c("cmap_name", "target")]
  #     # spearmanMOA <- merge(comp_cell_corres, MOA, by.x = "row.names", by.y = "cmap_name")# Add in compound targets to spearman correlation dataframe
  #     #
  #     # anno <- data.frame(compounds = spearmanMOA$Row.names, target = spearmanMOA$target)
  #     #
  #     # anno2 <- subset(anno, duplicated(anno$compounds) == F)
  #     # duplicated(anno2)
  #     # duplicated(anno2$compounds)
  #     # rownames(anno2) <- anno2$compounds
  #     # anno2$compounds <- NULL
  #
  #     comp_cell_corres <- corrMat()
  #     anno2 <- corrMatMOA()
  #
  #     targets <- c("Bromodomain")
  #     keep <- subset(anno2, anno2$target %in% targets | rownames(anno2) %in% c("alisertib", "temozolomide"))
  #
  #     anno3 <- subset(anno2, rownames(anno2) %in% rownames(keep))
  #
  #     comp_cell_corrMat <- pheatmap(comp_cell_corres,
  #                                   color = viridis(256, option = "C"), show_rownames = T, show_colnames = F,
  #                                   height = 4,
  #                                   width = 4,
  #                                   border_color = NA,
  #                                   # cellwidth = 12,
  #                                   # cellheight = 12,
  #                                   silent = T)#,
  #     # annotation_col = anno3,
  #     # annotation_row = anno3)
  #     # comp_cell_corrMat
  #
  #     # Can we add select labels?
  #
  #     # add.flag <- function(pheatmap,
  #     #                      kept.labels,
  #     #                      repel.degree) {
  #     #
  #     #   # repel.degree = number within [0, 1], which controls how much
  #     #   #                space to allocate for repelling labels.
  #     #   ## repel.degree = 0: spread out labels over existing range of kept labels
  #     #   ## repel.degree = 1: spread out labels over the full y-axis
  #     #
  #     #   heatmap <- pheatmap$gtable
  #     #   new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]
  #     #   # keep only labels in kept.labels, replace the rest with ""
  #     #   new.label$label <- ifelse(new.label$label %in% kept.labels,
  #     #                             new.label$label, "")
  #     #   # calculate evenly spaced out y-axis positions
  #     #   repelled.y <- function(d, d.select, k = repel.degree){
  #     #     # d = vector of distances for labels
  #     #     # d.select = vector of T/F for which labels are significant
  #     #
  #     #     # recursive function to get current label positions
  #     #     # (note the unit is "npc" for all components of each distance)
  #     #     strip.npc <- function(dd){
  #     #       if(!"unit.arithmetic" %in% class(dd)) {
  #     #         return(as.numeric(dd))
  #     #       }
  #     #
  #     #       d1 <- strip.npc(dd$arg1)
  #     #       d2 <- strip.npc(dd$arg2)
  #     #       fn <- dd$fname
  #     #       return(lazyeval::lazy_eval(paste(d1, fn, d2)))
  #     #     }
  #     #
  #     #     full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
  #     #     selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
  #     #
  #     #     return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
  #     #                     to = min(selected.range) - k*(min(selected.range) - min(full.range)),
  #     #                     length.out = sum(d.select)),
  #     #                 "npc"))
  #     #   }
  #     #   new.y.positions <- repelled.y(new.label$y,
  #     #                                 d.select = new.label$label != "")
  #     #   new.flag <- segmentsGrob(x0 = new.label$x,
  #     #                            x1 = new.label$x + unit(0.15, "npc"),
  #     #                            y0 = new.label$y[new.label$label != ""],
  #     #                            y1 = new.y.positions)
  #     #
  #     #   # shift position for selected labels
  #     #   new.label$x <- new.label$x + unit(0.2, "npc")
  #     #   new.label$y[new.label$label != ""] <- new.y.positions
  #     #
  #     #   # add flag to heatmap
  #     #   heatmap <- gtable::gtable_add_grob(x = heatmap,
  #     #                                      grobs = new.flag,
  #     #                                      t = 4,
  #     #                                      l = 4
  #     #   )
  #     #
  #     #   # replace label positions in heatmap
  #     #   heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  #     #
  #     #   # plot result
  #     #   grid.newpage()
  #     #   grid.draw(heatmap)
  #     #
  #     #   # return a copy of the heatmap invisibly
  #     #   invisible(heatmap)
  #     # }
  #     #
  #     # comp_cell_corrMat2 <- add.flag(pheatmap = comp_cell_corrMat,
  #     #                                kept.labels = rownames(keep),
  #     #                                repel.degree = 0.2)
  #
  #     comp_cell_corrMat2
  #   }
  # })


  ####################################################################################################
  #                                         scSynergySeq Tab                                         #
  ####################################################################################################

  # Based on single-cell correlations - set cut-off to remove cells from dataset,
  # and calculate potential compounds to pair with reference compound. (script exists).

  # Need to subset seurat object here, then calculate and display top compounds
  # targeting the remainder of cells. (triangle plot?)

  # use cowplot to knit multiple plots together... one = display cutoff cells.

  # output$cutOffRange <- renderText({
  #   c(paste("Sensitivity Threshold =", input$correlationCutoff[1], ";"),
  #     paste("Resistance Threshold =", input$correlationCutoff[2]))
  # })
  #
  # output$scSynergySeq1 <- renderPlot({
  #   tumorCells <- WhichCells(testdat, idents = input$cancerCellIdents)
  #
  #   compoundSpearmans <- testdat@meta.data[input$referenceCompound]
  #   living <- subset(compoundSpearmans,
  #                    subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff[2])
  #   dead <- subset(compoundSpearmans,
  #                  subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff[1])
  #   livingCells <- rownames(living)
  #   deadCells <- rownames(dead)
  #   p1 <- DimPlot(testdat, reduction = input$reductionUse,
  #                 cells.highlight = intersect(deadCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("blue"))
  #   p1 <- p1 + ggtitle("Sensitive Cells") + NoLegend()
  #   p2 <- DimPlot(testdat, reduction = input$reductionUse,
  #                 cells.highlight = intersect(livingCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("red"))
  #   p2 <- p2 + ggtitle("Resistant Cells") + NoLegend()
  #
  #
  #   # Make a plot to show percentage of live vs dead cells...
  #   liveDead <- data.frame(row.names = NULL, resistant = length(livingCells), sensitive = length(deadCells))
  #   liveDead <- as.data.frame(t(liveDead))
  #   liveDead$liveDead <- rownames(liveDead)
  #   colnames(liveDead) <- c("Freq", "liveDead")
  #   liveDead
  #   bp <- ggplot(liveDead, aes(x="", y = Freq, fill = liveDead))+
  #     geom_bar(width = 1, stat = "identity")
  #   pie <- bp + coord_polar("y", start = 0)
  #   pie <- pie + theme_minimal() + scale_fill_nejm()
  #   p3 <- plot_grid(plotlist = list(p1, p2, pie), ncol = 3, labels = "AUTO")
  #   p3
  # })

  output$cutOffRange_RDS <- renderText({
    c(paste("Sensitivity Threshold =", input$correlationCutoff_RDS[1], ";"),
      paste("Resistance Threshold =", input$correlationCutoff_RDS[2]))
  })

  output$referenceCompound <- renderText({
    print(paste0("Treat with " , c(input$referenceCompound)))
  })


  sensitiveCells <- reactiveVal()
  resistantCells <- reactiveVal()

  output$scSynergySeq1RDS <- renderPlot({ # need to update RDS section with dual slider selection
    RDSseurat <- seurat_corradded()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
    tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)

    compoundSpearmans <- RDSseurat@meta.data[input$referenceCompound] # How to do this for the RDS upload... where is the spearman matrix and do we add it as metadata to the seurat object?

    living <- subset(compoundSpearmans, subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff_RDS)
    dead <- subset(compoundSpearmans, subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff_RDS)
    livingCells <- rownames(living)
    deadCells <- rownames(dead)

    resistantCells(livingCells)
    sensitiveCells(deadCells)

    p1 <- DimPlot(RDSseurat, reduction = input$reductionUseRDS,
                  cells.highlight = intersect(deadCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("blue")) +
      theme_void()
    p1 <- p1 + ggtitle("Sensitive Cells") + NoLegend()
    p2 <- DimPlot(RDSseurat, reduction = input$reductionUseRDS,
                  cells.highlight = intersect(livingCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("red")) +
      theme_void()
    p2 <- p2 + ggtitle("Resistant Cells") + NoLegend()


    # Make a plot to show percentage of live vs dead cells...
    liveDead <- data.frame(row.names = NULL, resistant = length(livingCells), sensitive = length(deadCells))
    liveDead <- as.data.frame(t(liveDead))
    liveDead$liveDead <- rownames(liveDead)
    colnames(liveDead) <- c("Freq", "liveDead")
    liveDead
    bp <- ggplot(liveDead, aes(x="", y = Freq, fill = liveDead))+
      geom_bar(width = 1, stat = "identity")
    pie <- bp + coord_polar("y", start = 0)
    pie <- pie + theme_void() + scale_fill_nejm()
    p3 <- plot_grid(plotlist = list(p1, p2, pie), ncol = 3)
    p3
  })


  sensitiveCells_customTCS <- reactiveVal()
  resistantCells_customTCS <- reactiveVal()
  output$scSynergySeq_custTCS <- renderPlot({ # need to update RDS section with dual slider selection
    RDSseurat <- seurat_custom_corradded()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
    tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)

    compoundSpearmans <- RDSseurat@meta.data[input$Custom_TCS_nameInput] # How to do this for the RDS upload... where is the spearman matrix and do we add it as metadata to the seurat object?

    living <- subset(compoundSpearmans, subset = compoundSpearmans[input$Custom_TCS_nameInput] > input$correlationCutoff_customTCS)
    dead <- subset(compoundSpearmans, subset = compoundSpearmans[input$Custom_TCS_nameInput] < input$correlationCutoff_customTCS)
    livingCells <- rownames(living)
    deadCells <- rownames(dead)
    resistantCells_customTCS(livingCells)
    sensitiveCells_customTCS(deadCells)

    p1 <- DimPlot(RDSseurat, reduction = input$reductionUseRDS,
                  cells.highlight = intersect(deadCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("blue")) +
      theme_void()
    p1 <- p1 + ggtitle("Sensitive Cells") + NoLegend()
    p2 <- DimPlot(RDSseurat, reduction = input$reductionUseRDS,
                  cells.highlight = intersect(livingCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("red")) +
      theme_void()
    p2 <- p2 + ggtitle("Resistant Cells") + NoLegend()


    # Make a plot to show percentage of live vs dead cells...
    liveDead <- data.frame(row.names = NULL, resistant = length(livingCells), sensitive = length(deadCells))
    liveDead <- as.data.frame(t(liveDead))
    liveDead$liveDead <- rownames(liveDead)
    colnames(liveDead) <- c("Freq", "liveDead")
    liveDead
    bp <- ggplot(liveDead, aes(x="", y = Freq, fill = liveDead))+
      geom_bar(width = 1, stat = "identity")
    pie <- bp + coord_polar("y", start = 0)
    pie <- pie + theme_void() + scale_fill_nejm()
    p3 <- plot_grid(plotlist = list(p1, p2, pie), ncol = 3)
    p3
  })

  ##############################################################################
  #
  ##############################################################################

  output$scSynergySeq2_RDS <- renderPlot({ # need to update with dual slider
    if (input$perturbationButton >= 1){
      RDSseurat <-  seurat_corradded()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      if (input$uploadCorrelationMatrix == TRUE){
        print(class(corrMatUpload()))
        testSpearman <- as.data.frame(t(corrMatUpload()))
      } else {
        print(class(RDS_Final_CorrMat()))
        testSpearman <- as.data.frame(t(RDS_Final_CorrMat()))
      }
      tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)

      compoundSpearmans <- RDSseurat@meta.data[input$referenceCompound] # maybe need referenceCompoundRDS
      compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)

      living <- subset(compoundSpearmans,
                       subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff_RDS)
      dead <- subset(compoundSpearmans,
                     subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff_RDS)
      livingCells <- rownames(living)
      deadCells <- rownames(dead)

      testSpearman <- subset(testSpearman, rownames(testSpearman) %in% tumorCells) # added this line to calculate synergy predictions using only tumor cells (data is scaled with non-tumor cells so still representative of disease signatures)

      deadData <- testSpearman[rownames(dead),]
      liveData <- testSpearman[rownames(living),]
      sensitiveMeans <- as.data.frame(colMeans(deadData))
      resistantMeans <- as.data.frame(colMeans(liveData))
      sensitiveMeans$compound <- rownames(sensitiveMeans)
      resistantMeans$compound <- rownames(resistantMeans)

      colnames(resistantMeans) <- c("mean", "compound")
      colnames(sensitiveMeans) <- c("mean", "compound")

      # merge live and dead so we can do calculations across columns.

      newMeans <- merge(x = resistantMeans, y = sensitiveMeans, by.x = "compound", by.y = "compound")
      colnames(newMeans) <- c("compound", "resistantMean", "sensitiveMean")
      newMeans$deltaMean <- newMeans$resistantMean - newMeans$sensitiveMean


      # Which compounds are most predicted to work on this population?
      # And which compounds means changed most from dead population?

      resistantMeansMin <- resistantMeans %>% slice_min(order_by = mean, n = input$nCompoundSlider_RDS / 2)
      resistantMeansMax <- resistantMeans %>% slice_max(order_by = mean, n = input$nCompoundSlider_RDS / 2)

      sensitiveMeansMax <- sensitiveMeans %>% slice_max(order_by = mean, n = input$nCompoundSlider_RDS / 2)
      sensitiveMeansMin <- sensitiveMeans %>% slice_min(order_by = mean, n = input$nCompoundSlider_RDS / 2)

      deltaMeansMax <- newMeans %>% slice_max(order_by = deltaMean, n = input$nCompoundSlider_RDS / 2)
      deltaMeansMin <- newMeans %>% slice_min(order_by = deltaMean, n = input$nCompoundSlider_RDS / 2)

      # resistantMeansMin <- resistantMeans %>% slice_min(order_by = mean, n = 15)
      # resistantMeansMax <- resistantMeans %>% slice_max(order_by = mean, n = 15)
      resistantMeans <- rbind(resistantMeansMin, resistantMeansMax)
      sensitiveMeans <- rbind(sensitiveMeansMin, sensitiveMeansMax)
      deltaMeans <- rbind(deltaMeansMin, deltaMeansMax)


      # # Should just do this now because will want to include this in deltaMean part.
      for (i in 1:length(rownames(resistantMeans))){ # Assign positive or negative values so we can color by that.
        if (resistantMeans$mean[i] > 0){
          resistantMeans$PosNeg[i] <- "Pos"
        }
        if (resistantMeans$mean[i] < 0){
          resistantMeans$PosNeg[i] <- "Neg"
        }
      }

      for (i in 1:length(rownames(sensitiveMeans))){
        if (sensitiveMeans$mean[i] > 0){
          sensitiveMeans$PosNeg[i] <- "Pos"
        }
        if (sensitiveMeans$mean[i] < 0){
          sensitiveMeans$PosNeg[i] <- "Neg"
        }
      }

      for (i in 1:length(rownames(deltaMeans))){
        if (deltaMeans$deltaMean[i] > 0){
          deltaMeans$PosNeg[i] <- "Pos"
        }
        if (deltaMeans$deltaMean[i] < 0){
          deltaMeans$PosNeg[i] <- "Neg"
        }
      }
      #
      # resistantMeans <- as_tibble(resistantMeans, rownames = NA) # change class to tibble
      resistantMeans <- resistantMeans %>% arrange(desc(mean))
      sensitiveMeans <- sensitiveMeans %>% arrange(desc(mean))
      deltaMeans <- deltaMeans %>% arrange(desc(deltaMean))
      # resistantMeans$compound <- factor(resistantMeans$compound, levels = resistantMeans$compound[order(resistantMeans[resistantMeans$mean])]) ##################### !!!!!!!!!! THIS IS THE PROBLEM LINE
      #
      #

      # newMeans$pc <- predict(prcomp(~resistantMean+deltaMean, newMeans))[,1]
      # newMeans$density <- fields::interp.surface(MASS::kde2d(newMeans$resistantMean, newMeans$deltaMean),
      #                                            newMeans[,c("resistantMean", "deltaMean")])
      #
      # newMeansLabelDF <- rbind(slice_min(newMeans, order_by = pc, n = 10),
      #                          slice_max(newMeans, order_by = pc, n = 10))
      #
      # What if we plot the resistant correlation vs the delta correlation, the most negative should pair best with alisertib.



      p4 <- ggplot(
        resistantMeans,
        aes_string(x = "compound", y = "mean", fill = "PosNeg")) +
        geom_col(position = "identity", colour = "black", size = 0.25) +
        scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
        ggtitle(paste("Mean correlations with 'remaining' resistant population")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

      p5 <- ggplot(
        sensitiveMeans,
        aes_string(x = "compound", y = "mean", fill = "PosNeg")) +
        geom_col(position = "identity", colour = "black", size = 0.25) +
        scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
        ggtitle(paste("Mean correlations with 'killed' cell population")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


      p6 <- ggplot(
        deltaMeans,
        aes_string(x = "compound", y = "deltaMean", fill = "PosNeg")) +
        geom_col(position = "identity", colour = "black", size = 0.25) +
        scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
        ggtitle(paste("Differences in mean correlations to compound signatures between 'live' and 'dead' cells:")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      # p5

      # p7 <- ggplot(newMeans, aes(resistantMean, deltaMean, color = pc, alpha = 1/density)) +
      #   geom_point(shape = 16, size = 2, show.legend = F) +
      #   # geom_text_repel(data = newMeansLabelDF,
      #   #                 aes(x=resistantMean, y = deltaMean, label = compound),
      #   #                 color = "black") +
      #   theme_minimal() +
      #   scale_color_gradient(low = "#0091ff", high = "#f0650e")


      p8 <- plot_grid(plotlist = list(p4, p5), labels = "AUTO", ncol = 2)
      p9 <- plot_grid(plotlist = list(p8, p6), labels = c("", "C"), ncol = 1)
      p9
      #
      # p6
    }

  })

  tumorCellsRDS <- reactiveVal()
  livingCellsRDS <- reactiveVal()
  deadCellsRDS <- reactiveVal()
  output$scSynergySeq3_RDS <- renderPlot({ # need to update with dual slider
    if (input$perturbationButton >= 1){
      RDSseurat <-  seurat_corradded()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      if (input$uploadCorrelationMatrix == TRUE){
        print(class(corrMatUpload()))
        testSpearman <- as.data.frame(t(corrMatUpload()))
      } else {
        print(class(RDS_Final_CorrMat()))
        testSpearman <- as.data.frame(t(RDS_Final_CorrMat()))
      }
      tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
      compoundSpearmans <- RDSseurat@meta.data[input$referenceCompound]
      compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)
      living <- subset(compoundSpearmans,
                       subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff_RDS)
      dead <- subset(compoundSpearmans,
                     subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff_RDS)
      livingCells <- rownames(living)
      deadCells <- rownames(dead)
      tumorCellsRDS(tumorCells)
      livingCellsRDS(livingCells)
      deadCellsRDS(deadCells)
      testSpearman <- subset(testSpearman, rownames(testSpearman) %in% tumorCells) # added this line to calculate synergy predictions using only tumor cells (data is scaled with non-tumor cells so still representative of disease signatures)
      deadData <- testSpearman[rownames(dead),]
      liveData <- testSpearman[rownames(living),]
      sensitiveMeans <- as.data.frame(colMeans(deadData))
      resistantMeans <- as.data.frame(colMeans(liveData))
      sensitiveMeans$compound <- rownames(sensitiveMeans)
      resistantMeans$compound <- rownames(resistantMeans)
      colnames(resistantMeans) <- c("mean", "compound")
      colnames(sensitiveMeans) <- c("mean", "compound")

      # merge live and dead so we can do calculations across columns.
      newMeans <- merge(x = resistantMeans, y = sensitiveMeans, by.x = "compound", by.y = "compound")
      colnames(newMeans) <- c("compound", "resistantMean", "sensitiveMean")
      newMeans$deltaMean <- newMeans$resistantMean - newMeans$sensitiveMean
      # newMeans$pc <- predict(prcomp(~resistantMean+deltaMean, newMeans))[,1]
      newMeans$density <- fields::interp.surface(MASS::kde2d(newMeans$resistantMean, newMeans$deltaMean),
                                                 newMeans[,c("resistantMean", "deltaMean")])

      newMeansLabelDF <- rbind(slice_min(newMeans, order_by = pc, n = 100),
                               slice_max(newMeans, order_by = pc, n = 100))

      # create new to label df based on values below 0, select top n from slider? How to calculate most double negative? (size of triangle?)

      # 1) subset double negative compounds
      # 2) calculate absolute values of each axis
      # 3) Calculate hypotenuse
      # 4) Calculate area
      # 5) What's better way to visualize?


      # 1)
      newMeansLabelDF2 <- subset(newMeans, newMeans$resistantMean < 0 & newMeans$deltaMean < 0)
      print(dim(newMeansLabelDF2)) # 514 x 6

      # 2)
      newMeansLabelDF2$absResistantMean <- abs(newMeansLabelDF2$resistantMean)
      newMeansLabelDF2$absDeltaMean <- abs(newMeansLabelDF2$deltaMean)

      # 3) calculate hypotenuse
      newMeansLabelDF2$absHyp <- sqrt(newMeansLabelDF2$absResistantMean^2 + newMeansLabelDF2$absDeltaMean^2)

      newMeansLabelDF3 <- slice_max(newMeansLabelDF2, order_by = absHyp, n = input$synergyPredictorSlider)
      # newMeansLabelDF2 <- rbind(slice_min)

      p7 <- ggplot(newMeans, aes(resistantMean, deltaMean, color = deltaMean)) + # alpha = 1/density
        geom_point(shape = 16, size = 2, show.legend = F) +
        geom_text_repel(data = newMeansLabelDF3, max.overlaps = 15,
                        aes(x=resistantMean, y = deltaMean, label = compound),
                        color = "black") +
        # theme_minimal() +
        scale_color_gradient(low = "#0091ff", high = "#f0650e") + facet_zoom(xlim = c(min(newMeans$resistantMean) - 0.01, 0), ylim = c(min(newMeans$deltaMean) - 0.01, 0))

      p7
    }
  })

  ##############################################################################
  # scFOCAL Combination Indexing
  ##############################################################################

  ## Custom Input Combination Indexing
  tumorCells_customTCS <- reactiveVal()
  # resistantCells_customTCS <- reactiveVal()
  # sensitiveCells_customTCS <- reactiveVal()
  limmaRes_custTCS <- reactiveVal()
  zmat_custTCS <- reactiveVal()
  seurRDS_customTCSsensitivity <- reactiveVal()
  seurRDS_customTCS_sensitivityAssigner <- reactiveVal()
  observeEvent(eventExpr = input$perturbationButton_customTCS, {
    if (input$perturbationButton_customTCS >= 1){
      RDSseurat <- seurat_custom_corradded()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      if (input$uploadCorrelationMatrix == TRUE){
        testSpearman <- as.data.frame(t(corrMatUpload()))
      } else {
        testSpearman <- as.data.frame(t(RDS_Final_CorrMat()))
      }
      tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
      compoundSpearmans <- RDSseurat@meta.data[input$Custom_TCS_nameInput]
      compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)
      testSpearman <- testSpearman[which(rownames(testSpearman) %in% tumorCells),]
      z_matrix <- 0.5 * log((1 + testSpearman) / (1 - testSpearman)) # apply Fisher z-transformation
      resistantCells <- resistantCells_customTCS()
      sensitiveCells <- sensitiveCells_customTCS()
      #####
      sensitivityAssigner <- data.frame(cells = c(resistantCells, sensitiveCells),
                                        sensitivity = c(rep("resistant", length(resistantCells)),
                                                        rep("sensitive", length(sensitiveCells))))
      rownames(sensitivityAssigner) <- sensitivityAssigner$cells
      sensitivityAssigner$cells <- NULL

      RDSseurat <- AddMetaData(RDSseurat, metadata = sensitivityAssigner)
      seurRDS_customTCSsensitivity(RDSseurat)
      sensitivityAssigner <- subset(sensitivityAssigner, rownames(sensitivityAssigner) %in% tumorCells)
      seurRDS_customTCS_sensitivityAssigner(sensitivityAssigner)
      #####
      # set up groupings and construct design matrix...
      snames <- colnames(testSpearman)
      subjectID <- data.frame(subjectID = RDSseurat@meta.data[,input$subject_replicate_ident], row.names = rownames(RDSseurat@meta.data))
      subjectID <- subjectID[which(rownames(subjectID) %in% rownames(testSpearman)),]

      sensitivity <- data.frame(sensitivity = RDSseurat$sensitivity, row.names = colnames(RDSseurat))
      # sensitivity <- subset(sensitivity, rownames(sensitivity) %in% tumorCells)
      sensitivity <- sensitivity[which(rownames(sensitivity) %in% rownames(testSpearman)),]
      # subject_group <- interaction(sensitivity$sensitivity, subjectID$subjectID)
      sensitivity <- factor(sensitivity)
      subjectID <- factor(subjectID)
      design <- model.matrix(~ subjectID + sensitivity)
      z_matrix <- t(z_matrix)
      zmat_custTCS(z_matrix)
      fit <- lmFit(z_matrix, design)
      fit <- eBayes(fit)
      results <- topTable(fit, coef = 2, adjust.method = "fdr", number = Inf)
      limmaRes_custTCS(results)
      print(head(results))
      print(head(colnames(z_matrix)))
    }
  })

  ## L1000 Input Combination Indexing
  tumorCells <- reactiveVal()
  limmaRes <- reactiveVal()
  zmat <- reactiveVal()
  seurRDSsensitivity <- reactiveVal()
  seurRDS_sensitivityAssigner <- reactiveVal()
  observeEvent(eventExpr = input$perturbationButton, {
    if (input$perturbationButton >= 1){
      RDSseurat <- seurat_corradded()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      if (input$uploadCorrelationMatrix == TRUE){
        testSpearman <- as.data.frame(t(corrMatUpload()))
      } else {
        testSpearman <- as.data.frame(t(RDS_Final_CorrMat()))
      }
      tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
      compoundSpearmans <- RDSseurat@meta.data[input$referenceCompound]
      compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)
      testSpearman <- testSpearman[which(rownames(testSpearman) %in% tumorCells),]
      z_matrix <- 0.5 * log((1 + testSpearman) / (1 - testSpearman)) # apply Fisher z-transformation
      resistantCells <- resistantCells()
      sensitiveCells <- sensitiveCells()
      #####
      sensitivityAssigner <- data.frame(cells = c(resistantCells, sensitiveCells),
                                        sensitivity = c(rep("resistant", length(resistantCells)),
                                                        rep("sensitive", length(sensitiveCells))))
      rownames(sensitivityAssigner) <- sensitivityAssigner$cells
      sensitivityAssigner$cells <- NULL

      RDSseurat <- AddMetaData(RDSseurat, metadata = sensitivityAssigner)
      seurRDSsensitivity(RDSseurat)
      sensitivityAssigner <- subset(sensitivityAssigner, rownames(sensitivityAssigner) %in% tumorCells)
      seurRDS_sensitivityAssigner(sensitivityAssigner)
      #####
      # set up groupings and construct design matrix...
      snames <- colnames(testSpearman)
      subjectID <- data.frame(subjectID = RDSseurat@meta.data[,input$subject_replicate_ident], row.names = rownames(RDSseurat@meta.data))
      subjectID <- subjectID[which(rownames(subjectID) %in% rownames(testSpearman)),]

      sensitivity <- data.frame(sensitivity = RDSseurat$sensitivity, row.names = colnames(RDSseurat))
      # sensitivity <- subset(sensitivity, rownames(sensitivity) %in% tumorCells)
      sensitivity <- sensitivity[which(rownames(sensitivity) %in% rownames(testSpearman)),]
      # subject_group <- interaction(sensitivity$sensitivity, subjectID$subjectID)
      sensitivity <- factor(sensitivity)
      subjectID <- factor(subjectID)
      design <- model.matrix(~ subjectID + sensitivity)
      z_matrix <- t(z_matrix)
      zmat(z_matrix)
      fit <- lmFit(z_matrix, design)
      fit <- eBayes(fit)
      results <- topTable(fit, coef = 2, adjust.method = "fdr", number = Inf)
      limmaRes(results)
      print(head(results))
      print(head(colnames(z_matrix)))
    }
  })

  ##############################################################################

  output$customTCS_limmaVolc <- renderPlot({
    if (input$perturbationButton_customTCS >= 1){
      res <- limmaRes_custTCS()
      volc <- EnhancedVolcano(toptable = res, # title = paste0(input$Custom_TCS_nameInput, " induced alterations"),
                              # subtitle = "fisher's Z transformation of raw connectivity values + limma",
                              lab = rownames(res),
                              x = "logFC",   xlab = bquote("Connectivity" ~Log[2] ~ "fold change"),
                              y = "adj.P.Val", FCcutoff = 0.01, title = NULL, subtitle = NULL,
                              xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.1, max(res[["logFC"]], na.rm = TRUE) + 0.1)
                              # ylim = c(0, max(-log10(res[["adj.P.Val"]]), na.rm = TRUE) + 5)
      )
      volc_facet <- volc + ggforce::facet_zoom(xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.05,0))
      volc_facet <- volc_facet + theme(panel.spacing.y = unit(0.1, "cm"))
      volc_facet
    }
  })

  output$l1000in_limmaVolc <- renderPlot({
    if (input$perturbationButton >= 1){
      res <- limmaRes()
      volc <- EnhancedVolcano(toptable = res, # title = paste0(input$Custom_TCS_nameInput, " induced alterations"),
                              # subtitle = "fisher's Z transformation of raw connectivity values + limma",
                              lab = rownames(res),
                              x = "logFC",   xlab = bquote("Connectivity" ~Log[2] ~ "fold change"),
                              y = "adj.P.Val", FCcutoff = 0.01, title = NULL, subtitle = NULL,
                              xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.1, max(res[["logFC"]], na.rm = TRUE) + 0.1)
                              # ylim = c(0, max(-log10(res[["adj.P.Val"]]), na.rm = TRUE) + 5)
                              )
      volc_facet <- volc + ggforce::facet_zoom(xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.05,0))
      volc_facet <- volc_facet + theme(panel.spacing.y = unit(0.1, "cm"))
      volc_facet
    }
  })

  output$download_customTCS_sensitivityAssigner <- downloadHandler(
    filename = function() {
      paste(input$Custom_TCS_nameInput, "_sensitivityByTCS_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(seurRDS_customTCS_sensitivityAssigner(), file)

    }
  )

  output$download_sensitivityAssigner <- downloadHandler(
    filename = function() {
      paste(input$referenceCompound, "_sensitivityByTCS_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(seurRDS_sensitivityAssigner(), file)

    }
  )

  output$download_customTCS_limmaVolc <- downloadHandler(
    filename = function() {
      paste(input$Custom_TCS_nameInput, "diffConnectivityVolcano_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      res <- limmaRes_custTCS()
      pdf(file, height = 14)
      volc <- EnhancedVolcano(toptable = res, # title = paste0(input$Custom_TCS_nameInput, " induced alterations"),
                              # subtitle = "fisher's Z transformation of raw connectivity values + limma",
                              lab = rownames(res), title = NULL, subtitle = NULL,
                              x = "logFC",   xlab = bquote("Connectivity" ~Log[2] ~ "fold change"),
                              y = "adj.P.Val", FCcutoff = 0.01,
                              xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.1, max(res[["logFC"]], na.rm = TRUE) + 0.1),
                              ylim = c(0, max(-log10(res[["adj.P.Val"]]), na.rm = TRUE) + 5))
      volc_facet <- volc + ggforce::facet_zoom(xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.05,0))
      volc_facet <- volc_facet + theme(panel.spacing.y = unit(0.1, "cm"))
      print(volc_facet)
      dev.off()
    }
  )

  output$download_limmaVolc <- downloadHandler(
    filename = function() {
      paste(input$referenceCompound, "diffConnectivityVolcano_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      res <- limmaRes()
      pdf(file, height = 14)
      volc <- EnhancedVolcano(toptable = res, # title = paste0(input$Custom_TCS_nameInput, " induced alterations"),
                              # subtitle = "fisher's Z transformation of raw connectivity values + limma",
                              lab = rownames(res), title = NULL, subtitle = NULL,
                              x = "logFC",   xlab = bquote("Connectivity" ~Log[2] ~ "fold change"),
                              y = "adj.P.Val", FCcutoff = 0.01,
                              xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.1, max(res[["logFC"]], na.rm = TRUE) + 0.1),
                              ylim = c(0, max(-log10(res[["adj.P.Val"]]), na.rm = TRUE) + 5))
      volc_facet <- volc + ggforce::facet_zoom(xlim = c(min(res[["logFC"]], na.rm = TRUE) - 0.05,0))
      volc_facet <- volc_facet + theme(panel.spacing.y = unit(0.1, "cm"))
      print(volc_facet)
      dev.off()
    }
  )

  mrc_custTCS <- reactiveVal()
  resistMat_custTCS <- reactiveVal()
  output$customTCS_rcm_overall <- renderPlot({
    if (input$perturbationButton_customTCS >= 1){
      zmat <- zmat_custTCS()
      resistantCells <- resistantCells_customTCS()
      tumorCells <- colnames(zmat)
      resistMat <- zmat[,which(colnames(zmat) %in% resistantCells)]
      resistMat_custTCS(resistMat)
      mrc <- data.frame(mrc = rowMeans(resistMat), row.names = rownames(resistMat))
      mrc$Group <- rownames(mrc)
      mrc_custTCS(mrc)

      # Create a barplot mapping the mrc values to fill.
      ggplot(mrc, aes(x = reorder(Group, mrc), y = mrc, fill = mrc)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(option = "viridis", direction = -1) +
        labs(x = "Group", y = "mrc",
             title = paste0("Mean connectivities to ", input$Custom_TCS_nameInput, " resistant cells")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  mrc <- reactiveVal()
  resistMat <- reactiveVal()
  output$rcm_overall <- renderPlot({
    if (input$perturbationButton >= 1){
      zmat <- zmat()
      resistantCells <- resistantCells()
      tumorCells <- colnames(zmat)
      resistMat <- zmat[,which(colnames(zmat) %in% resistantCells)]
      resistMat(resistMat)
      mrc <- data.frame(mrc = rowMeans(resistMat), row.names = rownames(resistMat))
      mrc$Group <- rownames(mrc)
      mrc(mrc)

      # Create a barplot mapping the mrc values to fill.
      ggplot(mrc, aes(x = reorder(Group, mrc), y = mrc, fill = mrc)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(option = "viridis", direction = -1) +
        labs(x = "Group", y = "mrc",
             title = paste0("Mean connectivities to ", input$referenceCompound, " resistant cells")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  combinationScoreDF_customTCS <- reactiveVal()
  output$customTCS_combinationScoreBar <- renderPlot({
    if (input$perturbationButton_customTCS >= 1){
      mrc <- mrc_custTCS()
      res <- limmaRes_custTCS()
      mer <- merge(mrc, res, by.x = "Group", by.y = "row.names")
      mer <- subset(mer, mer$logFC < 0 & mer$mrc < 0 & adj.P.Val < 0.05)
      mer$combinationScore <- mer$logFC * mer$mrc
      mer_sub <- mer[,c("Group", "mrc", "logFC", "combinationScore")]
      combinationScoreDF_customTCS(mer_sub)

      ggplot(mer, aes(x = reorder(Group, combinationScore), y = combinationScore, fill = combinationScore)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(option = "viridis", direction = -1) +
        labs(x = "Group", y = "mrc",
             title = paste0("scFOCAL combinations score for ", input$Custom_TCS_nameInput)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  combinationScoreDF <- reactiveVal()
  output$combinationScoreBar <- renderPlot({
    if (input$perturbationButton >= 1){
      mrc <- mrc()
      res <- limmaRes()
      mer <- merge(mrc, res, by.x = "Group", by.y = "row.names")
      mer <- subset(mer, mer$logFC < 0 & mer$mrc < 0 & adj.P.Val < 0.05)
      mer$combinationScore <- mer$logFC * mer$mrc
      mer_sub <- mer[,c("Group", "mrc", "logFC", "combinationScore")]
      combinationScoreDF(mer_sub)

      ggplot(mer, aes(x = reorder(Group, combinationScore), y = combinationScore, fill = combinationScore)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(option = "viridis", direction = -1) +
        labs(x = "Group", y = "mrc",
             title = paste0("scFOCAL combinations score for ", input$referenceCompound)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  output$customTCS_CombinationHeatmap <- renderPlot({
    resistMat <- resistMat_custTCS()
    ranking <- combinationScoreDF_customTCS()
    resCells <- resistantCells_customTCS()
    rank_top <- ranking[order(ranking$combinationScore, decreasing = T)[1:30],]
    print(head(rank_top))
    compound_list <- rank_top$Group
    print(head(rownames(resistMat)))
    print(head(colnames(resistMat)))
    print(head(intersect(rownames(resistMat), compound_list)))
    print(head(intersect(rownames(resistMat), compound_list)))
    mat_list <- resistMat[which(rownames(resistMat) %in% compound_list),
                          which(colnames(resistMat) %in% resCells)]
    CombinationScore <- rowAnnotation(bar = anno_barplot(
      as.numeric(rank_top$combinationScore)),
      annotation_label = "Combination Score")
    p <- Heatmap(mat_list,
                 name = paste0(input$Custom_TCS_nameInput, " Combination Scoring"),
                 left_annotation = CombinationScore, show_column_names = FALSE,
                 row_names_rot = 30, row_names_gp = gpar(fontsize = 10))
    p
  })

  output$CombinationHeatmap <- renderPlot({
    resistMat <- resistMat()
    ranking <- combinationScoreDF()
    resCells <- resistantCells()
    rank_top <- ranking[order(ranking$combinationScore, decreasing = T)[1:30],]
    print(head(rank_top))
    compound_list <- rank_top$Group
    print(head(rownames(resistMat)))
    print(head(colnames(resistMat)))
    print(head(intersect(rownames(resistMat), compound_list)))
    print(head(intersect(rownames(resistMat), compound_list)))
    mat_list <- resistMat[which(rownames(resistMat) %in% compound_list),
                          which(colnames(resistMat) %in% resCells)]
    CombinationScore <- rowAnnotation(bar = anno_barplot(
      as.numeric(rank_top$combinationScore)),
      annotation_label = "Combination Score")
    p <- Heatmap(mat_list,
                 name = paste0(input$referenceCompound, " Combination Scoring"),
                 left_annotation = CombinationScore, show_column_names = FALSE,
                 row_names_rot = 30, row_names_gp = gpar(fontsize = 10))
    p
  })

  output$customTCS_CombinationHeatmap_2 <- renderPlot({
    resistMat <- resistMat_custTCS()
    ranking <- combinationScoreDF_customTCS()
    resCells <- resistantCells_customTCS()
    rank_top <- ranking[order(ranking$combinationScore, decreasing = T),]
    print(head(rank_top))
    compound_list <- rank_top$Group
    print(head(rownames(resistMat)))
    print(head(colnames(resistMat)))
    print(head(intersect(rownames(resistMat), compound_list)))
    print(head(intersect(rownames(resistMat), compound_list)))
    mat_list <- resistMat[which(rownames(resistMat) %in% compound_list),
                          which(colnames(resistMat) %in% resCells)]
    CombinationScore <- rowAnnotation(bar = anno_barplot(
      as.numeric(rank_top$combinationScore)),
      annotation_label = "Combination Score")
    p <- Heatmap(mat_list,
                 name = paste0(input$Custom_TCS_nameInput, " Combination Scoring"),
                 left_annotation = CombinationScore, show_column_names = FALSE,
                 show_row_names = F, cluster_rows = F)
    # row_names_rot = 30, row_names_gp = gpar(fontsize = 10))
    p
  })

  output$CombinationHeatmap_2 <- renderPlot({
    resistMat <- resistMat()
    ranking <- combinationScoreDF()
    resCells <- resistantCells()
    rank_top <- ranking[order(ranking$combinationScore, decreasing = T),]
    print(head(rank_top))
    compound_list <- rank_top$Group
    print(head(rownames(resistMat)))
    print(head(colnames(resistMat)))
    print(head(intersect(rownames(resistMat), compound_list)))
    print(head(intersect(rownames(resistMat), compound_list)))
    mat_list <- resistMat[which(rownames(resistMat) %in% compound_list),
                          which(colnames(resistMat) %in% resCells)]
    CombinationScore <- rowAnnotation(bar = anno_barplot(
      as.numeric(rank_top$combinationScore)),
      annotation_label = "Combination Score")
    p <- Heatmap(mat_list,
                 name = paste0(input$referenceCompound, " Combination Scoring"),
                 left_annotation = CombinationScore, show_column_names = FALSE,
                 show_row_names = F, cluster_rows = F)
    # row_names_rot = 30, row_names_gp = gpar(fontsize = 10))
    p
  })

  output$customTCS_combinationScoreTable <- DT::renderDataTable({
    DT::datatable(combinationScoreDF_customTCS())
  })

  output$combinationScoreTable <- DT::renderDataTable({
    DT::datatable(combinationScoreDF())
  })

  observeEvent(input$customTCS_combinationTableButton, {
    showModal(modalDialog(title = "Combination Scoring Data",
                          tags$head(tags$style(".modal-dialog{ width:85% !important;}")),
                          DT::dataTableOutput("customTCS_combinationScoreTable"),
                          br(),
                          "set up download option!!",
                          easyClose = TRUE,
                          footer = modalButton("Close")
    ))
  })

  observeEvent(input$combinationTableButton, {
    showModal(modalDialog(title = "Combination Scoring Data",
                          tags$head(tags$style(".modal-dialog{ width:85% !important;}")),
                          DT::dataTableOutput("combinationScoreTable"),
                          br(),
                          "set up download option!!",
                          easyClose = TRUE,
                          footer = modalButton("Close")
    ))
  })

  output$combinationScoreDownload_customTCS <- downloadHandler(
    filename = function() {
      paste(input$Custom_TCS_nameInput, "_asRef_combinationScoreTable.csv", sep = "")
    },
    content = function(file) {
      write.csv(combinationScoreDF_customTCS(), file, row.names = FALSE)
    }
  )

  output$combinationScoreDownload <- downloadHandler(
    filename = function() {
      paste(input$referenceCompound, "_asRef_combinationScoreTable.csv", sep = "")
    },
    content = function(file) {
      write.csv(combinationScoreDF(), file, row.names = FALSE)
    }
  )


  #################################################################################################
  # Calculate differential expression between sensitive and resistant cell populations
  #################################################################################################

  resVsSensDEGS <- reactiveVal()
  observeEvent(input$calcSensVsResistantDEGS, {
    RDSseurat <- seurat_corradded()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
    tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
    compoundSpearmans <- RDSseurat@meta.data[input$referenceCompound]
    compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)
    living <- subset(compoundSpearmans,
                     subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff_RDS)
    dead <- subset(compoundSpearmans,
                   subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff_RDS)
    livingCells <- rownames(living)
    deadCells <- rownames(dead)
    tumorObj <- subset(RDSseurat, cells = tumorCells)
    tumorObj <- SetIdent(tumorObj, cells = livingCells, value = 'resistant')
    tumorObj <- SetIdent(tumorObj, cells = deadCells, value = 'sensitive')
    print(head(Idents(tumorObj)))
    ResVsSensDEGS <- FindMarkers(tumorObj,
                                 ident.1 = 'resistant',
                                 ident.2 = 'sensitive',
                                 test.use = input$pertDEx_testUse,
                                 only.pos = F,
                                 max.cells.per.ident = input$cellSubsetMax_pert
    )
    resVsSensDEGS(ResVsSensDEGS)
  })

  # volcano plot
  output$resVsSensVolcano <- renderPlot({
    req(resVsSensDEGS())
    ResVsSensDEGS <- resVsSensDEGS()
    EnhancedVolcano(toptable = ResVsSensDEGS,
                    lab = rownames(ResVsSensDEGS),
                    x = 'avg_log2FC', y = "p_val_adj", title = NULL, subtitle = NULL,
                    ylim = c(0, max(-log10(ResVsSensDEGS[["p_val_adj"]]), na.rm=TRUE) + 0.2)
    )
  })

  # datatable
  output$resVsSensTable <- DT::renderDataTable({
    resVsSensDEGS()
  })

  # results download
  output$resVsSensDEGSdownload <- downloadHandler(
    filename = function() {
      paste(input$referenceCompound, "_resVsSensitive_DEGS.csv", sep = "")
    },
    content = function(file) {
      write.csv(resVsSensDEGS(), file, row.names = TRUE)
    }
  )

  # conditional events
  output$resVsSensCalced <- reactive({
    return(!is.null(resVsSensDEGS()))
  })
  outputOptions(output, 'resVsSensCalced', suspendWhenHidden = FALSE)

  # Custom TCS
  ##############################################################################

  resVsSensDEGS_customTCS <- reactiveVal()
  observeEvent(input$calcSensVsResistantDEGS_customTCS, {
    RDSseurat <- seurRDS_customTCSsensitivity()
    RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
    tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
    tumorObj <- subset(RDSseurat, cells = tumorCells)
    Idents(tumorObj) <- tumorObj$sensitivity
    print(head(Idents(tumorObj)))
    ResVsSensDEGS <- FindMarkers(tumorObj,
                                 ident.1 = 'resistant', latent.vars = c(input$subject_replicate_ident),
                                 ident.2 = 'sensitive',
                                 test.use = input$pertDEx_testUse_customTCS,
                                 only.pos = F,
                                 max.cells.per.ident = input$cellSubsetMax_pert_customTCS
    )
    resVsSensDEGS_customTCS(ResVsSensDEGS)
  })

  # volcano plot
  output$resVsSensVolcano_customTCS <- renderPlot({
    req(resVsSensDEGS_customTCS())
    ResVsSensDEGS <- resVsSensDEGS_customTCS()
    EnhancedVolcano(toptable = ResVsSensDEGS,
                    lab = rownames(ResVsSensDEGS),
                    x = 'avg_log2FC', y = "p_val_adj",
                    ylim = c(0, max(-log10(ResVsSensDEGS[["p_val_adj"]]), na.rm=TRUE) + 0.2)
    )
  })

  # datatable
  output$resVsSensTable_customTCS <- DT::renderDataTable({
    resVsSensDEGS_customTCS()
  })

  # results download
  output$resVsSensDEGSdownload_customTCS <- downloadHandler(
    filename = function() {
      paste(input$Custom_TCS_nameInput, "_resVsSensitive_DEGS.csv", sep = "")
    },
    content = function(file) {
      write.csv(resVsSensDEGS_customTCS(), file, row.names = TRUE)
    }
  )

  # conditional events
  output$resVsSensCalced_customTCS <- reactive({
    return(!is.null(resVsSensDEGS_customTCS()))
  })
  outputOptions(output, 'resVsSensCalced_customTCS', suspendWhenHidden = FALSE)


  #####

  #################################################################################################
  #
  #################################################################################################

  # output$ISOSCOLES_sensitivityUMAP <- renderPlot({
  #   if (input$RankSecondaryCompounds >=1){ # How do we make this re-analyze when we change cutoff?
  #
  #     ###########################################################################################
  #
  #     # Create a connected barplot comparing two groups of cells in terms of shift of proportions
  #
  #     ###########################################################################################
  #
  #     connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {
  #       b <- barplot(dat, col=color, space = space, ...)
  #
  #       for (i in seq_len(ncol(dat) - 1)) {
  #         lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line
  #
  #         for (j in seq_len(nrow(dat))) {
  #           if (j == 1) {
  #             lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))
  #             polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
  #                     c(0, dat[j,i], dat[j,i+1], 0),
  #                     col=adjustcolor(color[j], alpha.f=alpha))
  #           }
  #           if (j == 2) {
  #             lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
  #             polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
  #                     c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),
  #                     col=adjustcolor(color[j], alpha.f=alpha))
  #           }
  #           if (j > 2) {
  #             lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))
  #             polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
  #                     c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1],
  #                       colSums(dat[1:(j-1),])[i+1]),
  #                     col=adjustcolor(color[j], alpha.f=alpha))
  #           }
  #         }
  #       }
  #     }
  #
  #     # need to determine sensitive and resistant cells from previous calculation I think.
  #     tumorCells <- tumorCells()
  #     livingCells <- livingCells()
  #     deadCells <- deadCells()
  #     # p1 <- DimPlot(testdat, reduction = "umap",
  #     #               cells.highlight = intersect(deadCells, tumorCells), sizes.highlight = 0.001, cols.highlight = c("blue"))
  #     # p1 <- p1 + ggtitle("Sensitive Cells") + NoLegend()
  #     # p1
  #
  #     # Need to assign sensitivity or resistance as a new identity in @meta.data
  #     livingDF <- as.data.frame(livingCells)
  #     livingDF$sensitivity <- "Resistant"
  #     colnames(livingDF) <- c("cellID", "sensitivity")
  #
  #     deadDF <- as.data.frame(deadCells)
  #     deadDF$sensitivity <- "Sensitive"
  #     colnames(deadDF) <- c("cellID", "sensitivity")
  #
  #     comparisonDF <- rbind(livingDF, deadDF)
  #     rownames(comparisonDF) <- comparisonDF$cellID
  #     comparisonDF$cellID <- NULL
  #
  #
  #     testdat <- SetIdent(testdat, value = input$groupBy)
  #     # tumordat <- subset(testdat, idents = input$cancerCellIdents)
  #     testdat <- AddMetaData(testdat, comparisonDF)
  #     DimPlot(testdat, group.by = "sensitivity", reduction = input$reductionUse)
  #   }
  # })

  ###########################################################################################

  # Create a connected barplot comparing two groups of cells in terms of shift of proportions

  ###########################################################################################

  output$scFOCAL_StateShiftAlluvial <- renderPlot({
    if (input$perturbationButton >= 1){ # How do we make this re-analyze when we change cutoff?
      # need to determine sensitive and resistant cells from previous calculation I think.
      # tumorCells <- tumorCells()
      # livingCells <- livingCellsRDS()
      # deadCells <- deadCellsRDS()
      RDSseurat <- seurat_corradded()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
      compoundSpearmans <- RDSseurat@meta.data[input$referenceCompound]
      compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)

      living <- subset(compoundSpearmans,
                       subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff_RDS)
      dead <- subset(compoundSpearmans,
                     subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff_RDS)
      livingCells <- rownames(living)
      deadCells <- rownames(dead)

      # Need to assign sensitivity or resistance as a new identity in @meta.data
      livingDF <- as.data.frame(livingCells)
      print("livingDF")
      print(livingDF)
      livingDF$sensitivity <- "Resistant"
      colnames(livingDF) <- c("cellID", "sensitivity")

      deadDF <- as.data.frame(deadCells)
      deadDF$sensitivity <- "Sensitive"
      colnames(deadDF) <- c("cellID", "sensitivity")

      comparisonDF <- rbind(livingDF, deadDF)
      rownames(comparisonDF) <- comparisonDF$cellID
      comparisonDF$cellID <- NULL

      print(head(comparisonDF))


      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      RDSseurat <- AddMetaData(RDSseurat, comparisonDF)

      groupingCounts <- table(RDSseurat@meta.data[tumorCells,input$groupByRDS],
                              RDSseurat@meta.data[tumorCells,"sensitivity"])
      sensitiveSum <- sum(groupingCounts[,"Sensitive"])
      resistantSum <- sum(groupingCounts[,"Resistant"])
      groupingCounts2 <- as.data.frame(groupingCounts)

      colnames(groupingCounts2) <- c("groupBy", "sensitivity", "Freq")
      for (i in 1:length(groupingCounts2$Freq)){
        if (groupingCounts2$sensitivity[i] == "Sensitive"){
          groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/sensitiveSum*100
        }
        if (groupingCounts2$sensitivity[i] == "Resistant"){
          groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/resistantSum*100
        }
      }

      # groupingCounts2 <- groupingCounts2[which(is.na(groupingCounts2$sensitivity) != TRUE),]

      groupingCounts2$groupBy <- as.factor(groupingCounts2$groupBy)
      groupingCounts2$sensitivity <- factor(groupingCounts2$sensitivity, levels = c("Sensitive", "Resistant"))

      library(ggalluvial)
      # is_alluvia_form(groupingCounts2)
      ggplot(groupingCounts2,
             aes(x = sensitivity, stratum = groupBy, alluvium = groupBy,
                 y = Percent,
                 fill = groupBy, label = groupBy)) +
        scale_fill_brewer(type = "qual", palette = "Spectral") +
        # scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
        geom_flow(stat = "alluvium", lode.guidance = "frontback",
                  color = "darkgray") +
        geom_stratum() +
        theme(legend.position = "bottom") + theme_minimal() +
        ggtitle("Predicted cell type/state shift", subtitle = "Proportions of cell type/state split on drug connectivity")

    }
  })

  output$scFOCAL_StateShiftAlluvial_customTCS <- renderPlot({
    if (input$perturbationButton_customTCS >= 1){ # How do we make this re-analyze when we change cutoff?
      # need to determine sensitive and resistant cells from previous calculation I think.
      # tumorCells <- tumorCells()
      # livingCells <- livingCellsRDS()
      # deadCells <- deadCellsRDS()
      RDSseurat <- seurRDS_customTCSsensitivity()
      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
      compoundSpearmans <- RDSseurat@meta.data[input$Custom_TCS_nameInput]
      compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)

      # living <- subset(compoundSpearmans,
      #                  subset = compoundSpearmans[input$referenceCompound] > input$correlationCutoff_RDS)
      # dead <- subset(compoundSpearmans,
      #                subset = compoundSpearmans[input$referenceCompound] < input$correlationCutoff_RDS)
      # livingCells <- rownames(living)
      # deadCells <- rownames(dead)

      livingCells <- resistantCells_customTCS()
      deadCells <- sensitiveCells_customTCS()

      # Need to assign sensitivity or resistance as a new identity in @meta.data
      livingDF <- as.data.frame(livingCells)
      print("livingDF")
      print(livingDF)
      livingDF$sensitivity <- "Resistant"
      colnames(livingDF) <- c("cellID", "sensitivity")

      deadDF <- as.data.frame(deadCells)
      deadDF$sensitivity <- "Sensitive"
      colnames(deadDF) <- c("cellID", "sensitivity")

      comparisonDF <- rbind(livingDF, deadDF)
      rownames(comparisonDF) <- comparisonDF$cellID
      comparisonDF$cellID <- NULL

      print(head(comparisonDF))


      RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
      RDSseurat <- AddMetaData(RDSseurat, comparisonDF)

      groupingCounts <- table(RDSseurat@meta.data[tumorCells,input$groupByRDS],
                              RDSseurat@meta.data[tumorCells,"sensitivity"])
      sensitiveSum <- sum(groupingCounts[,"Sensitive"])
      resistantSum <- sum(groupingCounts[,"Resistant"])
      groupingCounts2 <- as.data.frame(groupingCounts)

      colnames(groupingCounts2) <- c("groupBy", "sensitivity", "Freq")
      for (i in 1:length(groupingCounts2$Freq)){
        if (groupingCounts2$sensitivity[i] == "Sensitive"){
          groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/sensitiveSum*100
        }
        if (groupingCounts2$sensitivity[i] == "Resistant"){
          groupingCounts2$Percent[i] <- groupingCounts2$Freq[i]/resistantSum*100
        }
      }

      # groupingCounts2 <- groupingCounts2[which(is.na(groupingCounts2$sensitivity) != TRUE),]

      groupingCounts2$groupBy <- as.factor(groupingCounts2$groupBy)
      groupingCounts2$sensitivity <- factor(groupingCounts2$sensitivity, levels = c("Sensitive", "Resistant"))

      library(ggalluvial)
      # is_alluvia_form(groupingCounts2)
      ggplot(groupingCounts2,
             aes(x = sensitivity, stratum = groupBy, alluvium = groupBy,
                 y = Percent,
                 fill = groupBy, label = groupBy)) +
        scale_fill_brewer(type = "qual", palette = "Spectral") +
        # scale_fill_discrete(guide = guide_legend(reverse = TRUE)) +
        geom_flow(stat = "alluvium", lode.guidance = "frontback",
                  color = "darkgray") +
        geom_stratum() +
        theme(legend.position = "bottom") + theme_minimal() +
        ggtitle("Predicted cell type/state shift", subtitle = "Proportions of cell type/state split on drug connectivity")

    }
  })


  # Create paired boxplot comparisons for stats....
  # output$customTCS_pertBoxPlots <- renderPlot({
  #   if (input$perturbationButton_customTCS >= 1){
  #     RDSseurat <- seurRDS_sensitivity()
  #     RDSseurat <- SetIdent(RDSseurat, value = input$groupByRDS)
  #     tumorCells <- WhichCells(RDSseurat, idents = input$cancerCellIdentsRDS)
  #     compoundSpearmans <- RDSseurat@meta.data[input$Custom_TCS_nameInput]
  #     compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)
  #
  #     livingCells <- resistantCells_customTCS()
  #     deadCells <- sensitiveCells_customTCS()
  #
  #     # Need to assign sensitivity or resistance as a new identity in @meta.data
  #     livingDF <- as.data.frame(livingCells)
  #     livingDF$sensitivity <- "Resistant"
  #     colnames(livingDF) <- c("cellID", "sensitivity")
  #     deadDF <- as.data.frame(deadCells)
  #     deadDF$sensitivity <- "Sensitive"
  #     colnames(deadDF) <- c("cellID", "sensitivity")
  #     comparisonDF <- rbind(livingDF, deadDF)
  #     rownames(comparisonDF) <- comparisonDF$cellID
  #     comparisonDF$cellID <- NULL
  #   }
  # })

  optimal_quantiles_calc <- reactive({
    req(
      seurat_corradded(), #
      input$referenceCompound, #
      input$subject_replicate_ident,
      input$groupByRDS
    )
    if (input$perturbationButton >= 1){
      obj <- seurRDSsensitivity()
      obj <- SetIdent(obj, value = input$groupByRDS)
      tumorCells <- WhichCells(obj, idents = input$cancerCellIdentsRDS)
      # filter for tumor cells only...

      # These are the column names you want to use, stored as strings
      score_col <- input$referenceCompound
      subject_col <- input$subject_replicate_ident
      group_col <- input$groupByRDS

      plot_data <- obj@meta.data[tumorCells, c(score_col, subject_col, group_col)]
      quantiles_to_test <- seq(0.05, 0.45, by = 0.05)
      iterative_results <- data.frame()
      non_zero_scores <- plot_data[[score_col]][plot_data[[score_col]] != 0]

      for (q in quantiles_to_test) {
        thresholds <- quantile(non_zero_scores, probs = c(q, 1 - q), na.rm = TRUE)
        temp_grouped_data <- plot_data %>%
          mutate(compound_group = case_when(
            .data[[score_col]] <= thresholds[1] ~ "Low",
            .data[[score_col]] >= thresholds[2] ~ "High",
            TRUE ~ "Middle"
          ))

        temp_proportions_df <- temp_grouped_data %>%
          filter(compound_group %in% c("Low", "High")) %>%
          # Use !!sym() to correctly evaluate the column names
          group_by(!!sym(subject_col), compound_group, !!sym(group_col)) %>%
          tally(name = "cell_count") %>%
          ungroup() %>%
          group_by(!!sym(subject_col), compound_group) %>%
          mutate(proportion = cell_count / sum(cell_count)) %>%
          ungroup() %>%
          # Also correct it in complete()
          complete(!!sym(subject_col), compound_group, !!sym(group_col), fill = list(proportion = 0, cell_count = 0))

        test_results <- temp_proportions_df %>%
          group_by(!!sym(group_col)) %>%
          pivot_wider(names_from = compound_group, values_from = proportion) %>%
          do({
            test_res <- tryCatch({
              # Check if there's enough data for a paired test
              if(sum(!is.na(.$Low) & !is.na(.$High)) > 1) {
                wilcox.test(.$Low, .$High, paired = TRUE)
              } else {
                list(p.value = 1.0)
              }
            }, error = function(e) {
              list(p.value = 1.0)
            })
            data.frame(p_value = test_res$p.value)
          }) %>%
          ungroup() %>%
          mutate(quantile = q)

        iterative_results <- bind_rows(iterative_results, test_results)
      }

      # Identify and return the optimal quantile for each cell state
      optimal_quantiles <- iterative_results %>%
        group_by(!!sym(group_col)) %>%
        slice_min(order_by = p_value, n = 1)

      return(optimal_quantiles)
    }
  })


  optimal_quantiles_calc_customTCS <- reactive({
    req(
      seurRDS_customTCSsensitivity(),
      input$Custom_TCS_nameInput,
      input$subject_replicate_ident,
      input$groupByRDS
    )
    if (input$perturbationButton_customTCS >= 1){
      obj <- seurRDS_customTCSsensitivity()
      obj <- SetIdent(obj, value = input$groupByRDS)
      tumorCells <- WhichCells(obj, idents = input$cancerCellIdentsRDS)
      # filter for tumor cells only...

      # These are the column names you want to use, stored as strings
      score_col <- input$Custom_TCS_nameInput
      subject_col <- input$subject_replicate_ident
      group_col <- input$groupByRDS

      plot_data <- obj@meta.data[tumorCells, c(score_col, subject_col, group_col)]
      quantiles_to_test <- seq(0.05, 0.45, by = 0.05)
      iterative_results <- data.frame()
      non_zero_scores <- plot_data[[score_col]][plot_data[[score_col]] != 0]

      for (q in quantiles_to_test) {
        thresholds <- quantile(non_zero_scores, probs = c(q, 1 - q), na.rm = TRUE)
        temp_grouped_data <- plot_data %>%
          mutate(compound_group = case_when(
            .data[[score_col]] <= thresholds[1] ~ "Low",
            .data[[score_col]] >= thresholds[2] ~ "High",
            TRUE ~ "Middle"
          ))

        temp_proportions_df <- temp_grouped_data %>%
          filter(compound_group %in% c("Low", "High")) %>%
          # Use !!sym() to correctly evaluate the column names
          group_by(!!sym(subject_col), compound_group, !!sym(group_col)) %>%
          tally(name = "cell_count") %>%
          ungroup() %>%
          group_by(!!sym(subject_col), compound_group) %>%
          mutate(proportion = cell_count / sum(cell_count)) %>%
          ungroup() %>%
          # Also correct it in complete()
          complete(!!sym(subject_col), compound_group, !!sym(group_col), fill = list(proportion = 0, cell_count = 0))

        test_results <- temp_proportions_df %>%
          group_by(!!sym(group_col)) %>%
          pivot_wider(names_from = compound_group, values_from = proportion) %>%
          do({
            test_res <- tryCatch({
              # Check if there's enough data for a paired test
              if(sum(!is.na(.$Low) & !is.na(.$High)) > 1) {
                wilcox.test(.$Low, .$High, paired = TRUE)
              } else {
                list(p.value = 1.0)
              }
            }, error = function(e) {
              list(p.value = 1.0)
            })
            data.frame(p_value = test_res$p.value)
          }) %>%
          ungroup() %>%
          mutate(quantile = q)

        iterative_results <- bind_rows(iterative_results, test_results)
      }

      # Identify and return the optimal quantile for each cell state
      optimal_quantiles <- iterative_results %>%
        group_by(!!sym(group_col)) %>%
        slice_min(order_by = p_value, n = 1)

      return(optimal_quantiles)
    }
  })

  # Add this new observe block for debugging purposes
  observe({
    # Use the same req() checks as your main reactive
    req(
      seurRDS_customTCSsensitivity(),
      input$Custom_TCS_nameInput,
      input$subject_replicate_ident,
      input$groupBy_RDS
    )

    # Trigger this observer when the button is clicked
    req(input$perturbationButton_customTCS >= 1)

    # Prepare the data frame just for printing
    obj <- seurRDS_customTCSsensitivity()

    score_col <- input$Custom_TCS_nameInput
    subject_col <- input$subject_replicate_ident
    group_col <- input$groupBy_RDS

    # Ensure columns exist before trying to select them
    req(all(c(score_col, subject_col, group_col) %in% colnames(obj@meta.data)))

    plot_data <- obj@meta.data[, c(score_col, subject_col, group_col)]

    # cat("--- Debugging plot_data ---\n")
    # dplyr::glimpse(plot_data)
    # cat("---------------------------\n\n")

  })

  # observe({
  #   # Use the same req() checks as your main reactive
  #   req(
  #     seurat_corradded(),
  #     input$referenceCompound,
  #     input$subject_replicate_ident,
  #     input$groupBy_RDS
  #   )
  #
  #   # Trigger this observer when the button is clicked
  #   req(input$perturbationButton >= 1)
  #
  #   # Prepare the data frame just for printing
  #   obj <- seurat_corradded()
  #
  #   score_col <- input$referenceCompound
  #   subject_col <- input$subject_replicate_ident
  #   group_col <- input$groupBy_RDS
  #
  #   # Ensure columns exist before trying to select them
  #   req(all(c(score_col, subject_col, group_col) %in% colnames(obj@meta.data)))
  #
  #   plot_data <- obj@meta.data[, c(score_col, subject_col, group_col)]
  #
  #   # This will now print to your R console every time the inputs change
  #   cat("--- Debugging plot_data ---\n")
  #   dplyr::glimpse(plot_data)
  #   cat("---------------------------\n\n")
  #
  # })

  ##############################################################################

  # output$optimal_quantiles_output <- DT::renderDT({
  #   optimal_quantiles_calc()
  # }, options = list(pageLength = 5), rownames = FALSE)

  # 4. Main renderPlot() function for visualization
  output$main_plot <- renderPlot({

    # --- DATA PREPARATION ---
    # Use the same req() checks to ensure all inputs are ready
    req(
      seurat_corradded(),
      input$referenceCompound,
      input$subject_replicate_ident,
      input$groupByRDS,
      input$quantile_slider # Also require the slider input
    )

    # Ensure the button has been clicked before trying to plot
    req(input$perturbationButton >= 1)

    obj <- seurat_corradded()
    obj <- SetIdent(obj, value = input$groupByRDS)
    tumorCells <- WhichCells(obj, idents = input$cancerCellIdentsRDS)

    # Get column names from the correct inputs
    score_col <- input$referenceCompound
    subject_col <- input$subject_replicate_ident
    group_col <- input$groupByRDS

    # Extract metadata, filtering for tumor cells
    plot_data <- obj@meta.data[tumorCells, c(score_col, subject_col, group_col)]

    # Get the quantile from the slider input
    q_selected <- input$quantile_slider / 100

    non_zero_scores <- plot_data[[score_col]][plot_data[[score_col]] != 0]

    # Calculate thresholds based on the selected quantile
    thresholds <- quantile(non_zero_scores, probs = c(q_selected, 1 - q_selected), na.rm = TRUE)
    q_low_threshold <- thresholds[1]
    q_high_threshold <- thresholds[2]

    # Create grouping column with descriptive names
    low_group_name <- paste0("Low (Bottom ", q_selected * 100, "%)")
    high_group_name <- paste0("High (Top ", q_selected * 100, "%)")

    grouped_data_for_plot <- plot_data %>%
      mutate(compound_group = case_when(
        .data[[score_col]] <= q_low_threshold ~ low_group_name,
        .data[[score_col]] >= q_high_threshold ~ high_group_name,
        TRUE ~ "Middle"
      ))

    proportions_for_plot_df <- grouped_data_for_plot %>%
      filter(compound_group != "Middle") %>%
      group_by(!!sym(subject_col), compound_group, !!sym(group_col)) %>%
      tally(name = "cell_count") %>%
      ungroup() %>%
      group_by(!!sym(subject_col), compound_group) %>%
      mutate(proportion = cell_count / sum(cell_count)) %>%
      ungroup() %>%
      complete(!!sym(subject_col), compound_group, !!sym(group_col), fill = list(proportion = 0, cell_count = 0)) %>%
      filter(compound_group != "Middle") %>%
      mutate(compound_group = factor(compound_group, levels = c(low_group_name, high_group_name)))

    # --- PLOT GENERATION ---

    # Plot 1: Histogram of scores
    color_map <- setNames(c("blue", "red", "gray"), c(low_group_name, high_group_name, "Middle"))
    h <- ggplot(grouped_data_for_plot, aes(x = .data[[score_col]], fill = compound_group)) +
      geom_histogram(bins = 100, alpha = 0.8) +
      scale_fill_manual(values = color_map, name = "Quantile Group") +
      labs(title = paste("Distribution of", score_col, "Scores"), x = paste(score_col, "Score"), y = "Cell Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

    # Plot 2: Boxplot of proportions
    p <- ggplot(proportions_for_plot_df, aes(x = compound_group, y = proportion, fill = compound_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_line(aes(group = .data[[subject_col]]), color = "gray50", alpha = 0.8) + # Corrected group aesthetic
      geom_point(size = 2.5, aes(color = compound_group)) +
      stat_compare_means(method = "wilcox.test", label = "p.format", paired = TRUE) +
      facet_wrap(vars(!!sym(group_col)), scales = "free_y", nrow = 1) + # Corrected facet_wrap
      labs(title = paste("Proportions at Quantile Cutoff:", q_selected * 100, "%"), x = NULL, y = "Proportion") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # Plot 3: Difference plot
    diff_df <- proportions_for_plot_df %>%
      select(!!sym(subject_col), compound_group, !!sym(group_col), proportion) %>%
      pivot_wider(names_from = compound_group, values_from = proportion) %>%
      mutate(prop_diff = (.data[[high_group_name]] - .data[[low_group_name]]) * 100) %>%
      filter(.data[[group_col]] != "Non-Neoplastic")

    d <- ggplot(diff_df, aes(x = .data[[group_col]], y = prop_diff, fill = .data[[group_col]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      labs(title = "Difference in Proportion", subtitle = "(High - Low)", x = NULL, y = "Difference (% points)") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # --- COMBINE PLOTS AND RETURN ---

    # Combine the plots into a single figure using cowplot
    c1 <- cowplot::plot_grid(h, d, ncol = 2, rel_widths = c(1,1))
    combined_row <- cowplot::plot_grid(c1, p, ncol = 1)

    # Return the final combined plot to be rendered
    return(combined_row)

  }, res = 110) # Set resolution for a sharper plot


  output$optimal_quantiles_output_customTCS <- DT::renderDT({
    optimal_quantiles_calc_customTCS()
  }, options = list(pageLength = 5), rownames = FALSE)

  # 4. Main renderPlot() function for visualization
  output$main_plot_customTCS <- renderPlot({

    # --- DATA PREPARATION ---
    # Use the same req() checks to ensure all inputs are ready
    req(
      seurRDS_customTCSsensitivity(),
      input$Custom_TCS_nameInput,
      input$subject_replicate_ident,
      input$groupByRDS,
      input$quantile_slider_customTCS # Also require the slider input
    )

    # Ensure the button has been clicked before trying to plot
    req(input$perturbationButton_customTCS >= 1)

    obj <- seurRDS_customTCSsensitivity()
    obj <- SetIdent(obj, value = input$groupByRDS)
    tumorCells <- WhichCells(obj, idents = input$cancerCellIdentsRDS)

    # Get column names from the correct inputs
    score_col <- input$Custom_TCS_nameInput
    subject_col <- input$subject_replicate_ident
    group_col <- input$groupByRDS

    # Extract metadata, filtering for tumor cells
    plot_data <- obj@meta.data[tumorCells, c(score_col, subject_col, group_col)]

    # Get the quantile from the slider input
    q_selected <- input$quantile_slider_customTCS / 100

    non_zero_scores <- plot_data[[score_col]][plot_data[[score_col]] != 0]

    # Calculate thresholds based on the selected quantile
    thresholds <- quantile(non_zero_scores, probs = c(q_selected, 1 - q_selected), na.rm = TRUE)
    q_low_threshold <- thresholds[1]
    q_high_threshold <- thresholds[2]

    # Create grouping column with descriptive names
    low_group_name <- paste0("Low (Bottom ", q_selected * 100, "%)")
    high_group_name <- paste0("High (Top ", q_selected * 100, "%)")

    grouped_data_for_plot <- plot_data %>%
      mutate(compound_group = case_when(
        .data[[score_col]] <= q_low_threshold ~ low_group_name,
        .data[[score_col]] >= q_high_threshold ~ high_group_name,
        TRUE ~ "Middle"
      ))

    proportions_for_plot_df <- grouped_data_for_plot %>%
      filter(compound_group != "Middle") %>%
      group_by(!!sym(subject_col), compound_group, !!sym(group_col)) %>%
      tally(name = "cell_count") %>%
      ungroup() %>%
      group_by(!!sym(subject_col), compound_group) %>%
      mutate(proportion = cell_count / sum(cell_count)) %>%
      ungroup() %>%
      complete(!!sym(subject_col), compound_group, !!sym(group_col), fill = list(proportion = 0, cell_count = 0)) %>%
      filter(compound_group != "Middle") %>%
      mutate(compound_group = factor(compound_group, levels = c(low_group_name, high_group_name)))

    # --- PLOT GENERATION ---

    # Plot 1: Histogram of scores
    color_map <- setNames(c("blue", "red", "gray"), c(low_group_name, high_group_name, "Middle"))
    h <- ggplot(grouped_data_for_plot, aes(x = .data[[score_col]], fill = compound_group)) +
      geom_histogram(bins = 100, alpha = 0.8) +
      scale_fill_manual(values = color_map, name = "Quantile Group") +
      labs(title = paste("Distribution of", score_col, "Scores"), x = paste(score_col, "Score"), y = "Cell Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

    # Plot 2: Boxplot of proportions
    p <- ggplot(proportions_for_plot_df, aes(x = compound_group, y = proportion, fill = compound_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_line(aes(group = .data[[subject_col]]), color = "gray50", alpha = 0.8) + # Corrected group aesthetic
      geom_point(size = 2.5, aes(color = compound_group)) +
      stat_compare_means(method = "wilcox.test", label = "p.format", paired = TRUE) +
      facet_wrap(vars(!!sym(group_col)), scales = "free_y", nrow = 1) + # Corrected facet_wrap
      labs(title = paste("Proportions at Quantile Cutoff:", q_selected * 100, "%"), x = NULL, y = "Proportion") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # Plot 3: Difference plot
    diff_df <- proportions_for_plot_df %>%
      select(!!sym(subject_col), compound_group, !!sym(group_col), proportion) %>%
      pivot_wider(names_from = compound_group, values_from = proportion) %>%
      mutate(prop_diff = (.data[[high_group_name]] - .data[[low_group_name]]) * 100) %>%
      filter(.data[[group_col]] != "Non-Neoplastic")

    d <- ggplot(diff_df, aes(x = .data[[group_col]], y = prop_diff, fill = .data[[group_col]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      labs(title = "Difference in Proportion", subtitle = "(High - Low)", x = NULL, y = "Difference (% points)") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # --- COMBINE PLOTS AND RETURN ---

    # Combine the plots into a single figure using cowplot
    c1 <- cowplot::plot_grid(h, d, ncol = 2, rel_widths = c(1,1))
    combined_row <- cowplot::plot_grid(c1, p, ncol = 1)

    # Return the final combined plot to be rendered
    return(combined_row)

  }, res = 110) # Set resolution for a sharper plot

  ##

  # 4. Main renderPlot() function for visualization
  output$main_plot_userSet <- renderPlot({

    # --- DATA PREPARATION ---
    # Use the same req() checks, but for the new value slider
    req(
      seurat_corradded(),
      input$referenceCompound,
      input$subject_replicate_ident,
      input$groupByRDS,
      input$value_slider
    )

    # Ensure the button has been clicked before trying to plot
    req(input$perturbationButton >= 1)

    obj <- seurat_corradded()
    obj <- SetIdent(obj, value = input$groupByRDS)
    tumorCells <- WhichCells(obj, idents = input$cancerCellIdentsRDS)

    # Get column names from the correct inputs
    score_col <- input$Custom_TCS_nameInput
    subject_col <- input$subject_replicate_ident
    group_col <- input$groupByRDS

    # Extract metadata, filtering for tumor cells
    plot_data <- obj@meta.data[tumorCells, c(score_col, subject_col, group_col)]

    q_low_threshold <- input$correlationCutoff[1]
    q_high_threshold <- input$correlationCutoff[2]

    # Create grouping column with descriptive names based on values
    low_group_name <- paste0("Score <= ", round(q_low_threshold, 2))
    high_group_name <- paste0("Score >= ", round(q_high_threshold, 2))
    # --------------------------------------------------------------------

    # Use .data pronoun for the score column
    grouped_data_for_plot <- plot_data %>%
      mutate(compound_group = case_when(
        .data[[score_col]] <= q_low_threshold ~ low_group_name,
        .data[[score_col]] >= q_high_threshold ~ high_group_name,
        TRUE ~ "Middle"
      ))

    proportions_for_plot_df <- grouped_data_for_plot %>%
      filter(compound_group != "Middle") %>%
      group_by(!!sym(subject_col), compound_group, !!sym(group_col)) %>%
      tally(name = "cell_count") %>%
      ungroup() %>%
      group_by(!!sym(subject_col), compound_group) %>%
      mutate(proportion = cell_count / sum(cell_count)) %>%
      ungroup() %>%
      complete(!!sym(subject_col), compound_group, !!sym(group_col), fill = list(proportion = 0, cell_count = 0)) %>%
      filter(compound_group != "Middle") %>%
      mutate(compound_group = factor(compound_group, levels = c(low_group_name, high_group_name)))

    # --- PLOT GENERATION ---

    # Plot 1: Histogram of scores
    color_map <- setNames(c("blue", "red", "gray"), c(low_group_name, high_group_name, "Middle"))
    h <- ggplot(grouped_data_for_plot, aes(x = .data[[score_col]], fill = compound_group)) +
      geom_histogram(bins = 100, alpha = 0.8) +
      scale_fill_manual(values = color_map, name = "Cutoff Group") +
      labs(title = paste("Distribution of", score_col, "Scores"), x = paste(score_col, "Score"), y = "Cell Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

    # Plot 2: Boxplot of proportions
    p <- ggplot(proportions_for_plot_df, aes(x = compound_group, y = proportion, fill = compound_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_line(aes(group = .data[[subject_col]]), color = "gray50", alpha = 0.8) +
      geom_point(size = 2.5, aes(color = compound_group)) +
      stat_compare_means(method = "wilcox.test", label = "p.format", paired = TRUE) +
      facet_wrap(vars(!!sym(group_col)), scales = "free_y", nrow = 1) +
      labs(title = paste("Proportions by Score Cutoff"), x = NULL, y = "Proportion") + # Updated title
      theme_bw() +
      theme(strip.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # Plot 3: Difference plot (no major changes needed here)
    diff_df <- proportions_for_plot_df %>%
      select(!!sym(subject_col), compound_group, !!sym(group_col), proportion) %>%
      pivot_wider(names_from = compound_group, values_from = proportion) %>%
      mutate(prop_diff = (.data[[high_group_name]] - .data[[low_group_name]]) * 100) %>%
      filter(.data[[group_col]] != "Non-Neoplastic")

    d <- ggplot(diff_df, aes(x = .data[[group_col]], y = prop_diff, fill = .data[[group_col]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      labs(title = "Difference in Proportion", subtitle = "(High - Low)", x = NULL, y = "Difference (% points)") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # --- COMBINE PLOTS AND RETURN ---
    c1 <- cowplot::plot_grid(h, d, ncol = 2, rel_widths = c(1,1))
    combined_row <- cowplot::plot_grid(c1, p, ncol = 1)
    return(combined_row)

  }, res = 110)

  output$main_plot_userSet_customTCS <- renderPlot({

    # --- DATA PREPARATION ---
    # Use the same req() checks, but for the new value slider
    req(
      seurRDS_customTCSsensitivity(),
      input$Custom_TCS_nameInput,
      input$subject_replicate_ident,
      input$groupByRDS,
      input$value_slider)

    # Ensure the button has been clicked before trying to plot
    req(input$perturbationButton_customTCS >= 1)

    obj <- seurRDS_customTCSsensitivity()
    obj <- SetIdent(obj, value = input$groupByRDS)
    tumorCells <- WhichCells(obj, idents = input$cancerCellIdentsRDS)

    # Get column names from the correct inputs
    score_col <- input$Custom_TCS_nameInput
    subject_col <- input$subject_replicate_ident
    group_col <- input$groupByRDS

    # Extract metadata, filtering for tumor cells
    plot_data <- obj@meta.data[tumorCells, c(score_col, subject_col, group_col)]

    q_low_threshold <- input$correlationCutoff_customTCS[1]
    q_high_threshold <- input$correlationCutoff_customTCS[2]

    # Create grouping column with descriptive names based on values
    low_group_name <- paste0("Score <= ", round(q_low_threshold, 2))
    high_group_name <- paste0("Score >= ", round(q_high_threshold, 2))
    # --------------------------------------------------------------------

    # Use .data pronoun for the score column
    grouped_data_for_plot <- plot_data %>%
      mutate(compound_group = case_when(
        .data[[score_col]] <= q_low_threshold ~ low_group_name,
        .data[[score_col]] >= q_high_threshold ~ high_group_name,
        TRUE ~ "Middle"
      ))

    proportions_for_plot_df <- grouped_data_for_plot %>%
      filter(compound_group != "Middle") %>%
      group_by(!!sym(subject_col), compound_group, !!sym(group_col)) %>%
      tally(name = "cell_count") %>%
      ungroup() %>%
      group_by(!!sym(subject_col), compound_group) %>%
      mutate(proportion = cell_count / sum(cell_count)) %>%
      ungroup() %>%
      complete(!!sym(subject_col), compound_group, !!sym(group_col), fill = list(proportion = 0, cell_count = 0)) %>%
      filter(compound_group != "Middle") %>%
      mutate(compound_group = factor(compound_group, levels = c(low_group_name, high_group_name)))

    # --- PLOT GENERATION ---

    # Plot 1: Histogram of scores
    color_map <- setNames(c("blue", "red", "gray"), c(low_group_name, high_group_name, "Middle"))
    h <- ggplot(grouped_data_for_plot, aes(x = .data[[score_col]], fill = compound_group)) +
      geom_histogram(bins = 100, alpha = 0.8) +
      scale_fill_manual(values = color_map, name = "Cutoff Group") +
      labs(title = paste("Distribution of", score_col, "Scores"), x = paste(score_col, "Score"), y = "Cell Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

    # Plot 2: Boxplot of proportions
    p <- ggplot(proportions_for_plot_df, aes(x = compound_group, y = proportion, fill = compound_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_line(aes(group = .data[[subject_col]]), color = "gray50", alpha = 0.8) +
      geom_point(size = 2.5, aes(color = compound_group)) +
      stat_compare_means(method = "wilcox.test", label = "p.format", paired = TRUE) +
      facet_wrap(vars(!!sym(group_col)), scales = "free_y", nrow = 1) +
      labs(title = paste("Proportions by Score Cutoff"), x = NULL, y = "Proportion") + # Updated title
      theme_bw() +
      theme(strip.text = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # Plot 3: Difference plot (no major changes needed here)
    diff_df <- proportions_for_plot_df %>%
      select(!!sym(subject_col), compound_group, !!sym(group_col), proportion) %>%
      pivot_wider(names_from = compound_group, values_from = proportion) %>%
      mutate(prop_diff = (.data[[high_group_name]] - .data[[low_group_name]]) * 100) %>%
      filter(.data[[group_col]] != "Non-Neoplastic")

    d <- ggplot(diff_df, aes(x = .data[[group_col]], y = prop_diff, fill = .data[[group_col]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      labs(title = "Difference in Proportion", subtitle = "(High - Low)", x = NULL, y = "Difference (% points)") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # --- COMBINE PLOTS AND RETURN ---
    combined_row <- cowplot::plot_grid(h, p, d, ncol = 3, rel_widths = c(1, 5, 1.2))
    return(combined_row)

  }, res = 110)

}
