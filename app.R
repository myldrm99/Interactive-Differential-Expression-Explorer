library(shiny)
library(GEOquery)
library(DT)
library(ggplot2)
library(plotly)
library(limma)
library(pheatmap)
library(grid)

ui <- fluidPage(
  tags$head(tags$style(HTML("
body { background-color: #f0f2f5; }
    .app-header { background-color: #4A90E2; color: white; padding: 20px 0; text-align: center; font-size: 28px; font-weight: bold; }
        .app-subtitle { background-color: #4A90E2; color: white; text-align: center; font-style: italic; padding-bottom: 15px; margin-bottom: 20px; }
            .sidebar-panel { background-color: #F5F5F5; padding: 15px; border-radius: 4px; }
                .group-btn { display: block; width: 100%; margin-bottom: 8px; }
              .main-panel { background-color: white; padding: 15px; border-radius: 4px; box-shadow: 0 0 5px rgba(0,0,0,0.1); }
              .app-footer { background-color: #E0E0E0; padding: 10px 0; text-align: center; margin-top: 20px; }
                  .app-footer .btn, .app-footer .shiny-link { margin: 0 8px; }
                "))),
  
  # Header
  div(class="app-header",   "Interactive Differential Expression Explorer"),
  div(class="app-subtitle","Analyze GEO or User Data with Advanced Parameter Tuning"),
  
  sidebarLayout(
    sidebarPanel(class="sidebar-panel",
                 # Data source choice
                 radioButtons("data_source","Data source:",
                              choices = c("GEO accession"="geo","Upload file"="upload")),
                 conditionalPanel("input.data_source=='geo'",
                                  textInput("geo","GEO Accession:",value="GSE10072")
                 ),
                 conditionalPanel("input.data_source=='upload'",
                                  fileInput("upload","Upload expression matrix (CSV/TSV)",
                                            accept=c(".csv",".tsv"))
                 ),
                 actionButton("load","Load Data",class="btn-primary group-btn"),
                 hr(),
                 
                 # Only for GEO, allow selecting metadata columns
                 conditionalPanel("input.data_source=='geo'",
                                  uiOutput("columnSelector")
                 ),
                 
                 hr(),
                 actionButton("runAnalysis","Run Analysis",class="btn-success group-btn"),
                 actionButton("resetGroups","Reset Groups",class="btn-warning group-btn"),
                 hr(),
                 actionButton("assignG1","► Add Selected to Group 1",class="btn-info group-btn"),
                 actionButton("assignG2","► Add Selected to Group 2",class="btn-info group-btn"),
                 actionButton("clearSel","✕ Clear Selected",class="btn-danger group-btn"),
                 hr(),
                 
                 # 1) Threshold sliders
                 sliderInput("fc_thr","Log₂ Fold-Change threshold:",
                             min=0, max=5, value=1, step=0.1),
                 sliderInput("pval_thr","Adjusted P-value threshold:",
                             min=0, max=0.1, value=0.05, step=0.005),
                 
                 # Side-panel summary of DE Genes
                 hr(),
                 htmlOutput("deSummary")
    ),
    
    mainPanel(class="main-panel",
              tabsetPanel(
                tabPanel("Data Info",   verbatimTextOutput("dataInfo")),
                tabPanel("Samples",     DTOutput("samplesTable")),
                tabPanel("Volcano",     plotlyOutput("volcanoPlot", height="600px")),
                tabPanel("Scatter",     plotlyOutput("scatterPlot", height="600px")),
                tabPanel("Heatmap",     plotOutput("heatmapPlot", height="600px")),
                # ---- Added PCA tab ----
                tabPanel("PCA",         plotlyOutput("pcaPlot",    height="600px"))
              )
    )
  ),
  
  # Footer
  div(class="app-footer",
      downloadButton("downloadVolcano", "Download Volcano (PNG)"),
      downloadButton("downloadScatter","Download Scatter (PNG)"),
      downloadButton("downloadHeatmap","Download Heatmap (PNG)"),
      downloadButton("downloadPCA",    "Download PCA (PNG)"),
      downloadButton("downloadResults","Download Results (CSV)"),
      actionLink("help","Help / About", class="shiny-link")
  )
)

server <- function(input, output, session) {
  
  # — 1) Load either GEO or uploaded data on button click
  rawData <- eventReactive(input$load, {
    if (input$data_source == "upload") {
      req(input$upload)
      f <- input$upload
      ext <- tools::file_ext(f$name)
      df  <- if (ext=="csv") {
        read.csv(f$datapath, row.names=1, check.names=FALSE)
      } else {
        read.delim(f$datapath, row.names=1, check.names=FALSE)
      }
      validate(need(ncol(df)>1, "Uploaded file needs at least 2 samples"))
      return(as.matrix(df))
    } else {
      req(input$geo)
      g <- getGEO(input$geo, GSEMatrix=TRUE)
      validate(need(length(g)>0, "No series-matrix found"))
      return(g[[1]])  # ExpressionSet
    }
  })
  
  # — 2) Expose general info about the expression matrix
  exprMat <- reactive({
    obj <- rawData()
    if (inherits(obj,"ExpressionSet")) {
      mat <- exprs(obj)
    } else {
      mat <- obj
    }
    validate(need(is.matrix(mat) && nrow(mat)>0, "No numeric data available"))
    mat
  })
  
  output$dataInfo <- renderPrint({
    mat <- exprMat()
    cat("Expression matrix:\n")
    cat(" • Features (rows):", nrow(mat), "\n")
    cat(" • Samples  (cols):", ncol(mat), "\n\n")
    cat("Data type:", typeof(mat), "\n")
    cat("Contains NA?        ", any(is.na(mat)), "\n")
    cat("Contains negatives? ", any(mat<0), "\n")
    cat("Contains non-integers?", any(mat!=round(mat)), "\n\n")
    cat("Value summary:\n")
    print(summary(as.numeric(mat)))
  })
  
  # — 3) Sample metadata / names
  sampleDF <- reactive({
    obj <- rawData()
    if (inherits(obj,"ExpressionSet")) {
      df <- pData(obj)
      df$Sample <- rownames(df)
    } else {
      mat <- exprMat()
      df  <- data.frame(Sample=colnames(mat), row.names=NULL, stringsAsFactors=FALSE)
    }
    df
  })
  
  # — 4) Column selector (only for GEO)
  output$columnSelector <- renderUI({
    df <- sampleDF(); req(df)
    cols <- setdiff(names(df),"Sample")
    selectInput("show_cols","Metadata columns:",
                choices = cols,
                selected= head(cols,5),
                multiple=TRUE, width="100%")
  })
  
  # — 5) Group state
  groups <- reactiveValues(g1=character(0), g2=character(0))
  observeEvent(input$resetGroups, {
    groups$g1 <- character(0); groups$g2 <- character(0)
  })
  
  # — 6) Annotate sample table
  annotatedSamples <- reactive({
    df <- sampleDF(); req(df)
    if (!is.null(input$show_cols)) {
      df <- df[, c("Sample", input$show_cols), drop=FALSE]
    } else {
      df <- df[, "Sample", drop=FALSE]
    }
    df$Group <- ifelse(df$Sample %in% groups$g1, "Group 1",
                       ifelse(df$Sample %in% groups$g2, "Group 2", "—"))
    df
  })
  
  # — 7) Render sample table with proxy for group coloring
  output$samplesTable <- renderDT({
    datatable(
      annotatedSamples(),
      rownames   = FALSE,
      selection  = "multiple",
      extensions = c("Buttons","FixedHeader"),
      options    = list(
        dom         = "Bfrtip",
        buttons     = c("copy","csv","excel","pdf","print"),
        deferRender = TRUE,
        scrollY     = "400px",
        scrollX     = TRUE,
        paging      = FALSE,
        fixedHeader = TRUE
      ),
      class="stripe hover compact"
    ) %>% formatStyle(
      "Group", target="row",
      backgroundColor = styleEqual(
        c("Group 1","Group 2","—"),
        c("#e1f5fe","#ffebee","white")
      )
    )
  }, server=TRUE)
  
  proxy <- dataTableProxy("samplesTable")
  observe({
    replaceData(proxy, annotatedSamples(),
                resetPaging=FALSE, rownames=FALSE)
  })
  
  # — 8) Table‐driven grouping buttons
  observeEvent(input$assignG1, {
    sel <- input$samplesTable_rows_selected
    df  <- annotatedSamples()
    new <- df$Sample[sel]
    groups$g1 <- union(groups$g1, setdiff(new,groups$g2))
    groups$g2 <- setdiff(groups$g2, groups$g1)
  })
  observeEvent(input$assignG2, {
    sel <- input$samplesTable_rows_selected
    df  <- annotatedSamples()
    new <- df$Sample[sel]
    groups$g2 <- union(groups$g2, setdiff(new,groups$g1))
    groups$g1 <- setdiff(groups$g1, groups$g2)
  })
  observeEvent(input$clearSel, {
    sel <- input$samplesTable_rows_selected
    df  <- annotatedSamples()
    clr <- df$Sample[sel]
    groups$g1 <- setdiff(groups$g1, clr)
    groups$g2 <- setdiff(groups$g2, clr)
  })
  
  # — 9) Differential expression via limma
  deResults <- eventReactive(input$runAnalysis, {
    mat <- exprMat(); req(mat)
    sel <- c(groups$g1, groups$g2)
    validate(need(length(sel)>=2, "Need at least one sample in each group"))
    validate(need(all(sel%in%colnames(mat)), "Selections missing in data"))
    
    m2  <- mat[, sel, drop=FALSE]
    grp <- factor(c(rep("G1", length(groups$g1)),
                    rep("G2", length(groups$g2))))
    design <- model.matrix(~0 + grp)
    colnames(design) <- levels(grp)
    
    fit  <- lmFit(m2, design)
    ct   <- makeContrasts(G2-G1, levels=design)
    fit2 <- contrasts.fit(fit, ct) |> eBayes()
    topTable(fit2, number=Inf, sort.by="P")
  })
  
  # Side-panel summary of DE Genes
  deCounts <- reactive({
    df <- deResults()
    sig <- df$adj.P.Val <= input$pval_thr & abs(df$logFC) >= input$fc_thr
    up    <- sum(sig & df$logFC >  0)
    down  <- sum(sig & df$logFC <  0)
    total <- up + down
    list(total=total, up=up, down=down)
  })
  
  output$deSummary <- renderUI({
    req(deResults())
    cnt <- deCounts()
    tags$div(
      style = "background:#ffffff; padding:10px; border-radius:4px; box-shadow:0 0 3px rgba(0,0,0,0.1);",
      tags$strong("DE Genes Summary"), tags$br(),
      sprintf("Total: %d", cnt$total), tags$br(),
      sprintf("  Up: %d", cnt$up), tags$br(),
      sprintf("  Down: %d", cnt$down)
    )
  })
  
  # — 10) Interactive Volcano (Plotly + threshold sliders)
  output$volcanoPlot <- renderPlotly({
    df <- deResults(); req(nrow(df)>0)
    df$logP <- -log10(df$adj.P.Val)
    df$sig  <- factor(
      ifelse(df$adj.P.Val <= input$pval_thr &
               abs(df$logFC) >= input$fc_thr,
             ifelse(df$logFC>0,"Up","Down"),"NS"),
      levels=c("Up","Down","NS")
    )
    
    plot_ly(df,
            x     = ~logFC,
            y     = ~logP,
            color = ~sig,
            colors= c("Up"="red","Down"="blue","NS"="grey50"),
            text  = ~paste(
              "Gene:", rownames(df),
              "<br>log2FC:", round(logFC,2),
              "<br>adj.P.Val:", signif(adj.P.Val,3)
            ),
            hoverinfo="text",
            type="scatter", mode="markers"
    ) %>%
      layout(
        title = "Volcano Plot",
        xaxis = list(title="log₂ Fold Change"),
        yaxis = list(title="-log₁₀ Adjusted P-value")
      )
  })
  
  # — 11) Interactive MA / Scatter Plot
  output$scatterPlot <- renderPlotly({
    df <- deResults(); req(nrow(df)>0)
    df$sig  <- factor(
      ifelse(df$adj.P.Val <= input$pval_thr &
               abs(df$logFC) >= input$fc_thr,
             ifelse(df$logFC>0,"Up","Down"),"NS"),
      levels=c("Up","Down","NS")
    )
    
    plot_ly(df,
            x     = ~AveExpr,
            y     = ~logFC,
            color = ~sig,
            colors= c("Up"="red","Down"="blue","NS"="grey50"),
            text  = ~paste(
              "Gene:", rownames(df),
              "<br>AveExpr:", round(AveExpr,2),
              "<br>log2FC:", round(logFC,2),
              "<br>adj.P.Val:", signif(adj.P.Val,3)
            ),
            hoverinfo="text",
            type="scatter", mode="markers"
    ) %>%
      layout(
        title = "MA Plot",
        xaxis = list(title="Average Expression"),
        yaxis = list(title="log₂ Fold Change")
      )
  })
  
  # — 12) Heatmap
  output$heatmapPlot <- renderPlot({
    df <- deResults()
    req(nrow(df) > 0, length(groups$g1) > 0, length(groups$g2) > 0)
    
    sel <- c(groups$g1, groups$g2)
    topN <- head(rownames(df), min(50, nrow(df)))
    mat <- log2(exprMat()[topN, sel, drop = FALSE] + 1)
    
    annotation_col <- data.frame(
      Group = factor(
        c(rep("Group 1", length(groups$g1)),
          rep("Group 2", length(groups$g2)))
      ),
      row.names = sel
    )
    
    grid::grid.newpage()
    p <- pheatmap::pheatmap(
      mat,
      annotation_col = annotation_col,
      scale          = "row",
      main           = "Top 50 DE Genes",
      fontsize_col   = 10,
      silent         = TRUE
    )
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })
  
  # — 13) PCA
  pcaResults <- eventReactive(input$runAnalysis, {
    sel <- c(groups$g1, groups$g2)
    mat <- exprMat()[, sel, drop = FALSE]
    vsd <- log2(mat + 1)
    pca <- prcomp(t(vsd), scale. = TRUE)
    varPerc <- pca$sdev^2 / sum(pca$sdev^2) * 100
    dfPCA <- as.data.frame(pca$x[, 1:2])
    dfPCA$Sample <- rownames(dfPCA)
    dfPCA$Group  <- ifelse(dfPCA$Sample %in% groups$g1, "Group 1", "Group 2")
    list(df = dfPCA, var = varPerc)
  })
  
  output$pcaPlot <- renderPlotly({
    res <- pcaResults()
    df  <- res$df
    v   <- res$var
    plot_ly(df,
            x     = ~PC1, y     = ~PC2,
            color = ~Group,
            colors= c("Group 1" = "#1f77b4", "Group 2" = "#ff7f0e"),
            text  = ~paste(
              "Sample:", Sample,
              "<br>PC1:", round(PC1, 2),
              "<br>PC2:", round(PC2, 2)
            ),
            hoverinfo = "text",
            mode      = "markers"
    ) %>%
      layout(
        title = "PCA: PC1 vs PC2",
        xaxis = list(title = sprintf("PC1 (%.1f%%)", v[1])),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", v[2]))
      )
  })
  
  # — 14) Downloads
  output$downloadResults <- downloadHandler(
    filename = "differential_expression_results.csv",
    content  = function(file) write.csv(deResults(), file, row.names=TRUE)
  )
  
  output$downloadVolcano <- downloadHandler(
    filename = "volcano.png",
    content  = function(file) {
      png(file, width = 800, height = 600)
      df <- deResults(); req(nrow(df) > 0)
      df$logP <- -log10(df$P.Value)
      df$sig  <- factor(ifelse(df$P.Value < 0.05 & abs(df$logFC) > 1,
                               ifelse(df$logFC > 0, "Up","Down"), "NS"),
                        levels = c("Up","Down","NS"))
      print(
        ggplot(df, aes(logFC, logP, color = sig)) +
          geom_point(alpha = 0.6) +
          scale_color_manual(values = c(Up="red",Down="blue",NS="grey50")) +
          geom_vline(xintercept = c(-1,1), linetype="dashed") +
          geom_hline(yintercept = -log10(0.05), linetype="dashed") +
          labs(x="Log2 Fold Change", y="-Log10 P-value", color="Status") +
          theme_minimal() + theme(legend.position="bottom")
      )
      dev.off()
    }
  )
  
  output$downloadScatter <- downloadHandler(
    filename = "scatter.png",
    content  = function(file) {
      png(file, width = 800, height = 600)
      df <- deResults(); req(nrow(df) > 0)
      df$sig <- factor(ifelse(df$P.Value < 0.05 & abs(df$logFC) > 1,
                              ifelse(df$logFC > 0, "Up","Down"), "NS"),
                       levels = c("Up","Down","NS"))
      print(
        ggplot(df, aes(AveExpr, logFC, color = sig)) +
          geom_point(alpha = 0.6) +
          scale_color_manual(values = c(Up="red",Down="blue",NS="grey50")) +
          geom_hline(yintercept = 0, linetype="dashed") +
          labs(x="Average Expression", y="Log2 Fold Change", color="Status") +
          theme_minimal() + theme(legend.position="bottom")
      )
      dev.off()
    }
  )
  
  output$downloadHeatmap <- downloadHandler(
    filename = "heatmap.png",
    content  = function(file) {
      png(file, width = 800, height = 800)
      df <- deResults()
      sel <- c(groups$g1, groups$g2)
      topN  <- head(rownames(df), 50)
      mat   <- exprMat()[ topN, sel, drop = FALSE ]
      matLog<- log2(mat + 1)
      anno <- data.frame(
        Group = factor(
          c(rep("Group 1", length(groups$g1)),
            rep("Group 2", length(groups$g2)))
        )
      )
      rownames(anno) <- sel
      phm <- pheatmap(
        matLog,
        annotation_col = anno,
        scale           = "row",
        main            = "Top 50 DE Genes",
        silent          = TRUE
      )
      grid.newpage()
      grid.draw(phm$gtable)
      dev.off()
    }
  )
  
  output$downloadPCA <- downloadHandler(
    filename = "pca_plot.png",
    content = function(file) {
      png(file, width = 800, height = 600)
      res <- pcaResults()
      df  <- res$df
      v   <- res$var
      p <- ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 4, alpha = 0.8) +
        labs(
          title = "PCA: PC1 vs PC2",
          x     = sprintf("PC1 (%.1f%%)", v[1]),
          y     = sprintf("PC2 (%.1f%%)", v[2])
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
      print(p)
      dev.off()
    }
  )
  
  # — 15) Help modal
  observeEvent(input$help, {
    showModal(modalDialog(
      title = "Help / About",
      HTML("
          <ul>
            <li>Choose GEO or upload your own matrix.</li>
            <li>Select samples and assign to Group 1/Group 2.</li>
            <li>Adjust log₂FC & adjusted p-value thresholds with sliders.</li>
            <li>Explore interactive Volcano and MA plots (hover & zoom).</li>
            <li>View a static heatmap of top DE genes.</li>
            <li>Explore PCA of your selected samples.</li>
            <li>Download any plot or the full DE results.</li>
          </ul>
        "),
      easyClose = TRUE, footer = NULL
    ))
  })
}

shinyApp(ui, server)
