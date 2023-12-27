# User Interface ####
ui <- dashboardPage(
 dashboardHeader(title = "Multi-Omics Integration", titleWidth = 400,
                 ... = tags$li(class = "dropdown", actionLink(inputId = "testButton", label = "Stop Application"))),
 # Sidebar ####
 dashboardSidebar(width = 400, collapsed = FALSE,
                  sidebarMenu(
                   # Datasets Menu ####
                   menuItem("Datasets", tabName = "datasets", icon = icon("table")),
                   # Clustering Menu ####
                   menuItem("Clustering", icon = icon("project-diagram"),
                            menuItem("Bayesian Consensus Clustering", tabName = "clusterSub1"),
                            menuItem("iClusterPlus", tabName = "clusterSub2"),
                            menuItem("LRAcluster", tabName = "clusterSub3"),
                            menuItem("PINSPlus", tabName = "clusterSub4"),
                            menuItem("Similarity Network Fusion", tabName = "clusterSub5")),
                   # Analysis Menu ####
                   menuItem("Analysis", icon = icon("dna"),
                            menuItem("MAF Summary",
                                     menuItem("Bayesian Consensus Clustering", tabName = "mafSummaryBcc"),
                                     menuItem("iClusterPlus", tabName = "mafSummaryicp"),
                                     menuItem("LRAcluster", tabName = "mafSummaryLra"),
                                     menuItem("PINSPlus", tabName = "mafSummaryPins"),
                                     menuItem("Similarity Network Fusion", tabName = "mafSummarySnf")),
                            menuItem("Oncoplots", 
                                     menuItem("Bayesian Consensus Clustering", tabName = "oncoSub1"),
                                     menuItem("iClusterPlus", tabName = "oncoSub2"),
                                     menuItem("LRAcluster", tabName = "oncoSub3"),
                                     menuItem("PINSPlus", tabName = "oncoSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "oncoSub5")),
                            menuItem("Recurrent CNV Plots",
                                     menuItem("Bayesian Consensus Clustering", tabName = "recCnvSub1"),
                                     menuItem("iClusterPlus", tabName = "recCnvSub2"),
                                     menuItem("LRAcluster", tabName = "recCnvSub3"),
                                     menuItem("PINSPlus", tabName = "recCnvSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "recCnvSub5")),
                            menuItem("Permutation Tests",
                                     menuItem("Bayesian Consensus Clustering", tabName = "pTestsSub1"),
                                     menuItem("iClusterPlus", tabName = "pTestsSub2"),
                                     menuItem("LRAcluster", tabName = "pTestsSub3"),
                                     menuItem("PINSPlus", tabName = "pTestsSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "pTestsSub5")),
                            # menuItem("Differentially Expressed Genes",
                            #          menuItem("Bayesian Consensus Clustering", tabName = "degSub1"),
                            #          menuItem("iClusterPlus", tabName = "degSub2"),
                            #          menuItem("LRAcluster", tabName = "degSub3"),
                            #          menuItem("PINSPlus", tabName = "degSub4"),
                            #          menuItem("Similarity Network Fusion", tabName = "degSub5")),
                            menuItem("Enrichment Analysis",
                                     menuItem("Bayesian Consensus Clustering", tabName = "enrAnalysisSub1"),
                                     menuItem("iClusterPlus", tabName = "enrAnalysisSub2"),
                                     menuItem("LRAcluster", tabName = "enrAnalysisSub3"),
                                     menuItem("PINSPlus", tabName = "enrAnalysisSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "enrAnalysisSub5")),
                            menuItem("Gene Expression Heatmaps",
                                     menuItem("Bayesian Consensus Clustering", tabName = "geHeatmapsSub1"),
                                     menuItem("iClusterPlus", tabName = "geHeatmapsSub2"),
                                     menuItem("LRAcluster", tabName = "geHeatmapsSub3"),
                                     menuItem("PINSPlus", tabName = "geHeatmapsSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "geHeatmapsSub5")),
                            menuItem("Survival Analysis",
                                     menuItem("Bayesian Consensus Clustering", tabName = "survAnalysisSub1"),
                                     menuItem("iClusterPlus", tabName = "survAnalysisSub2"),
                                     menuItem("LRAcluster", tabName = "survAnalysisSub3"),
                                     menuItem("PINSPlus", tabName = "survAnalysisSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "survAnalysisSub5")),
                            menuItem("Clinical Statistics",
                                     menuItem("Bayesian Consensus Clustering", tabName = "clinicStatsSub1"),
                                     menuItem("iClusterPlus", tabName = "clinicStatsSub2"),
                                     menuItem("LRAcluster", tabName = "clinicStatsSub3"),
                                     menuItem("PINSPlus", tabName = "clinicStatsSub4"),
                                     menuItem("Similarity Network Fusion", tabName = "clinicStatsSub5")))
                  )),
 # Body Pages ####
 dashboardBody(
  tabItems(
   # 1. Dataset query page ####
   tabItem(tabName = "datasets",
           h1("Query Datasets"),
           fluidRow(
            box(title = "Input Box", width = 12,
                fluidRow(
                 column(12, h4("Copy Number Variation")),
                 column(4, selectInput(inputId = "catCnv", label = "Data Category", choices = c("Copy Number Variation"))),
                 column(4, selectInput(inputId = "dType", label = "Data Type", choices = c("Masked Copy Number Segment"))),
                 column(4, selectInput(inputId = "sType", label = "Sample Type", choices = c("Primary Tumor")))
                ), hr(),
                fluidRow(
                 column(12, h4("RNA")),
                 column(4, selectInput(inputId = "catRna", label = "Data Category", choices = c("Transcriptome Profiling"))),
                 column(4, selectInput(inputId = "expStrat", label = "Exp. Str.", choices = c("RNA-Seq"))),
                 column(4, selectInput(inputId = "workType", label = "Workflow Type", choices = c("HTSeq - Counts")))
                ), hr(),
                fluidRow(
                 column(12, h4("Protein")),
                 column(4, selectInput(inputId = "catProtein", label = "Data Category", choices = c("Protein expression"))),
                ), hr(),
                fluidRow(
                 column(4,
                        actionButton(inputId = "button", "Get Dataset")
                 )
                )
            )
           ),
           fluidRow(
            tabBox(title = "Preview Data", width = 12,
                   tabPanel(title = "CNV",
                            dataTableOutput("cnvTable")
                   ),
                   tabPanel(title = "RNA",
                            dataTableOutput("rnaTable")
                   ),
                   tabPanel(title = "Protein",
                            dataTableOutput("proteinTable")
                   ),
                   tabPanel(title = "Clinical",
                            dataTableOutput("clinicalTable")
                   )
            )
           )
   ),
   # 2. Clustering Pages ####
   # 2.1. BCC ####
   tabItem(tabName = "clusterSub1",
           fluidRow(box(width = 12, status = "primary", h1("Bayesian Consensus Clustering"))),
           fluidRow(
            box(width = 4,
                h4("Preparing Datasets"),
                fluidRow(
                 column(width = 6, numericInput(inputId = "bayesEpsilon", label = "Epsilon", value = 0.0015, min = 0, max = 0.1, step = 0.0001)),
                 column(width = 6, numericInput(inputId = "bayesQuantileCut", label = "Quantile Cut-Off", value = 0.78, min = 0, max = 1))
                ),
                hr(),
                h4("Optimal K Value"),
                fluidRow(
                 column(width = 6, checkboxInput(inputId = "bayesCheckBox1", label = "Individual Alpha", value = TRUE)),
                 column(width = 6, sliderInput(inputId = "maxIter1", label = "Max Iteration", min = 1000, value = 1000, max = 10000, step = 1000))
                ),
                hr(),
                fluidRow(
                 column(width = 12,
                        actionButton(inputId = "bayesOptButton", label = "Calculate")
                 ))),
            box(title = "Boxplot of Mean-Adjusted Adherence Vectors", width = 8,
                fluidRow(
                 column(width = 12, plotOutput(outputId = "bayesOptimalK")),
                ),
                hr(),
                fluidRow(
                 column(width = 12, verbatimTextOutput(outputId = "bayesText"))
                ))),
           fluidRow(
            box(title = "", width = 4,
                fluidRow(
                 column(width = 12, sliderInput(inputId = "numberOfK", label = "# of Clusters", min = 2, value = 2, max = 5)),
                ),
                fluidRow(
                 column(width = 12, sliderInput(inputId = "maxIter2", label = "Max Iteration", min = 1000, value = 1000, max = 10000, step = 1000)),
                ),
                fluidRow(
                 column(width = 12, actionButton(inputId = "bayesRunButton", label = "Run"))
                ),
                fluidRow(
                   column(width = 6, sliderInput(inputId = "bayesPhi", label = "Phi", min = 0, max = 90, value = 40)),
                   column(width = 6, sliderInput(inputId = "bayesTheta", label = "Theta", min = 0, max = 90, value = 40))
                )
            ),
            tabBox(width = 8, title = "3D Clustering Plots",
                   tabPanel(title = "CNV", plotOutput(outputId = "bayesCnv", height = 700)),
                   tabPanel(title = "RNA", plotOutput(outputId = "bayesRna", height = 700)),
                   tabPanel(title = "Protein", plotOutput(outputId = "bayesProtein", height = 700))
            ))),
   # 2.2. iClusterPlus ####
   tabItem(tabName = "clusterSub2",
           fluidRow(box(width = 12, status = "primary", h1("iClusterPlus"))),
           fluidRow(
            box(title = "Optimal K Value", width = 4,
                h4("Preparing Datasets"),
                fluidRow(
                 column(width = 12, numericInput(inputId = "iClusterInput3", label = "Epsilon", value = 0.0015, min = 0, max = 0.1, step = 0.0001))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "iClusterInput4", label = "Quantile Cut-Off", value = 0.78, min = 0, max = 1))),
                fluidRow(
                 column(width = 12,
                        actionButton(inputId = "iClusterOptButton", label = "Calculate")))),
            box(title = "Deviance Ratio (Percentage of Explained Variation)", width = 8,
                fluidRow(
                 column(width = 12, plotOutput(outputId = "iClusterOut1"))))),
           fluidRow(
            box(title = "Configure Method", width = 4,
                fluidRow(
                 column(width = 12, numericInput(inputId = "iClusterInput2", label = "# of Clusters (K+1)", value = 2, max = 10)),
                 column(width = 12, actionButton(inputId = "iClusterRunButton", label = "Run")))),
            tabBox(title = "3D Clustering Plot", width = 8,
                   tabPanel(title = "CNV", plotOutput(outputId = "iClusterTabCnv", height = 700)),
                   tabPanel(title = "RNA", plotOutput(outputId = "iClusterTabRna", height = 700)),
                   tabPanel(title = "Protein", plotOutput(outputId = "iClusterTabProtein", height = 700))
            ))),
   # 2.3. LRAcluster ####
   tabItem(tabName = "clusterSub3",
           fluidRow(box(width = 12, status = "primary", h1("LRAcluster"))),
           fluidRow(
            box(title = "Preparing Datasets", width = 4,
                fluidRow(
                 column(width = 12, numericInput(inputId = "lraInput3", label = "Epsilon", value = 0.0015, min = 0, max = 0.1, step = 0.0001))
                ),
                fluidRow(
                 column(width = 12, numericInput(inputId = "lraInput4", label = "Quantile Cut-Off", value = 0.78, min = 0, max = 1))
                ),
                fluidRow(
                 column(width = 12,
                        actionButton(inputId = "lraOptButton", label = "Calculate")
                 ))),
            box(title = "Silhouette Scores - 1", width = 4,
                fluidRow(
                 column(width = 12, plotOutput(outputId = "lraOut1"))
                )),
            box(title = "Silhouette Scores - 2", width = 4,
                fluidRow(
                 column(width = 12, plotOutput(outputId = "lraOut2")),
                ))),
           fluidRow(
            box(title = "Configure Method", width = 4,
                fluidRow(
                 column(width = 12, sliderInput(inputId = "lraInput2", label = "# of Clusters", min = 2, value = 2, max = 5)),
                ),
                fluidRow(
                 column(width = 12, actionButton(inputId = "lraRunButton", label = "Run"))
                )),
            tabBox(width = 8, title = "3D Clustering Plots", 
                   tabPanel(title = "CNV", plotOutput(outputId = "lraOut4", height = 700)),
                   tabPanel(title = "RNA", plotOutput(outputId = "lraOut5", height = 700)),
                   tabPanel(title = "Protein", plotOutput(outputId = "lraOut6", height = 700))
            ))),
   # 2.4. PINSPlus ####
   tabItem(tabName = "clusterSub4",
           fluidRow(box(width = 12, status = "primary", h1("PINSPlus"))),
           fluidRow(
            box(title = "Preparing Datasets", width = 4,
                fluidRow(
                 column(width = 12, numericInput(inputId = "pinsplusInput1", label = "Epsilon", value = 0.0015, min = 0, max = 0.1, step = 0.0001))
                ),
                fluidRow(
                 column(width = 12, numericInput(inputId = "pinsplusInput2", label = "Quantile Cut-Off", value = 0.78, min = 0, max = 1))
                ),
                hr(),
                h4("Method Configuration"),
                fluidRow(
                 column(width = 12, selectInput(inputId = "pinsplusCpus", label = "# of CPUs", choices = c(1, 2, 3, 4)))
                ),
                fluidRow(
                 column(width = 12, actionButton(inputId = "pinsplusRun", label = "Run"))
                )),
            # box(title = "Plot Clustering", width = 8, plotOutput(outputId = "pinsplusPlot"))
            tabBox(width = 8, title = "3D Clustering Plots", 
                   tabPanel(title = "CNV", plotOutput(outputId = "pinsplusOut4", height = 700)),
                   tabPanel(title = "RNA", plotOutput(outputId = "pinsplusOut5", height = 700)),
                   tabPanel(title = "Protein", plotOutput(outputId = "pinsplusOut6", height = 700))
            ))),
   # 2.5. SNF ####
   tabItem(tabName = "clusterSub5",
           fluidRow(box(width = 12, status = "primary", h1("Similarity Network Fusion"))),
           fluidRow(
            box(title = "Preparing Datasets", width = 4,
                # fluidRow(
                #   column(width = 12, sliderInput(inputId = "lraInput1", label = "Select Dataset Size", value = 361, min = 5, max = 361)),
                # ),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput1", label = "Epsilon", value = 0.0015, min = 0, max = 0.1, step = 0.0001))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput2", label = "Quantile Cut-Off", value = 0.78, min = 0, max = 1))),
                hr(),
                h4("Optimal K"),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput3", label = "Number of Neighbours(K)", value = 20, min = 10, max = 30, step = 1))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput4", label = "Hyperparameter(Alpha)", min = 0.3, value = 0.5, max = 0.8, step = 0.1))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput5", label = "Number of Iterations(T)", min = 10, value = 20, max = 20, step = 1))),
                fluidRow(
                 column(width = 12, actionButton(inputId = "snfOptButton", label = "Calculate")))),
            column(width = 4, valueBoxOutput(outputId = "snfBestEigenvalue", width = 12)),
            column(width = 4, valueBoxOutput(outputId = "snfRotCost", width = 12)),
            column(width = 4, valueBoxOutput(outputId = "snf2Eigenvalue", width = 12)),
            column(width = 4, valueBoxOutput(outputId = "snf2RotCost", width = 12))
           ),
           fluidRow(
            box(title = "Configure Method", width = 4,
                fluidRow(
                 column(width = 12, sliderInput(inputId = "snfInput6", label = "# of Clusters", min = 2, value = 3, max = 5))),
                fluidRow(
                 column(width = 12, actionButton(inputId = "snfRunButton", label = "Run"))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput7", label = "Number of Neighbours(K)", value = 20, min = 10, max = 30, step = 1))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput8", label = "Hyperparameter(Alpha)", min = 0.3, value = 0.5, max = 0.8, step = 0.1))),
                fluidRow(
                 column(width = 12, numericInput(inputId = "snfInput9", label = "Number of Iterations(T)", min = 10, value = 20, max = 20, step = 1)))
            ),
            # box(title = "Clustering Plot", width = 8,
            #     fluidRow(
            #       column(width = 12, plotOutput(outputId = "snfOut3"))
            #     ))
            tabBox(width = 8, title = "3D Clustering Plots", 
                   tabPanel(title = "CNV", plotOutput(outputId = "snfOut4", height = 700)),
                   tabPanel(title = "RNA", plotOutput(outputId = "snfOut5", height = 700)),
                   tabPanel(title = "Protein", plotOutput(outputId = "snfOut6", height = 700))
            ))),
   # 3. Analysis Pages ####
   # 3.1. MAF Summary ####
   tabItem(tabName = "mafSummaryBcc",
           fluidRow(box(width = 12, status = "primary", h3("MAF Summary \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "mafSummaryBccCluster1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "mafSummaryBccCluster2", height = "auto"))
                           )))),
   tabItem(tabName = "mafSummaryicp",
           fluidRow(box(width = 12, status = "primary", h3("MAF Summary \\ iClusterPlus"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", height = "auto", width = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "mafSummaryicpCluster1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "mafSummaryicpCluster2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "mafSummaryicpCluster3", height = "auto"))
                           )))),
   tabItem(tabName = "mafSummaryLra",
           fluidRow(box(width = 12, status = "primary", h3("MAF Summary \\ LRAluster"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", height = "auto", width = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "mafSummaryLraCluster1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "mafSummaryLraCluster2", height = "auto"))
                           )))),
   tabItem(tabName = "mafSummaryPins",
           fluidRow(box(width = 12, status = "primary", h3("MAF Summary \\ PINSPlus"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", height = "auto", width = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "mafSummaryPinsCluster1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "mafSummaryPinsCluster2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "mafSummaryPinsCluster3", height = "auto")),
                                  tabPanel(title = "Cluster 4", imageOutput(outputId = "mafSummaryPinsCluster4", height = "auto"))
                           )))),
   tabItem(tabName = "mafSummarySnf",
           fluidRow(box(width = 12, status = "primary", h3("MAF Summary \\ Similarity Network Fusion"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", height = "auto", width = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "mafSummarySnfCluster1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "mafSummarySnfCluster2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "mafSummarySnfCluster3", height = "auto"))
                           )))),
   # 3.2. Oncoplots ####
   tabItem(tabName = "oncoSub1",
           fluidRow(box(width = 12, status = "primary", h3("Oncoplots \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "oncoplotBcc1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "oncoplotBcc2", height = "auto"))
                           )))),
   tabItem(tabName = "oncoSub2",
           fluidRow(box(width = 12, status = "primary", h3("Oncoplots \\ iClusterPlus"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "oncoploticp1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "oncoploticp2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "oncoploticp3", height = "auto"))
                            )))),
   tabItem(tabName = "oncoSub3",
           fluidRow(box(width = 12, status = "primary", h3("Oncoplots \\ LRAluster"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "oncoplotLra1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "oncoplotLra2", height = "auto"))
                           )))),
   tabItem(tabName = "oncoSub4",
           fluidRow(box(width = 12, status = "primary", h3("Oncoplots \\ PINSPlus"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "oncoplotPins1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "oncoplotPins2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "oncoplotPins3", height = "auto")),
                                  tabPanel(title = "Cluster 4", imageOutput(outputId = "oncoplotPins4", height = "auto"))
                           )))),
   tabItem(tabName = "oncoSub5",
           fluidRow(box(width = 12, status = "primary", h3("Oncoplots \\ Similarity Network Fusion"))),
           fluidRow(column(width = 7, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "oncoplotSnf1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "oncoplotSnf2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "oncoplotSnf3", height = "auto"))
                           )))),
   # 3.3. Recurrent CNV ####
   tabItem(tabName = "recCnvSub1",
           fluidRow(box(width = 12, status = "primary", h3("Recurrent CNV Plots \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "recCnvBcc1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "recCnvBcc2", height = "auto"))
                                  )))),
   tabItem(tabName = "recCnvSub2",
           fluidRow(box(width = 12, status = "primary", h3("Recurrent CNV Plots \\ iClusterPlus"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "recCnvicp1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "recCnvicp2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "recCnvicp3", height = "auto"))
                                  )))),
   tabItem(tabName = "recCnvSub3",
           fluidRow(box(width = 12, status = "primary", h3("Recurrent CNV Plots \\ LRAluster"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "recCnvLra1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "recCnvLra2", height = "auto"))
                                  )))),
   tabItem(tabName = "recCnvSub4",
           fluidRow(box(width = 12, status = "primary", h3("Recurrent CNV Plots \\ PINSPlus"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "recCnvPins1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "recCnvPins2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "recCnvPins3", height = "auto")),
                                  tabPanel(title = "Cluster 4", imageOutput(outputId = "recCnvPins4", height = "auto"))
                                  )))),
   tabItem(tabName = "recCnvSub5",
           fluidRow(box(width = 12, status = "primary", h3("Recurrent CNV Plots \\ Similarity Network Fusion"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "recCnvSnf1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "recCnvSnf2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "recCnvSnf3", height = "auto"))
                                  )))),
   # 3.4. Permutation Tests ####
   tabItem(tabName = "pTestsSub1",
           fluidRow(box(width = 12, status = "primary", h3("Permutation Tests \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 6, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "pTestsBcc1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "pTestsBcc2", height = "auto"))
                           )))),
   tabItem(tabName = "pTestsSub2",
           fluidRow(box(width = 12, status = "primary", h3("Permutation Tests \\ iClusterPlus"))),
           fluidRow(column(width = 6, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "pTestsicp1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "pTestsicp2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "pTestsicp3", height = "auto"))
                           )))),
   tabItem(tabName = "pTestsSub3",
           fluidRow(box(width = 12, status = "primary", h3("Permutation Tests \\ LRAluster"))),
           fluidRow(column(width = 6, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "pTestsLra1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "pTestsLra2", height = "auto"))
                           )))),
   tabItem(tabName = "pTestsSub4",
           fluidRow(box(width = 12, status = "primary", h3("Permutation Tests \\ PINSPlus"))),
           fluidRow(column(width = 6, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "pTestsPins1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "pTestsPins2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "pTestsPins3", height = "auto")),
                                  tabPanel(title = "Cluster 4", imageOutput(outputId = "pTestsPins4", height = "auto"))
                           )))),
   tabItem(tabName = "pTestsSub5",
           fluidRow(box(width = 12, status = "primary", h3("Permutation Tests \\ Similarity Network Fusion"))),
           fluidRow(column(width = 6, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "pTestsSnf1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "pTestsSnf2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "pTestsSnf3", height = "auto"))
                           )))),
   # 3.5. Differentially Expressed Genes ####
   # 3.6. Enrichment Analysis ####
   tabItem(tabName = "enrAnalysisSub1",
           fluidRow(box(width = 12, status = "primary", h3("Enrichment Analysis \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 12,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "enrAnalysisBcc1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "enrAnalysisBcc2", height = "auto"))
                           )))),
   tabItem(tabName = "enrAnalysisSub2",
           fluidRow(box(width = 12, status = "primary", h3("Enrichment Analysis \\ iClusterPlus"))),
           fluidRow(column(width = 12,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "enrAnalysisicp1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "enrAnalysisicp2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "enrAnalysisicp3", height = "auto"))
                           )))),
   tabItem(tabName = "enrAnalysisSub3",
           fluidRow(box(width = 12, status = "primary", h3("Enrichment Analysis \\ LRAluster"))),
           fluidRow(column(width = 12,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "enrAnalysisLra1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "enrAnalysisLra2", height = "auto"))
                           )))),
   tabItem(tabName = "enrAnalysisSub4",
           fluidRow(box(width = 12, status = "primary", h3("Enrichment Analysis \\ PINSPlus"))),
           fluidRow(column(width = 12,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "enrAnalysisPins1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "enrAnalysisPins2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "enrAnalysisPins3", height = "auto")),
                                  tabPanel(title = "Cluster 4", imageOutput(outputId = "enrAnalysisPins4", height = "auto"))
                           )))),
   tabItem(tabName = "enrAnalysisSub5",
           fluidRow(box(width = 12, status = "primary", h3("Enrichment Analysis \\ Similarity Network Fusion"))),
           fluidRow(column(width = 12,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "enrAnalysisSnf1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "enrAnalysisSnf2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "enrAnalysisSnf3", height = "auto"))
                           )))),
   # 3.7. Gene Expression Heatmaps ####
   tabItem(tabName = "geHeatmapsSub1",
           fluidRow(box(width = 12, status = "primary", h3("Gene Expression Heatmaps \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 12,
                           box(id = "", title = "(2 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "geHeatmapsBcc", height = "auto"),
                           )))),
   tabItem(tabName = "geHeatmapsSub2",
           fluidRow(box(width = 12, status = "primary", h3("Gene Expression Heatmaps \\ iClusterPlus"))),
           fluidRow(column(width = 12,
                           box(id = "", title = "(3 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "geHeatmapsicp", height = "auto"),
                           )))),
   tabItem(tabName = "geHeatmapsSub3",
           fluidRow(box(width = 12, status = "primary", h3("Gene Expression Heatmaps \\ LRAluster"))),
           fluidRow(column(width = 12,
                           box(id = "", title = "(2 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "geHeatmapsLra", height = "auto"),
                           )))),
   tabItem(tabName = "geHeatmapsSub4",
           fluidRow(box(width = 12, status = "primary", h3("Gene Expression Heatmaps \\ PINSPlus"))),
           fluidRow(column(width = 12,
                           box(id = "", title = "(4 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "geHeatmapsPins", height = "auto"),
                           )))),
   tabItem(tabName = "geHeatmapsSub5",
           fluidRow(box(width = 12, status = "primary", h3("Gene Expression Heatmaps \\ Similarity Network Fusion"))),
           fluidRow(column(width = 12,
                           box(id = "", title = "(3 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "geHeatmapsSnf", height = "auto"),
                           )))),
   # 3.8. Survival Analysis ####
   tabItem(tabName = "survAnalysisSub1",
           fluidRow(box(width = 12, status = "primary", h3("Survival Analysis \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 8, offset = 2,
                           box(id = "", title = "(2 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "survAnalysisBcc", height = "auto"),
                           )))),
   tabItem(tabName = "survAnalysisSub2",
           fluidRow(box(width = 12, status = "primary", h3("Survival Analysis \\ iClusterPlus"))),
           fluidRow(column(width = 8, offset = 2,
                           box(id = "", title = "(3 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "survAnalysisicp", height = "auto"),
                           )))),
   tabItem(tabName = "survAnalysisSub3",
           fluidRow(box(width = 12, status = "primary", h3("Survival Analysis \\ LRAluster"))),
           fluidRow(column(width = 8, offset = 2,
                           box(id = "", title = "(2 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "survAnalysisLra", height = "auto"),
                           )))),
   tabItem(tabName = "survAnalysisSub4",
           fluidRow(box(width = 12, status = "primary", h3("Survival Analysis \\ PINSPlus"))),
           fluidRow(column(width = 8, offset = 2,
                           box(id = "", title = "(4 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "survAnalysisPins", height = "auto"),
                           )))),
   tabItem(tabName = "survAnalysisSub5",
           fluidRow(box(width = 12, status = "primary", h3("Survival Analysis \\ Similarity Network Fusion"))),
           fluidRow(column(width = 8, offset = 2,
                           box(id = "", title = "(3 Clusters)", width = "auto", height = "auto",
                               imageOutput(outputId = "survAnalysisSnf", height = "auto"),
                           )))),
   # 3.9. Clinical Statistics ####
   tabItem(tabName = "clinicStatsSub1",
           fluidRow(box(width = 12, status = "primary", h3("Clinical Statistics \\ Bayesian Consensus Cluster"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1",
                                  fluidRow(column(width = 6,
                                          imageOutput(outputId = "clinicStatsBcc1", height = "auto"),
                                          imageOutput(outputId = "clinicStatsBcc2", height = "auto")),
                                  column(width = 6,
                                         imageOutput(outputId = "clinicStatsBcc3", height = "auto"),
                                         imageOutput(outputId = "clinicStatsBcc4", height = "auto")))),
                                   tabPanel(title = "Cluster 2",
                                   fluidRow(column(width = 6,
                                          imageOutput(outputId = "clinicStatsBcc5", height = "auto"),
                                          imageOutput(outputId = "clinicStatsBcc6", height = "auto"),),
                                          column(width = 6,
                                          imageOutput(outputId = "clinicStatsBcc7", height = "auto"))))
                                   )))),
   tabItem(tabName = "clinicStatsSub2",
           fluidRow(box(width = 12, status = "primary", h3("Clinical Statistics \\ iClusterPlus"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "clinicStatsicp1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "clinicStatsicp2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "clinicStatsicp3", height = "auto"))
                           )))),
   tabItem(tabName = "clinicStatsSub3",
           fluidRow(box(width = 12, status = "primary", h3("Clinical Statistics \\ LRAluster"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "clinicStatsLra1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "clinicStatsLra2", height = "auto"))
                           )))),
   tabItem(tabName = "clinicStatsSub4",
           fluidRow(box(width = 12, status = "primary", h3("Clinical Statistics \\ PINSPlus"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "clinicStatsPins1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "clinicStatsPins2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "clinicStatsPins3", height = "auto")),
                                  tabPanel(title = "Cluster 4", imageOutput(outputId = "clinicStatsPins4", height = "auto"))
                           )))),
   tabItem(tabName = "clinicStatsSub5",
           fluidRow(box(width = 12, status = "primary", h3("Clinical Statistics \\ Similarity Network Fusion"))),
           fluidRow(column(width = 8, offset = 2,
                           tabBox(id = "", title = "Clusters", width = "auto", height = "auto",
                                  tabPanel(title = "Cluster 1", imageOutput(outputId = "clinicStatsSnf1", height = "auto")),
                                  tabPanel(title = "Cluster 2", imageOutput(outputId = "clinicStatsSnf2", height = "auto")),
                                  tabPanel(title = "Cluster 3", imageOutput(outputId = "clinicStatsSnf3", height = "auto"))
                           ))))
  )
 )
)