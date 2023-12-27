# Define server ####
server <- function(input, output, session) {
 cat("\nServer is running...\n")
 cat("\n\tSonuçların küme sayıları kontrol edilecek!..\n\n")
 bin2dec <- function(cluster) {
  vec <- integer(nrow(cluster))
  for (i in seq(ncol(cluster))) {
   vec[which(cluster[,i] == 1)] <- i
  }
  return(vec)
 }
 play.sound <- function() {
  beepr::beep("/usr/share/sounds/Yaru/stereo/complete.oga")
  Sys.sleep(0.67)
  beepr::beep("/usr/share/sounds/Yaru/stereo/complete.oga")
 }

 values.query <- reactiveValues(out = NULL)
 values.bcc <- reactiveValues(alphaStarD = NULL, data = NULL, out = NULL)
 values.lra <- reactiveValues(rlist = NULL, avg.sil = NULL, data = NULL)
 values.icluster <- reactiveValues(data = NULL, best.fit = NULL)
 values.snf <- reactiveValues(data = NULL, out = NULL, opt.k = NULL)
 values.pinsplus <- reactiveValues(data = NULL)

 observeEvent(input$testButton, {
  # ifelse(is.null(values$react), values$react <<- 0, values$react <<- values$react + 1)
  stopApp()
 })
 
 # Render data frames ####
 observeEvent(input$button, {
  cnv <- c(d.category = input$catCnv, d.type = input$dType, s.type = input$sType)
  rna <- c(d.category = input$catRna, exp.strategy = input$expStrat, w.type = input$workType)
  protein <- c(d.category = input$catProtein)
  values.query$out <<- data.analysis$run(cnv = cnv, rna = rna, protein = protein)
 })
 
 query.event <- eventReactive(values.query$out, {
  y <- values.query$out
  return(y)
 })
 
 output$cnvTable <- renderDataTable({
  x <- query.event()[["CNV"]]
  x <- as.data.frame(x)
  return(x)
 }, list(scrollX = TRUE))
 output$rnaTable <- renderDataTable({
  x <- query.event()[["RNA"]]
  x <- as.data.frame(SummarizedExperiment::assay(x))
  return(x)
 }, list(scrollX = TRUE))
 output$proteinTable <- renderDataTable({
  x <- query.event()[["Protein"]]
  x <- as.data.frame(x)
  return(x)
 }, list(scrollX = TRUE))
 output$clinicalTable <- renderDataTable({
  x <- query.event()[["Clinical"]]
  x <- as.data.frame(x)
  return(x)
 }, list(scrollX = TRUE))
 
 # Bayesian CC ####
 observeEvent(input$bayesOptButton, {
  # if (is.null(values.query$out)) return(message("No dataset queried!..\n"))
  library(bayesCC)
  library(cluster)
  library(parallel)
  # out <- values.query$out
  out <- queried.data
  values.bcc$data <<- bcc$prepare.data(data = out, eps = input$bayesEpsilon, quantile.cut = input$bayesQuantileCut)
  # values.bcc$alphaStarD <<- bcc$calc.opt.k(data = values.bcc$data, .maxiter = input$maxIter1)
  results <- readRDS("bcc_results.rds")
  alpha.s <- lapply(results, bayesCC::alphaStar)
  values.bcc$alphaStarD <<- data.frame(alpha.s)
 })
 bcc.opt.k.event <- eventReactive(values.bcc$alphaStarD, {
  x <- values.bcc$alphaStarD
  return(x)
 })
 output$bayesOptimalK <- renderPlot({
  x <- bcc.opt.k.event()
  beepr::beep(10)
  boxplot(x, main = "Mean-adjusted adherence by K")
 })
 output$bayesText <- renderPrint({
  alpha.sd <- bcc.opt.k.event()
  Ks <- 2:5
  names(Ks) <- paste0("K", Ks)
  opt.Ks <- lapply(alpha.sd, max)
  opt.K <- Ks[which.max(unlist(opt.Ks))]
  paste("Optimal clustering value K:", opt.K)
 })
 observeEvent(input$bayesRunButton, {
  data <- values.bcc$data
  if(is.null(data)) return(cat("No input datasets provided!\n"))
  cat("Running Bayes CC algorithm...\n")
  values.bcc$out <<- bayesCC::bayesCC(X = data[-4],
                                      K = input$numberOfK,
                                      IndivAlpha = TRUE,
                                      maxiter = input$maxIter2)
 })
 bcc.out.event <- eventReactive(values.bcc$out, {
  cluster <- bin2dec(values.bcc$out$Cbest)
  return(cluster)
 })
 output$bayesCnv <- renderPlot({
  cluster <- bcc.out.event()
  cat("Plotting 3D CNV clusters...\n")
  plot.3D.clusters(values.bcc$data[["CNV"]], cluster, .phi = input$bayesPhi, .theta = input$bayesTheta)
 })
 output$bayesRna <- renderPlot({
  cluster <- bcc.out.event()
  cat("Plotting 3D RNA clusters...\n")
  plot.3D.clusters(values.bcc$data[["RNA"]], cluster, .phi = input$bayesPhi, .theta = input$bayesTheta)
 })
 output$bayesProtein <- renderPlot({
  cluster <- bcc.out.event()
  cat("Plotting 3D protein clusters...\n")
  plot.3D.clusters(values.bcc$data[["Protein"]], cluster, .phi = input$bayesPhi, .theta = input$bayesTheta)
 })
 
 # LRAcluster ####
 observeEvent(input$lraOptButton, {
  library(LRAcluster)
  values.lra$data <<- lra$prepare.data(queried.data)
  values.lra$rlist <<- lra$lra.list(values.lra$data, dim = 2)
 })
 lra.opt.k.event <- eventReactive(values.lra$rlist, {
  coor <- t(values.lra$rlist$coordinate)
  avg.sil <- lra$optimize_k(coor)
 })
 output$lraOut1 <- renderPlot({
  avg_sil <- lra.opt.k.event()
  k <- 2:5
  plot(x = k, type = 'b', y = avg_sil, xlab = 'Number of clusters', ylab = 'Average Silhouette Scores', frame = FALSE)
 })
 output$lraOut2 <- renderPlot({
  lra.opt.k.event()
  data <- t(values.lra$rlist$coordinate)
  play.sound()
  factoextra::fviz_nbclust(data, FUNcluster = kmeans, method = "silhouette", k.max = 5)
 })
 run.lra.event <- eventReactive(input$lraRunButton, {
  library(plot3D)
  if(is.null(values.lra$rlist)) return(cat("No coordinate values calculated!\n"))
  x <- t(values.lra$rlist$coordinate)
  centers <- input$lraInput2
  rclust <- kmeans(x, centers)
 })
 output$lraOut3 <- renderPlot({
  # library(plot3D)
  rclust <- run.lra.event()
  lra$plot.clusters(rclust$cluster, values.lra$rlist)
 })
 output$lraOut4 <- renderPlot({
  rclust <- run.lra.event()
  .name <- c(method = "LRAcluster", data = "CNV")
  plot.3D.clusters(values.lra$data[[1]], rclust$cluster, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(values.lra$data[[1]], rclust$cluster, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(values.lra$data[[1]], rclust$cluster, .phi = 0, .theta = 90, name = .name)
 })
 output$lraOut5 <- renderPlot({
  rclust <- run.lra.event()
  .name <- c(method = "LRAcluster", data = "RNA")
  plot.3D.clusters(values.lra$data[[2]], rclust$cluster, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(values.lra$data[[2]], rclust$cluster, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(values.lra$data[[2]], rclust$cluster, .phi = 0, .theta = 90, name = .name)
 })
 output$lraOut6 <- renderPlot({
  rclust <- run.lra.event()
  .name <- c(method = "LRAcluster", data = "Protein")
  plot.3D.clusters(values.lra$data[[3]], rclust$cluster, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(values.lra$data[[3]], rclust$cluster, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(values.lra$data[[3]], rclust$cluster, .phi = 0, .theta = 90, name = .name)
 })

 # iClusterPlus ####
 observeEvent(input$iClusterOptButton, {
  library(iClusterPlus)
  inp <- queried.data
  values.icluster$data <- icluster$prepare.data(data = inp,
                                                eps = input$iClusterInput3,
                                                quantile.cut = input$iClusterInput4)
 })
 icluster.opt.event <- eventReactive(values.icluster$data, {
  # cv.fit <- icluster$tune.model(values.icluster$data)
  cat("Reading from file: icluster_cv.fit.rds...\t")
  cv.fit <- readRDS("icluster_cv.fit.rds")
  cat("DONE.\n")
  out <- icluster$select.model(cv.fit, values.icluster$data)
  # values.icluster$best.fit <<- out$best.fit
  return(out)
 })
 output$iClusterOut1 <- renderPlot({
  out.list <- icluster.opt.event()
  play.sound()
  plot(1:(out.list$nK + 1),
       c(0, out.list$devRatMinBIC),
       type = "b",
       xlab = "Number of clusters (K+1)",
       ylab = "%Explained Variation")
 })
 icluster.run.event <- eventReactive(input$iClusterRunButton, {
  if(is.null(values.icluster$data))
   return(cat("\tNo input dataset selected!\n"))
  out <- icluster$run.iClusterPlus(values.icluster$data, k = input$iClusterInput2)
  play.sound()
  out
 })
 output$iClusterTabCnv <- renderPlot({
  out <- icluster.run.event()
  .name <- c(method = "iClusterPlus", data = "CNV")
  plot.3D.clusters(t(out$data$CNV), out$fit$clusters, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(out$data$CNV), out$fit$clusters, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(out$data$CNV), out$fit$clusters, .phi = 0, .theta = 90, name = .name)
 })
 output$iClusterTabRna <- renderPlot({
  out <- icluster.run.event()
  .name <- c(method = "iClusterPlus", data = "RNA")
  plot.3D.clusters(t(out$data$RNA), out$fit$clusters, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(out$data$RNA), out$fit$clusters, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(out$data$RNA), out$fit$clusters, .phi = 0, .theta = 90, name = .name)
 })
 output$iClusterTabProtein <- renderPlot({
  out <- icluster.run.event()
  .name <- c(method = "iClusterPlus", data = "Protein")
  plot.3D.clusters(t(out$data$Protein), out$fit$clusters, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(out$data$Protein), out$fit$clusters, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(out$data$Protein), out$fit$clusters, .phi = 0, .theta = 90, name = .name)
 })
 
 # PINSPlus ####
 observeEvent(input$pinsplusRun, {
  library(PINSPlus)
  values.pinsplus$data <<- pinsplus$prepare.data(queried.data, eps = input$pinsplusInput1, quantile.cut = input$pinsplusInput2)
 })
 pinsplus.event <- eventReactive(values.pinsplus$data, {
  n <- as.integer(input$pinsplusCpus)
  out <- pinsplus$run(data = values.pinsplus$data, cpus = n)
  play.sound()
  out
 })
 # output$pinsplusPlot <- renderPlot({
 #   out <- pinsplus.event()
 #   pinsplus$plot.survival(res = r)
 # })
 output$pinsplusOut4 <- renderPlot({
  out <- pinsplus.event()
  .name <- c(method = "PINSPlus", data = "CNV")
  plot.3D.clusters(t(values.pinsplus$data$CNV), out$cluster1, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(values.pinsplus$data$CNV), out$cluster1, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(values.pinsplus$data$CNV), out$cluster1, .phi = 0, .theta = 90, name = .name)
 })
 output$pinsplusOut5 <- renderPlot({
  out <- pinsplus.event()
  .name <- c(method = "PINSPlus", data = "RNA")
  plot.3D.clusters(t(values.pinsplus$data$RNA), out$cluster1, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(values.pinsplus$data$RNA), out$cluster1, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(values.pinsplus$data$RNA), out$cluster1, .phi = 0, .theta = 90, name = .name)
 })
 output$pinsplusOut6 <- renderPlot({
  out <- pinsplus.event()
  .name <- c(method = "PINSPlus", data = "Protein")
  plot.3D.clusters(t(values.pinsplus$data$Protein), out$cluster1, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(values.pinsplus$data$Protein), out$cluster1, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(values.pinsplus$data$Protein), out$cluster1, .phi = 0, .theta = 90, name = .name)
 })

 # SNF ####
 observeEvent(input$snfOptButton, {
  library(SNFtool)
  values.snf$data <<- snf$prepare.data(data = queried.data, eps = input$snfInput1, quantile.cut = input$snfInput2)
 })
 observeEvent(values.snf$data, {
  values.snf$opt.k <<- snf$calc.opt.k(D1 = values.snf$data[[1]], D2 = values.snf$data[[2]], D3 = values.snf$data[[3]], K = input$snfInput3, alpha = input$snfInput4, T = input$snfInput5)
  play.sound()
 })
 output$snfBestEigenvalue <- renderValueBox({
  sub <- names(values.snf$opt.k[[2]])[1]
  valueBox(value = values.snf$opt.k$estimationResult[[1]], subtitle = sub)
 })
 output$snf2Eigenvalue <- renderValueBox({
  sub <- names(values.snf$opt.k[[2]])[2]
  valueBox(values.snf$opt.k$estimationResult[[2]], subtitle = sub)
 })
 output$snfRotCost <- renderValueBox({
  sub <- names(values.snf$opt.k[[2]])[3]
  valueBox(values.snf$opt.k$estimationResult[[3]], subtitle = sub)
 })
 output$snf2RotCost <- renderValueBox({
  sub <- names(values.snf$opt.k[[2]])[4]
  valueBox(values.snf$opt.k$estimationResult[[4]], subtitle = sub)
 })
 snf.run.event <- eventReactive(input$snfRunButton, {
  out <- snf$run.SNF(D1 = values.snf$data[[1]], D2 = values.snf$data[[2]], D3 = values.snf$data[[3]], K = input$snfInput7, alpha = input$snfInput8, T = input$snfInput9)
  return(out)
 })
 output$snfOut4 <- renderPlot({
  out <- snf.run.event()
  .name <- c(method = "SNF", data = "CNV")
  plot.3D.clusters(t(values.snf$data$CNV), out, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(values.snf$data$CNV), out, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(values.snf$data$CNV), out, .phi = 0, .theta = 90, name = .name)
 })
 output$snfOut5 <- renderPlot({
  out <- snf.run.event()
  .name <- c(method = "SNF", data = "RNA")
  plot.3D.clusters(t(values.snf$data$RNA), out, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(values.snf$data$RNA), out, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(values.snf$data$RNA), out, .phi = 0, .theta = 90, name = .name)
 })
 output$snfOut6 <- renderPlot({
  out <- snf.run.event()
  .name <- c(method = "SNF", data = "Protein")
  plot.3D.clusters(t(values.snf$data$Protein), out, .phi = 0, .theta = 0, name = .name)
  plot.3D.clusters(t(values.snf$data$Protein), out, .phi = 0, .theta = 45, name = .name)
  plot.3D.clusters(t(values.snf$data$Protein), out, .phi = 0, .theta = 90, name = .name)
 })

 # 3. Analysis ####
 # 3.1. MAF Summary ####
 library(pdftools)
 output$mafSummaryBccCluster1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/bayesCC_g1.pdf", dpi = 100)
  list(src = "bayesCC_g1_1.png", width = "100%")
 })
 output$mafSummaryBccCluster2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/bayesCC_g2.pdf", dpi = 100)
  list(src = "bayesCC_g2_1.png", width = "100%")
 })
 output$mafSummaryicpCluster1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/iclusterplus_g1.pdf", dpi = 100)
  list(src = "iclusterplus_g1_1.png", width = "100%")
 })
 output$mafSummaryicpCluster2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/iclusterplus_g2.pdf", dpi = 100)
  list(src = "iclusterplus_g2_1.png", width = "100%")
 })
 output$mafSummaryicpCluster3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/iclusterplus_g3.pdf", dpi = 100)
  list(src = "iclusterplus_g3_1.png", width = "100%")
 })
 output$mafSummaryLraCluster1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/lra_g1.pdf", dpi = 100)
  list(src = "lra_g1_1.png", width = "100%")
 })
 output$mafSummaryLraCluster2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/lra_g2.pdf", dpi = 100)
  list(src = "lra_g2_1.png", width = "100%")
 })
 output$mafSummaryPinsCluster1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/pinsplus_g1.pdf", dpi = 100)
  list(src = "pinsplus_g1_1.png", width = "100%")
 })
 output$mafSummaryPinsCluster2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/pinsplus_g2.pdf", dpi = 100)
  list(src = "pinsplus_g2_1.png", width = "100%")
 })
 output$mafSummaryPinsCluster3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/pinsplus_g3.pdf", dpi = 100)
  list(src = "pinsplus_g3_1.png", width = "100%")
 })
 output$mafSummaryPinsCluster4 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/pinsplus_g4.pdf", dpi = 100)
  list(src = "pinsplus_g4_1.png", width = "100%")
 })
 output$mafSummarySnfCluster1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/snf_g1.pdf", dpi = 100)
  list(src = "snf_g1_1.png", width = "100%")
 })
 output$mafSummarySnfCluster2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/snf_g2.pdf", dpi = 100)
  list(src = "snf_g2_1.png", width = "100%")
 })
 output$mafSummarySnfCluster3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/MAF Summary/snf_g3.pdf", dpi = 100)
  list(src = "snf_g3_1.png", width = "100%")
 })
 # 3.2. Oncoplots ####
 output$oncoplotBcc1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/bayesCC_g1.pdf", dpi = 100)
  list(src = "bayesCC_g1_1.png", width = "100%")
 })
 output$oncoplotBcc2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/bayesCC_g2.pdf", dpi = 100)
  list(src = "bayesCC_g2_1.png", width = "100%")
 })
 output$oncoploticp1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/iclusterplus_g1.pdf", dpi = 100)
  list(src = "iclusterplus_g1_1.png", width = "100%")
 })
 output$oncoploticp2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/iclusterplus_g2.pdf", dpi = 100)
  list(src = "iclusterplus_g2_1.png", width = "100%")
 })
 output$oncoploticp3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/iclusterplus_g3.pdf", dpi = 100)
  list(src = "iclusterplus_g3_1.png", width = "100%")
 })
 output$oncoplotLra1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/lra_g1.pdf", dpi = 100)
  list(src = "lra_g1_1.png", width = "100%")
 })
 output$oncoplotLra2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/lra_g2.pdf", dpi = 100)
  list(src = "lra_g2_1.png", width = "100%")
 })
 output$oncoplotPins1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/pinsplus_g1.pdf", dpi = 100)
  list(src = "pinsplus_g1_1.png", width = "100%")
 })
 output$oncoplotPins2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/pinsplus_g2.pdf", dpi = 100)
  list(src = "pinsplus_g2_1.png", width = "100%")
 })
 output$oncoplotPins3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/pinsplus_g3.pdf", dpi = 100)
  list(src = "pinsplus_g3_1.png", width = "100%")
 })
 output$oncoplotPins4 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/pinsplus_g4.pdf", dpi = 100)
  list(src = "pinsplus_g4_1.png", width = "100%")
 })
 output$oncoplotSnf1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/snf_g1.pdf", dpi = 100)
  list(src = "snf_g1_1.png", width = "100%")
 })
 output$oncoplotSnf2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/snf_g2.pdf", dpi = 100)
  list(src = "snf_g2_1.png", width = "100%")
 })
 output$oncoplotSnf3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Oncoplot/snf_g3.pdf", dpi = 100)
  list(src = "snf_g3_1.png", width = "100%")
 })
 # 3.3. Rec. CNV ####
 output$recCnvBcc1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/BayesCC/BayesCC_g1.pdf", dpi = 100)
  list(src = "BayesCC_g1_1.png", width = "100%")
 })
 output$recCnvBcc2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/BayesCC/BayesCC_g2.pdf", dpi = 100)
  list(src = "BayesCC_g2_1.png", width = "100%")
 })
 output$recCnvicp1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Iclusterplus/Iclusterplus_g1.pdf", dpi = 100)
  list(src = "Iclusterplus_g1_1.png", width = "100%")
 })
 output$recCnvicp2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Iclusterplus/Iclusterplus_g2.pdf", dpi = 100)
  list(src = "Iclusterplus_g2_1.png", width = "100%")
 })
 output$recCnvicp3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Iclusterplus/Iclusterplus_g3.pdf", dpi = 100)
  list(src = "Iclusterplus_g3_1.png", width = "100%")
 })
 output$recCnvLra1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Lra/Lra_g1.pdf", dpi = 100)
  list(src = "Lra_g1_1.png", width = "100%")
 })
 output$recCnvLra2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Lra/Lra_g2.pdf", dpi = 100)
  list(src = "Lra_g2_1.png", width = "100%")
 })
 output$recCnvPins1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Pinsplus/Pinsplus_g1.pdf", dpi = 100)
  list(src = "Pinsplus_g1_1.png", width = "100%")
 })
 output$recCnvPins2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Pinsplus/Pinsplus_g2.pdf", dpi = 100)
  list(src = "Pinsplus_g2_1.png", width = "100%")
 })
 output$recCnvPins3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Pinsplus/Pinsplus_g3.pdf", dpi = 100)
  list(src = "Pinsplus_g3_1.png", width = "100%")
 })
 output$recCnvPins4 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Pinsplus/Pinsplus_g4.pdf", dpi = 100)
  list(src = "Pinsplus_g4_1.png", width = "100%")
 })
 output$recCnvSnf1 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Snf/Snf_g1.pdf", dpi = 100)
  list(src = "Snf_g1_1.png", width = "100%")
 })
 output$recCnvSnf2 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Snf/Snf_g2.pdf", dpi = 100)
  list(src = "Snf_g2_1.png", width = "100%")
 })
 output$recCnvSnf3 <- renderImage(deleteFile = TRUE, {
  pdf_convert("../../../../Desktop/Outputs/CNV/Recurrent CNV Plots/Snf/Snf_g3.pdf", dpi = 100)
  list(src = "Snf_g3_1.png", width = "100%")
 })
 # 3.4. Permutation Tests ####
 output$pTestsBcc1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/BayesCC/Permutation Test/Group 1/1.pdf", dpi = 75)
   list(src = "1_1.png")
 })
 output$pTestsBcc2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/BayesCC/Permutation Test/Group 2/2.pdf", dpi = 75)
   list(src = "2_1.png")
 })
 output$pTestsicp1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/iClusterPlus/Permutation Test/Group 1/1.pdf", dpi = 75)
   list(src = "1_1.png")
 })
 output$pTestsicp2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/iClusterPlus/Permutation Test/Group 2/2.pdf", dpi = 75)
   list(src = "2_1.png")
 })
 output$pTestsicp3 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/iClusterPlus/Permutation Test/Group 3/3.pdf", dpi = 75)
   list(src = "3_1.png")
 })
 output$pTestsLra1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/LRAcluster/Permutation Test/Group 1/1.pdf", dpi = 75)
   list(src = "1_1.png")
 })
 output$pTestsLra2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/LRAcluster/Permutation Test/Group 2/2.pdf", dpi = 75)
   list(src = "2_1.png")
 })
 output$pTestsPins1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/PINSPlus/Group 1/1.pdf", dpi = 75)
   list(src = "1_1.png")
 })
 output$pTestsPins2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/PINSPlus/Group 2/2.pdf", dpi = 75)
   list(src = "2_1.png")
 })
 output$pTestsPins3 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/PINSPlus/Group 3/3.pdf", dpi = 75)
   list(src = "3_1.png")
 })
 output$pTestsPins4 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/PINSPlus/Group 4/4.pdf", dpi = 75)
   list(src = "4_1.png")
 })
 output$pTestsSnf1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/SNF/Group 1/1.pdf", dpi = 75)
   list(src = "1_1.png")
 })
 output$pTestsSnf2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/SNF/Group 2/2.pdf", dpi = 75)
   list(src = "2_1.png")
 })
 output$pTestsSnf3 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../CNV Analysis/SNF/Group 3/3.pdf", dpi = 75)
   list(src = "3_1.png")
 })
 # 3.5. Enrichment Analysis ####
 output$enrAnalysisBcc1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_bayesCC_g1.pdf", dpi = 100)
   list(src = "EA_bayesCC_g1_1.png", width = "100%")
 })
 output$enrAnalysisBcc2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_bayesCC_g2.pdf", dpi = 100)
   list(src = "EA_bayesCC_g2_1.png", width = "100%")
 })
 output$enrAnalysisicp1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_iclusterplus_g1.pdf", dpi = 100)
   list(src = "EA_iclusterplus_g1_1.png", width = "100%")
 })
 output$enrAnalysisicp2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_iclusterplus_g2.pdf", dpi = 100)
   list(src = "EA_iclusterplus_g2_1.png", width = "100%")
 })
 output$enrAnalysisicp3 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_iclusterplus_g3.pdf", dpi = 100)
   list(src = "EA_iclusterplus_g3_1.png", width = "100%")
 })
 output$enrAnalysisLra1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_lra_g1.pdf", dpi = 100)
   list(src = "EA_lra_g1_1.png", width = "100%")
 })
 output$enrAnalysisLra2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_lra_g2.pdf", dpi = 100)
   list(src = "EA_lra_g2_1.png", width = "100%")
 })
 output$enrAnalysisPins1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_pinsplus_g1.pdf", dpi = 100)
   list(src = "EA_pinsplus_g1_1.png", width = "100%")
 })
 output$enrAnalysisPins2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_pinsplus_g2.pdf", dpi = 100)
   list(src = "EA_pinsplus_g2_1.png", width = "100%")
 })
 output$enrAnalysisPins3 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_pinsplus_g3.pdf", dpi = 100)
   list(src = "EA_pinsplus_g3_1.png", width = "100%")
 })
 output$enrAnalysisPins4 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_pinsplus_g4.pdf", dpi = 100)
   list(src = "EA_pinsplus_g4_1.png", width = "100%")
 })
 output$enrAnalysisSnf1 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_snf_g1.pdf", dpi = 100)
   list(src = "EA_snf_g1_1.png", width = "100%")
 })
 output$enrAnalysisSnf2 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_snf_g2.pdf", dpi = 100)
   list(src = "EA_snf_g2_1.png", width = "100%")
 })
 output$enrAnalysisSnf3 <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Enrichment Analysis/EA_snf_g3.pdf", dpi = 100)
   list(src = "EA_snf_g3_1.png", width = "100%")
 })
 # 3.7. Gene Expression Heatmaps ####
 output$geHeatmapsBcc <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Heatmaps/bayesCC_2groups.pdf", dpi = 45)
   list(src = "bayesCC_2groups_1.png", width = "99%")
 })
 output$geHeatmapsicp <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Heatmaps/iclusterplus_3groups.pdf", dpi = 45)
   list(src = "iclusterplus_3groups_1.png", width = "99%")
 })
 output$geHeatmapsLra <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Heatmaps/lra_2groups.pdf", dpi = 45)
   list(src = "lra_2groups_1.png", width = "99%")
 })
 output$geHeatmapsPins <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Heatmaps/pinsplus_4groups.pdf", dpi = 45)
   list(src = "pinsplus_4groups_1.png", width = "99%")
 })
 output$geHeatmapsSnf <- renderImage(deleteFile = TRUE, {
   pdf_convert("../../../../Desktop/Outputs/RNA/Heatmaps/snf_3groups.pdf", dpi = 45)
   list(src = "snf_3groups_1.png", width = "99%")
 })
 # 3.8. Survival Analysis ####
 output$survAnalysisBcc <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/RNA/Survival/bayesCC.png", width = "100%")
 })
 output$survAnalysisicp <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/RNA/Survival/iclusterplus.png", width = "100%")
 })
 output$survAnalysisLra <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/RNA/Survival/lra.png", width = "100%")
 })
 output$survAnalysisPins <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/RNA/Survival/pinsplus.png", width = "100%")
 })
 output$survAnalysisSnf <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/RNA/Survival/snf.png", width = "100%")
 })
 # 3.9. Clinical Statistics ####
  output$clinicStatsBcc1 <- renderImage(deleteFile = FALSE, {
    list(src = "../../../../Desktop/Outputs/Statistical/Group 1 - Path. Stage.png", width = "100%")
  })
  output$clinicStatsBcc2 <- renderImage(deleteFile = FALSE, {
    list(src = "../../../../Desktop/Outputs/Statistical/Group 1 - Prior treatment.png", width = "100%")
  })
 output$clinicStatsBcc3 <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/Statistical/Group 1 - Site of Resection.png", width = "100%")
 })
 output$clinicStatsBcc4 <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/Statistical/Group 1 - Vital Status.png", width = "100%")
 })
 output$clinicStatsBcc5 <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/Statistical/Group 2 - Pathologic Stage.png", width = "100%")
 })
  output$clinicStatsBcc6 <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/Statistical/Group 2 - Prior Treatment.png", width = "100%")
 })
 output$clinicStatsBcc7 <- renderImage(deleteFile = FALSE, {
   list(src = "../../../../Desktop/Outputs/Statistical/Group 2 - Site of Resection.png", width = "100%")
 })
}