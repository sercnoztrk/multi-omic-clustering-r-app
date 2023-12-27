# Init ####
if (Sys.info()["sysname"] == "Linux") {
  setwd("/media/sercan/908CBCDF8CBCC0D0/Users/Sercan/Documents/YÃ¼ksek Lisans/Tez/Shiny UI/")
  load(".RData")
  library(shiny, quietly = TRUE)
  library(shinydashboard, quietly = TRUE)
  library(TCGAbiolinks, quietly = TRUE)
} else if (Sys.info()["sysname"] == "Windows") {
  setwd("Yüksek Lisans/Tez/Shiny UI/")
  load(".RData")
  library(shiny, quietly = TRUE)
  library(shinydashboard, quietly = TRUE)
  library(TCGAbiolinks, quietly = TRUE)
}

source("ui.R")
source("server.R")

if(!exists("bcc")) {
  source("../BayesCC/BayesCC.R")
  bcc <- list(prepare.data = prepare.data, calc.opt.k = calc.opt.k)
}
if(!exists("data.analysis")) {
  source("../Data Analysis/Data_Analysis.R")
  data.analysis <- list(barcode.filtering = barcode.filtering,
                        common.patients = common.patients,
                        dl.data = dl.data,
                        get.patients.data = get.patients.data,
                        load.data = load.data,
                        query.data = query.data,
                        run = run,
                        tcga_replicateFilter = tcga_replicateFilter)
}
if(!exists("lra")) {
  source("../LRAcluster/LRAcluster.R")
  lra <- list(lra.list = lra.list,
              optimize_k = optimize_k,
              plot.clusters = plot.clusters,
              prepare.data = prepare.data)
  remove(lra.list, optimize_k, plot.clusters, prepare.data)
}
if(!exists("icluster")) {
  source("../iCluster+/run.iCluster+.R")
  icluster <- list(generate.heatmap = generate.heatmap,
                   prepare.data = prepare.data,
                   run.iClusterPlus = run.iClusterPlus,
                   select.model = select.model,
                   tune.model = tune.model)
  remove("generate.heatmap", "prepare.data", "run.iClusterPlus", "select.model", "tune.model")
}
if(!exists("snf")) {
  snf <- list(prepare.data = readRDS("prepare.data.rds"),
              run.SNF = readRDS("run.SNF.rds"))
  file.remove("prepare.data.rds", "run.SNF.rds")
}
if(!exists("pinsplus")) {
  pinsplus <- list(prepare.data = readRDS("prepare.data.rds"), 
                   run = readRDS("run.rds"),
                   survival = readRDS("../PINSPlus/survival.rds"))
  file.remove("prepare.data.rds", "run.rds")
}

queried.data <- data.analysis$run(param$cnv, param$rna, param$protein)

plot.3D.clusters <- function(data, cluster, .phi = 40, .theta = 40, name) {
  library(plot3D, quietly = TRUE)
  pca <- prcomp(x = data, rank. = 3, scale. = TRUE)
  x <- pca$rotation[, 1]
  y <- pca$rotation[, 2]
  z <- pca$rotation[, 3]
  # data(iris)
  # x <- sep.l <- iris$Sepal.Length
  # y <- pet.l <- iris$Petal.Length
  # z <- sep.w <- iris$Sepal.Width
  len <- length(unique(cluster))
  .color <- c("dark red", "dark blue", "dark green", "dark orange")[1:len]
  .labels <- c("1", "2", "3", "4")[1:len]
  if(len == 2) { .at <- c(1.25, 1.75) }
  else if(len == 3) { .at <- c(1.3333333, 2, 2.6666667) }
  else if(len == 4) { .at <- c(1.375, 2.125, 2.875, 3.625) }
  panelfirst <- function(pmat) {
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = cluster, pch = ".", cex = 2, add = TRUE, colkey = FALSE)
    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = cluster, pch = ".", cex = 2, add = TRUE, colkey = FALSE)
  }
  file <- paste0("/home/sercan/Desktop/", name[["method"]], "_", name[["data"]], "_", "Theta", .theta, ".pdf")
  cat(file, "\n")
  # pdf(file)
  plot3D::scatter3D(x = x,
                    y = y,
                    z = z,
                    colvar = cluster,
                    phi = .phi,
                    theta = .theta,
                    labels = cluster,
                    col = .color,
                    pch = 19,
                    bty = "g",
                    # type = "h",
                    clab = "Cluster",
                    colkey = list(addlines = TRUE, side = 4, labels = .labels, at = .at),
                                  panel.first = panelfirst)
  # dev.off()
  # plot3D::text3D(x = x,
  #                y = y,
  #                z = z,
  #                labels = cluster,
  #                add = TRUE)
}

# Return ####
runApp(shinyApp(ui, server), launch.browser = FALSE)