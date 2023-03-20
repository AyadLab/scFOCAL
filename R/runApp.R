#' Launch ISOSCELES
#' @examples
#' runISOSCELES()
#' @import shiny
#' @import Seurat
#' @import ggplot2
#' @import tibble
#' @import cowplot
#' @import viridis
#' @import dplyr
#' @import ggsci
#' @import ggrepel
#' @import tidyverse
#' @import plotly
#' @import htmlwidgets
#' @import reshape2
#' @import Hmisc
#' @import corrplot
#' @import pheatmap
#' @import grid
#' @import MAST
#' @import shinydashboard
#' @import shinythemes
#' @import scales
#' @import ggforce
#' @import EnhancedVolcano
#' @import DT
#' @import shiny
#' @export
runISOSCELES <- function() {
  shiny::runApp(appDir = system.file('shiny', package = "ISOSCELES"))
  }

