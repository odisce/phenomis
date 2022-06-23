#' Postprocessing and univariate analysis of omics data
#'
#' The 'phenomis' package provides methods to perform post-processing (i.e. quality control
#' and normalization) as well as univariate statistical analysis of single and multi-omics data sets.
#' These methods include quality control metrics, signal drift and batch effect correction,
#' intensity transformation, univariate hypothesis testing, but also clustering (as well as
#' annotation of metabolomics data). The data are handled in the standard Bioconductor formats
#' (i.e. SummarizedExperiment and MultiAssayExperiment for single and multi-omics datasets,
#' respectively; the alternative ExpressionSet and MultiDataSet formats are also supported for
#' convenience). As a result, all methods can be readily chained as workflows. The pipeline can
#' be further enriched by multivariate analysis and feature selection, by using the 'ropls' and
#' 'biosigner' packages, which support the same formats. Data can be conveniently imported from
#' and exported to text files. Although the methods were initially targeted to metabolomics data,
#' most of the methods can be applied to other types of omics data (e.g., transcriptomics, proteomics).
#'
#' @import Biobase biodb biodbChebi ggplot2 limma MultiAssayExperiment SummarizedExperiment
#' @importFrom data.table fread
#' @importFrom futile.logger flog.threshold
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices boxplot.stats colorRampPalette dev.off pdf rainbow
#' @importFrom graphics abline arrows axis box image layout mtext par plot rect text title
#' @importFrom grid grid.draw grid.newpage
#' @importFrom htmlwidgets saveWidget
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom methods as is validObject
#' @importFrom MultiDataSet add_eset as.list createMultiDataSet
#' @importFrom plotly ggplotly
#' @importFrom PMCMRplus kwAllPairsNemenyiTest
#' @importFrom ranger ranger
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ropls getOpls getPcaVarVn getScoreMN opls view
#' @importFrom stats aggregate anova aov as.dendrogram as.dist cor cor.test cutree dist hclust IQR kruskal.test lm loess loess.control median model.matrix na.omit p.adjust pf pnorm predict qf quantile rnorm runif sd t.test TukeyHSD var wilcox.test
#' @importFrom tibble as_tibble is_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom utils head read.table tail write.table
#' @importFrom VennDiagram venn.diagram
#' @name phenomis-package
#' @aliases phenomis phenomis-package
#' @docType package
#' @author E. A. Thévenot (CEA)
#'
#' Maintainer: Etienne Thevenot <etienne.thevenot@@cea.fr>
#' @keywords package
#' @examples
#' # See the package vignette
#'
NULL