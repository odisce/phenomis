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
#' @import biodbChebi MultiAssayExperiment SummarizedExperiment
#' @importClassesFrom MultiDataSet MultiDataSet
#' @importFrom grDevices boxplot.stats rainbow
#' @importFrom graphics title
#' @importFrom methods as
#' @importFrom stats IQR cor median predict quantile rnorm runif sd
#' @importFrom utils head read.table tail
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

#' Searching patterns in code
#'
#' This function is useful for code inspection, development, and debugging
#' 
#' @param pattern.c Character: pattern to be searched in the files
#' @param dir.c Character: parent directory of the files
#' @param fixed.l Logical: should the pattern be searched 'as is'
#' @param filename_pattern.c Character: optional filter of file names to be searched
#' @return nothing is returned
#' @export
#' @examples
#' \dontrun{
#' search_code("aes", filename_pattern.c = ".R")
#' }
search_code <- function(pattern.c,
                        dir.c = getwd(),
                        fixed.l = TRUE,
                        filename_pattern.c = NULL) {
  
  warn_option.n <- getOption("warn")
  options(warn = -1)
  
  if (is.na(file.info(dir.c)[["isdir"]]))
    stop(dir.c, " is not a readable directory")
  
  file.vc <- list.files(dir.c,
                        full.names = TRUE,
                        pattern = filename_pattern.c,
                        recursive = TRUE)
  
  found.l <- FALSE
  
  for (file.c in file.vc) {
    
    lines.vc <- readLines(file.c)
    
    pattern.vi <- grep(pattern.c,
                       x = lines.vc,
                       fixed = fixed.l)
    
    if (length(pattern.vi) > 0) {
      
      message("\nPattern found in '",
              gsub(dir.c, "", file.c, fixed = TRUE),
              "' in the following lines:")
 
        for (pattern.i in pattern.vi)
          message("line ", pattern.i, ": ", lines.vc[pattern.i])
      
      found.l <- TRUE
      
    }
  }
  
  if (!found.l) {
    message("The pattern was not found in any of the searched files:")
    message(paste(file.vc, collapse = "\n"))
  }
  
  options(warn = warn_option.n)
  
  invisible(NA)
  
}


.identi_check <- function(x, y) {
  
  colors.vc <- c("TRUE" = 32, "FALSE" = 31)
  
  .check <- function(test.c) {
    test_result.c <- as.character(eval(parse(text = test.c)))
    cat(paste0("\033[0;",
               colors.vc[test_result.c],
               "m",
               test.c, ": ", test_result.c,
               "\033[0m","\n"))
  }
  
  .check("identical(x, y)")
  
  .check("identical(unname(x), unname(y))")
  
  .check("identical(length(x), length(y))")
  
  .check("identical(class(x), class(y))")
  
  .check("identical(mode(x), mode(y))")
  
  .check("!any(duplicated(x))")
  
  .check("!any(duplicated(y))")
  
  .check("identical(sort(x), sort(y))")
  
  bind.m <- cbind(sort(x), sort(y))
  
  ident.vl <- apply(bind.m, 1,
                    function(bind.v) {
                      identical(unname(bind.v[1]), unname(bind.v[2]))
                    })
  
  diff.vi <- union(which(is.na(ident.vl)), which(!ident.vl))
  
  if(length(diff.vi))
    print(bind.m[diff.vi, , drop = FALSE])
  
}