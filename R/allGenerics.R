#### annotating ####

#' MS annotation
#'
#' Annotation with chemical and biological databases by using the 'biodb' 
#' package suite. The present implementation currently enables to query the 
#' ChEBI database or a local database.
#' 
#' @param x An S4 object of class \code{SummarizedExperiment} or
#' \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param database.c character(1): database to be used for annotation; either
#' the ChEBI distant database ('chebi'), or a local database ('local.ms')
#' @param param.ls list: parameters for database query; the database can be
#' queried by either the mass to charge ratio (mz) or the chebi ID; other query
#' parameters include the ionization mode (ms.mode), the mz tolerance (mz.tol;
#' e.g. 5 ppm for Orbitrap Mass Spectrometers), the fields to retrieve (fields),
#' the maximum number of items to retrieve when a field contains more than one 
#' value (fieldsLimit), the maximum number of results to provide for each 
#' query (max.results), prefix of the new columns providing the queried
#' information in the feature metadata (prefix), separator in case of multiple
#' retrieved values (sep), local data base to be queried (local.ms.db); 
#' additional information is provided by the vignettes from the biodb and 
#' biodbChebi packages on Bioconductor
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) 
#' including the appended rowData data frame(s)
#' @rdname annotating
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' # see the (default) parameters (e.g. for ChEBI query)
#' annotating_parameters("chebi")
#' # mz annotation with ChEBI
#' \donttest{
#' sacurine.se <- annotating(sacurine.se, database.c = "chebi",
#' param.ls = list(query.type = "mz", query.col = "mass_to_charge",
#' ms.mode = "neg", prefix = "chebiMZ."))
#' }
#' # mz annotation with local database
#' msdbDF <- read.table(system.file("extdata/local_ms_db.tsv",
#' package = "phenomis"),
#' header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' sacurine.se <- annotating(sacurine.se, database.c = "local.ms",
#' param.ls = list(query.type = "mz", query.col = "mass_to_charge",
#' ms.mode = "neg",
#' mz.tol = 5, mz.tol.unit = "ppm", local.ms.db = msdbDF, prefix = "localMS."))
#' rowData(sacurine.se)[!is.na(rowData(sacurine.se)[, "localMS.accession"]), ]
#' # annotation from ChEBI identifiers
#' \donttest{
#' sacurine.se <- annotating(sacurine.se, database.c = "chebi",
#' param.ls = list(query.type = "chebi.id", query.col = "database_identifier",
#' prefix = "chebiID."))
#' head(rowData(sacurine.se))
#' }
setGeneric("annotating",
           function(x,
                    database.c = c("chebi", "local.ms")[1],
                    param.ls = list(query.type = c("mz",
                                                   "chebi.id")[1],
                                    query.col = "mz",
                                    ms.mode = "pos",
                                    mz.tol = 10,
                                    mz.tol.unit = "ppm",
                                    fields = c("chebi.id",
                                               "name",
                                               "formula",
                                               "molecular.mass",
                                               "monoisotopic.mass"),
                                    fieldsLimit = 1,
                                    max.results = 3,
                                    local.ms.db = data.frame(),
                                    prefix = paste0(database.c, "."),
                                    sep = "|"),
                    report.c = c("none", "interactive", "myfile.txt")[2])
           standardGeneric("annotating"))


#### clustering ####

#' clustering
#'
#' Hierarchical clustering of both samples and variables
#'
#' @param x An S4 object of class \code{SummarizedExperiment}
#' or \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param dissym.c character(1): dissymilarity to be used in the hierarchical
#' clustering (as provided by the hclust package)
#' @param correl.c character(1): correlation coefficient (in case
#' '1-cor' or '1-abs(cor)' are selected as dissymilarity)
#' @param agglo.c character(1): agglomeration method
#' @param clusters.vi integer(2): number of sample and variable
#' clusters, respectively; the default values (2) are only provided as starting
#' guess (e.g. in case of two groups of samples)
#' @param cex.vn numeric(2) [Plot parameter]; size of the sample and
#' variable labels
#' @param palette.c character(1) [Plot parameter]: color palette
#' @param scale_plot.l logical(1) [Plot parameter]: scaling (mean-centering and
#' unit variance scaling) to enhance contrast (for plotting only)
#' @param title.c character(1) [Plot parameter]: Graphic the subtitle
#' @param figure.c character(1): File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively;
#' if 'none', no figure will be generated
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including columns 
#' indicating the clusters in rowData and colData if clusters.vi' has been 
#' specified
#' @rdname clustering
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.se <- correcting(sacurine.se)
#' sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
#' sacurine.se <- transforming(sacurine.se)
#' sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
#' sacurine.se <- clustering(sacurine.se)
#' utils::head(rowData(sacurine.se))
#' 
#' # MultiAssayExperiment
#' 
#' prometis.mae <- reading(system.file("extdata/prometis", package="phenomis"))
#' prometis.mae <- clustering(prometis.mae)
setGeneric("clustering",
           function(x,
                    dissym.c = c("euclidean",
                                 "maximum",
                                 "manhattan",
                                 "canberra",
                                 "binary",
                                 "minkowski",
                                 "1-cor",
                                 "1-abs(cor)")[7],
                    correl.c = c("pearson",
                                 "kendall",
                                 "spearman")[1],
                    agglo.c = c("ward.D",
                                "ward.D2",
                                "single",
                                "complete",
                                "average",
                                "mcquitty",
                                "median",
                                "centroid")[2],
                    clusters.vi = c(2, 2),
                    cex.vn = c(1, 1),
                    palette.c = c("blueOrangeRed",
                                  "redBlackGreen")[1],
                    scale_plot.l = TRUE,
                    title.c = NA,
                    figure.c = c("none", "interactive", "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
           standardGeneric("clustering"))


#### correcting ####

#' correcting
#'
#' Signal drift and batch effect correction. The normalization strategy relies 
#' on the measurements of a pooled (or QC) sample injected periodically: 
#' for each variable, a regression model is fitted to the values of the pool and
#' subsequently used to adjust the intensities of the samples of interest (van
#' der Kloet et al, 2009; Dunn et al, 2011). In case the number of pool
#' observations is below 5, the linear method is used (for all variables) and a
#' warning is generated. In case no pool is available, the samples themselves
#' can be used to computed the regression model (Thevenot et al., 2015). The
#' sample metadata of each datasets (e.g. colData Data Frames) must contain
#' 3 columns: 1) 'sampleType' (character): either 'sample', 'blank', or 'pool',
#' 2) 'injectionOrder' (integer): order of injection in the instrument, and
#' 3) 'batch' (character): batch name.
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or
#' \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param method.vc character of length 1 or the total number of datasets:
#' method(s) to be used for each dataset (either 'serrf' or 'loess');
#' for the 'serrf' approach, the seed is internally set to 123 for
#' reproducibility; in case the parameter is of length 1 and x contains
#' multiple datasets, the same method will be used for all datasets
#' @param reference.vc character of length 1 or the total number of datasets:
#' sample type to be used as reference for the correction (as indicated in the
#' 'sampleType' column from the colData(x); e.g. 'pool' [default]); should be
#' set to 'pool' for the 'serrf' method; in case the parameter is of length 1
#' and x contains multiple datasets, the same reference sample type will be
#' used for all datasets
#' @param loess_span.vn character of length 1 or the total number of datasets:
#' smoothing parameter for the loess regression; between 0 and 1; (default set
#' to 1); in case the parameter is of length 1 and x contains multiple datasets,
#' the same span value will be used for all datasets
#' @param serrf_corvar.vi character of length 1 or the total number of datasets:
#' number of correlated features for the random forest regression;
#' (default set to 10); in case the parameter is of length 1 and x contains
#' multiple datasets, the same value will be used for all datasets
#' @param sample_intensity.c character(1): metric to be used when displaying 
#' the sample intensities
#' @param title.c character(1): Graphic title: if NA [default] the
#' 'title' slot from the experimentData will be used (metadata)
#' @param figure.c character(1): File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively;
#' if 'none', no figure will be generated
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including the corrected 
#' intensities in the assay matrix (matrices)
#' @rdname correcting
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.se <- correcting(sacurine.se)
#' 
#' # MultiDataSet (to be done)
setGeneric("correcting",
           function(x,
                    method.vc = c("loess", "serrf")[1],
                    reference.vc = c("pool", "sample")[1],
                    loess_span.vn = 1,
                    serrf_corvar.vi = 10,
                    sample_intensity.c = c("median", "mean", "sum")[2],
                    title.c = NA,
                    figure.c = c("none", "interactive", "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("correcting"))

#### filtering ####

#' Filtering of the features (and/or samples) with a high proportion of NAs or
#' a low variance
#'
#' Filtering of the features (and/or samples) with a high proportion of NAs or
#' a low variance
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or
#' \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param class.c character(1): name of the column of the sample metadata giving
#' the classification groups: the filtering will be applied on each class
#' (default: "" meaning that there are no specific classes to consider)
#' @param max_na_prop.n numeric(1): maximum proportion of NAs for a feature (or
#' sample) to be kept (e.g. the default 20% maximum proportion of missing 
#' values); in case 'class.c' is provided, the maximum proportion of NAs for a 
#' feature must be achieved in at least one sample class)
#' @param min_variance.n numeric(1): minimum variance for a feature (or sample)
#' to be kept (e.g. the default 0 value to discard constant features
#' (or samples); in case 'class.c' is provided, the minimum variance for a 
#' feature must be achieved in all sample classes
#' @param dims.vc Vector of one or two characters: dimension(s) to which the
#' filtering should be applied; either 'features', 'samples', c('features',
#' 'samples'), or c('samples', 'features'); in the two latter cases, the 
#' dimensions indicated in the dims.vc are filtered sequentially
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including the filtered 
#' data and metadata
#' @rdname filtering
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' assay.mn <- assay(sacurine.se)
#' ropls::view(assay.mn)
#' filtering(sacurine.se)
#' assay.mn[assay.mn < 1e5] <- NA
#' ropls::view(assay.mn)
#' assay(sacurine.se) <- assay.mn
#' filtering(sacurine.se)
#' filtering(sacurine.se, class.c = "gender")
#' filtering(sacurine.se, class.c = "sampleType")
#' 
#' # MultiAssayExperiment
#' 
#' prometis.mae <- reading(system.file("extdata/prometis", package="phenomis"))
#' filtering(prometis.mae)
#' for (set.c in names(prometis.mae)) {
#' set.se <- prometis.mae[[set.c]]
#' assay.mn <- assay(set.se)
#' assay.mn[assay.mn < quantile(c(assay.mn), 0.2)] <- NA
#' assay(set.se) <- assay.mn
#' prometis.mae[[set.c]] <- set.se
#' }
#' filtering(prometis.mae)
#' 
#' # MultiDataSet
#' 
#' prometis.mset <- reading(system.file("extdata/prometis", package="phenomis"),
#'                          output.c = "set")
#' filtering(prometis.mset)
#' for (set.c in names(prometis.mset)) {
#' eset <- prometis.mset[[set.c]]
#' exprs.mn <- Biobase::exprs(eset)
#' exprs.mn[exprs.mn < quantile(c(exprs.mn), 0.2)] <- NA
#' Biobase::exprs(eset) <- exprs.mn
#' prometis.mset <- MultiDataSet::add_eset(prometis.mset, eset,
#'                                         dataset.type = set.c,
#'                                         GRanges = NA, overwrite = TRUE,
#'                                         warnings = FALSE)
#' }
#' filtering(prometis.mset)
setGeneric("filtering",
           function(x,
                    class.c = "",
                    max_na_prop.n = 0.2,
                    min_variance.n = .Machine$double.eps,
                    dims.vc = c("features", "samples"),
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("filtering"))


#### hypotesting ####

#' Univariate hypothesis testing
#'
#' The hypotesting method is a wrapper of the main R functions for
#' hypothesis testing and corrections for multiple testing. The list of
#' available tests includes two sample tests (t-test and Wilcoxon rank test,
#' but also the limma test), analysis of variance (for one and two factors)
#' and Kruskal-Wallis rank test, and correlation tests (by using either the
#' pearson or the spearman correlation). 
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or
#' \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param test.c character(1): One of the 9 available hypothesis tests can be
#' selected (either 'ttest', 'limma', 'wilcoxon', 'anova', 'kruskal', 'pearson',
#' 'spearman', limma2ways', 'limma2waysInter', 'anova2ways', 'anova2waysInter')
#' @param factor_names.vc (Vector of) character(s): Factor(s) of interest
#' (up to two), i.e. name(s) of a column from the pData(x)
#' @param factor_levels.ls List: for each factor of interest (up to two),
#' the levels of the factor can be specified (i.e. re-ordered) by including a
#' character vector with those levels in the list; by default
#' (no specification), the two vectors are set to "default".
#' @param adjust.c character(1): Name of the method for correction of multiple
#' testing (the p.adjust function is used)
#' @param adjust_thresh.n numeric(1): Threshold for (corrected) p-values
#' @param signif_maxprint.i integer(1): Maximum number of significant feature to
#' display on the screen (by default, 'NA', all significant features are
#' displayed)
#' @param title.c character(1): Title of the graphics
#' @param display_signif.l logical(1): In case of two sample tests (or
#' correlation test), should individual boxplots (or scatterplots) of 
#' significant features be shown?
#' @param prefix.c character(1): prefix to be added to the supplementary columns 
#' from the variableMetadata to prevent overwriting of pre-existing columns with
#' identical names [default: ""]
#' @param figure.c character(1): File name with '.pdf' extension for the figure
#' (for venn diagrams, e.g. in the 'anova2ways' test, the extension will be
#' internally changed to '.tiff' for compatibility with the 
#' VennDiagram package); if interactive' (default), figures will be displayed 
#' interactively; if 'none', no figure will be generated
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including the difference 
#' in means/medians or correlations and the adjusted p-values in feature 
#' metadata
#' @rdname hypotesting
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.se <- correcting(sacurine.se, figure.c = 'none')
#' sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
#' sacurine.se <- transforming(sacurine.se)
#' sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
#' # Student's T test
#' sacurine.se <- hypotesting(sacurine.se, "ttest", "gender")
#' # Pearson correlation test
#' sacurine.se <- hypotesting(sacurine.se, "pearson", "age")
#' # ANOVA
#' colData(sacurine.se)[, "ageGroup"] <- vapply(colData(sacurine.se)[, "age"],
#'                                                  function(x) {
#'                                                    if (x < 35) {
#'                                                      return("thirty")
#'                                                    } else if (x < 50) {
#'                                                      return("fourty")
#'                                                    } else {
#'                                                      return("fifty")}},
#'                                                  FUN.VALUE = character(1))
#' sacurine.se <- hypotesting(sacurine.se, "anova", "ageGroup")
#' 
#' # MultiAssayExperiment
#' 
#' prometis.mae <- reading(system.file("extdata/prometis", package="phenomis"))
#' prometis.mae <- hypotesting(prometis.mae, "limma", "gene")
#' 
#' # MultiDataSet
#' 
#' prometis.mset <- reading(system.file("extdata/prometis", package="phenomis"),
#'                          output.c = "set")
#' prometis.mset <- hypotesting(prometis.mset, "limma", "gene")
setGeneric("hypotesting",
           function(x,
                    test.c = c("ttest", "limma", "wilcoxon",
                               "anova", "kruskal",
                               "pearson", "spearman",
                               "limma2ways", "limma2waysInter",
                               "anova2ways", "anova2waysInter")[2],
                    factor_names.vc,
                    factor_levels.ls = list(factor1Vc = "default",
                                            factor2Vc = "default"),
                    adjust.c = c("holm", "hochberg", "hommel", "bonferroni",
                                 "BH", "BY", "fdr", "none")[5],
                    adjust_thresh.n = 0.05,
                    signif_maxprint.i = NA,
                    title.c = NA,
                    display_signif.l = FALSE,
                    prefix.c = "",
                    figure.c = c("none",
                                 "interactive",
                                 "interactive_plotly",
                                 "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
           standardGeneric("hypotesting"))


#### inspecting ####

#' Inspecting
#'
#' Provides numerical metrics and graphical overview of SummarizedExperiment,
#' MultiAssayExperiment, ExpressionSet, or MultiDataSet instance
#' Please note that all variables with a proportion of 
#' missing values > 'max_na_prop.n'
#' or a variance of 0 will be filtered out at the beginning of the method and
#' therefore in the output object
#'
#' @param x An S4 object of class \code{SummarizedExperiment}
#' or \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param pool_as_pool1.l logical(1): should pool be included (as pool1)
#' in the correlation with the dilution factor?
#' @param pool_cv.n numeric(1): threshold for the coefficient of variation 
#' of the pools; the default value (30\%) is often used in metabolomics
#' @param loess_span.n numeric(1): span parameter used in the loess trend 
#' estimation; the default value is set to 1 to prevent overfitting
#' @param sample_intensity.c Character: function to be used to display the 
#' global sample intensity; default: 'mean'
#' @param title.c character(1): MultiAssayExperiment: title of the barplot 
#' showing the number of samples and variables in each dataset; ExpressionSet: 
#' title of the multipanel graphic displaying the metrics (if NA -default- the 
#' title slot from the experimentData will be used)
#' @param plot_dims.l (MultiAssayExperiment) logical(1): should an overview of
#' the number of samples and variables in all datasets be barplotted?
#' @param figure.c character(1): File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; 
#' if 'none', no figure will be generated
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including the computed 
#' sample and variable metrics in the rowData and colData metadata.
#' @rdname inspecting
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.se <- inspecting(sacurine.se)
#' sacurine.se <- correcting(sacurine.se)
#' sacurine.se <- inspecting(sacurine.se)
#' sacurine.se <- transforming(sacurine.se)
#' sacurine.se <- inspecting(sacurine.se)
#' 
#' # MultiAssayExperiment
#' prometis.mae <- reading(system.file("extdata/prometis", 
#'                                     package = "phenomis"))
#' prometis.mae <- inspecting(prometis.mae)
setGeneric("inspecting",
           function(x,
                    pool_as_pool1.l = FALSE,
                    pool_cv.n = 0.3,
                    loess_span.n = 1,
                    sample_intensity.c = c("median", "mean", "sum")[2],
                    title.c = NA,
                    plot_dims.l = TRUE,
                    figure.c = c("none", "interactive", "myfile.pdf")[2],
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("inspecting"))


#### normalizing ####

#' Normalization of the data matrix intensities
#'
#' The matrix intensities may be normalized by using the Probabilistic Quotient
#' Normalization to scale the spectra to the same virtual overall concentration
#'
#' @param x An S4 object of class \code{SummarizedExperiment} 
#' or \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param method.vc character of length 1 or the total number of datasets: 
#' method(s) to be used for each dataset (default is 'pqn'); in case the 
#' parameter is of length 1 and x contains multiple datasets, the same method 
#' will be used for all datasets
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment} 
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including 
#' the (list of) matrix with normalized intensities
#' @rdname normalizing
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.se <- sacurine.se[, colnames(sacurine.se) != 'HU_neg_096_b2']
#' sacurine.se <- transforming(sacurine.se, method.vc = "log10")
#' norm.se <- normalizing(sacurine.se, method.vc = "pqn")
#' 
#' # MultiDataSet

setGeneric("normalizing",
           function(x,
                    method.vc = "pqn",
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("normalizing"))


#### reducing ####

#' Grouping chemically redundant MS1 features
#'
#' This method groups chemically redundant features from a peak table, based on
#' 1) correlation of sample profiles, 2) retention time window, 3) referenced
#' m/z differences. The initial algorithm is named 'Analytic Correlation 
#' Filtration' (Monnerie et al., 2019; DOI:10.3390/metabo9110250) and is 
#' available in Perl and on the Workflow4Metabolomics platform.
#' Here, the algorithm described in the paper was implemented in R as follows:
#' An adjacency matrix of all pairs of features is built, containing a 1 when 
#' the features have a (Pearson) correlation above the (0.9) threshold, a 
#' retention time difference between the (6) seconds threshold, and an m/z 
#' difference belonging to referenced adducts, isotopes and 
#' fragments m/z difference, and containing a 0 otherwise. The connex components
#' of this adjacency matrix are extracted ('igraph' package).
#' Within each component, the features are ranked 
#' by decreasing average intensity in samples; all features except the first one
#' are flagged as 'redundant'. Note: the algorithm relies on the 'mzdiff_db.tsv'
#' file referencing the known adducts, isotopes, and fragments.
#'
#' @param x An S4 object of class \code{SummarizedExperiment} 
#' or \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported): the dataset(s) must contain the dataMatrix 
#' and the variableMetadata (with the mz' and 'rt' columns)
#' @param cor_method.c character(1): correlation method (default: 'pearson')
#' @param cor_threshold.n numeric(1): correlation threshold (default: 0.9)
#' @param rt_tol.n numeric(1): retention time width in seconds (default: 6 s);
#' the time window may be increased when using hydrophilic interaction (HILIC)
#' chromatography
#' @param rt_colname.c character(1): column name for the retention time
#' in the rowData/fData (default: 'rt')
#' @param mzdiff_tol.n numeric(1): tolerance in Da for the matching of 
#' m/z differences and referenced adducts, isotopes, 
#' and fragments (default: 0.005 Da)
#' @param mz_colname.c character(1): column name for the m/z 
#' in the rowData/fData (default: 'mz')
#' @param return_adjacency.l logical(1): should the adjacency matrix be returned
#' (in addition to the updated SummarizedExperiment/ExpressionSet)?
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return updated \code{SummarizedExperiment} or \code{MultiAssayExperiment}
#' (or \code{ExpressionSet} and \code{MultiDataSet}): the 
#' SummarizedExperiment(s) (resp. ExpressionSet(s)) now include(s) 5 new columns 
#' in the rowData (resp. fData): redund_samp_mean', 'redund_is', redund_group', 
#' redund_iso_add_frag', redund_repres' and 'redund_relative' containing, 
#' respectively, the redundant features (coded by 1; i.e. features 
#' with a relative annotation distinct from '' and 'M'), the connected 
#' components, the m/z diff. chemical annotations, the representative ion 
#' of each group, and the annotations relative to this representative ion 
#' within each group
#' @export
#' @examples
#' metabo.se <- reading(system.file("extdata/prometis/metabo",
#'                                  package = "phenomis"),
#'                       report.c = "none")
#' metabo.se <- reducing(metabo.se,
#'                       rt_tol.n = 15)
#' # Note: in the 'prometis' example data set from this package, the chemical
#' # redundancy has already been filtered out
setGeneric("reducing",
           function(x,
                    cor_method.c = "pearson",
                    cor_threshold.n = 0.9,
                    rt_tol.n = 6,
                    rt_colname.c = "rt",
                    mzdiff_tol.n = 0.005,
                    mz_colname.c = "mz",
                    return_adjacency.l = FALSE,
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("reducing"))


#### transforming ####

#' Transformation of the data matrix intensities
#'
#' A logarithmic or square root transformation may be applied to the data matrix
#' intensities in (each of) the data set (e.g. to stabilize the variance)
#'
#' @param x An S4 object of class \code{SummarizedExperiment} 
#' or \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param method.vc character of length 1 or the total number of datasets:
#' transformation to be used for each dataset (either 'log2', 
#' 'log10', 'sqrt', or 'none')
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return \code{SummarizedExperiment} or \code{MultiAssayExperiment}
#' (or \code{ExpressionSet} and \code{MultiDataSet}) including the (list of) 
#' matrix with transformed intensities
#' @rdname transforming
#' @export
#' @examples
#' sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
#' sacurine.se <- correcting(sacurine.se)
#' sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
#' sacurine.se <- transforming(sacurine.se)
#' # MultiAssayExperiment
#' prometis.mae <- reading(system.file("extdata/prometis", 
#'                                     package = "phenomis"))
#' prometis.mae <- transforming(prometis.mae, method.vc = c("log2", "none"))
#' # Note: in the 'prometis' example data set from the package, the data are
#' # already log2 transformed
setGeneric("transforming",
           function(x,
                    method.vc = c("log2", "log10", "sqrt", "none")[1],
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("transforming"))


#### writing ####

#' Exporting a SummarizedExperiment (or MultiAssayExperiment) instance into 
#' (subfolders with) the 3 tabulated files 'dataMatrix.tsv', 
#' sampleMetadata.tsv', 'variableMetadata.tsv'
#'
#' Note that the \code{dataMatrix} is transposed before export (e.g., the 
#' samples are written column wise in the 'dataMatrix.tsv' exported file).
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or 
#' \code{MultiAssayExperiment} (\code{ExpressionSet} and \code{MultiDataSet}
#' are still supported)
#' @param dir.c character(1): directory where each dataset should be written
#' @param prefix.c character(1): prefix to be used (followed by '_') in the
#' 'dataMatrix.tsv', 'sampleMetadata.tsv', and 'variableMetadata.tsv' file names
#' @param files.ls list: alternatively to the dir.c argument, the full names of 
#' the files can be provided as a list
#' @param overwrite.l logical(1): should existing files be overwritten?
#' @param report.c character(1): File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return No object returned.
#' @rdname writing
#' @export
#' @examples
#' metabo.se <- reading(system.file("extdata/prometis/metabo", 
#'                      package="phenomis"))
#'\donttest{
#' writing(metabo.se, dir.c = file.path(getwd(), "metabo"))
#'}
#'# MultiAssayExperiment
#' prometis.mae <- reading(system.file("extdata/prometis",package="phenomis"))
#'\donttest{
#' writing(prometis.mae, dir.c = file.path(getwd(), "prometis"))
#' # alternatively
#' writing(prometis.mae,
#'          dir.c = NA,
#'          files.ls = list(metabo = list(dataMatrix = file.path(getwd(),
#'                                              "met_dataMatrix.tsv"),
#'                                       sampleMetadata = file.path(getwd(),
#'                                       "met_sampleMetadata.tsv"),
#'                                       variableMetadata = file.path(getwd(),
#'                                       "met_variableMetadata.tsv")),
#'                         proteo = list(dataMatrix = file.path(getwd(),
#'                                       "pro_dataMatrix.tsv"),
#'                                       sampleMetadata = file.path(getwd(),
#'                                       "pro_sampleMetadata.tsv"),
#'                                       variableMetadata = file.path(getwd(),
#'                                       "pro_variableMetadata.tsv"))))
#'}
setGeneric("writing",
           function(x,
                    dir.c,
                    prefix.c = "",
                    files.ls = NULL,
                    overwrite.l = FALSE,
                    report.c = c("none", "interactive", "myfile.txt")[2])
             standardGeneric("writing"))
