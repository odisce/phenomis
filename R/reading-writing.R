#' reading
#'
#' Reading dataset(s) in the 3 tables 'dataMatrix' (or 'DM'), sampleMetadata'
#' (or 'SM') and variableMetadata' (or 'VM') tabular format. In case of a single
#' dataset (3 tables in the specified directory), a SummarizedExperiment 
#' instance is returned. In case of a multiple dataset (several subfolders 
#' containing 3 tables), a MultiAssayExperiment instance is created.
#'
#' @param dir.c character(1): directory containing the 3 tabular files (single
#' dataset), or containing several subdirectories with 3 tabular files
#' (multiple datasets)
#' @param subsets.vc character(): specifying a subset of the
#' subdirectories to be included in the MultiAssayExperiment (by default, all
#' subdirectories containing the 3 tables will be considered as datasets)
#' @param files.ls list: if dir.c is set to NA, the full names of the
#' individual files can be provided; in case of a SummarizedExperiment, the
#' names of the list must be 'dataMatrix', 'sampleMetadata',
#' and 'variableMetadata' with the corresponding file full names;
#' in case of a MultiAssayExperiment, the list must consists of one such sublist
#' per dataset
#' @param output.c character(1): Either 'exp' for SummarizedExperiment 
#' (or MultiAssayExperiment), or 'set' for ExpressionSet (or MultiDataSet) 
#' output formats (the latter are supported for convenience)
#' @param report.c character(1): File name for the printed results (call to
#' 'sink()'); if NA (default), messages will be printed on the screen; if NULL,
#' no verbose will be generated
#' @return \code{SummarizedExperiment} (one dataset) or 
#' \code{MultiAssayExperiment} (multiple datasets) instance containing the 
#' dataset(s)
#' @examples
#' data_dir.c <- system.file("extdata", package = "phenomis")
#' ## 1) Single set
#' sacurine_dir.c <- file.path(data_dir.c, "sacurine")
#' sacurine.se <- reading(sacurine_dir.c)
#' # or
#' sacurine.se <- reading(NA,
#'                   files.ls = list(dataMatrix = file.path(sacurine_dir.c,
#'                                   "Galaxy1_dataMatrix.tabular"),
#'                                   sampleMetadata = file.path(sacurine_dir.c,
#'                                   "Galaxy2_sampleMetadata.tabular"),
#'                                   variableMetadata = file.path(sacurine_dir.c,
#'                                   "Galaxy3_variableMetadata.tabular")))
#' ## 2) Multiple sets
#' prometis_dir.c <- file.path(data_dir.c, "prometis")
#' prometis.mae <- reading(prometis_dir.c)
#' metabo.mae <- reading(prometis_dir.c, subsets.vc = "metabo")
#' # or
#' prometis.mae <- reading(NA,
#'                        files.ls = list(metabo = list(dataMatrix = file.path(prometis_dir.c,
#'                                                            "metabo", "dataMatrix.tsv"),
#'                                                      sampleMetadata = file.path(prometis_dir.c,
#'                                                            "metabo", "sampleMetadata.tsv"),
#'                                                      variableMetadata = file.path(prometis_dir.c,
#'                                                            "metabo", "variableMetadata.tsv")),
#'                                        proteo = list(dataMatrix = file.path(prometis_dir.c,
#'                                                            "proteo", "dataMatrix.tsv"),
#'                                                      sampleMetadata = file.path(prometis_dir.c,
#'                                                            "proteo", "sampleMetadata.tsv"),
#'                                                      variableMetadata = file.path(prometis_dir.c,
#'                                                            "proteo", "variableMetadata.tsv"))))
#' @rdname reading
#' @export
reading <- function(dir.c,
                    files.ls = NULL,
                    subsets.vc = NA,
                    output.c = c("exp", "set")[1],
                    report.c = c("none", "interactive", "myfile.txt")[2]) {
  
  if (!(report.c %in% c("none", "interactive")))
    sink(report.c, append = TRUE)
  
  x <- NULL
  
  ## Parameter check
  
  stopifnot(output.c %in% c("exp", "set"))
  
  stopifnot(!is.na(dir.c) || !is.null(files.ls))
  stopifnot(!(!is.na(dir.c) && !is.null(files.ls)))
  
  ## Creating the ExpressionSet/SummarizedExperiment
  ## or building the list for the MultiDataSet/MultiAssayExperiment
  
  if (!is.na(dir.c)) {
    
    if (!file.exists(dir.c))
      stop("Directory '", dir.c, "' was not found.")
    
    if (!file.info(dir.c)[, "isdir"])
      stop(dir.c, "' is not a directory.")
    
    dir.vc <- dir(dir.c, full.names = TRUE)
    
    dir.vl <- file.info(dir.vc)[, "isdir"]
    
    subdir.vc <- dir.vc[dir.vl]
    
    if (length(subdir.vc) == 0) { ## ExpressionSet/SummarizedExperiment (x)
      
      x <- .read_se(dir.c,
                    dataMatrix = NA,
                    sampleMetadata = NA,
                    variableMetadata = NA,
                    output.c = output.c,
                    report.c = report.c)
      
    } else { ## MultiDataSet/MultiAssayExperiment (files.ls)
      
      files.ls <- .get_mae_files(subdir.vc = subdir.vc,
                                 report.c = report.c)
      
    }
    
  } else if (is.na(dir.c)) {
    
    subNamVc <- names(files.ls)
    
    if (sum(vapply(files.ls, is.list,
                   FUN.VALUE = logical(1))) == 0) { ## ExpressionSet/SummarizedExperiment (x)
      
      if (length(subNamVc) == 3 &&
          identical(subNamVc, c("dataMatrix",
                                "sampleMetadata",
                                "variableMetadata"))) {
        
        x <- .read_se(NA,
                      dataMatrix = files.ls[["dataMatrix"]],
                      sampleMetadata = files.ls[["sampleMetadata"]],
                      variableMetadata = files.ls[["variableMetadata"]],
                      output.c = output.c)
        
      } else {
        
        stop("'files.ls does not contain any sublist nor is a list 
             with names 'dataMatrix', 'sampleMetadata' and 'variableMetadata' 
             giving the corresponding file full names.")
        
      }
      
    } else { ## MultiDataSet/MultiAssayExperiment (files.ls)
      
      for (set.c in names(files.ls)) {
        
        set.ls <- files.ls[[set.c]]
        
        if (!identical(names(set.ls),
                       c("dataMatrix",
                         "sampleMetadata",
                         "variableMetadata"))) {
          
          if (report.c != "none")
            message("The names of the following sublist are not 'dataMatrix', 
                    'sampleMetadata' and 'variableMetadata':\n",
                    set.c,
                    "\nThe corresponding datasets will be skipped.")
          
          files.ls[[set.c]] <- NULL
          
        } else if (!file.exists(set.ls[["dataMatrix"]])) {
          
          if (report.c != "none")
            message("No 'dataMatrix' file was found in the following sublist:\n",
                    set.c,
                    "\nThe corresponding datasets will be skipped.")
          
          files.ls[[set.c]] <- NULL
          
        }
        
      }
      
      if (length(files.ls) == 0)
        stop("None of the provided sublists meets the requirements
             for the creation of a dataset.")
      
    }
    
  }
  
  ## Creating the MultiDataSet/MultiAssayExperiment
  
  if (is.null(x)) {
    
    mae.ls <- .read_mae_tables(files.ls = files.ls,
                               subsets.vc = subsets.vc,
                               output.c = output.c,
                               report.c = report.c)
    x <- mae.ls[["x"]]
    disagree.mc <- mae.ls[["disagree.mc"]]
    
  }
  
  # Printing
  
  if (report.c != "none") {
    
    if ((is(x, "MultiAssayExperiment") || is(x, "MultiDataSet")) &&
        length(x) > 1 && any(is.na(disagree.mc))) {
      warning("Discrepancies between the sampleMetadata from the datasets:")
      print(disagree.mc)
    }
    
    print(x)
  }
  
  if (!(report.c %in% c("none", "interactive")))
    sink()
  
  
  # Returning
  
  return(invisible(x))
  
}


.read_se <- function(dir.c,
                     dataMatrix = NA,
                     sampleMetadata = NA,
                     variableMetadata = NA,
                     output.c = "exp",
                     report.c = "interactive") {
  
  if (!is.na(dir.c) && !is.na(dataMatrix))
    stop("Either 'dir.c' or 'dataMatrix' argument must be set to NA")
  
  ## getting file paths
  
  if (!is.na(dir.c)) {
    
    tab_file.vc <- .get_se_files(dir.c = dir.c,
                                 report.c = report.c)
    
  } else {
    
    tab_file.vc <- c(dataMatrix = dataMatrix,
                     sampleMetadata = sampleMetadata,
                     variableMetadata = variableMetadata)
    
    for (tab.c in names(tab_file.vc)) {
      
      tab_file.c <- tab_file.vc[tab.c]
      
      if (!file.exists(tab_file.c)) {
        
        if (tab.c == "dataMatrix") {
          
          stop("The provided dataMatrix file was not found:\n",
               tab_file.c)
          
        } else if (!file.exists(tab_file.c)) {
          
          if (report.c != "none")
            warning("The following '", tab.c,
                    "' file was not found:\n", tab_file.c,
                    "\nThe corresponding '",
                    ifelse(tab.c == "sampleMetadata",
                           'phenoData',
                           'featureData'),
                    "' slot will be empty.",
                    immediate. = TRUE)
          
          tab_file.vc[tab.c] <- NA
          
        }
      }
    }
    
  }
  
  ## loading tables
  
  x_tables.ls <- .read_se_tables(tab_file.vc = tab_file.vc)
  
  tdat.mn <- x_tables.ls[["tdat.mn"]]
  sam.df <- x_tables.ls[["sam.df"]]
  var.df <- x_tables.ls[["var.df"]]
  
  chkLs <- .checkW4Mformat(t(tdat.mn), sam.df, var.df,
                           infCw = report.c)
  
  if (!chkLs[["chkL"]]) {
    stop("Sample and/or variable names do not match between your tables. 
         Use the 'report.c = NA' argument to get more feedback.")
  } else if (chkLs[["ordL"]]) {
    tdat.mn <- t(chkLs[["dat.mn"]])
  }
  
  # Building the ExpressionSet
  
  x <- Biobase::ExpressionSet(assayData = tdat.mn,
                              phenoData = Biobase::AnnotatedDataFrame(data = sam.df),
                              featureData = Biobase::AnnotatedDataFrame(data = var.df),
                              experimentData = Biobase::MIAME(title = ifelse(!is.na(dir.c),
                                                                             basename(dir.c),
                                                                             "")))
  # Converting to SummarizedExperiment
  
  if (output.c == "exp") {
    
    x <- as(x, "SummarizedExperiment")
    
  }
  
  methods::validObject(x)
  
  return(x)
  
}


.get_se_files <- function(dir.c,
                          report.c) {
  
  dm_file.c <- list.files(dir.c, pattern = "(dataMatrix|DM)",
                          full.names = TRUE)
  
  if (length(dm_file.c) == 0) {
    stop("No 'dataMatrix' (or 'DM') file was found in the directory:\n",
         dir.c)
  } else if (length(dm_file.c) > 1) {
    stop("Multiple 'dataMatrix' (or 'DM') files were found 
           in the directory:\n",
         dir.c)
  }
  
  sm_file.c <- list.files(dir.c, pattern = "(sampleMetadata|SM)",
                          full.names = TRUE)
  
  if (length(sm_file.c) == 0) {
    if (report.c != "none")
      warning("No 'sampleMetada' (or 'SM') file was found 
                in the directory:\n",
              dir.c,
              "\nThe corresponding 'pData' slot will be empty.")
    sm_file.c <- NA_character_
  } else if (length(sm_file.c) > 1) {
    stop("Multiple 'sampleMetadata' files were found in the directory:\n",
         dir.c)
  }
  
  vm_file.c <- list.files(dir.c,
                          pattern = "(variableMetadata|VM)",
                          full.names = TRUE)
  
  if (length(vm_file.c) == 0) {
    if (report.c != "none")
      warning("No 'variableMetadata' (or 'VM') file was found 
                in the directory:\n",
              dir.c,
              "\nThe corresponding 'fData' slot will be empty.")
    vm_file.c <- NA_character_
  } else if (length(vm_file.c) > 1) {
    stop("Multiple 'variableMetadata' files were found in the directory:\n",
         dir.c)
  }
  
  c(dataMatrix = dm_file.c,
    sampleMetadata = sm_file.c,
    variableMetadata = vm_file.c)
  
}


.read_se_tables <- function(tab_file.vc) {
  
  for (tab.c in names(tab_file.vc)) {
    
    tab_file.c <- tab_file.vc[tab.c]
    
    if (!is.na(tab_file.c)) {
      
      ## R standards for row and column names in matrices and data frames
      .checkRformat(tab_file.c)
      
      tab.df <- data.frame(data.table::fread(tab_file.c,
                                             header = TRUE,
                                             sep = "\t"),
                           check.names = FALSE,
                           row.names = 1,
                           stringsAsFactors = FALSE)
      
      # looking for duplicates in column names
      colname_dup.vi <- which(duplicated(colnames(tab.df)))
      if (length(colname_dup.vi)) {
        stop_coldup.c <- paste(colnames(tab.df[colname_dup.vi]),
                               collapse = "', '")
        stop("The '", tab.c, "' file has several columns named: '",
             stop_coldup.c, "'")
      }
      
      switch(tab.c,
             dataMatrix = {
               tdat.mn <- as.matrix(tab.df)
             },
             sampleMetadata = {
               sam.df <- tab.df
             },
             variableMetadata = {
               var.df <- tab.df
             })
      
    } else {
      
      switch(tab.c,
             sampleMetadata = {
               sam.df <- data.frame(row.names = colnames(tdat.mn))
             },
             variableMetadata = {
               var.df <- data.frame(row.names = rownames(tdat.mn))
             })
      
    }
  }
  
  list(tdat.mn = tdat.mn,
       sam.df = sam.df,
       var.df = var.df)
  
}


.get_mae_files <- function(subdir.vc,
                           report.c) {
  
  names(subdir.vc) <- basename(subdir.vc)
  
  subdir.vl <- vapply(subdir.vc,
                      function(sub_dir.c) {
                        fileC <- list.files(sub_dir.c,
                                            pattern = "(dataMatrix|DM)",
                                            full.names = TRUE)
                        length(fileC) == 1 && file.exists(fileC)
                      }, FUN.VALUE = logical(1))
  
  if (sum(subdir.vl) == 0) {
    
    stop_subdir.c <- paste(subdir.vc, collapse = "\n")
    stop("All subfolders have none or multiple 'dataMatrix' or 
               'DM' file(s):\n",
         stop_subdir.c,
         "\n", call. = FALSE)
    
  } else if (sum(!subdir.vl) > 0) {
    
    if (report.c != "none") {
      stop_subdir.c <- paste(subdir.vc[!subdir.vl], collapse = "\n")
      message("No or multiple 'dataMatrix' or 'DM' file(s) was/were found 
                    in the following subfolders:\n",
              stop_subdir.c,
              "\nThe corresponding datasets will be skipped.")
    }
    
    subdir.vc <- subdir.vc[subdir.vl]
    
  }
  
  files.ls <- vector(mode = "list", length = length(subdir.vc))
  names(files.ls) <- names(subdir.vc)
  
  for (set.c in names(files.ls)) {
    
    files.ls[[set.c]] <- list(list.files(subdir.vc[set.c],
                                         pattern = "(dataMatrix|DM)",
                                         full.names = TRUE),
                              list.files(subdir.vc[set.c],
                                         pattern = "(sampleMetadata|SM)",
                                         full.names = TRUE),
                              list.files(subdir.vc[set.c],
                                         pattern = "(variableMetadata|VM)",
                                         full.names = TRUE))
    
    names(files.ls[[set.c]]) <- c("dataMatrix",
                                  "sampleMetadata",
                                  "variableMetadata")
    
  }
  
  files.ls
  
}


.read_mae_tables <- function(files.ls,
                             subsets.vc,
                             output.c,
                             report.c) {
  
  # in case subsets have been specified
  
  if (!any(is.na(subsets.vc))) {
    
    missing_sets.vc <- subsets.vc[!(subsets.vc %in% names(files.ls))]
    
    if (length(missing_sets.vc)) {
      stop_missing.c <- paste(missing_sets.vc, collapse = "\n")
      stop("The following selected subsets were not found 
             in the subfolder(s) or sublist(s):\n",
           stop_missing.c)
    }
    
    files.ls <- files.ls[subsets.vc]
    
  }
  
  # uploading the sets as a list of ExpressionSet
  
  read_eset.ls <- .read_mae_eset_ls(files.ls = files.ls,
                                    report.c)
  
  eset.ls <- read_eset.ls[["eset.ls"]]
  disagree.mc <- read_eset.ls[["disagree.mc"]]
  
  # building the MultiDataSet or MultiAssayExperiment
  
  x <- .read_mae_build(eset.ls = eset.ls,
                       output.c = output.c)
  
  list(x = x,
       disagree.mc = disagree.mc)
  
} ## end .read_mae_tables


.read_mae_eset_ls <- function(files.ls,
                              report.c) {
  
  set.vc <- names(files.ls)
  set.n <- length(set.vc)
  
  # uploading the sets as ExpressionSet
  
  eset.ls <- lapply(set.vc,
                    function(set.c) {
                      if (report.c != "none")
                        message("Reading the '", set.c, "' dataset...")
                      .read_se(NA,
                               files.ls[[set.c]][["dataMatrix"]],
                               files.ls[[set.c]][["sampleMetadata"]],
                               files.ls[[set.c]][["variableMetadata"]],
                               output.c = "set")
                    })
  names(eset.ls) <- set.vc
  
  # checking the agreement between the sample metadata
  
  disagree.mc <- matrix("",
                        nrow = set.n,
                        ncol = set.n,
                        dimnames = list(set.vc, set.vc))
  
  if (length(eset.ls) > 1) {
    
    for (i in seq(1, set.n - 1, by = 1)) {
      for (j in seq(i + 1, set.n, by = 1)) {
        sam_i.df <- Biobase::pData(eset.ls[[i]])
        sam_j.df <- Biobase::pData(eset.ls[[j]])
        sam_com.vc <- intersect(rownames(sam_i.df), rownames(sam_j.df))
        sam_i.df <- sam_i.df[sam_com.vc, ]
        sam_j.df <- sam_j.df[sam_com.vc, ]
        sam_var.vc <- intersect(Biobase::varLabels(eset.ls[[i]]),
                                Biobase::varLabels(eset.ls[[j]]))
        disagree.vc <- ""
        if (length(sam_var.vc)) {
          for (sam_var.c in sam_var.vc) {
            if (!identical(sam_i.df[, sam_var.c],
                           sam_j.df[, sam_var.c])) {
              disagree.vc <- c(disagree.vc, sam_var.c)
            }
          }
        }
        disagree.mc[i, j] <- paste(disagree.vc, collapse = ", ")
      }
    }
    
  }
  
  list(eset.ls = eset.ls,
       disagree.mc = disagree.mc)
  
}

.read_mae_build <- function(eset.ls,
                            output.c) {
  
  if (output.c == "set") { # Creating the MultiDataSet
    
    x <- MultiDataSet::createMultiDataSet()
    
    for (set.c in names(eset.ls)) {
      
      eset <- eset.ls[[set.c]]
      
      eset$id <- rownames(Biobase::pData(eset))
      
      x <- MultiDataSet::add_eset(x, eset,
                                  dataset.type = set.c,
                                  GRanges = NA)
      
    }
    
  } else if (output.c == "exp") {  # Creating the MultiAssayExperiment
    
    sam_all.vc <- sort(unique(Reduce("union", lapply(eset.ls,
                                                     Biobase::sampleNames))))
    sam_var.vc <- setdiff(sort(unique(Reduce("intersect",
                                             lapply(eset.ls,
                                                    Biobase::varLabels)))),
                          c("sampleType", "injectionOrder", "batch"))
    
    first.eset <- eset.ls[[1]]
    sam.df <- Biobase::pData(first.eset)[, sam_var.vc, drop = FALSE]
    
    # removing the common sample metadata from the individual eset
    Biobase::pData(first.eset)[, ".id"] <- Biobase::sampleNames(first.eset)
    Biobase::pData(first.eset) <- Biobase::pData(first.eset)[,
                                                             setdiff(Biobase::varLabels(first.eset),
                                                                     sam_var.vc), drop = FALSE]
    eset.ls[[1]] <- first.eset
    
    map.ls <- vector(mode = "list", length = length(eset.ls))
    names(map.ls) <- names(eset.ls)
    map.ls[[1]] <- data.frame(primary = Biobase::sampleNames(first.eset),
                              colname = Biobase::sampleNames(first.eset))
    
    if (length(eset.ls) > 1) {
      for (set.c in names(eset.ls)[-1]) {
        
        set.eset <- eset.ls[[set.c]]
        sam_add.df <- Biobase::pData(set.eset)[, sam_var.vc, drop = FALSE]
        Biobase::pData(set.eset)[, ".id"] <- Biobase::sampleNames(set.eset)
        Biobase::pData(set.eset) <- Biobase::pData(set.eset)[,
                                                             setdiff(Biobase::varLabels(set.eset),
                                                                     sam_var.vc), drop = FALSE]
        
        sam_add.vc <- setdiff(rownames(sam_add.df), rownames(sam.df))
        
        if (length(sam_add.vc)) {
          sam.df <- rbind.data.frame(sam.df,
                                     sam_add.df[sam_add.vc, , drop = FALSE])
        }
        
        map.ls[[set.c]] <- data.frame(primary = Biobase::sampleNames(set.eset),
                                      colname = Biobase::sampleNames(set.eset))
        
        
        eset.ls[[set.c]] <- set.eset
      }
    }
    
    sam.df <- sam.df[sort(rownames(sam.df)), ]
    map.df <- MultiAssayExperiment::listToMap(map.ls)
    
    x <- MultiAssayExperiment::MultiAssayExperiment(experiments = lapply(eset.ls,
                                                                         function(eset)
                                                                           as(eset, "SummarizedExperiment")),
                                                    colData = sam.df,
                                                    sampleMap = map.df)
    
  }
  
  methods::validObject(x)
  
  x
  
}

#### writing (MultiAssayExperiment) ####

#' @rdname writing
#' @export
setMethod("writing", "MultiAssayExperiment",
          function(x,
                   dir.c,
                   prefix.c = "",
                   files.ls = NULL,
                   overwrite.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            sets.vc <- names(x)
            
            if (!is.na(dir.c)) {
              
              set_dir.vc <- file.path(dir.c, sets.vc)
              names(set_dir.vc) <- sets.vc
              
              if (file.exists(dir.c) && file.info(dir.c)[, "isdir"]) {
                
                dir.vc <- dir(dir.c, full.names = TRUE)
                
                dir.vl <- file.info(dir.vc)[, "isdir"]
                
                subdir.vc <- dir.vc[dir.vl]
                
                subDupVc <- intersect(set_dir.vc,
                                      subdir.vc)
                
                if (length(subDupVc) && !overwrite.l) {
                  stop_subdup.c <- paste(subDupVc, collapse = "\n")
                  stop("The following subfolder(s) were detected 
                       in your directory, please remove them or specify another 
                       parent directory to avoid overwriting:\n",
                       stop_subdup.c)
                }
                
              } else {
                
                dir.create(dir.c,
                           showWarnings = report.c != "none")
                
              }
              
              for (set.c in names(set_dir.vc)) {
                
                if (report.c != "none")
                  message("Writing the '", set.c, "' dataset...")
                
                set_dir.c <- set_dir.vc[set.c]
                
                if (!(file.exists(set_dir.c) &&
                      file.info(set_dir.c)[, "isdir"]))
                  dir.create(set_dir.c,
                             showWarnings = report.c != "none")
                
                set.se <- x[[set.c]]
                ## including the common sample metadata in each data set
                colData(set.se) <- cbind(colData(x)[colnames(set.se), ],
                                         colData(set.se))
                colData(set.se)[, ".id"] <- NULL
                
                writing(set.se,
                        set_dir.c,
                        prefix.c = prefix.c,
                        overwrite.l = overwrite.l,
                        report.c = report_set.c)
                
              }
              
              if (report.c != "none")
                message("The subfolders have been written in the directory:\n",
                        dir.c)
              
            } else if (is.na(dir.c)) {
              
              if (is.null(files.ls))
                stop("'files.ls' must be provided when 'dir.c' is set to NA")
              
              if (is.null(names(files.ls)) || any(is.na(names(files.ls))))
                stop("All names of the sublists must be provided
                     (they should match the names of the MultiDataSet
                     datasets)")
              
              filLisVl <- vapply(files.ls, is.list, FUN.VALUE = logical(1))
              
              if (!all(filLisVl)) {
                stop_sublist.c <- paste(names(filLisVl)[!filLisVl], collapse = "\n")
                stop("The following element(s) of 'files.ls' 
                     is/are not sublist(s):\n",
                     stop_sublist.c)
              }
              
              if (!identical(sets.vc, names(filLisVl))) {
                stop_set.c <- paste(sets.vc, collapse = ", ")
                stop_sublist.c <- paste(names(filLisVl), collapse = ", ")
                stop("The name(s) of the 'x' MultiAssayExperiment:\n",
                     stop_set.c,
                     "\ndo(es) not match the names of the sublists:\n",
                     stop_sublist.c)
              }
              
              for (set.c in sets.vc) {
                
                filLs <- files.ls[[set.c]]
                
                if (length(filLs) != 3 ||
                    !identical(names(filLs), c("dataMatrix",
                                               "sampleMetadata",
                                               "variableMetadata")))
                  stop("The names of the '", set.c,
                       "' sublist of 'files.ls are not identical to
                       'dataMatrix', 'sampleMetadata' and 'variableMetadata'.")
                
                if (is.na(filLs[["dataMatrix"]]))
                  stop("The 'dataMatrix' file name from the '", set.c,
                       "' sublist is missing (ie set to NA).")
                
                for (file.c in names(filLs)) {
                  
                  filFulNamC <- filLs[[file.c]]
                  
                  if (!is.na(filFulNamC)
                      && file.exists(filFulNamC)
                      && !overwrite.l)
                    stop("The following file from the '",
                         set.c, "' sublist already exists:\n",
                         filFulNamC,
                         "\nPlease remove it 
                         or choose another name to avoid overwriting.")
                  
                }
                
              }
              
              for (set.c in sets.vc) {
                
                filLs <- files.ls[[set.c]]
                
                if (report.c != "none")
                  message("Writing the '", set.c, "' dataset")
                
                set.se <- x[[set.c]]
                colData(set.se) <- cbind(colData(x)[colnames(set.se), ],
                                         colData(set.se))
                
                writing(set.se,
                        NA,
                        files.ls = list(dataMatrix = filLs[["dataMatrix"]],
                                        sampleMetadata = filLs[["sampleMetadata"]],
                                        variableMetadata = filLs[["variableMetadata"]]),
                        overwrite.l = overwrite.l,
                        report.c = report_set.c)
                
              }
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
          })



#### writing (SummarizedExperiment) ####

#' @rdname writing
#' @export
setMethod("writing", "SummarizedExperiment",
          function(x,
                   dir.c,
                   prefix.c = "",
                   files.ls = NULL,
                   overwrite.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]){
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (!is.na(dir.c)) {
              
              if (!(file.exists(dir.c) && file.info(dir.c)[, "isdir"]))
                dir.create(dir.c,
                           showWarnings = report.c != "none")
              
              if (prefix.c != "")
                prefix.c <- paste0(prefix.c, "_")
              
              tab_file.vc <- c(dataMatrix = file.path(dir.c,
                                                      paste0(prefix.c,
                                                             "dataMatrix.tsv")),
                               sampleMetadata = file.path(dir.c,
                                                          paste0(prefix.c,
                                                                 "sampleMetadata.tsv")),
                               variableMetadata = file.path(dir.c,
                                                            paste0(prefix.c,
                                                                   "variableMetadata.tsv")))
              
            } else if (is.na(dir.c)) {
              
              if (is.null(files.ls))
                stop("'files.ls' must be provided when 'dir.c' is set to NA",
                     call. = FALSE)
              
              if (length(files.ls) != 3 ||
                  !identical(names(files.ls), c("dataMatrix",
                                                "sampleMetadata",
                                                "variableMetadata")))
                stop("The names of the 'files.ls' list are not identical to
                     'dataMatrix', 'sampleMetadata' and 'variableMetadata'.")
              
              if (is.na(files.ls[["dataMatrix"]]))
                stop("The 'dataMatrix' file name from the 'files.ls' list 
                     is missing (ie set to NA).")
              
              tab_file.vc <- c(dataMatrix = files.ls[["dataMatrix"]],
                               sampleMetadata = files.ls[["sampleMetadata"]],
                               variableMetadata = files.ls[["variableMetadata"]])
              
            }
            
            for (tab.c in names(tab_file.vc)) {
              
              tab_file.c <- tab_file.vc[tab.c]
              
              if (!is.na(tab_file.c) && file.exists(tab_file.c) && !overwrite.l)
                stop("The following file already exists:\n", tab_file.c,
                     "\nPlease choose another file name.")
              
            }
            
            ## Writing
            
            tdat.mn <- SummarizedExperiment::assay(x)
            sam.df <- as.data.frame(SummarizedExperiment::colData(x))
            var.df <- as.data.frame(SummarizedExperiment::rowData(x))
            chkLs <- .checkW4Mformat(t(tdat.mn), sam.df, var.df, infCw = report.c)
            
            if (!chkLs[["chkL"]]) {
              stop("Sample and/or variable names do not match between 
                   your tables. Use the report.c = 'interactive' argument 
                   to get more feedback.")
            } else if (chkLs[["ordL"]]) {
              tdat.mn <- t(chkLs[["dat.mn"]])
            }
            
            datDF <- cbind.data.frame(dataMatrix = rownames(tdat.mn),
                                      as.data.frame(tdat.mn))
            
            utils::write.table(datDF,
                               file = tab_file.vc['dataMatrix'],
                               quote = FALSE,
                               row.names = FALSE,
                               sep = "\t")
            
            if (!is.na(dir.c) || !is.na(tab_file.vc["sampleMetadata"])) {
              
              sam.df <- cbind.data.frame(sampleMetadata = rownames(sam.df),
                                         sam.df)
              utils::write.table(sam.df,
                                 file = tab_file.vc['sampleMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
              
            }
            
            if (!is.na(dir.c) || !is.na(tab_file.vc["variableMetadata"])) {
              
              var.df <- cbind.data.frame(variableMetadata = rownames(var.df),
                                         var.df)
              utils::write.table(var.df,
                                 file = tab_file.vc['variableMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
              
            }
            
            if (report.c != "none") {
              message_file.c <- paste(tab_file.vc[!is.na(basename(tab_file.vc))],
                                      collapse = "\n")
              message("The following file(s) have been written:\n",
                      message_file.c)
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
          })


#### writing (MultiDataSet) ####

#' @rdname writing
#' @export
setMethod("writing", "MultiDataSet",
          function(x,
                   dir.c,
                   prefix.c = "",
                   files.ls = NULL,
                   overwrite.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            sets.vc <- names(x)
            
            if (!is.na(dir.c)) {
              
              set_dir.vc <- file.path(dir.c, sets.vc)
              names(set_dir.vc) <- sets.vc
              
              if (file.exists(dir.c) && file.info(dir.c)[, "isdir"]) {
                
                dir.vc <- dir(dir.c, full.names = TRUE)
                
                dir.vl <- file.info(dir.vc)[, "isdir"]
                
                subdir.vc <- dir.vc[dir.vl]
                
                subDupVc <- intersect(set_dir.vc,
                                      subdir.vc)
                
                if (length(subDupVc) && !overwrite.l) {
                  stop_subdup.c <- paste(subDupVc, collapse = "\n")
                  stop("The following subfolder(s) were detected 
                       in your directory, please remove them 
                       or specify another parent directory 
                       to avoid overwriting:\n",
                       stop_subdup.c)
                }
                
              } else {
                
                dir.create(dir.c,
                           showWarnings = report.c != "none")
                
              }
              
              for (set.c in names(set_dir.vc)) {
                
                if (report.c != "none")
                  message("Writing the '", set.c, "' dataset...")
                
                set_dir.c <- set_dir.vc[set.c]
                
                if (!(file.exists(set_dir.c) &&
                      file.info(set_dir.c)[, "isdir"]))
                  dir.create(set_dir.c,
                             showWarnings = report.c != "none")
                
                writing(x[[set.c]],
                        set_dir.c,
                        prefix.c = prefix.c,
                        overwrite.l = overwrite.l,
                        report.c = report_set.c)
                
              }
              
              if (report.c != "none")
                message("The subfolders have been written in the directory:\n",
                        dir.c)
              
            } else if (is.na(dir.c)) {
              
              if (is.null(files.ls))
                stop("'files.ls' must be provided when 'dir.c' is set to NA")
              
              if (is.null(names(files.ls)) || any(is.na(names(files.ls))))
                stop("All names of the sublists must be provided (they should 
                     match the names of the MultiDataSet datasets)")
              
              filLisVl <- vapply(files.ls, is.list, FUN.VALUE = logical(1))
              
              if (!all(filLisVl)) {
                stop_sublist.c <- paste(names(filLisVl)[!filLisVl],
                                        collapse = "\n")
                stop("The following element(s) of 'files.ls' 
                     is/are not sublist(s):\n",
                     stop_sublist.c)
              }
              
              if (!identical(sets.vc, names(filLisVl))) {
                stop_set.c <- paste(sets.vc, collapse = ", ")
                stop_sublist.c <- paste(names(filLisVl), collapse = ", ")
                stop("The name(s) of the 'x' MultiDataSet:\n",
                     stop_set.c,
                     "\ndo(es) not match the names of the sublists:\n",
                     stop_sublist.c)
              }
              
              for (set.c in sets.vc) {
                
                filLs <- files.ls[[set.c]]
                
                if (length(filLs) != 3 ||
                    !identical(names(filLs), c("dataMatrix",
                                               "sampleMetadata",
                                               "variableMetadata")))
                  stop("The names of the '", set.c,
                       "' sublist of 'files.ls are not identical to
                       'dataMatrix', 'sampleMetadata' and 'variableMetadata'.")
                
                if (is.na(filLs[["dataMatrix"]]))
                  stop("The 'dataMatrix' file name from the '", set.c,
                       "' sublist is missing (ie set to NA).")
                
                for (file.c in names(filLs)) {
                  
                  filFulNamC <- filLs[[file.c]]
                  
                  if (!is.na(filFulNamC)
                      && file.exists(filFulNamC)
                      && !overwrite.l)
                    stop("The following file from the '",
                         set.c, "' sublist already exists:\n",
                         filFulNamC,
                         "\nPlease remove it or choose another name 
                         to avoid overwriting.")
                  
                }
                
              }
              
              for (set.c in sets.vc) {
                
                filLs <- files.ls[[set.c]]
                
                if (report.c != "none")
                  message("Writing the '", set.c, "' dataset")
                
                writing(x[[set.c]],
                        NA,
                        files.ls = list(dataMatrix = filLs[["dataMatrix"]],
                                        sampleMetadata = filLs[["sampleMetadata"]],
                                        variableMetadata = filLs[["variableMetadata"]]),
                        overwrite.l = overwrite.l,
                        report.c = report_set.c)
                
              }
            }
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
          })


#### writing (ExpressionSet) ####

#' @rdname writing
#' @export
setMethod("writing", "ExpressionSet",
          function(x,
                   dir.c,
                   prefix.c = "",
                   files.ls = NULL,
                   overwrite.l = FALSE,
                   report.c = c("none", "interactive", "myfile.txt")[2]){
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            if (!is.na(dir.c)) {
              
              if (!(file.exists(dir.c) && file.info(dir.c)[, "isdir"]))
                dir.create(dir.c,
                           showWarnings = report.c != "none")
              
              if (prefix.c != "")
                prefix.c <- paste0(prefix.c, "_")
              
              tab_file.vc <- c(dataMatrix = file.path(dir.c,
                                                      paste0(prefix.c,
                                                             "dataMatrix.tsv")),
                               sampleMetadata = file.path(dir.c,
                                                          paste0(prefix.c,
                                                                 "sampleMetadata.tsv")),
                               variableMetadata = file.path(dir.c,
                                                            paste0(prefix.c,
                                                                   "variableMetadata.tsv")))
              
            } else if (is.na(dir.c)) {
              
              if (is.null(files.ls))
                stop("'files.ls' must be provided when 'dir.c' is set to NA")
              
              if (length(files.ls) != 3 ||
                  !identical(names(files.ls), c("dataMatrix",
                                                "sampleMetadata",
                                                "variableMetadata")))
                stop("The names of the 'files.ls' list are not identical to
                     'dataMatrix', 'sampleMetadata' and 'variableMetadata'.")
              
              if (is.na(files.ls[["dataMatrix"]]))
                stop("The 'dataMatrix' file name from the 'files.ls' list 
                     is missing (ie set to NA).")
              
              tab_file.vc <- c(dataMatrix = files.ls[["dataMatrix"]],
                               sampleMetadata = files.ls[["sampleMetadata"]],
                               variableMetadata = files.ls[["variableMetadata"]])
              
            }
            
            for (tab.c in names(tab_file.vc)) {
              
              tab_file.c <- tab_file.vc[tab.c]
              
              if (!is.na(tab_file.c) && file.exists(tab_file.c) && !overwrite.l)
                stop("The following file already exists:\n", tab_file.c,
                     "\nPlease choose another file name.")
              
            }
            
            ## Writing
            
            tdat.mn <- Biobase::exprs(x)
            sam.df <- Biobase::pData(x)
            var.df <- Biobase::fData(x)
            chkLs <- .checkW4Mformat(t(tdat.mn), sam.df, var.df, 
                                     infCw = report.c)
            
            if (!chkLs[["chkL"]]) {
              stop("Sample and/or variable names do not match between 
                   your tables. Use the report.c = 'interactive' argument 
                   to get more feedback.")
            } else if (chkLs[["ordL"]]) {
              tdat.mn <- t(chkLs[["dat.mn"]])
            }
            
            datDF <- cbind.data.frame(dataMatrix = rownames(tdat.mn),
                                      as.data.frame(tdat.mn))
            
            utils::write.table(datDF,
                               file = tab_file.vc['dataMatrix'],
                               quote = FALSE,
                               row.names = FALSE,
                               sep = "\t")
            
            if (!is.na(dir.c) || !is.na(tab_file.vc["sampleMetadata"])) {
              
              sam.df <- cbind.data.frame(sampleMetadata = rownames(sam.df),
                                         as.data.frame(sam.df))
              utils::write.table(sam.df,
                                 file = tab_file.vc['sampleMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
              
            }
            
            if (!is.na(dir.c) || !is.na(tab_file.vc["variableMetadata"])) {
              
              var.df <- cbind.data.frame(variableMetadata = rownames(var.df),
                                         as.data.frame(var.df))
              utils::write.table(var.df,
                                 file = tab_file.vc['variableMetadata'],
                                 quote = FALSE,
                                 row.names = FALSE,
                                 sep = "\t")
              
            }
            
            if (report.c != "none")
              message("The following file(s) have been written:\n",
                      paste(tab_file.vc[!is.na(basename(tab_file.vc))], 
                            collapse = "\n"))
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
          })


.checkRformat <- function(filCa) {
  
  rowVc <- data.table::fread(filCa,
                             header = TRUE,
                             sep = "\t")[[1]]
  
  colVc <- unlist(data.table::fread(filCa,
                                    header = FALSE,
                                    nrows = 1,
                                    sep = "\t")[1])[-1]
  
  if (any(duplicated(rowVc))) {
    stop_rowdup.c <- paste(rowVc[duplicated(rowVc)], collapse = "', '")
    stop("The following ",
         ifelse(names(filCa) == 'sampleMetadata', 'sample', 'variable'),
         " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         stop_rowdup.c, "'")
  }
  
  if (any(duplicated(colVc))) {
    stop_dupcol.c <- paste(colVc[duplicated(colVc)], collapse = "', '")
    stop("The following ", ifelse(names(filCa) == 'sampleMetadata', 'variable', 
                                  'sample'), 
         " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         stop_dupcol.c, "'")
  }
  
}


.checkW4Mformat <- function(dat.mnw, sam.dfw, var.dfw,
                            infCw = "interactive") {
  
  chkL <- TRUE
  ordL <- FALSE
  
  if (mode(dat.mnw) != "numeric") {
    message("The dataMatrix is not of the 'numeric' type.")
    chkL <- FALSE
  }
  
  if (!identical(rownames(dat.mnw), rownames(sam.dfw))) {
    ## checking sample names
    
    datSamDifVc <- setdiff(rownames(dat.mnw), rownames(sam.dfw))
    
    if (length(datSamDifVc)) {
      if (infCw != "none") {
        message("The following samples were found in the dataMatrix column names
                but not in the sampleMetadata row names:\n")
        print(cbind.data.frame(col = as.numeric(vapply(datSamDifVc,
                                                       function(samC)
                                                         which(rownames(dat.mnw) == samC),
                                                       integer(1))),
                               name = datSamDifVc))
      }
      chkL <- FALSE
    }
    
    samDatDifVc <- setdiff(rownames(sam.dfw), rownames(dat.mnw))
    
    if (length(samDatDifVc)) {
      if (infCw != "none") {
        message("The following samples were found in the sampleMetadata row 
                names but not in the dataMatrix column names:\n")
        print(cbind.data.frame(row = as.numeric(vapply(samDatDifVc,
                                                       function(samC)
                                                         which(rownames(sam.dfw) == samC), FUN.VALUE = integer(1))),
                               name = samDatDifVc))
      }
      chkL <- FALSE
    }
    
    if (nrow(dat.mnw) != nrow(sam.dfw)) {
      if (infCw != "none")
        message("The dataMatrix has ", nrow(dat.mnw), " columns (ie samples) 
                whereas the sampleMetadata has ", nrow(sam.dfw), " rows.")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(dat.mnw)),
                         rownames(sam.dfw))) {
      if (infCw != "none")
        message("The dataMatrix column names start with an 'X' but not 
                the sampleMetadata row names.")
      chkL <- FALSE
    } else if (identical(gsub("^X", "", rownames(sam.dfw)), 
                         rownames(dat.mnw))) {
      if (infCw != "none")
        message("The sampleMetadata row names start with an 'X' but not 
                the dataMatrix column names.")
      chkL <- FALSE
    } else if (identical(sort(rownames(dat.mnw)), sort(rownames(sam.dfw)))) {
      if (infCw != "none")
        message("Re-ordering dataMatrix sample names to match sampleMetadata.")
      dat.mnw <- dat.mnw[rownames(sam.dfw), , drop = FALSE]
      stopifnot(identical(sort(rownames(dat.mnw)), sort(rownames(sam.dfw))))
      ordL <- TRUE
    } else {
      if (infCw != "none") {
        message("The dataMatrix column names and the sampleMetadata row names 
                are not identical:\n")
        print(cbind.data.frame(indice = seq_len(nrow(dat.mnw)),
                               dataMatrix_columnnames = rownames(dat.mnw),
                               sampleMetadata_rownames = rownames(sam.dfw))[rownames(dat.mnw) != rownames(sam.dfw), , drop = FALSE])
      }
      chkL <- FALSE
    }
    
  }
  
  if (!identical(colnames(dat.mnw), rownames(var.dfw))) {
    ## checking variable names
    
    datVarDifVc <- setdiff(colnames(dat.mnw), rownames(var.dfw))
    
    if (length(datVarDifVc)) {
      if (infCw != "none") {
        message("The following variables were found in the dataMatrix row names
                but not in the variableMetadata row names:\n")
        print(cbind.data.frame(row = as.numeric(vapply(datVarDifVc,
                                                       function(varC)
                                                         which(colnames(dat.mnw) == varC),
                                                       FUN.VALUE = integer(1))),
                               name = datVarDifVc))
      }
      chkL <- FALSE
    }
    
    varDatDifVc <- setdiff(rownames(var.dfw), colnames(dat.mnw))
    
    if (length(varDatDifVc)) {
      if (infCw != "none") {
        message("The following variables were found in the variableMetadata 
                row names but not in the dataMatrix row names:\n")
        print(cbind.data.frame(row = as.numeric(vapply(varDatDifVc,
                                                       function(varC)
                                                         which(rownames(var.dfw) == varC), FUN.VALUE = integer(1))),
                               name = varDatDifVc))
      }
      chkL <- FALSE
    }
    
    if (ncol(dat.mnw) != nrow(var.dfw)) {
      if (infCw != "none")
        message("The dataMatrix has ",
                nrow(dat.mnw),
                " rows (ie variables) whereas the variableMetadata has ",
                nrow(var.dfw),
                " rows.")
      chkL <- FALSE
    } else if (identical(sort(colnames(dat.mnw)), sort(rownames(var.dfw)))) {
      if (infCw != "none")
        message("Re-ordering dataMatrix variable names to match 
                variableMetadata.")
      dat.mnw <- dat.mnw[, rownames(var.dfw), drop = FALSE]
      stopifnot(identical(sort(colnames(dat.mnw)), sort(rownames(var.dfw))))
      ordL <- TRUE
    } else {
      if (infCw != "none") {
        message("The dataMatrix row names and the variableMetadata row names 
                are not identical:\n")
        print(cbind.data.frame(row = seq_len(ncol(dat.mnw)),
                               dataMatrix_rownames = colnames(dat.mnw),
                               variableMetadata_rownames = rownames(var.dfw))[colnames(dat.mnw) != rownames(var.dfw), , drop = FALSE])
      }
      chkL <- FALSE
    }
  }
  
  return(list(chkL = chkL,
              ordL = ordL,
              dat.mn = dat.mnw))
  
}















