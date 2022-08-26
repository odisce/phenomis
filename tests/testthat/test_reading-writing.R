testthat::context("Testing 'read-write'")

#### .reading ####

testthat::test_that(".reading", {
  
  sacdir.c <- system.file("extdata/sacurine", package = "phenomis")
  
  sacurine.se1 <- .read_se(sacdir.c)
  
  testthat::expect_true(class(sacurine.se1) == "SummarizedExperiment")
  
  testthat::expect_equal(assay(sacurine.se1)[1, 1],
                         477491,
                         tolerance = 1e-3)
  # alternatively
  sacurine.se2 <- .read_se(NA,
                           file.path(sacdir.c, "Galaxy1_dataMatrix.tabular"),
                           file.path(sacdir.c, 
                                     "Galaxy2_sampleMetadata.tabular"),
                           file.path(sacdir.c, 
                                     "Galaxy3_variableMetadata.tabular"))
  
  testthat::expect_true(class(sacurine.se2) == "SummarizedExperiment")
  
  testthat::expect_equal(assay(sacurine.se2)[1, 1],
                         477491,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(.read_se(NA,
                                  file.path(sacdir.c, "Galaxy1_dataMatrix.tsv"),
                                  file.path(sacdir.c, 
                                            "Galaxy2_sampleMetadata.tabular"),
                                  file.path(sacdir.c, 
                                            "Galaxy3_variableMetadata.tabular")))
  
})

#### reading_se ####

testthat::test_that("reading_SummarizedExperiment", {
  
  sacdir.c <- system.file("extdata/sacurine", package = "phenomis")
  
  sacurine.se2 <- reading(sacdir.c)
  
  testthat::expect_true(class(sacurine.se2) == "SummarizedExperiment")
  
  testthat::expect_equal(assay(sacurine.se2)[1, 1],
                         477491,
                         tolerance = 1e-3)
  # alternatively
  sacurine.se3 <- reading(NA,
                          files.ls = list(dataMatrix = file.path(sacdir.c, "Galaxy1_dataMatrix.tabular"),
                                          sampleMetadata = file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
                                          variableMetadata = file.path(sacdir.c, "Galaxy3_variableMetadata.tabular")))
  
  testthat::expect_true(class(sacurine.se3) == "SummarizedExperiment")
  
  testthat::expect_equal(assay(sacurine.se3)[1, 1],
                         477491,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(reading(NA,
                                 files.ls = list(dataMatrix = file.path(sacdir.c, "Galaxy1_dataMatrix.tsv"),
                                                 sampleMetadata = file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
                                                 variableMetadata = file.path(sacdir.c, "Galaxy3_variableMetadata.tabular"))))
  
})

#### reading_mae ####

testthat::test_that("reading_MultiAssayExperiment", {
  
  prometis_dir.c <- system.file("extdata/prometis", package = "phenomis")
  
  prometis.mae <- reading(prometis_dir.c)
  
  testthat::expect_true(class(prometis.mae) == "MultiAssayExperiment")
  
  testthat::expect_identical(names(prometis.mae), c("metabo", "proteo"))
  
  testthat::expect_equal(MultiAssayExperiment::assays(prometis.mae)[["metabo"]]["isobutyric acid", "L818f"],
                         24.64327,
                         tolerance = 1e-3)
  
  
  met.mae <- reading(prometis_dir.c, subsets.vc = "metabo")
  
  testthat::expect_true(class(met.mae) == "MultiAssayExperiment")
  
  testthat::expect_identical(names(met.mae), "metabo")
  
  testthat::expect_equal(MultiAssayExperiment::assays(met.mae)[["metabo"]]["isobutyric acid", "L818f"],
                         24.64327,
                         tolerance = 1e-3)
  
})

#### reading_mds ####

testthat::test_that("reading_MultiDataSet", {
  
  prometis_dir.c <- system.file("extdata/prometis", package = "phenomis")
  
  prometis.mset1 <- reading(prometis_dir.c, output.c = "set")
  
  testthat::expect_true(class(prometis.mset1) == "MultiDataSet")
  
  testthat::expect_identical(names(prometis.mset1), c("metabo", "proteo"))
  
  testthat::expect_equal(Biobase::exprs(prometis.mset1[["metabo"]])["isobutyric acid", "L818f"],
                         24.64327,
                         tolerance = 1e-3)
  
  
  prometis.mset2 <- reading(NA,
                            files.ls = list(metabo = list(dataMatrix = file.path(prometis_dir.c, "metabo", "dataMatrix.tsv"),
                                                          sampleMetadata = file.path(prometis_dir.c, "metabo", "sampleMetadata.tsv"),
                                                          variableMetadata = file.path(prometis_dir.c, "metabo", "variableMetadata.tsv")),
                                            proteo = list(dataMatrix = file.path(prometis_dir.c, "proteo", "dataMatrix.tsv"),
                                                          sampleMetadata = file.path(prometis_dir.c, "proteo", "sampleMetadata.tsv"),
                                                          variableMetadata = file.path(prometis_dir.c, "proteo", "variableMetadata.tsv"))),
                            output.c = "set")
  
  testthat::expect_true(class(prometis.mset2) == "MultiDataSet")
  
  testthat::expect_identical(names(prometis.mset2), c("metabo", "proteo"))
  
  testthat::expect_equal(Biobase::exprs(prometis.mset2[["metabo"]])["isobutyric acid", "L818f"],
                         24.64327,
                         tolerance = 1e-3)
  
  testthat::expect_identical(colnames(Biobase::pData(prometis.mset2[["metabo"]])),
                             c("gene", "mouse_id", "sex", "id"))
  
  
  
  metMset1 <- reading(prometis_dir.c, subsets.vc = "metabo", output.c = "set")
  
  testthat::expect_true(class(metMset1) == "MultiDataSet")
  
  testthat::expect_identical(names(metMset1), "metabo")
  
  testthat::expect_equal(Biobase::exprs(metMset1[["metabo"]])["isobutyric acid",
                                                              "L818f"],
                         24.64327,
                         tolerance = 1e-3)
  
  
  metMset2 <- reading(NA,
                      files.ls = list(metabo = list(dataMatrix = file.path(prometis_dir.c, "metabo", "dataMatrix.tsv"),
                                                    sampleMetadata = file.path(prometis_dir.c, "metabo", "sampleMetadata.tsv"),
                                                    variableMetadata = file.path(prometis_dir.c, "metabo", "variableMetadata.tsv")),
                                      proteo = list(dataMatrix = file.path(prometis_dir.c, "proteo", "dataMatrix.tsv"),
                                                    sampleMetadata = file.path(prometis_dir.c, "proteo", "sampleMetadata.tsv"),
                                                    variableMetadata = file.path(prometis_dir.c, "proteo", "variableMetadata.tsv"))),
                      subsets.vc = "metabo", output.c = "set")
  
  testthat::expect_true(class(metMset2) == "MultiDataSet")
  
  testthat::expect_identical(names(metMset2), "metabo")
  
  testthat::expect_equal(Biobase::exprs(metMset2[["metabo"]])["isobutyric acid",
                                                              "L818f"],
                         24.64327,
                         tolerance = 1e-3)
  
  
  testthat::expect_error(reading(NA,
                                 list(metabo = list(dataMatrix = file.path(prometis_dir.c, "metabo", "dataMatrix.tsv_XXX"),
                                                    sampleMetadata = file.path(prometis_dir.c, "metabo", "sampleMetadata.tsv"),
                                                    variableMetadata = file.path(prometis_dir.c, "metabo", "variableMetadata.tsv")),
                                      proteo = list(dataMatrix = file.path(prometis_dir.c, "proteo", "dataMatrix.tsv_XXX"),
                                                    sampleMetadata = file.path(prometis_dir.c, "proteo", "sampleMetadata.tsv"),
                                                    variableMetadata = file.path(prometis_dir.c, "proteo", "variableMetadata.tsv"))),
                                 output.c = "set"))
  
  testthat::expect_warning(reading(NA,
                                   list(metabo = list(dataMatrix = file.path(prometis_dir.c, "metabo", "dataMatrix.tsv"),
                                                      sampleMetadata = file.path(prometis_dir.c, "metabo", "sampleMetadata.tsv_XXX"),
                                                      variableMetadata = file.path(prometis_dir.c, "metabo", "variableMetadata.tsv")),
                                        proteo = list(dataMatrix = file.path(prometis_dir.c, "proteo", "dataMatrix.tsv"),
                                                      sampleMetadata = file.path(prometis_dir.c, "proteo", "sampleMetadata.tsv"),
                                                      variableMetadata = file.path(prometis_dir.c, "proteo", "variableMetadata.tsv"))),
                                   output.c = "set"))
  
  testthat::expect_true(identical(colnames(Biobase::pData(prometis.mset2[["metabo"]])),
                                  c("gene", "mouse_id", "sex", "id")))
  
  
})

#### writing ####

testthat::test_that("writing_MultiDataSet", {
  
  prometis_dir.c <- system.file("extdata/prometis", package = "phenomis")
  
  prometis.mset4 <- reading(prometis_dir.c, output.c = "set")
  
  testthat::expect_error(writing(prometis.mset4,
                                 dir.c = NA,
                                 files.ls = list(metabo = list(dataMatrix = NA,
                                                               sampleMetadata = file.path(getwd(), "metabo_sampleMetadata.tsv"),
                                                               variableMetadata = file.path(getwd(), "metabo_variableMetadata.tsv")),
                                                 proteo = list(dataMatrix = file.path(getwd(), "proteo_dataMatrix.tsv"),
                                                               sampleMetadata = file.path(getwd(), "proteo_sampleMetadata.tsv"),
                                                               variableMetadata = file.path(getwd(), "proteo_variableMetadata.tsv")))))
  
})
