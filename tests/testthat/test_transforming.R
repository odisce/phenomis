testthat::context("Testing 'transforming'")

testthat::test_that("transforming-se", {
  
  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
  sacurine.se <- correcting(sacurine.se, figure.c = "none")
  sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
  sacurine.se <- transforming(sacurine.se)
  
  testthat::expect_equal(assay(sacurine.se)["Testosterone glucuronide", "HU_neg_020"],
                         15.4297,
                         tolerance = 1e-6)
  
})

testthat::test_that("transforming-mae", {
  
  prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  
  prometis.mae1 <- transforming(prometis.mae)
  testthat::expect_equivalent(assays(prometis.mae1)[["metabo"]]["isovaleric acid", "L824m"],
                              4.348297,
                              tolerance = 1e-6)
  
  prometis.mae2 <- transforming(prometis.mae, "sqrt")
  testthat::expect_equivalent(assays(prometis.mae2)[["metabo"]]["isovaleric acid", "L824m"],
                              4.513194,
                              tolerance = 1e-6)
  
})


testthat::test_that("transforming-mset", {
  
  prometis.mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
  
  prometis.mset1 <- transforming(prometis.mset)
  testthat::expect_equal(Biobase::exprs(prometis.mset1[["metabo"]])["isovaleric acid", "L824m"],
                         4.348297,
                         tolerance = 1e-6)
  
  prometis.mset2 <- transforming(prometis.mset, "sqrt")
  testthat::expect_equal(Biobase::exprs(prometis.mset2[["metabo"]])["isovaleric acid", "L824m"],
                         4.513194,
                         tolerance = 1e-6)
  
})
