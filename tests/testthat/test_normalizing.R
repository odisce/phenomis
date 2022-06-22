testthat::context("Testing 'normalizing'")

testthat::test_that("normalizing-se", {
  
  sacurine.se <- reading(file.path(system.file(package = "phenomis"),
                                   "extdata/sacurine"))
  sacurine.se <- sacurine.se[, colnames(sacurine.se) != 'HU_neg_096_b2']
  norm.se <- normalizing(sacurine.se, method.vc = "pqn")
  
  testthat::expect_equivalent(assay(norm.se)["Testosterone glucuronide", "HU_neg_020"],
                              148872.7,
                              tolerance = 1e-6)
  
})

testthat::test_that("normalizing-mae", {
  
  mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  
  norm.mae <- normalizing(mae, method.vc = "pqn")
  
  testthat::expect_equivalent(assays(norm.mae)[["metabo"]]["isovaleric acid", "L824m"],
                              20.25144,
                              tolerance = 1e-6)
  
})


testthat::test_that("normalizing-mset", {
  
  mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
  
  norm.mset <- normalizing(mset, method.vc = "pqn")
  
  testthat::expect_equivalent(Biobase::exprs(norm.mset[["metabo"]])["oxalic acid", "L819f"],
                              23.48427,
                              tolerance = 1e-6)
  
})
