testthat::context("Testing 'reducing'")

testthat::test_that("reducing-se", {
  
  metabo.se <- reading(system.file("extdata/prometis/metabo",
                                   package = "phenomis"),
                       report.c = "none")
  metabo.se <- reducing(metabo.se)
  testthat::expect_identical(unname(rowData(metabo.se)["glycolic acid", 
                                                       "redund_relative"]),
                             "")
  testthat::expect_identical(as.numeric(table(rowData(metabo.se)[, 
                                                                 "redund_group"])),
                             numeric())
  
  
})
