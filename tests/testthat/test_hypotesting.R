testthat::context("Testing 'hypotesting'")

#### ttest ####

testthat::test_that("ttest", {
  
  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
  sacurine.se <- correcting(sacurine.se, figure.c = "none")
  sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
  sacurine.se <- transforming(sacurine.se)
  sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
  
  # .twoSampCorTests
  
  ttestLs <- .twoSampCorTests(data.mn = t(assay(sacurine.se)),
                              samp.df = colData(sacurine.se),
                              feat.df = rowData(sacurine.se),
                              test.c = "ttest",
                              factorNameC = "gender",
                              factorLevelsVc = "default",
                              adjust.c = "BH",
                              adjust_thresh.n = 0.05,
                              prefix.c = "",
                              figure.c = "none")
  testthat::expect_equivalent(ttestLs[["feat.df"]]["Testosterone glucuronide", "ttest_gender_Female.Male_diff"],
                              2.42603,
                              tolerance = 1e-6)
  testthat::expect_equivalent(ttestLs[["feat.df"]]["Testosterone glucuronide", "ttest_gender_Female.Male_BH"],
                              9.054552e-10,
                              tolerance = 1e-6)
  testthat::expect_equivalent(ttestLs[["feat.df"]]["1,7-Dimethyluric acid", "ttest_gender_Female.Male_BH"],
                              0.5868704,
                              tolerance = 1e-6)
  
  
  # se
  
  sacurine.se <- hypotesting(sacurine.se,
                             test.c = "ttest",
                             factor_names.vc = "gender",
                             figure.c = "none",
                             report.c = "none")
  
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "ttest_gender_Female.Male_diff"],
                              2.42603,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "ttest_gender_Female.Male_BH"],
                              9.054552e-10,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["1,7-Dimethyluric acid", "ttest_gender_Female.Male_BH"],
                              0.5868704,
                              tolerance = 1e-6)
  
  sacurine.se <- hypotesting(sacurine.se,
                             test.c = "ttest",
                             factor_names.vc = "gender",
                             factor_levels.ls = list(factor1 = c("Male", "Female")),
                             figure.c = "none",
                             report.c = "none")
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "ttest_gender_Male.Female_diff"],
                              -2.42603,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "ttest_gender_Male.Female_BH"],
                              9.054552e-10,
                              tolerance = 1e-6)
  
  sacurine.se <- hypotesting(sacurine.se,
                             test.c = "wilcoxon",
                             factor_names.vc = "gender",
                             figure.c = "none",
                             report.c = "none")
  
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "wilcoxon_gender_Female.Male_diff"],
                              1.919634,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "wilcoxon_gender_Female.Male_BH"],
                              3.957732e-12,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["1-Methylxanthine", "wilcoxon_gender_Female.Male_BH"],
                              0.03904366,
                              tolerance = 1e-6)
  
  sacurine.se <- hypotesting(sacurine.se,
                             test.c = "limma",
                             factor_names.vc = "gender",
                             figure.c = "none",
                             report.c = "none")
  
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "limma_gender_Female.Male_diff"],
                              2.42603,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["Testosterone glucuronide", "limma_gender_Female.Male_BH"],
                              1.221723e-11,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(sacurine.se)["1-Methylxanthine", "limma_gender_Female.Male_BH"],
                              0.07010866,
                              tolerance = 1e-6)
  
  sacurine.se <- hypotesting(sacurine.se, "pearson", "age")
  
  # mset
  
  prometis.mset <- reading(system.file("extdata/prometis", package = "phenomis"), output.c = "set")
  
  prometis.mset <- hypotesting(prometis.mset,
                               test.c = "ttest",
                               factor_names.vc = "gene",
                               figure.c = "none",
                               report.c = "none")
  
  testthat::expect_identical(vapply(Biobase::fData(prometis.mset),
                                    function(fdaDF) {
                                      sum(fdaDF[, "ttest_gene_LAT.WT_signif"])
                                    }, FUN.VALUE = double(1)),
                             c(metabo = 33, proteo = 3))
  
})

#### anova ####

testthat::test_that("anova", {
  
  sacurine.se <- reading(system.file("extdata/sacurine",
                                     package = "phenomis"))
  sacurine.se <- correcting(sacurine.se, figure.c = "none")
  sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
  sacurine.se <- transforming(sacurine.se)
  sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
  colData(sacurine.se)[, "ageGroup"] <- vapply(colData(sacurine.se)[, "age"],
                                               function(x) {
                                                 if (x < 35) {
                                                   return("thirty")
                                                 } else if (x < 50) {
                                                   return("fourty")
                                                 } else {
                                                   return("fifty")}},
                                               FUN.VALUE = character(1))
  feat.df <- .anovas(data.mn = t(assay(sacurine.se)),
                     samp.df = colData(sacurine.se),
                     feat.df = rowData(sacurine.se),
                     test.c = "anova",
                     factorNameC = "ageGroup",
                     factorLevelsVc = "default",
                     adjust.c = "BH",
                     adjust_thresh.n = 0.05,
                     prefix.c = "",
                     figure.c = "none")[["feat.df"]]
  testthat::expect_equivalent(feat.df["1-Methylxanthine", "anova_ageGroup_BH"],
                              0.003023911,
                              tolerance = 1e-6)
  testthat::expect_equivalent(feat.df["1-Methylxanthine", "anova_ageGroup_fourty.thirty_diff"],
                              -0.8841594,
                              tolerance = 1e-6)
  testthat::expect_equivalent(feat.df["1-Methylxanthine", "anova_ageGroup_fifty.thirty_BH"],
                              0.01544621,
                              tolerance = 1e-6)
  
  feat.df <- .anovas(data.mn = t(assay(sacurine.se)),
                     samp.df = colData(sacurine.se),
                     feat.df = rowData(sacurine.se),
                     test.c = "anova",
                     factorNameC = "ageGroup",
                     factorLevelsVc = c("thirty", "fourty", "fifty"),
                     adjust.c = "BH",
                     adjust_thresh.n = 0.05,
                     prefix.c = "",
                     figure.c = "none")[["feat.df"]]
  testthat::expect_equivalent(feat.df["1-Methylxanthine", "anova_ageGroup_thirty.fourty_diff"],
                              0.8841594,
                              tolerance = 1e-6)
  testthat::expect_equivalent(feat.df["1-Methylxanthine", "anova_ageGroup_thirty.fourty_BH"],
                              0.01602643,
                              tolerance = 1e-6)
  testthat::expect_equivalent(feat.df["Testosterone glucuronide", "anova_ageGroup_thirty.fifty_diff"],
                              -2.314257,
                              tolerance = 1e-6)
  testthat::expect_equivalent(feat.df["Testosterone glucuronide", "anova_ageGroup_thirty.fifty_BH"],
                              9.217707e-05,
                              tolerance = 1e-6)
})

#### kruskal ####

testthat::test_that("kruskal", {
  
  sacurine.se <- reading(system.file("extdata/sacurine",
                                     package = "phenomis"))
  
  sacurine.se <- correcting(sacurine.se, figure.c = "none")
  sacurine.se <- sacurine.se[, colData(sacurine.se)[, "sampleType"] != "pool"]
  sacurine.se <- transforming(sacurine.se)
  sacurine.se <- sacurine.se[, colnames(sacurine.se) != "HU_neg_096_b2"]
  colData(sacurine.se)[, "ageGroup"] <- vapply(colData(sacurine.se)[, "age"],
                                               function(x) {
                                                 if (x < 35) {
                                                   return("thirty")
                                                 } else if (x < 50) {
                                                   return("fourty")
                                                 } else {
                                                   return("fifty")}},
                                               FUN.VALUE = character(1))
  
  feat.df <- .anovas(data.mn = t(assay(sacurine.se)),
                     samp.df = colData(sacurine.se),
                     feat.df = rowData(sacurine.se),
                     test.c = "kruskal",
                     factorNameC = "ageGroup",
                     factorLevelsVc = "default",
                     adjust.c = "BH",
                     adjust_thresh.n = 0.05,
                     prefix.c = "testthat_",
                     figure.c = "none")[["feat.df"]]
  testthat::expect_equivalent(feat.df["1-Methylxanthine", "testthat_kruskal_ageGroup_fourty.thirty_BH"],
                              0.03497385,
                              tolerance = 1e-6)
  
})

#### anova2ways

testthat::test_that("anova2ways", {
  
  metabo.se <- reading(system.file("extdata/prometis/metabo",
                                   package = "phenomis"))
  
  ## .anova2ways
  anova2ways.ls <- .anovas2ways(data.mn = t(assay(metabo.se)),
                                samp.df = colData(metabo.se),
                                feat.df = rowData(metabo.se),
                                test.c = "anova2ways",
                                factor_names.vc = c("gene", "sex"),
                                factor_levels.ls = list(factor1Vc = c("WT", "LAT"),
                                                        factor2Vc = c("M", "F")),
                                adjust.c = "BH",
                                adjust_thresh.n = 0.05,
                                prefix.c = "prefix_",
                                figure.c = "none")
  
  testthat::expect_equivalent(anova2ways.ls[["metric.mn"]]["guanosine", "prefix_anova2ways_sex_M.F_diff"],
                              -0.8415528,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(anova2ways.ls[["feat.df"]]["isobutyric acid", "prefix_anova2ways_gene_WT.LAT_BH"],
                              0.5926903,
                              tolerance = 1e-6)
  
  anova2ways.ls <- .anovas2ways(data.mn = t(assay(metabo.se)),
                                samp.df = colData(metabo.se),
                                feat.df = rowData(metabo.se),
                                test.c = "anova2ways",
                                factor_names.vc = c("gene", "sex"),
                                factor_levels.ls = list(factor1 = "default",
                                                        factor2 = "default"),
                                adjust.c = "BH",
                                adjust_thresh.n = 0.05,
                                prefix.c = "",
                                figure.c = "none")
  
  testthat::expect_equivalent(anova2ways.ls[["metric.mn"]]["L-methionine", "anova2ways_sex_F.M_diff"],
                              -0.2421317773,
                              tolerance = 1e-6)
  
  testthat::expect_error(.anovas2ways(data.mn = t(assay(metabo.se)),
                                      samp.df = colData(metabo.se),
                                      feat.df = rowData(metabo.se),
                                      test.c = "anova2ways",
                                      factor_names.vc = c("gene", "sex"),
                                      factor_levels.ls = list(factor1 = c("WT", "LAT2"),
                                                              factor2 = c("M", "F")),
                                      adjust.c = "BH",
                                      adjust_thresh.n = 0.05,
                                      prefix.c = "",
                                      figure.c = "none"))
  
  metabo.se <- hypotesting(metabo.se,
                           test.c = "anova2ways",
                           factor_names.vc = c("gene", "sex"),
                           figure.c = "interactive",
                           report.c = "none")
  
  testthat::expect_equivalent(rowData(metabo.se)["citric acid", "anova2ways_sex_F.M_diff"],
                              0.630733,
                              tolerance = 1e-6)
  
})

#### anova2waysInter ####

testthat::test_that("anova2waysInter", {
  
  metabo.se <- reading(system.file("extdata/prometis/metabo",
                                   package = "phenomis"))
  
  ## .anova2waysInter
  anova2waysInterLs <- .anovas2ways(data.mn = t(assay(metabo.se)),
                                    samp.df = colData(metabo.se),
                                    feat.df = rowData(metabo.se),
                                    test.c = "anova2waysInter",
                                    factor_names.vc = c("gene", "sex"),
                                    factor_levels.ls = list(factor1Vc = c("WT", "LAT"),
                                                            factor2Vc = c("M", "F")),
                                    adjust.c = "BH",
                                    adjust_thresh.n = 0.05,
                                    prefix.c = "",
                                    figure.c = "none")
  
  testthat::expect_equivalent(anova2waysInterLs[["metric.mn"]]["adenosine 5'-monophosphate", "anova2waysInter_sex_M.F_diff"],
                              -1.1272963726,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(anova2waysInterLs[["feat.df"]]["L-serine", "anova2waysInter_gene_WT.LAT_BH"],
                              5.17659e-05,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(anova2waysInterLs[["feat.df"]]["isovaleric acid", "anova2waysInter_gene:sex_BH"],
                              0.763541,
                              tolerance = 1e-6)
  
  anova2waysInterLs <- .anovas2ways(data.mn = t(assay(metabo.se)),
                                    samp.df = colData(metabo.se),
                                    feat.df = rowData(metabo.se),
                                    test.c = "anova2waysInter",
                                    factor_names.vc = c("gene", "sex"),
                                    factor_levels.ls = list(factor1 = "default",
                                                            factor2 = "default"),
                                    adjust.c = "BH",
                                    adjust_thresh.n = 0.05,
                                    prefix.c = "",
                                    figure.c = "none")
  
  testthat::expect_equivalent(anova2waysInterLs[["metric.mn"]]["ADP", "anova2waysInter_sex_F.M_diff"],
                              0.5175954271,
                              tolerance = 1e-6)
  
  testthat::expect_error(.anovas2ways(data.mn = t(assay(metabo.se)),
                                      samp.df = colData(metabo.se),
                                      feat.df = rowData(metabo.se),
                                      test.c = "anova2waysInter",
                                      factor_names.vc = c("gene", "sex"),
                                      factor_levels.ls = list(factor1 = c("WT", "LAT2"),
                                                              factor2 = c("M", "F")),
                                      adjust.c = "BH",
                                      adjust_thresh.n = 0.05,
                                      prefix.c = "",
                                      figure.c = "none"))
  
  metabo.se <- hypotesting(metabo.se,
                           test.c = "anova2waysInter",
                           factor_names.vc = c("gene", "sex"),
                           figure.c = "interactive",
                           report.c = "none")
  
  testthat::expect_equivalent(rowData(metabo.se)["UDP-alpha-D-galactose", "anova2waysInter_sex_F.M_diff"],
                              0.892837,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(metabo.se)["glycocholate", "anova2waysInter_sex_F.M_BH"],
                              0.0120354,
                              tolerance = 1e-6)
  
})

#### limma2ways ####

testthat::test_that("limma2ways", {
  
  metabo.se <- reading(system.file("extdata/prometis/metabo",
                                   package = "phenomis"))
  
  ## .limma2ways
  limma2waysLs <- .anovas2ways(data.mn = t(assay(metabo.se)),
                               samp.df = colData(metabo.se),
                               feat.df = rowData(metabo.se),
                               test.c = "limma2ways",
                               factor_names.vc = c("gene", "sex"),
                               factor_levels.ls = list(factor1Vc = c("WT", "LAT"),
                                                       factor2Vc = c("M", "F")),
                               adjust.c = "BH",
                               adjust_thresh.n = 0.05,
                               prefix.c = "",
                               figure.c = "none")
  
  testthat::expect_equivalent(limma2waysLs[["metric.mn"]]["L-threonine", "limma2ways_sex_M.F_diff"],
                              0.1144511178,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(limma2waysLs[["feat.df"]]["ADP-D-ribose", "limma2ways_gene_WT.LAT_BH"],
                              5.50691e-02,
                              tolerance = 1e-6)
  
    limma2waysLs <- .anovas2ways(data.mn = t(assay(metabo.se)),
                               samp.df = colData(metabo.se),
                               feat.df = rowData(metabo.se),
                               test.c = "limma2ways",
                               factor_names.vc = c("gene", "sex"),
                               factor_levels.ls = list(factor1 = "default",
                                                       factor2 = "default"),
                               adjust.c = "BH",
                               adjust_thresh.n = 0.05,
                               prefix.c = "testing_",
                               figure.c = "none")
  
  testthat::expect_equivalent(limma2waysLs[["metric.mn"]]["Tyr-Phe-Met", "testing_limma2ways_sex_F.M_diff"],
                              0.6177232690,
                              tolerance = 1e-6)
  
  testthat::expect_error(.anovas2ways(data.mn = t(assay(metabo.se)),
                                      samp.df = colData(metabo.se),
                                      feat.df = rowData(metabo.se),
                                      test.c = "limma2ways",
                                      factor_names.vc = c("gene", "sex"),
                                      factor_levels.ls = list(factor1 = c("WT", "LAT2"),
                                                              factor2 = c("M", "F")),
                                      adjust.c = "BH",
                                      adjust_thresh.n = 0.05,
                                      prefix.c = "",
                                      figure.c = "none"))
  
  metabo.se <- hypotesting(metabo.se,
                           test.c = "limma2ways",
                           factor_names.vc = c("gene", "sex"),
                           figure.c = "interactive",
                           report.c = "none")
  
  testthat::expect_equivalent(rowData(metabo.se)["GDP-L-fucose", "limma2ways_sex_F.M_diff"],
                              0.176800,
                              tolerance = 1e-6)
  testthat::expect_equivalent(rowData(metabo.se)["ADP-D-ribose", "limma2ways_sex_F.M_BH"],
                              0.1906116,
                              tolerance = 1e-6)
  
})

#### limma2waysInter ####

testthat::test_that("limma2waysInter", {
  
  metabo.se <- reading(system.file("extdata/prometis/metabo",
                                   package = "phenomis"))
  
  ## .limma2waysInter
  limma2waysInterLs <- .anovas2ways(data.mn = t(assay(metabo.se)),
                                    samp.df = colData(metabo.se),
                                    feat.df = rowData(metabo.se),
                                    test.c = "limma2waysInter",
                                    factor_names.vc = c("gene", "sex"),
                                    factor_levels.ls = list(factor1Vc = c("WT", "LAT"),
                                                            factor2Vc = c("M", "F")),
                                    adjust.c = "BH",
                                    adjust_thresh.n = 0.05,
                                    prefix.c = "",
                                    figure.c = "none")
  
  testthat::expect_equivalent(limma2waysInterLs[["metric.mn"]]["pantothenic acid", "limma2waysInter_gene_WT.LAT_diff"],
                              -0.93364925,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(limma2waysInterLs[["feat.df"]]["GDP", "limma2waysInter_gene_WT.LAT_BH"],
                              8.53441e-04,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(limma2waysInterLs[["feat.df"]]["GDP-L-fucose", "limma2waysInter_gene:sex_BH"],
                              0.764520,
                              tolerance = 1e-6)
  
  limma2waysInterLs <- .anovas2ways(data.mn = t(assay(metabo.se)),
                                    samp.df = colData(metabo.se),
                                    feat.df = rowData(metabo.se),
                                    test.c = "limma2waysInter",
                                    factor_names.vc = c("gene", "sex"),
                                    factor_levels.ls = list(factor1 = "default",
                                                            factor2 = "default"),
                                    adjust.c = "BH",
                                    adjust_thresh.n = 0.05,
                                    prefix.c = "",
                                    figure.c = "none")
  
  testthat::expect_equivalent(limma2waysInterLs[["metric.mn"]]["9H-xanthine", "limma2waysInter_gene_LAT.WT_diff"],
                              0.47207710,
                              tolerance = 1e-6)
  
  testthat::expect_error(.anovas2ways(data.mn = t(assay(metabo.se)),
                                               samp.df = colData(metabo.se),
                                               feat.df = rowData(metabo.se),
                                               test.c = "limma2waysInter",
                                               factor_names.vc = c("gene", "sex"),
                                               factor_levels.ls = list(factor1 = c("WT", "LAT2"),
                                                                       factor2 = c("M", "F")),
                                               adjust.c = "BH",
                                               adjust_thresh.n = 0.05,
                                               prefix.c = "",
                                               figure.c = "none"))
  
  metabo.se <- hypotesting(metabo.se,
                           test.c = "limma2waysInter",
                           factor_names.vc = c("gene", "sex"),
                           figure.c = "interactive",
                           report.c = "none")
  
  testthat::expect_equivalent(rowData(metabo.se)["ADP-D-ribose", "limma2waysInter_gene_LAT.WT_diff"],
                              -1.25982,
                              tolerance = 1e-5)
  testthat::expect_equivalent(rowData(metabo.se)["isobutyric acid", "limma2waysInter_gene_LAT.WT_BH"],
                              0.577243199,
                              tolerance = 1e-6)
  
})
