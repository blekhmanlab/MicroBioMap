test_that("getting compendium with defaults works as expected", {
    skip_on_ci()
    cpd <- getCompendium(bfc = testbfc)
    expect_s4_class(cpd, "TreeSummarizedExperiment")
    expect_gt(nrow(cpd), 1000)
    expect_gt(ncol(cpd), 168000)
    expect_contains(assayNames(cpd), "counts")
    # v1.1.0 changed the name of this taxon:
    expect_equal(max(counts(cpd)['Bacteria.Bacillota.Clostridia.Eubacteriales.Alkalibacteraceae.Alkalibaculum',]), 64)
    expect_error(max(counts(cpd)['Bacteria.Firmicutes.Clostridia.Eubacteriales.(unclassified).Alkalibaculum',]))
})

test_that("getting compendium with specified version works as expected", {
    skip_on_ci()
    cpd <- getCompendium('1.0.1', bfc = testbfc)
    expect_s4_class(cpd, "TreeSummarizedExperiment")
    expect_gt(nrow(cpd), 1000)
    expect_gt(ncol(cpd), 1000)
    expect_contains(assayNames(cpd), "counts")
    expect_equal(max(counts(cpd)['Bacteria.Firmicutes.Clostridia.Eubacteriales.(unclassified).Alkalibaculum',]), 16)
    expect_error(max(counts(cpd)['Bacteria.Firmicutes.Clostridia.Eubacteriales.Alkalibaculum.NA',]))
})
