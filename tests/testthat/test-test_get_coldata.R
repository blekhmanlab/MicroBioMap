test_that("sample metadata download works", {
    coldat <- .getCompendiumColdata('1.1.0', testbfc)
    expect_equal(ncol(coldat), 11)
    expect_contains(colnames(coldat), c(
        "srs", "project", "srr", "library_strategy",
        "library_source", "pubdate", "total_bases",
        "instrument", "geo_loc_name",
        "region"
    ))
})
