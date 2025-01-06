test_that("sample metadata download works", {
    bfc <- BiocFileCache::BiocFileCache()
    coldat <- .getCompendiumColdata('1.1.0', bfc)
    expect_equal(ncol(coldat), 11)
    expect_contains(colnames(coldat), c(
        "srs", "project", "srr", "library_strategy",
        "library_source", "pubdate", "total_bases",
        "instrument", "geo_loc_name",
        "region"
    ))
})
