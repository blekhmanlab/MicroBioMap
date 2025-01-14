# Run before any test
testbfc <- BiocFileCache::BiocFileCache(cache = testthat::test_path())

# Remove the cache created during tests
withr::defer(BiocFileCache::cleanbfc(
  #testbfc,
  BiocFileCache::BiocFileCache(cache = testthat::test_path()),
  days=-Inf, ask=FALSE), envir=testthat::teardown_env(), priority='last')
