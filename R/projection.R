#' @import dplyr
#' @importFrom vegan rrarefy
#' @importFrom tibble rownames_to_column column_to_rownames

#library(tidyr)
.gm_mean = function(x){
  exp(mean(log(x[x > 0])))
}

.rclr <- function(a) {
  answer <- log(a/.gm_mean(a))
  answer[] <- lapply(answer, function(i) if(is.numeric(i)) ifelse(is.infinite(i), 0, i) else i)
  answer
}

#' Projection retrieval and formatting
#'
#' Filters a data frame containing coordinates for
#' multiple projections and formats it for further use.
#' @param pname string indicating the name of the projection to use
#' @param load_all a matrix of weights (with a "projection" field on which to filter)
#' @returns A data frame of weights for a single projection,
#'              plus a "taxon" field containing a string unique to that taxon.
.loadProjection <- function(pname, load_all) {
  load_all |>
    filter(projection==pname) |>
    mutate(taxon = paste(kingdom,phylum,class,order,family, sep='||')) # TODO: use numbers for this instead
}

#' Projection library builder
#'
#' Utility for working with remotely hosted projection files.
#' @param pname string indicating the name of the projection to use
#' @param bfc BiocFileCache object to use
#' @returns A **function** that allows a user to load a projection using
#'            its name.
#' @examples
#' lib <- projection_library()
#' projected <- project_it(my_data, projection_library('full'))
#'
#' @export
projection_library <- function(version=NA, bfc = BiocFileCache::BiocFileCache()) {
  versions <- .getVersions(bfc, entry='projection')
  print(versions)

  if(is.na(version)) {
    # If the user has not specified a version, grab whichever
    # is indicated in the manifest as the default (i.e. most recent)
    version <- versions[versions$default,]$version[1]
  }
  print(paste('Retrieving projections version',version))
  load_all <- .getProjectData(version, 'projection', bfc)

  curried <- function(pname) {
    .loadProjection(pname, load_all)
  }
  curried
}

#' User data formatting
#'
#' @description
#' Formats user input data from two input files into a single data frame for further analysis.
#' Consolidates read counts at the family level.
#' @details
#' Taxonomic information is provided across two files: the matrix of read counts uses
#' id-based column names ("TAX1", "TAX2", etc). The second file associates those column
#' names with inferred taxonomic ranks.
#' This function modifies the count table by annotating the column names with the
#' taxonomic information.
#' @param userdata data frame countaining samples in rows and taxa in columns
#' @param usertaxa data frame associating the column names from userdata with kingdom, phylum, class,
#'                 order, and family classifications.
#' @returns A data frame with the same number of rows as userdata and one column for each distinct
#'            family in the usertaxa file that was observed at least once in the dataset.
#'
#' @export
loadUserData <- function(userdata, usertaxa) {
  test <- userdata |>
        pivot_longer(!sample
            , names_to = 'col_id'
            , values_to = 'countnum'
        ) |> # one row for each sample/taxon count
        filter(countnum > 0) |>
        left_join(usertaxa, by='col_id') |> # add taxon names
        mutate(famlevel = paste(kingdom,phylum,class,order,family, sep='||')) |>
        select(!c('kingdom','phylum','class','order','family')) |>
        group_by(sample, famlevel) |> # summarize at family level
        summarise(
            aggcount = sum(countnum)
        ) |>
        ungroup() |>
        pivot_wider(names_from='famlevel' # make taxon table
            , values_from='aggcount'
            , values_fill=0
        )
}

#' Sample rarefaction
#'
#' @description
#' Helper function for rarefying samples in a taxonomic table to a set read count.
#' @details
#' Throwing away data is generally not helpful, but this step is done to ensure
#' samples being projected into an existing latent space are of the same size as the samples
#' that were included in the original ordination.
#' @param userdata data frame countaining samples in rows and taxa in columns. (Generally here the output of `loadUserData()`)
#' @param level The desired number of reads per sample. **Samples with fewer reads than this will be filtered out.**
#' @returns A data frame of the same dimensions as userdata, but with rarefied read counts that result in
#'            all rows summing to `level`.
#'
#' @export
rarefy_dataset <- function(userdata, level=3000, seed=NA) {
  if(!is.na(seed)) set.seed(seed)

  userdata |>
    unique(by=!c('sample')) |> # ensure no duplicate entries
    rowwise() |>
    mutate(
      totalreads = sum(c_across(where(is.numeric)))
    ) |>
    filter(
      totalreads >= 3000 # make sure each has enough reads
    ) |> # TODO log the results of this
    select(!totalreads) |>
    column_to_rownames('sample') |>
    rrarefy(3000) |>
    as.data.frame() |>
    rownames_to_column('sample')
}

#' Robust centered log-ratio transformation
#'
#' @description
#' Helper function for applying rCLR to a data frame with a "sample" column.
#' @details
#' This step should happen *after* rarefaction, and converts a matrix of non-negative read counts into
#' one containing unbounded floating-point numbers.
#' @param userdata data frame countaining samples in rows and taxa in columns.
#' @returns A data frame of the same dimensions as userdata, but with rCLR-transformed values
#'
#' @export
transform <- function(userdata) {
  userdata |>
    column_to_rownames('sample') |>
    .rclr() |>
    rownames_to_column('sample')
}

#' Projection operation
#'
#' @description
#' Projects a taxonomic table into an existing ordination using a set of pre-calculated loadings.
#' @details
#' This is the final step in the projection process. The input describes samples with taxa, and the output
#' describes the same taxa with principal components. Note that this step doesn't require rCLR-transformed
#' values *per se* -- input data should be transformed using the same process as the data used in the original
#' ordination. For the compendium ordinations currently available, this means rCLR.
#' @param indata data frame countaining samples in rows and taxa in columns.
#' @param loadings data frame containing taxa in each row (i.e. the indata column names) and
#'                  a principal component ("PC1", "PC2", etc) in each column. The values in each
#'                  cell indicate the weight that taxon's read count should be given when calculating
#'                   a sample's value for that PC.
#' @returns A data frame of the same rows as data and the same columns as loadings
#'
#' @export
project_it <- function(indata, loadings) {
    indata |>
        pivot_longer(cols=!c(sample)
            , names_to='taxon'
            , values_to='val'
            , values_drop_na=TRUE
        ) |>
        left_join(loadings, by='taxon') |>
        mutate(
          across(starts_with('PC')
          , \(d) val * d)
        ) |>
        ungroup() |>
        group_by(sample) |>
        summarise(
          across(starts_with('PC'), sum)
        )
}
