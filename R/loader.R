.getVersions <- function(bfc, verbose=FALSE) {
    # Determines the most recent version of the compendium
    # and retrieves the manifest that describes all available releases.
    # Returns a data.table listing all versions and the necessary URLs
    # This requires the canonical_doi configuration value stored in
    # constants.R, which always resolves to the most recent version.

    # Check if we've already cached the manifest.
    # If not, make an HTTP call to get the URL we need
    rpath <- BiocFileCache::bfcquery(bfc, 'manifest')$rpath
    if(length(rpath) == 0) {
      if(verbose) {
        print('Retrieving version information.')
      }
      resolve <- curl::curl_fetch_memory(canonical_doi)
      if(resolve$status_code != 200) {
        stop(paste0(
          'Could not resolve canonical DOI. Status code: ',
          resolve$status_code
        ))
      }

      if(verbose) {
        print('Determined data address:')
        print(resolve$url)
      }
      manifest <- paste0(resolve$url, '/files/manifest.csv')

      rpath <- tryCatch(
        {
          bfcrpath(bfc, manifest)
        },
        error = function(msg){
          print('Could not retrieve manifest file. Falling back to manifest as of v1.1.0')
          towrite <- data.table::data.table(
            version = c('1.1.0', '1.0.1'),
            zenodo_id = c('13733642', '10452633'),
            default = c(TRUE, FALSE)
          )
          # we save this to the cache so the app remembers not to keep looking online
          # for a manifest every time the version information is needed
          savepath <- BiocFileCache::bfcnew(bfc, 'manifest', ext='.csv')
          data.table::fwrite(towrite, file=savepath)
          savepath
        }
      )
    }
    else {
      if(verbose) {
        print('Cached version information found.')
      }
    }
    results <- data.table::fread(rpath)

    colnames(results) <- c('version','zenodo_id','default')
    results$data_url <- paste0('https://zenodo.org/record/', results$zenodo_id, '/files/taxonomic_table.csv.gz')
    results$coldata_url <- paste0('https://zenodo.org/record/', results$zenodo_id, '/files/sample_metadata.tsv')
    data.table::setkey(results, version)
    results
}

.getCompendiumData <- function(version, bfc) {
    versions <- .getVersions(bfc)
    rpath <- bfcrpath(bfc, versions[version]$data_url)
    data.table::fread(rpath)
}

.getCompendiumColdata <- function(version, bfc) {
    versions <- .getVersions(bfc)
    rpath <- bfcrpath(bfc, versions[version]$coldata_url)
    sampdat <- as.data.frame(data.table::fread(rpath))
    rownames(sampdat) <- paste(sampdat[[2]], sampdat[[3]], sep = "_")
    sampdat
}

#' load all compendium data into a TreeSummarizedExperiment
#'
#' @param bfc BiocFileCache object to use
#'
#' @returns a `TreeSummarizedExperiment`
#'
#' @importFrom data.table fread setkey
#' @importClassesFrom Matrix TsparseMatrix
#' @import TreeSummarizedExperiment
#' @import R.utils
#' @import ape
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcquery bfcnew
#'
#' @export
#'
#' @examples
#' cpd <- getCompendium()
#'
#' dim(cpd)
#' cpd
#' assayNames(cpd)
#' head(colData(cpd))
#'

getCompendium <- function(version=NA, bfc = BiocFileCache::BiocFileCache()) {
    versions <- .getVersions(bfc)

    if(is.na(version)) {
        # If the user has not specified a version, grab whichever
        # is indicated in the manifest as the default (i.e. most recent)
        version <- versions[versions$default,]$version[1]
    }
    print(paste('Retrieving compendium version',version))
    dat <-.getCompendiumData(version, bfc)
    coldat <- .getCompendiumColdata(version, bfc)

    sampnames <- dat[[2]]

    coldat <- coldat[match(sampnames, rownames(coldat)), ]

    taxa <- colnames(dat)[3:ncol(dat)]
    requireNamespace("Matrix")
    # mat = as(as.matrix(dat[,3:ncol(dat)]), 'TsparseMatrix')
    mat <- as.matrix(dat[, 3:ncol(dat)])
    rownames(mat) <- sampnames
    colnames(mat) <- taxa
    sampinfo <- do.call(rbind, strsplit(sampnames, "_"))
    colnames(sampinfo) <- c("project", "sample")
    coldata <- data.frame(sampinfo)
    rownames(coldata) <- sampnames
    splittaxa <- do.call(rbind, lapply(
        strsplit(taxa, "\\."),
        function(x) {
            c(x, rep(NA, 8 - length(x)))
        }
    ))
    colnames(splittaxa) <- c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain"
    )
    rowdata <- data.frame(splittaxa)
    rownames(rowdata) <- taxa
    td <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        colData = coldat,
        rowData = rowdata,
        assays = list(counts = t(mat))
    )
    td
}

taxonname2edgelist <- function(taxon) {
    v <- strsplit(taxon, "\\.")[[1]]
    v <- v[!v == "NA"]
    if (length(v) > 1) {
        lv <- length(v)
        df <- data.frame(from = v[seq_len(lv - 1)], to = v[1+seq_len(lv-1)])
    } else {
        df <- data.frame()
    }
    df
}

taxa2edgelist <- function(taxa) {
    taxa_edgelist <- lapply(taxa, taxonname2edgelist)
    df <- unique(do.call(rbind, taxa_edgelist))
    return(df)
    unique_names <- unique(c(df$from, df$to))
    l <- seq_along(unique_names)
    names(l) <- unique_names
    parents <- l[df$from]
    nodes <- l[df$to]
    df$parent <- parents
    df$node <- nodes
    df$label <- df$to
    df
}

taxa2phylo <- function(taxa) {
    edgelist <- taxa2edgelist(taxa)
    edgelist <- as.matrix(edgelist)

    edgelist <- edgelist[!is.na(edgelist[, 1]) & !is.na(edgelist[, 2]), ]

    from <- edgelist[, 1]
    to <- edgelist[, 2]
    ids <- unique(c(edgelist[, 1], edgelist[, 2]))

    tip.label <- setdiff(ids, from)
    node.label <- unique(from)

    # make a map from taxonomy ID to internal 1:n ids
    idmap <- seq_along(c(tip.label, node.label))
    names(idmap) <- c(tip.label, node.label)

    # make a phylo object
    tree <- list(
        edge       = matrix(c(idmap[as.character(from)], idmap[as.character(to)]), ncol = 2),
        tip.label  = unname(tip.label),
        # node.label = unname(node.label),
        Nnode      = length(node.label)
    )
    class(tree) <- "phylo"

    tree
}
