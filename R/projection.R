library(dplyr)
library(tidyr)
library(tibble)
library(vegan)

.gm_mean = function(x){
  exp(mean(log(x[x > 0])))
}

.rclr <- function(a) {
  answer <- log(a/.gm_mean(a))
  answer[] <- lapply(answer, function(i) if(is.numeric(i)) ifelse(is.infinite(i), 0, i) else i)
  return(answer)
}

.loadUserData <- function(userdata, usertaxa, seed=NA) {
    if(!is.na(seed)) set.seed(seed)

    userdata |>
        pivot_longer(!sample
            , names_to = 'col_id'
            , values_to = 'count'
        ) |> # one row for each sample/taxon count
        left_join(usertaxa, by='col_id') |> # add taxon names
        mutate(famlevel = paste(kingdom,phylum,class,order,family, sep='||')) |>
        select(!c('kingdom','phylum','class','order','family','genus')) |>
        group_by(sample, famlevel) |> # summarize at family level
        summarise(
            aggcount = sum(count)
        ) |>
        ungroup() |>
        pivot_wider(names_from='famlevel' # make taxon table
            , values_from='aggcount'
        ) |>
        column_to_rownames('sample') |>
        (\(x) {
            x[,colSums(x) > 0] # grab non-empty taxa
        })() |>
        rownames_to_column('sample') |>
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
        .rclr() |>
        rownames_to_column('sample')
}

.loadProjection <- function(pname, load_all) {
  load_all |>
    filter(projection==pname) |>
    mutate(taxon = paste(kingdom,phylum,class,order,family, sep='||')) # TODO: use numbers for this instead
}

project_it <- function(indata, loadings) {
    indata |>
        pivot_longer(cols=!c(sample)
            , names_to='taxon'
            , values_to='val'
            , values_drop_na=TRUE
        ) |>
        left_join(loadings, by='taxon') |>
        mutate(
          across(starts_with('PC'), \(d) val * d)
        ) |>
        ungroup() |>
        group_by(sample) |>
        summarise(
          across(starts_with('PC'), sum)
        )
}
