read_wbd <- function(hu = "2", wbd_dir = "rawdata/geodata/WBD_17_HU2_Shape") {
    fn <- paste0("WBDHU", hu, ".shp")
    path <- file.path(wbd_dir, "Shape", fn)
    read_sf(path)
}

check_nhd_feature <- function(feature) {
    nhd_features <- c(
            "Area",
            "AreaEventFC",
            "Flowline",
            "Line",
            "LineEventFC",
            "Point",
            "PointEventFC",
            "Waterbody")
    if (!(feature %in% nhd_features)) {
            stop(paste("feature must be one of",
             paste(nhd_features, sep = ", ")))
    }
    invisible(TRUE)
}

read_nhd <- function(feature, nhd_dir = "rawdata/geodata/NHD_H_17100104_HU8_Shape") {
    check_nhd_feature(feature)
    fn <- paste0("NHD", feature, ".shp")
    path <- file.path(nhd_dir, "Shape", fn)
    read_sf(path)
}

#' Function to iteratively add rows from `bigset` that intersect with `init`.
#' Careful to remove e.g. some downstream segment if you don't want to add those
#' as well.
#'
#' Iteratively select rows that intersect with existing rows.
#'
#' @param init spatial feather to build from
#' @param bigset set of all spatial features to iteratively add to a spatial
#'  data frame using \code{\link[sf]{st_intersect}}
#'
#' @return a data frame that includes all rows that intersect `init` directly or
#'  indirectly.
#' @export
iter_intersect <- function(init, bigset) {
  n_row <- c(0, 1)
  sfdf <- init
  while (diff(n_row) > 0) {
    sfdf <- st_filter(bigset, sfdf, .predicate = st_intersects)
    n_row[1] <- n_row[2]
    n_row[2] <- nrow(sfdf)
  }
  sfdf
}
