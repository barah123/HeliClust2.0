

#' Internal helper: coerce all columns to numeric
#'
#' @param x Data frame or matrix.
#' @keywords internal
heliclust_as_numeric <- function(x) {
  x <- as.data.frame(x)
  x[] <- lapply(x, function(col) as.numeric(as.character(col)))
  x
}

#' Compute distance matrix for HeliClust2.0
#'
#' This function scales the data (optional) and computes a distance matrix
#' suitable for hierarchical clustering.
#'
#' @param x Numeric data frame or matrix (rows = samples, columns = variables).
#' @param dist_method Character; distance method passed to \code{vegan::vegdist}
#'   (e.g. "euclidean", "bray").
#' @param scale_data Logical; if TRUE, standardizes variables using \code{scale()}.
#'
#' @return A numeric matrix representing pairwise distances.
#' @export
#' @examples
#' \dontrun{
#'   dist_mat <- heliclust_dist(my_data)
#' }
heliclust_dist <- function(x,
                           dist_method = "euclidean",
                           scale_data = TRUE) {

  x <- heliclust_as_numeric(x)

  if (scale_data) {
    x <- scale(x)
  }

  dist_obj <- vegan::vegdist(x, method = dist_method, na.rm = TRUE)
  as.matrix(dist_obj)
}

#' Hierarchical clustering for metabolite data
#'
#' Performs hierarchical clustering on metabolite / biomarker data.
#'
#' @inheritParams heliclust_dist
#' @param hclust_method Linkage method for \code{stats::hclust}
#'   (e.g. "ward.D2", "average", "complete").
#'
#' @return An object of class \code{hclust}.
#' @export
#' @examples
#' \dontrun{
#'   hc <- heliclust_metabolites(my_metabolite_data)
#'   plot(hc)
#' }
heliclust_metabolites <- function(x,
                                  dist_method = "euclidean",
                                  hclust_method = "ward.D2",
                                  scale_data = TRUE) {

  dist_mat <- heliclust_dist(x,
                             dist_method = dist_method,
                             scale_data = scale_data)

  stats::hclust(stats::as.dist(dist_mat), method = hclust_method)
}

#' Hierarchical clustering for sequence / OTU data
#'
#' Wrapper around \code{heliclust_dist()} and \code{stats::hclust()} for
#' sequence data (e.g. OTU/ASV tables).
#'
#' @inheritParams heliclust_metabolites
#'
#' @return An object of class \code{hclust}.
#' @export
#' @examples
#' \dontrun{
#'   hc_seq <- heliclust_sequences(my_otu_table,
#'                                 dist_method = "bray",
#'                                 hclust_method = "average",
#'                                 scale_data = FALSE)
#'   plot(hc_seq)
#' }
heliclust_sequences <- function(x,
                                dist_method = "bray",
                                hclust_method = "average",
                                scale_data = FALSE) {

  dist_mat <- heliclust_dist(x,
                             dist_method = dist_method,
                             scale_data = scale_data)

  stats::hclust(stats::as.dist(dist_mat), method = hclust_method)
}
