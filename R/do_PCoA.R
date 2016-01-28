#' @title Perform PCoA on a distance matrix
#' @export
#'
#' @description
#'
#' This is a convenience wrapper for the ‘cmdscale’ function from the vegan
#' package. It takes a distance matrix (an object of type "dist") and performs a
#' PCoA ordination with ‘cmdscale’, for two dimensions (k = 2). Returns a list
#' containing ‘Coordinates’, a tbl with the PCoA coordinates; and ‘Variance1’
#' and ‘Variance2’, the percentage variance explained by the first and second
#' PCoA axes respectively.
#'
#' The method for calculating variance explained comes from
#' \url{http://r-sig-ecology.471788.n2.nabble.com/Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-td6429547.html}.
#'
#' @param DistanceMatrix a dist object containing the distances between samples
#'
#' @seealso dist for info on the dist object
#' @seealso read_dist to read a distance matrix from a file (e.g. reading in
#' UniFrac distances)
#' @seealso cmdscale the vegan function that performs the actual ordination
do_PCoA <- function(DistanceMatrix) {

  # Return error if distance matrix is not a distance matrix
  if (class(DistanceMatrix) != "dist") {
    stop("Distance matrix must be an object of type dist. See ?read_dist for a
         way to read a distance matrix from a file into a dist object.")
  }
  
  # Run the PCoA
  PCoA <- cmdscale(DistanceMatrix, k = 2, eig = TRUE)

  # Extract PCoA coordinates
  Coordinates <- data.frame(
    Sample = row.names(PCoA$points),
    PCoA1 = PCoA$points[,1],
    PCoA2 = PCoA$points[,2],
    row.names = NULL
  ) %>%
    as.tbl()

  # Calculate variance explained
  # Method is from http://r-sig-ecology.471788.n2.nabble.com/\
  # Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-\
  # td6429547.html
  Eigenvalues <- eigenvals(PCoA) 
  Variance <- Eigenvalues / sum(Eigenvalues) 
  Variance1 <- 100 * signif(Variance[1], 2)
  Variance2 <- 100 * signif(Variance[2], 2)

  # Return
  Result <- list(
    Coordinates = Coordinates,
    Variance1 = Variance1,
    Variance2 = Variance2
  )
  return(Result)
}
