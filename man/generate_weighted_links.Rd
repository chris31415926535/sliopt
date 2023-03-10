% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{generate_weighted_links}
\alias{generate_weighted_links}
\title{Generate one-to-many weighted links between sets of spatial regions}
\usage{
generate_weighted_links(from_shp, from_idcol, to_shp, to_idcol)
}
\arguments{
\item{from_shp}{An `sf` object with the geometries we are linking from.}

\item{from_idcol}{A length-one character vector with the name of the `from_shp` column with unique region identifiers.}

\item{to_shp}{An `sf` object with the geometries we are linking to.}

\item{to_idcol}{A length-one character vector with the name of the `to_shp` column with unique region identifiers.}
}
\value{
A `tbl_df` with three columns, showing the linked-from region ids, the linked-to region ids, and the weights for each link. Each from-id set of weights should sum to 1.
}
\description{
Generate one-to-many weighted links between sets of spatial regions
}
\examples{
 \dontrun{
 generate_weighted_links(da_ott, "DAUID", ons_shp, "ONS_ID")
 }
}
