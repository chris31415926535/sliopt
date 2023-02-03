
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sliopt: Single Link Indicator OPTimization for geospatial analysis

<!-- badges: start -->
<!-- badges: end -->

Sliopt provides an algorithm to create optimized single-link indicators
from one set of (non-self-overlapping) geometries to a different
overlapping set of (non-self-overlapping) geometries. It requires you to
have one “source of truth” for each set of geometries.

Sometimes you have data measured over one set of geometries, like census
tracts, but you want to work with a different larger set of geometries,
like school catchment areas. One approach is to map each smaller
geometry to one and only one larger geometry using a “single link
indicator,” or “SLI.” This will allow you to aggregate the values from
the smaller regions up to the bigger regions.

If the two sets of geometries nest into each other perfectly, then it’s
trivial to create an SLI. But if the two geometry sets have borders that
don’t align, there is no single way to map each member of the first set
to a member of the second set.

This package assumes that you are able to get at least one set of values
that is reliable across all geometries. It then uses that set of “ground
truths” to create an optimized SLI linking the two sets of geometries,
minimizing a user-specified error function (mean absolute error, mean
absolute percentage error, or mean squared error). This SLI can then be
used to estimate distributions for other values which are only known at
one of the geometry sets.

## Installation

You can install the development version of sliopt like so:

``` r
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
library(sf)
#> Linking to GEOS 3.8.0, GDAL 3.0.4, PROJ 6.3.1; sf_use_s2() is TRUE
library(sliopt)


layer1b <- list(list(matrix(c(0,0, 0,1, 1,1, 1,0, 0,0 ), byrow = TRUE, ncol = 2)),
                list(matrix(c(0,0, 0,1, -1,1, -1,0, 0,0 ), byrow = TRUE, ncol = 2))
)


layer2b <- list(list( matrix(c( -1,1, -1,.4, -.6,.4, -.6, 1, -1,1 ), byrow = TRUE, ncol = 2)),
                list( matrix(c( -1,0, .4, 0,  .4,.4,  -1,.4, -1,0 ), byrow = TRUE, ncol = 2)),
                list( matrix(c( .4,0,  1,0,   1, 1,  .4, 1, .4, 0), byrow = TRUE, ncol = 2)),
                list( matrix(c(-.6,1, -.6,.4 ,  .4, .4, .4, 1, -.6, 1), byrow = TRUE, ncol = 2))
)



layer1_sf <- lapply(layer1b, sf::st_polygon) %>%
  sf::st_sfc(crs="WGS84") %>%
  sf::st_sf() %>%
  dplyr::mutate(layer1_id = paste0("A", 1:nrow(.)))

layer2_sf <- lapply(layer2b, sf::st_polygon) %>%
  sf::st_sfc(crs="WGS84") %>%
  sf::st_sf() %>%
  dplyr::mutate(layer2_id = paste0("B", 1:nrow(.)))



# create two population clusters
pop_clusters <- create_pop_cluster(x=0.25, y=0.25, dx = 0.15, dy = 0.15, n=100, crs = "WGS84") %>%
  dplyr::bind_rows(create_pop_cluster(x=-0.75, y=0.75, dx = 0.25, dy = 0.25, n=100, crs = "WGS84"))

# create a plot
ggplot() +
  geom_sf(data = layer1_sf, colour = "red", size = 2) +
  geom_sf(data = layer2_sf, fill="lightblue", alpha = 0.6) +
  geom_sf(data = pop_clusters) +
  labs(title = "Point values over two overlapping sets of geometries")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r


layer1_count  <- sf::st_join(pop_clusters, layer1_sf) %>%
  dplyr::group_by(layer1_id) %>%
  dplyr::count(name = "layer1_n") %>%
  dplyr::filter(!is.na(layer1_id)) %>%
  sf::st_set_geometry(NULL) %>%
  dplyr::left_join(layer1_sf, ., by = "layer1_id")


layer2_count  <- sf::st_join(pop_clusters, layer2_sf) %>% 
  dplyr::group_by(layer2_id) %>%
  dplyr::count(name = "layer2_n") %>%
  dplyr::filter(!is.na(layer2_id)) %>%
  sf::st_set_geometry(NULL) %>%
  dplyr::left_join(layer2_sf, ., by = "layer2_id")

sli <- greedy_sli_search(to_shp = layer1_count, to_idcol = "layer1_id", to_valuecol = "layer1_n",
                  from_shp = layer2_count, from_idcol = "layer2_id", from_valuecol = "layer2_n", verbose = FALSE )
#> Precomputing intersections...
#> 1/2: Destination Region A2 1/2: Destination Region A2: Best mape: 0.07212
#> 1/2: Destination Region A2: Best mape: 0.07212 2/2: Destination Region A1 2/2:
#> Destination Region A1: Best mape: 0.07212 2/2: Destination Region A1: Best
#> mape: 0.07212 1/2: Destination Region A2 1/2: Destination Region A2: Best mape:
#> 0.07212 1/2: Destination Region A2: Best mape: 0.07212 2/2: Destination Region
#> A1 2/2: Destination Region A1: Best mape: 0.07212 2/2: Destination Region A1:
#> Best mape: 0.07212 1/2: Destination Region A1 1/2: Destination Region A1: Best
#> mape: 0.07212 1/2: Destination Region A1: Best mape: 0.07212 2/2: Destination
#> Region A2 2/2: Destination Region A2: Best mape: 0.07212 2/2: Destination Region
#> A2: Best mape: 0.07212

layer2_sf %>%
  left_join(sli, by = "layer2_id") %>%
  ggplot() +
  geom_sf(aes(fill = layer1_id)) +
  labs(title="Geometries shaded by SLI values")
```

<img src="man/figures/README-example-2.png" width="100%" />