# library(tidyverse)
# library(sf)
#
#
# # creating test data for the optimization algorithm... not sure it's going so well
#
#
# # create a point and we'll only take hoods within 10km
# ottawa_downtown <- dplyr::tibble(lon=-75.75, lat=45.35) %>%
#     sf::st_as_sf(coords = c("lon", "lat"), crs = "WGS84") %>%
#   sf::st_transform(crs=32189)
#
#
# ons_shp <- neighbourhoodstudy::ons_shp_gen3 %>%
#   # dplyr::mutate(bbox = purrr::map(geometry, sf::st_bbox)) %>%
#   # tidyr::unnest(bbox)
# #  dplyr::filter(bbox["xmin"] > -76)
#   sf::st_transform(crs = 32189) %>%
#   dplyr::filter(ONS_Region == "OTTAWA") %>%
#   dplyr::mutate(test = sf::st_is_within_distance(geometry, ottawa_downtown, 5000) %>% as.numeric()) %>%
#   dplyr::filter(test == 1) %>%
#   dplyr::filter(!stringr::str_detect(ONS_Name, "GREENBELT"))
#
# #ggplot(z) + geom_sf()
#
# dbs_shp <- neighbourhoodstudy::ottawa_shp_dbs %>%
#   sf::st_transform(crs = 32189) %>%
#   dplyr::mutate(test = sf::st_is_within_distance(geometry, ottawa_downtown, 7500) %>% as.numeric()) %>%
#   dplyr::filter(test == 1)
#
# synthetic_population <- create_synthetic_population(ons_shp, cluster_pop = 1000, num_clusters = 100, cluster_size = 400)
#
# ggplot() +
#   geom_sf(data = ons_shp) +
#   geom_sf(data = synthetic_population)
#   #geom_density_2d(data = clusters, aes)
#
# ons_shp_counted <- sf::st_join(y = synthetic_population, x=ons_shp) %>%
#   sf::st_set_geometry(NULL) %>%
#   dplyr::group_by(ONS_Name) %>%
#   dplyr::count() %>% {
#     dplyr::left_join(ons_shp, ., by = "ONS_Name")
#   }
#
# tictoc::tic()
# dbs_shp_counted <- sf::st_join(x=dbs_shp, y=synthetic_population) %>%
#   sf::st_set_geometry(NULL) %>%
#   dplyr::group_by(DBUID) %>%
#   dplyr::count() %>% {
#     dplyr::left_join(dbs_shp, ., by = "DBUID")
#   }
# tictoc::toc()
#
# dbs_shp_counted <- sf::st_make_valid(dbs_shp_counted)
#
# ons_shp_counted <- sf::st_make_valid(ons_shp_counted)
#
# ons_shp_counted %>%
#   ggplot() +
#   geom_sf(aes(fill = n)) +
#   ggplot2::scale_fill_viridis_c()
#
# dbs_shp_counted %>%
#   ggplot() +
#   geom_sf(aes(fill = n)) +
#   ggplot2::scale_fill_viridis_c()
#
#
#
#
# test <- greedy_sli_search3 (from_shp = dbs_shp_counted, from_idcol = "DBUID", from_valuecol = "n", to_shp = ons_shp_counted, to_idcol = "ONS_ID", to_valuecol = "n", optimize_for = "mae", tolerance = 0.05, iterations = 1, input_sli = NA, shuffle_inputs = TRUE, verbose = FALSE)
# #from_shp = dbs_shp_counted; from_idcol = "DBUID"; from_valuecol = "n"; to_shp = ons_shp_counted; to_idcol = "ONS_Name"; to_valuecol = "n"; optimize_for = "mae"; tolerance = 0.05; iterations = 3; input_sli = NA; shuffle_inputs = TRUE; verbose = FALSE
#
# test
#
# plot_sli_progress(test)
#
# output <- measure_sli_performance4(sli_for_test = test, from_values = dbs_shp_counted, from_idcol = "DBUID", from_valuecol = "n", to_values = ons_shp_counted, to_idcol = "ONS_ID", to_valuecol = "n")
#
# output
#
#
#
#
# ######## TESTING WITH IMAGINARY SHAPES
# library(tidyverse)
# library(sf)
#
#
# layer1 <- list(list(matrix(c(0,0, 0,1, 1,1, 1,0, 0,0 ), byrow = TRUE, ncol = 2),
#                     matrix(c(0,0, 0,1, -1,1, -1,0, 0,0 ), byrow = TRUE, ncol = 2)
#                   )
#                )
#
# layer1a <- list(matrix(c(0,0, 0,1, 1,1, 1,0, 0,0 ), byrow = TRUE, ncol = 2),
#                     matrix(c(0,0, 0,1, -1,1, -1,0, 0,0 ), byrow = TRUE, ncol = 2)
# )
#
# layer1b <- list(list(matrix(c(0,0, 0,1, 1,1, 1,0, 0,0 ), byrow = TRUE, ncol = 2)),
#                 list(matrix(c(0,0, 0,1, -1,1, -1,0, 0,0 ), byrow = TRUE, ncol = 2))
# )
#
# layer2 <- list(list( matrix(c( -1,1, -1,.4, -.6,.4, -.6, 1, -1,1 ), byrow = TRUE, ncol = 2)
#                     ,matrix(c( -1,0, .4, 0,  .4,.4,  -1,.4, -1,0 ), byrow = TRUE, ncol = 2)
#                     ,matrix(c( .4,0,  1,0,   1, 1,  .4, 1, .4, 0), byrow = TRUE, ncol = 2)
#                     ,matrix(c(-.6,1, -.6,.4 ,  .4, .4, .4, 1, -.6, 1), byrow = TRUE, ncol = 2)
#                     ))
#
# layer2b <- list(list( matrix(c( -1,1, -1,.4, -.6,.4, -.6, 1, -1,1 ), byrow = TRUE, ncol = 2)),
#                 list( matrix(c( -1,0, .4, 0,  .4,.4,  -1,.4, -1,0 ), byrow = TRUE, ncol = 2)),
#                 list( matrix(c( .4,0,  1,0,   1, 1,  .4, 1, .4, 0), byrow = TRUE, ncol = 2)),
#                 list( matrix(c(-.6,1, -.6,.4 ,  .4, .4, .4, 1, -.6, 1), byrow = TRUE, ncol = 2))
# )
#
# layer1_sf <- sf::st_multipolygon(layer1) %>%
#   sf::st_sfc(crs = "WGS84")
#
# layer2_sf <- sf::st_multipolygon(layer2) %>%
#   sf::st_sfc(crs = "WGS84")
#
# layer1_sf <- lapply(layer1b, sf::st_polygon) %>%
#   sf::st_sfc(crs="WGS84") %>%
#   sf::st_sf() %>%
#   dplyr::mutate(layer1_id = paste0("A", 1:nrow(.)))
#
# layer2_sf <- lapply(layer2b, sf::st_polygon) %>%
#   sf::st_sfc(crs="WGS84") %>%
#   sf::st_sf() %>%
#   dplyr::mutate(layer2_id = paste0("B", 1:nrow(.)))
#
#
# #purrr::map(layer1, sf::st_polygon) #sf::st_polygon(layer1a) %>%
# #  sf::st_sf(crs="wGS84")
# # create two population clusters
# pop_clusters <- create_pop_cluster(x=0.25, y=0.25, dx = 0.15, dy = 0.15, n=100, crs = "WGS84") %>%
#   dplyr::bind_rows(create_pop_cluster(x=-0.75, y=0.75, dx = 0.25, dy = 0.25, n=100, crs = "WGS84"))
#
# # create a plot
# ggplot() +
#   geom_sf(data = layer1_sf, colour = "red", size = 2) +
#   geom_sf(data = layer2_sf, fill="lightblue", alpha = 0.6) +
#   geom_sf(data = pop_clusters)
#
#
# layer1_count  <- onsr::get_pts_neighbourhood(pts = pop_clusters, pgon = layer1_sf) %>%
#   group_by(layer1_id) %>%
#   count(name = "layer1_n") %>%
#   drop_na() %>%
#   select(-geometry) %>%
#   left_join(layer1_sf, ., by = "layer1_id")
#
#
# layer2_count  <- onsr::get_pts_neighbourhood(pts = pop_clusters, pgon = layer2_sf) %>%
#   group_by(layer2_id) %>%
#   count(name = "layer2_n") %>%
#   drop_na() %>%
#   select(-geometry) %>%
#   left_join(layer2_sf, ., by = "layer2_id")
#
#
#
#
# sli <- greedy_sli_search(to_shp = layer1_count, to_idcol = "layer1_id", to_valuecol = "layer1_n",
#                   from_shp = layer2_count, from_idcol = "layer2_id", from_valuecol = "layer2_n" )
#
# layer2_sf %>%
#   left_join(sli) %>%
#   ggplot() +
#   geom_sf(aes(fill = layer1_id))
#
