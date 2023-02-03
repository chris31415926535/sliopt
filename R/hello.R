
measure_sli_performance4 <- function(sli_for_test, from_values, from_idcol, from_valuecol, to_values, to_idcol, to_valuecol){


  .to_idcol <- .from_value <- sli_values<- .to_value<- diff_count<- diff_prop<- NULL

  # rename da input for easier processing
  from_values <- from_values %>%
    dplyr::rename(.from_value = from_valuecol,
                  .from_idcol = from_idcol) %>%
    sf::st_set_geometry(NULL)


  to_values <- to_values %>%
    dplyr::rename(.to_value = to_valuecol,
                  .to_idcol = to_idcol ) %>%
    #dplyr::filter(ONS_ID != "0") %>%
    sf::st_set_geometry(NULL)

  sli_for_test <- dplyr::rename(sli_for_test,
                                .from_idcol = from_idcol,
                                .to_idcol = to_idcol)

  # calculate sli results
  sli_created_results <-  from_values %>%
    dplyr::left_join(sli_for_test, by = ".from_idcol") %>%
    dplyr::group_by(.to_idcol) %>%
    dplyr::summarise(sli_values = sum(.from_value, na.rm = TRUE))  %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), tidyr::replace_na, 0))

  # combine da and goldstandard results, compare
  compare_values <- to_values %>%
    dplyr::left_join(sli_created_results, by = ".to_idcol") %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), tidyr::replace_na, 0))%>%
    dplyr::mutate(diff_count = sli_values - .to_value,
                  diff_prop = diff_count/.to_value)

  compare_values %>%
    dplyr::summarise(mae = mean(abs(diff_count)),
                     #  mae_sd = mean(abs(diff_count)),
                     mape = mean(abs(diff_prop), na.rm = TRUE),
                     #   mape_sd = sd(abs(diff_prop)),
                     mse = mean(diff_count^2)
    ) %>%
    dplyr::mutate(raw_comparison = list(compare_values))
}




# FIXME TODO
plot_sli_performance <- function(comparison, title = NA, plot_type = c("absolute", "proportional")) {
  stop("Function not implemented")
  # plot_type <- match.arg(plot_type, plot_type)
  #
  # if (plot_type == "absolute"){
  #   result <- comparison %>%
  #     ggplot2::ggplot() +
  #     ggplot2::geom_point(ggplot2::aes(x=dwellings_ons, y=dwellings_sli)) +
  #     ggplot2::geom_abline(slope = 1, intercept = 0) +
  #     # scale_y_continuous(limits = c(-2000, 2000)) +
  #     ggplot2::ggtitle(title)
  # }
  #
  # if (plot_type == "proportional"){
  #   result <- comparison %>%
  #     ggplot2::ggplot() +
  #     ggplot2::geom_point(ggplot2::aes(x=ONS_ID, y=diff_prop)) +
  #     ggplot2::geom_abline(slope = 0, intercept = 0) +
  #     ggplot2::scale_y_continuous(limits = c(-1, 1.1)) +
  #     ggplot2::ggtitle(title)
  # }
  #
  # return(result)
}


# FIXME TODO
leaflet_sli_performance <- function(comparison, ons_shp, var = c("diff_count", "diff_prop")){
  stop("Function not implemented")
#
#   `:=` <- rlang::`:=`()
#
#   # set plot colour ranges based on what we know about inputs
#   # so they'll be constant across maps
#   if (var == "diff_count")  pal_domain <- c(-2000,2000)
#   if (var == "diff_prop")  pal_domain <- c(-1.2,1.2)
#
#
#   var <- match.arg(var, var)
#
#   forplot <- ons_shp %>%
#     sf::st_transform(crs = "WGS84") %>%
#     dplyr::mutate(ONS_ID = as.character(ONS_ID)) %>%
#     dplyr::left_join(comparison, by = "ONS_ID") %>%
#     dplyr::rename(plot_var := {{var}})
#
#   pal <- leaflet::colorNumeric(palette = "viridis",
#                                domain = pal_domain #c(min(forplot$plot_var, na.rm = TRUE),#max(forplot$plot_var, na.rm = TRUE)))
#   )
#
#   if (var == "diff_prop"){
#     forplot <- forplot %>%
#       dplyr::mutate(labels = sprintf("%s (%s)<br>%% Difference: %.1f%%", Name, ONS_ID,  100*plot_var) %>%
#                       purrr::map(htmltools::HTML))
#   }
#
#
#   if (var == "diff_count"){
#     forplot <- forplot %>%
#       dplyr::mutate(labels = sprintf("%s (%s)<br>Count Difference: %s", Name, ONS_ID,  plot_var) %>%
#                       purrr::map(htmltools::HTML))
#   }
#
#   forplot %>%
#     leaflet::leaflet() %>%
#     leaflet::addTiles() %>%
#     leaflet::addPolygons(weight = 1,
#                          fillColor = ~ pal(plot_var),
#                          fillOpacity = 0.9,
#                          label = ~ labels
#                          #label = ~ purrr::map(sprintf("%s (%s)<br>%s: %s", Name, ONS_ID, var, plot_var), htmltools::HTML)
#     ) %>%
#     leaflet::addScaleBar(position = "bottomright") %>%
#     leaflet::addLegend(position = "bottomleft", pal = pal, values = pal_domain, # ~plot_var,
#                        title = var)


}





# # measure performance of different SLIs
# # FIXME TODO! Generalize this, right now it's very specific
# create_sli_performance3 <- function(da_ons_intersect, da_values, da_values_column, goldstandard_values, goldstandard_values_column, da_ons_sli, da_ons_sli_opt_mae, da_ons_sli_opt_mape, da_ons_sli_opt_mse) {
#
#   measure_weighted_performance2(da_ons_intersect, da_values, da_values_column = da_values_column, goldstandard_values, goldstandard_values_column = goldstandard_values_column) %>%
#     dplyr::mutate(method = "Weighted Link", .before = 1) %>%
#     dplyr::bind_rows(
#       measure_sli_performance3(da_ons_sli, da_values, da_values_column = da_values_column, goldstandard_values, goldstandard_values_column = goldstandard_values_column) %>%
#         dplyr::mutate(method = "SLI unoptimized", .before = 1)
#     ) %>%
#     dplyr::bind_rows(
#       measure_sli_performance3(da_ons_sli_opt_mae, da_values, da_values_column = da_values_column, goldstandard_values, goldstandard_values_column = goldstandard_values_column) %>%
#         dplyr::mutate(method = "SLI mean absolute error", .before = 1)
#     ) %>%
#     dplyr::bind_rows(
#       measure_sli_performance3(da_ons_sli_opt_mape, da_values, da_values_column = da_values_column, goldstandard_values, goldstandard_values_column = goldstandard_values_column) %>%
#         dplyr::mutate(method = "SLI mean absolute percentage error", .before = 1)
#     ) %>%
#     dplyr::bind_rows(
#       measure_sli_performance3(da_ons_sli_opt_mse, da_values, da_values_column = da_values_column, goldstandard_values, goldstandard_values_column = goldstandard_values_column) %>%
#         dplyr::mutate(method = "SLI mean squared error", .before = 1)
#     ) %>%
#     dplyr::bind_cols(dplyr::tibble(
#       sli = list(da_ons_intersect, da_ons_sli, da_ons_sli_opt_mae, da_ons_sli_opt_mape, da_ons_sli_opt_mse )
#     ))
#
# }


#' Generate one-to-many weighted links between sets of spatial regions
#'
#' @param from_shp An `sf` object with the geometries we are linking from.
#' @param from_idcol A length-one character vector with the name of the `from_shp` column with unique region identifiers.
#' @param to_shp An `sf` object with the geometries we are linking to.
#' @param to_idcol A length-one character vector with the name of the `to_shp` column with unique region identifiers.
#'
#' @return A `tbl_df` with three columns, showing the linked-from region ids, the linked-to region ids, and the weights for each link. Each from-id set of weights should sum to 1.
#' @export
#'
#' @examples
#'  \dontrun{
#'  generate_weighted_links(da_ott, "DAUID", ons_shp, "ONS_ID")
#'  }
generate_weighted_links <- function(from_shp, from_idcol, to_shp, to_idcol){

  `:=` <- rlang::`:=`
  .from_idcol <- .to_idcol <- .area <- geometry <- NULL

  # Basic input validation
  # Are our input shapefiles really shapefiles?
  if (!"sf" %in% class(from_shp) | !"sf" %in% class(to_shp)) {
    stop("Input variables from_shp and to_shp must both be simple feature objects with class sf.")
  }

  # Do the id columns exist in the input files?
  if (!from_idcol %in% names(from_shp) | !to_idcol %in% names(to_shp)) {
    stop ("Input variables from_idcol and to_idcol must be names of columns in from_shp and to_shp.")
  }

  # Create our working data frames. Rename columns to make it clearer.
  from_shp <- from_shp %>%
    dplyr::rename(.from_idcol := {{from_idcol}}) %>%
    dplyr::select(.from_idcol)

  to_shp <- to_shp %>%
    dplyr::rename(.to_idcol := {{to_idcol}}) %>%
    dplyr::select(.to_idcol)

  # Set constant attributes so we don't get a warning
  sf::st_agr(from_shp) = "constant"
  sf::st_agr(to_shp) = "constant"

  # Get intersection area in temporary column
  result <- sf::st_intersection(from_shp, to_shp) %>%
    dplyr::mutate(.area = sf::st_area(geometry)) %>%
    sf::st_set_geometry(NULL) %>%

    # Get weights from area proportions
    dplyr::group_by(.from_idcol) %>%
    #dplyr::filter(dplyr::n() > 1) %>%
    dplyr::arrange(.from_idcol) %>%
    dplyr::mutate(weight = as.numeric(.area/sum(.area))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.area) %>%

    # Set names back to originals
    dplyr::rename(!!rlang::sym(from_idcol) := .from_idcol,
                  !!rlang::sym(to_idcol) := .to_idcol)

  result
}



#' Created single-link indicator from weighted links using maximum weights
#'
#' @param weighted_links A `tbl_df` containing weighted links from one set of ids to another, such as created by `generate_weighted_links()`.
#' @param from_idcol Character. Name of column containing ids.
#' @param weights_col Character (default "weight"). Name of column containing weights.
#'
#' @return A `tbl_df` with two columns and one row for each unique value in input column `from_idcol`.
#' @export
create_sli <- function(weighted_links, from_idcol, weights_col = "weight"){

  `:=` <- rlang::`:=`()
  .weight <- NULL

  weighted_links %>%
    dplyr::rename(.weight := {{weights_col}}) %>%
    dplyr::group_by(!!rlang::sym(from_idcol)) %>%
    dplyr::filter(.weight == max(.weight)) %>%
    dplyr::select(-.weight) %>%
    dplyr::ungroup()

}



#' Create optimized single-link indicators using a greedy search algorithm
#'
#' Sometimes you have data measured over one set of geographical regions (e.g.
#' census tracts) and you want to estimate how that data would look in an
#' overlapping but different set of geographical regions (e.g. school boards).
#' If the one set of geometries fits inside the other perfectly, then aggregating
#' your data is straightforward. But if the two sets of geometries overlap
#' messily, there is no straightforward "correct" way to link origin regions to
#' destination regions.
#'
#' This function implements an algorithm to create a single-link indicator that
#' matches each "from" region to one and only one "to" region. It requires a
#' a source of truth"--no small requirement--that is, a set of known-to-be-
#' accurate measurements of the same quantity over both sets of regions. It then
#' uses a greedy algorithm to attempt to link origin regions to destination regions
#' so as to minimize some score (mean average error, mean average percentage
#' error, or mean squared error).
#'
#' Origin and destination regions are only candidates for linking if they overlap
#' by more than an optional threshold value. This is to prevent spurious linkages
#' due to small data discrepancies (e.g. if two regions share a real-world
#' boundary like a road, but their geospatial data representations overlap in a
#' tiny sliver because of small differences in geocoding).
#'
#' The algorithm can start with an arbitrary input SLI, or else creates one
#' by assigning each origin to the destination it overlaps the most.
#'
#' For each origin region, the algorithm considers each possible destination region,
#' assigns it to that region, computes the new overall score, and accepts that
#' assignment if it results in a lower score.
#'
#' Algorithm progress is returned as an attribute called "optimization_progress".
#'
#' @param from_shp A tibble of origins that is also a simplefeature object, where each row represents one region.
#' @param from_idcol Character. The name of the `from_shp` column containing a unique identifier.
#' @param from_valuecol Character. The name of the `from_shp` column containing origin-region value of interest.
#' @param to_shp A tibble of destinations that is also a simplefeature object, where each row represents one region.
#' @param to_idcol Character. The name of the `to_shp` column containing a unique identifier.
#' @param to_valuecol Character. The name of the `to_shp` column containing the destination-region value of interest. This is the "gold standard" that is required to optimize the linkage.
#' @param optimize_for Character. Which measure are we optimizing? Accepts "mse", "mae", or "mape".
#' @param tolerance Numeric. What is the minimum fractional overlap allowed? Defaults to 0.05 (i.e. 5\%).
#' @param iterations Integer. How many times should we repeat the entire process? Default is 3.
#' @param input_sli Optional. A tibble containing an SLI from the origin to destination regions.
#' @param shuffle_inputs Boolean. Should we shuffle the order in which we consider regions? Default TRUE.
#' @param verbose Boolean. For debugging, gives too much info. Default FALSE.
#'
#' @return An SLI
#' @export
#'
greedy_sli_search <- function(from_shp, from_idcol, from_valuecol, to_shp, to_idcol, to_valuecol, optimize_for = c("mape", "mae", "mse"), tolerance = 0.05, iterations = 3, input_sli = NA, shuffle_inputs = TRUE, verbose = FALSE){

  `:=` <- rlang::`:=`

  # for R CMD CHECK
  .from_idcol <- weight<-  .to_idcol<- .from_valuecol <- .to_valuecol <- NULL
  optimize_for <- match.arg(optimize_for, optimize_for)


  # set up working SLI using weighted intersections
  # get from/to intersections and weights
  from_to_weighted <- generate_weighted_links(from_shp = from_shp, from_idcol = from_idcol, to_shp = to_shp, to_idcol)
  from_to_weighted <- dplyr::rename(from_to_weighted, .from_idcol = !!rlang::sym(from_idcol), .to_idcol = !!rlang::sym(to_idcol))


  # get to/from intersections and weights
  to_from_weighted <- generate_weighted_links(from_shp = to_shp, from_idcol = to_idcol, to_shp = from_shp, from_idcol)
  to_from_weighted <- dplyr::rename(to_from_weighted, .from_idcol = !!rlang::sym(to_idcol), .to_idcol = !!rlang::sym(from_idcol))


  # If no input SLI provided, create one
  if (is.na(input_sli)) {

    sli_new <- from_to_weighted %>%
      dplyr::group_by(.from_idcol) %>%
      dplyr::arrange(dplyr::desc(weight)) %>%
      dplyr::slice_head(n=1)
  }

  # If input SLI provided, name its columns appropriately
  if (!is.na(input_sli)){
    sli_new <- input_sli %>%
      dplyr::rename(.to_idcol := {{to_idcol}},
                    .from_idcol := {{from_idcol}}) %>%
      dplyr::mutate(.to_idcol = as.character(.to_idcol))
  }


  # Create our working data frames. Rename columns to make it clearer.
  # We filted both the origin and destination shapes so that they only
  # include rows with corresponding SLI entries (don't optimize things that aren't
  # in the input SLI!!!)
  from_shp <- from_shp %>%
    dplyr::rename(.from_idcol := {{from_idcol}},
                  .from_valuecol := {{from_valuecol}}) %>%
    dplyr::select(.from_idcol, .from_valuecol) %>%
    dplyr::mutate(.from_idcol = as.character(.from_idcol)) %>%
    dplyr::filter(.from_idcol %in% sli_new$.from_idcol)

  to_shp <- to_shp %>%
    dplyr::rename(.to_idcol := {{to_idcol}},
                  .to_valuecol := {{to_valuecol}}) %>%
    dplyr::select(.to_idcol, .to_valuecol) %>%
    dplyr::mutate(.to_idcol = as.character(.to_idcol)) #%>%    dplyr::filter(.to_idcol %in% sli_new$.to_idcol)

  # By default we randomize the order in case that makes some difference
  if (shuffle_inputs){
    to_shp <- to_shp[sample(1:nrow(to_shp)), ]
  }

  from_ids <- from_shp$.from_idcol
  to_ids <- to_shp$.to_idcol




  # pre-compute intersections in both directions. this takes ~250ms up front and saves ~30s per DA..
  # we can then access results via indexing a named list
  # FIXME precompute_intersections() now takes a tolerance argument to remove
  # small overlaps that we should ignore. default value is 0.001, 0.1% overlap.
  message("   Precomputing intersections...                             \r", appendLF = FALSE)
  # from_to_intersections <- precompute_intersections(from_shp = from_shp, from_idcol = ".from_idcol",
  #                                                   to_shp = to_shp, to_idcol = ".to_idcol",
  #                                                   tolerance = tolerance)
  # to_from_intersections <- precompute_intersections(from_shp= to_shp,
  #                                                   from_idcol = ".to_idcol",
  #                                                   to_shp = from_shp,
  #                                                   to_idcol = ".from_idcol",
  #                                                   tolerance = tolerance)
  from_to_intersections <- precompute_intersections(from_to_weights = from_to_weighted, from_idcol = ".from_idcol", to_idcol = ".to_idcol")
  to_from_intersections <- precompute_intersections(from_to_weights = to_from_weighted, from_idcol = ".from_idcol", to_idcol = ".to_idcol")

  # keep track of iterations and results
  optimization_progress <- dplyr::tibble()


  # set up just the origin and destination values without geometry for faster operations
  from_values <- sf::st_set_geometry(from_shp , NULL)
  to_values <- sf::st_set_geometry(to_shp , NULL)


  # to quiet warning as per https://github.com/r-spatial/sf/issues/406
  sf::st_agr(from_shp) <- "constant"
  sf::st_agr(to_shp) <- "constant"

  for (iteration in 1:iterations){

    if (verbose) message("   Iteration ", iteration, "...                                 \r", appendLF = FALSE)

    # if we're shuffling the inputs, do it now with the ids
    if (shuffle_inputs) to_ids <- base::sample(x = to_ids, size = length(to_ids), replace = FALSE)

    # loop through all destination regions
    for (i in 1:length(to_ids)){
      to_id <- to_ids[[i]]

      message(sprintf("%s/%s: Destination Region %s\r", i, length(to_ids), to_id), appendLF = FALSE)

      # get all from_regions that intersect this to_region
      froms_in_to <- to_from_intersections[to_id][[1]]

      # for each originating region that's contained within the destination region...
      for (from in froms_in_to){

        # check to see if DA has valid data, if it does not, skip it
        # if (is.na(da_values[da_values$DAUID == da,]$value)) next

        # find ONS_IDs of neighbourhoods that intersect the DA
        from_possible_tos <- from_to_intersections[from][[1]]

        # set up our results tibble
        to_results <- dplyr::tibble()

        # if it can only go to one possible region, skip it
        if (length(from_possible_tos) == 1) next

        if(verbose) message(paste0("  Checking all assignments for originating region ", from))
        # message(paste0("    ", from_possible_tos))
        #from_results <- dplyr::tibble()


        # for all possible destination regions this origin could possibly be assigned to...
        for (from_to in from_possible_tos){

          # Create temporary working copy of SLI
          sli_temp <- sli_new

          # update the SLI for analysis to assign this origin to the current destination
          sli_temp[sli_temp$.from_idcol == from,]$.to_idcol <- from_to

          ## NEXT: compute differences between SLI predictions and true values

          result <- get_sli_diffs2(sli_temp, from_values, to_values)


          if(verbose)  if (optimize_for == "mape" ) message(sprintf("    %s: mean proportional difference %.2f%%", from_to, 100*result$sli_diff_mape))
          if(verbose)  if (optimize_for == "mae"  ) message(sprintf("    %s: mean absolute error difference %.1f", from_to, result$sli_diff_mae))
          if(verbose)  if (optimize_for == "mse"  ) message(sprintf("    %s: mean squared error difference %.1f", from_to, result$sli_diff_mse))

          # collect all results for possible destination assignments
          to_results <- dplyr::bind_rows(to_results,
                                         dplyr::tibble(.to_idcol = from_to) %>%
                                           dplyr::bind_cols(result)
          )
        } # end for from_to in from_possible_tos


        # assign origin to destination that gives the best results
        # NOTE! Replacing ~3ms dplyr operation with ~222us base operation
        # then with ~ 32us base operation with drop=TRUE. pointless optimization!!
        #to_results$abs_result <- unlist(abs(to_results[paste0("sli_diff_", optimize_for)]))
        to_results$abs_result <- to_results[, paste0("sli_diff_", optimize_for) ,  drop = TRUE]

        best_to <- to_results[to_results$abs_result == min(to_results$abs_result),]$.to_idcol[[1]]


        if(verbose) message(paste0("     Best destination assignment for origin region ", from, ": ", best_to, ": ", min(to_results$abs_result)))
        message(sprintf("%s/%s: Destination Region %s: Best %s: %s                     \r", i, length(to_ids), to_id, optimize_for, signif(min(to_results$abs_result), digits = 4)), appendLF = FALSE)
        #message(sprintf("    Best  %s: %s", optimize_for, min(to_results$abs_result)))

        # if we get one best unique answer, we use that one.
        # otherwise don't update, something went wrong (this shouldn't happen!!)
        if (length(best_to) == 1) sli_new[sli_new$.from_idcol == from, ]$.to_idcol <- best_to

      } # end for from in froms_in_to

      # keep track of progress
      # if we found results
      if (nrow(to_results) > 0){
        optimization_progress <- dplyr::bind_rows(optimization_progress,
                                                  dplyr::tibble(iteration = iteration, step = i, metric = optimize_for, value =  min(abs(to_results[paste0("sli_diff_", optimize_for)]))))
      }

      # if we found no results but have past results, no change.
      if (nrow(to_results) == 0 & nrow(optimization_progress) > 0) {
        optimization_progress <- dplyr::bind_rows(optimization_progress,
                                                  dplyr::tibble(iteration = iteration, step = i, metric = optimize_for, value = optimization_progress$value[[nrow(optimization_progress)]] ))
      }

    } # end for to_id in to_ids (i.e. all destinations we are analyzing)



  } # end for iteration in iterations (i.e. the number of times we are repeating the entire process)

  # Create a nice clean output set with proper column names
  result <- sli_new %>%
    dplyr::select(.from_idcol, .to_idcol) %>%
    dplyr::rename(!!rlang::sym(from_idcol) := .from_idcol,
                  !!rlang::sym(to_idcol) := .to_idcol)

  # add optimization progress to result metadata
  # each trial will be added, note that iterations are across trials
  opt_metadata <- dplyr::bind_rows(attr(result, "optimization_progress"), optimization_progress)
  #opt_metadata$iteration <- 1:nrow(opt_metadata)
  attr(result, "optimization_progress") <- opt_metadata

  # return the updated sli
  return (result)

} # end of function



# Internal function?
# pre-compute the to_shp regions that each from_shp region intersects, and return
# these as a named list for quick & easy indexing
# from_to_weights = to_from_weighted; from_shp = to_shp; from_idcol = ".to_idcol"; to_shp = from_shp; to_idcol = ".from_idcol"
precompute_intersections <- function(from_to_weights, from_idcol, to_idcol, tolerance = 0.001){

  weight <- NULL
  # nest the results to create our result list
  step1 <- from_to_weights %>%
    dplyr::filter(weight > tolerance) %>%
    dplyr::select(-weight) %>%
    dplyr::group_by(!!rlang::sym(from_idcol)) %>%
    tidyr::nest(result = !!rlang::sym(to_idcol)) %>%
    dplyr::mutate(result = purrr::map(result, dplyr::pull, to_idcol))
  #dplyr::mutate(result = purrr::map(result, function(x) x[,!!rlang::sym(.to_idcol), drop=TRUE]))

  # extract the result into a simple list, then set the names for indexing
  result <- step1$result
  names(result) <- dplyr::pull(step1, !!rlang::sym(from_idcol))

  # done!
  return(result)

} # precompute_intersections()

# function to compare SLI results with gold-standard results
get_sli_diffs2 <- function(sli_temp, from_values, to_values, do_not_summarize = FALSE) {

  .to_idcol <- .from_valuecol <- value_sli <- .to_valuecol <- sli_diff_count <- sli_diff_prop <- NULL

  # using single-link indicator, estimate destination values
  result_sli <- sli_temp %>%
    dplyr::left_join(from_values, by = ".from_idcol") %>%
    dplyr::group_by(.to_idcol) %>%
    dplyr::summarise(value_sli = sum(.from_valuecol, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.to_idcol) %>%
    dplyr::mutate(.to_idcol = as.character(.to_idcol))

  ### PUT RESULTS TOGETHER, replace NAs with 0
  result <- to_values %>%
    dplyr::left_join(result_sli, by = ".to_idcol") %>%
    dplyr::mutate(value_sli = dplyr::if_else(is.na(value_sli), 0, as.numeric(value_sli))) %>%

    # calculate differences
    dplyr::mutate(sli_diff_count = (value_sli - .to_valuecol),
                  sli_diff_prop = sli_diff_count/.to_valuecol) %>%
    dplyr::mutate(sli_diff_prop = dplyr::if_else(.to_valuecol == value_sli, 0, sli_diff_prop)) %>% # to handle 0/0 = NaN
    dplyr::mutate(sli_diff_prop = dplyr::if_else(sli_diff_prop > 5, 5, sli_diff_prop)) %>%      # to handle 0/n = Inf
    dplyr::mutate(sli_diff_prop = dplyr::if_else(sli_diff_prop < -5, -5, sli_diff_prop)) %>%      # to handle -0/n = -Inf
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), tidyr::replace_na, 0))

  if (!do_not_summarize) {
    result <- result %>%
      # summarise
      dplyr::summarise(
        sli_diff_mae = mean(abs(sli_diff_count)),
        sli_diff_mape = mean(abs(sli_diff_prop)),
        sli_diff_mse = mean(sli_diff_count^2)
      )
  }

  return(result)
}







plot_sli_progress <- function(sli) {

  # for clean R CMD CHECK
  iteration <- step <- value <- NULL

  optimization_progress <- attr(sli, "optimization_progress")

  ggplot2::ggplot(optimization_progress) +
    ggplot2::geom_point(ggplot2::aes(x=((iteration-1) * max(step)) + step,
                                     y = value)) +
    ggplot2::labs(title = sprintf("Optimization progress")) +
    ggplot2::ylab(optimization_progress$metric[[1]]) +
    ggplot2::xlab("Step")


}




#' Create point clusters for synthetic data testing
#'
#' @param x Centre of distribution
#' @param y Centre of distribution
#' @param dx Standard deviation in x dimension
#' @param dy Standard deviation in y dimension
#' @param n Number of points
#' @param crs CRS, optional. Defaults to 32189 for Ottawa, ON.
#'
#' @return A simple feature collection of points.
#' @export
create_pop_cluster <- function(x, y, dx, dy, n, crs = 32189) {

  xs <- stats::rnorm(n = n, mean=x, sd = dx)
  ys <- stats::rnorm(n = n, mean=y, sd = dy)

  dplyr::tibble(x=xs, y=ys) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = crs)
}



#' Create synthetic population data for testing
#'
#' @param source_shape Input shape whose bounding box we'll use
#' @param num_clusters Number of clusters to create.
#' @param cluster_size SD of cluster size, in units of source_shape.
#' @param cluster_pop Population of clusters.
#'
#' @return Simple feature tibble with population clusters.
#' @export
create_synthetic_population <- function(source_shape, num_clusters = 100, cluster_size = 2500, cluster_pop = 1000){

  bbox <- sf::st_bbox(source_shape)
  clusters <- NA
  result_list <- list()

  for (i in 1:num_clusters){
    message(i)
    cluster_x <- stats::runif(n=1, min=bbox["xmin"], max=bbox["xmax"])
    cluster_y <- stats::runif(n=1, min=bbox["ymin"], max=bbox["ymax"])

    cluster <- create_pop_cluster(cluster_x, cluster_y, cluster_size, cluster_size, cluster_pop)
    result_list[[i]] <- cluster

  } # for (i in 1:num_clusters)

  clusters <- dplyr::tibble(temp = result_list) %>%
    tidyr::unnest(cols = c("temp")) %>%
    sf::st_as_sf() %>%
    sf::st_set_crs(sf::st_crs(source_shape))

  return(clusters)

} # create_synthetic_population()
