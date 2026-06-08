set.seed(42)

# --- Profiling (script-side profvis; no Dispersal.R hooks) ---
profile_timestep <- FALSE
profile_replicate <- 1L
profile_timestep_n <- 25L
profile_output <- ""  # empty -> <outputdir>/profvis_timestep_r{r}_tm{tm}.html

# --- Checkpoint (warm path cache for profiling later timesteps) ---
# Save at end of tm N (writes CHECKPOINT_FILE_OUT or <outputdir>/checkpoint_tm{N}_r{r}.rds):
#   CHECKPOINT_TM=24 CHECKPOINT_STOP=1 Rscript bsspread-compare.R
# Resume from tm N+1 (reads CHECKPOINT_FILE_IN or legacy CHECKPOINT_FILE):
#   CHECKPOINT_RESUME=1 CHECKPOINT_FILE_IN=/path/checkpoint_tm24_r1.rds Rscript ...
# Example: tm=24 in, run tm=25 on main, save baseline tm=25 out:
#   CHECKPOINT_RESUME=1 CHECKPOINT_FILE_IN=.../checkpoint_tm24_r1.rds \
#     CHECKPOINT_TM=25 CHECKPOINT_FILE_OUT=.../checkpoint_tm25_r1.rds CHECKPOINT_STOP=1 \
#     Rscript bsspread-compare.R
checkpoint_tm <- as.integer(Sys.getenv("CHECKPOINT_TM", unset = "0"))
checkpoint_replicate <- as.integer(Sys.getenv(
    "CHECKPOINT_REPLICATE",
    unset = "2"))
checkpoint_resume <- identical(Sys.getenv("CHECKPOINT_RESUME"), "1")
checkpoint_stop <- identical(Sys.getenv("CHECKPOINT_STOP"), "1")
# Re-seed before each replicate's initialise (same session, warm path cache):
#   BENCHMARK_RESEED_PER_REPLICATE=1  (uses params$random_seed)
benchmark_reseed_per_replicate <- identical(
    Sys.getenv("BENCHMARK_RESEED_PER_REPLICATE"), "1")
# Run replicates in parallel (mclapply); inner region/dispersal stay serial:
#   BENCHMARK_PARALLEL_REPLICATES=1 TASK_CPUS=4
benchmark_parallel_replicates <- identical(
    Sys.getenv("BENCHMARK_PARALLEL_REPLICATES"), "1")
# Skip results$finalize() and save_rasters/csv/plots (timing-only benchmark runs):
#   BENCHMARK_SKIP_SAVE=1
benchmark_skip_save <- identical(Sys.getenv("BENCHMARK_SKIP_SAVE"), "1")

# data_utils.R
# Platform wrapper utils
# -----------------------------------------------

# Return TRUE if 'source' defines 'name'
name_exists <- function(source, name) {
    if (name %in% names(source)) { 
        return(TRUE)
    }
    return(FALSE)
}

#' Return value of source[name] if defined, otherwise return the specified 'default' value.
apply_default <- function(source, name, default) {
    if (name %in% names(source)) { 
        return(source[name])
    }
    return(default)
}

#' Return a \code{terra::SpatRaster} 
#' by looking up a filename with the convention:
#'    params[name].filename 
#' @param params A JSON object containing an object property that describes a raster file { "filename": "" }
#' @param name Name of the JSON object property
#' @param conform_rast Optional raster to conform the output using bsrmap::conform_layer
get_rast <- function(params, name, conform_rast = NULL, na_strategy = "zero") {
    
    # params[name].filename or params[name].raster.filename
    if (name %in% names(params) && 
        "filename" %in% c(names(params[[name]]), names(params[[name]]$raster))
    ) {
        if ("filename" %in% names(params[[name]])) {
            filename = params[[name]]$filename
        } else if ("filename" %in% names(params[[name]]$raster)) {
            filename = params[[name]]$raster$filename
        }
        print(paste("get_rast: Found (", name, ")", filename))
        rast <- terra::rast(filename)
        if (!is.null(conform_rast)){
            rast <- terra::rast(lapply(1:terra::nlyr(rast), function(i) {
                bsrmap::conform_layer(rast[[i]], conform_rast, 
                                      na_strategy = na_strategy)
            }))
        }
        return (rast)

    } else {
        print(paste("get_rast: Not found (", name, ")"))
        rast = NULL
    }
}

#' Conform an input layer to the resolution and extent of a template
conform_layer <- function(rast, template_rast, na_strategy = "zero", use_aggr_fun = "mean") {
    if (!is.null(rast) && !is.null(template_rast)){
        ret = terra::rast(lapply(1:terra::nlyr(rast), function(i) {
            bsrmap::conform_layer(
                rast[[i]], 
                template_rast,
                normalize = FALSE,
                binarize = FALSE,
                na_strategy = na_strategy,
                use_aggr_fun = use_aggr_fun)
        }))
        names(ret) <- names(rast)
    } else {
        ret = NULL
    }
    return(ret)
}

#' Return a dataframe read from a csv using platform parameter convention:
#' params[name].filename
#' or
#' params[name].data.filename
get_csv <- function(params, name) {
    # params[name].filename
    if (name %in% names(params) && 
        "filename" %in% names(params[[name]])
    ) {
        filename = params[[name]]$filename
    }
    # params[name].data.filename
    else if (name %in% names(params) && 
        "data" %in% names(params[[name]]) && 
        "filename" %in% names(params[[name]]$data)
    ) {
        filename = params[[name]]$data$filename
    } else {
        print(paste("get_csv: Not found (", name, ")"))
        return(NULL)
    }
    print(paste("get_csv: Found (", name, ")", filename))
    return(utils::read.csv(filename))
}

#' Get a named column from a data frame
get_df_col <- function(df, name, series = FALSE) {
    if (is.null(df)) { 
        return(NULL)
    }
    if (!series && name %in% names(df)) { 
        print(paste("get_df_col: Column found (", name, ")"))
        return(df[[name]])
    } else {
        name_series <- paste0(name, "_", 1:ncol(df))
        name_idx <- match(name_series, names(df))
        name_idx <- name_idx[which(!is.na(name_idx))]
        if (length(name_idx) > 0) {
            name_series <- name_series[1:length(name_idx)]
            print(paste("get_df_col: Columns found (", 
                        paste(name_series, collapse = ", "), ")"))
            return(as.list(df[name_series]))
        }
    }
    print(paste("get_df_col: Column not found (", name, ")"))
    return(NULL)
}

# bsspread_factory.R
# Platform wrapper for bsspread package functions
# https://github.com/cebra-analytics/bsspread/
# -----------------------------------------------

library(bsspread)
library(data.table)

INITIALIZER_CSV_COLNAME="population"
INITIALIZER_RANDOM_CSV_COLNAME="weights"

POPULATION_GROWTH_CSV_COLNAME="growth"
POPULATION_GROWTH_MULT_CSV_COLNAME="growth_mult"
POPULATION_CAPACITY_CSV_COLNAME="capacity"
POPULATION_ESTABLISH_PR_CSV_COLNAME="suitability"

#' bsspread::Region object factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsspread/bsspread_region.json
bsspread_region_factory <- function(params) {
    if (is.null(params$region_type)){
        stop("bsspread_region_factory: region_type is not set")
    }

    print(paste("bsspread_region_factory:", params$region_type))

    switch(params$region_type,
        raster = {
                region_rast <- terra::rast(params$region_rast$filename)
                region <- bsspread::Region(region_rast)
                if ('na_strategy' %in% names(params)) {
                  region$na_strategy <- input.params$region$na_strategy
                } else {
                  region$na_strategy <- "zero"
                }
                if ('two_tier' %in% names(params) && params$two_tier == TRUE) {
                    region$set_aggr(params$aggr_factor, params$inner_radius)
                }
                return(region)
            },
        network = {
                region_point <- get_csv(params, 'region_point')
                if (is.null(region_point)){
                    stop("bsspread_region_factory: type network, 'region_point' is NULL")
                }
                region <- bsspread::Region(region_point)
                return(region)
            },
        none = {
                region <- bsspread::Region()
                if (is.numeric(params$region_max_area) && 
                    !is.null(params$region_max_area_units)) {
                  mult <- ifelse(params$region_max_area_units == "km",
                                 1000000, 1)
                  region$set_max_implicit_area(params$region_max_area*mult)
                }
                return(region)
            }
    )
}

#' bsspread::Initializer object factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsspread/bsspread_initializer.json
#' @param region A \code{bsspread::Region} or inherited class object
#' @param population_model A \code{bsspread::PopulationModel} or inherited class object
#' @param outputdir Dir to use when saving rasterized initial layer.
bsspread_initializer_factory <- function(params, region, population_model, outputdir) {
    print(paste("bsspread_initializer_factory:", params$initializer_type))

    if (is.null(region)){
        stop("bsspread_initializer_factory: region is NULL")
    }
    if (is.null(population_model)){
        stop("bsspread_initializer_factory: population_model is NULL")
    }

    # Region: Spatially implicit
    # The presence-only (binary) population is simply initialised with a true value.
    # Otherwise the population is initialised with 'population_size'.
    if (region$spatially_implicit()){
        if (population_model$get_type() == 'presence_only'){
            initial_data = TRUE
        } else {
            initial_data = switch(params$initializer_type,
                initial_layer = {
                    if (params$population_type == "stage_structured") {
                        attr(params$population_size, "stages") <- unlist(params$initial_stages)
                    }
                    params$population_size
                },
                random = {
                    bsspread::Incursions(
                        1, region,
                        incursion_mean = params$incursion_mean,
                        incursion_stages = unlist(params$incursion_stages)
                    )
                }
            )
        }
        return(bsspread::Initializer(
            initial_data,
            region = region,
            population_model = population_model)
        )
    }

    # Region: Grid, Patch
    if (is.null(params$initializer_type)){
        stop("bsspread_initializer_factory: initializer_type not set")
    }

    initial_data = switch(region$get_type(),
        grid = {
            if (name_exists(params, 'initial_geojson') 
                && !is.null(params$initial_geojson)
                    && !is.null(params$initial_geojson$filename)) {
              print(paste('Rasterizing initial_geojson', params$initial_geojson$filename, 
                          'using', params$initial_rast$filename, 'as the region template'))
              if (params$population_type %in% c("unstructured", "stage_structured") &&
                  is.null(params$population_size)){
                stop("Cannot use initial_geojson without a population_size value")
              }
              
              initial_vect <- terra::vect(params$initial_geojson$filename)
              
              # Project to template
              template_rast <- region$get_rast(1)
              initial_vect <- terra::project(
                initial_vect, 
                template_rast
              )
              # Rasterize geometry
              initial_rast <- terra::rasterize(
                initial_vect, 
                template_rast, 
                touches = TRUE
              )
              # Conform to template
              initial_rast <- bsrmap::conform_layer(
                initial_rast,
                template_rast
              )
              # Write to file
              if (!is.null(outputdir)){
                terra::writeRaster(
                  initial_rast, 
                  file.path(outputdir, 'initial_rast.tif'), 
                  gdal = c("COMPRESS=LZW", "TILED=YES"),
                  overwrite = TRUE
                )
              }
              if (params$population_type %in% c("unstructured", "stage_structured")) {
                attr(initial_rast, "size") <- params$population_size
              }
              initial_rast 
              
            } else {
                get_rast(
                    params, 'initial_rast',
                    region$get_rast(1),
                    na_strategy = region$na_strategy
                )
            }
        },
        patch = {
            df = get_csv(params, 'initial_point')
            if (params$initializer_type == 'random'){
                get_df_col(df, INITIALIZER_RANDOM_CSV_COLNAME)
            } else {
                get_df_col(df, INITIALIZER_CSV_COLNAME)
            }
        }
    )

    if (is.null(initial_data)){
        stop("bsspread_initializer_factory: Could not determine initial data from available parameters.")
    }
    
    # Attach initial age and/or stages attributes
    if (params$initializer_type == "initial_layer") {
        if (!is.null(params$initial_age_type) 
            && "initial_age" %in% names(params) 
            && !is.null(params$initial_age)
        ) {
            attr(initial_data, "age") <- switch(params$initial_age_type,
                single = params$initial_age,
                spatial = switch(params$region_type,
                    raster = 
                      {
                          rast = get_rast(
                              params, 'initial_age',
                              region$get_rast(1),
                              na_strategy = region$na_strategy
                          )
                          if (!is.null(rast)){
                              indices <- rast[region$get_indices()][,1]
                          }
                      },
                    network = 
                        {
                            get_df_col(get_csv(params, 'initial_age'), "age")
                        }
                    )
            )
        }
        if (params$population_type == "stage_structured" 
            && !is.null(params$initial_stages)
        ) {
            attr(initial_data, "stages") <- unlist(params$initial_stages)
        }
    }

    switch(params$initializer_type,
        initial_layer = bsspread::Initializer(
            initial_data,
            region = region,
            population_model = population_model
        ),
        random = bsspread::Initializer(
            bsspread::Incursions(initial_data,
                region = region,
                type = "weight",
                #continued = params$continued,
                #time_steps = params$time_steps
                incursion_mean = params$incursion_mean,
                incursion_stages = unlist(params$incursion_stages)
                ),
            region = region,
            population_model = population_model
        )
    )
}

#' bsspread::Population object factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsspread/bsspread_population.json
#' @param region A \code{bsspread::Region} or inherited class object
bsspread_populaton_model_factory <- function(params, region) {
    print(paste("bsspread_populaton_model_factory:", params$population_type))
    capacity        = NULL
    establish_pr    = NULL

    if (is.null(region)){
        stop("bsspread_populaton_model_factory: region is NULL")
    }

    capacity <- switch(params$region_type,
        raster = 
            {
                rast = get_rast(
                    params, 'capacity_rast', 
                    region$get_rast(1), 
                    na_strategy = region$na_strategy
                )
                if (!is.null(rast)){
                    as.matrix(rast[region$get_indices()])
                }
            },
        network = 
            {
                capacity_data <- get_df_col(
                    get_csv(params, 'capacity_point'), 
                    POPULATION_CAPACITY_CSV_COLNAME,
                    series = (
                        !is.null(params$capacity_type) 
                        && params$capacity_type == "spatiotemporal"
                    )
                )
                if (is.list(capacity_data)) {
                    as.matrix(as.data.frame(capacity_data))
                } else if (!is.null(capacity_data)) {
                    as.matrix(capacity_data)
                }
            },
        none = 
            {
                if (!is.null(params$capacity_type)) {
                    switch(params$capacity_type,
                           single = params[["capacity"]],
                           temporal = t(as.matrix(
                                get_df_col(get_csv(params, 'capacity_temp'), 
                                    POPULATION_CAPACITY_CSV_COLNAME
                                )
                            )
                        )
                    )
                }
            }
    )
    establish_pr <- switch(params$region_type,
        raster = 
            {
                rast = get_rast(
                    params, 'establish_pr_rast', 
                    region$get_rast(1), 
                    na_strategy = region$na_strategy
                )
                if (!is.null(rast)){
                    as.matrix(rast[region$get_indices()])
                }
            },
        network = 
            {
                establish_pr_data <- get_df_col(
                    get_csv(params, 'establish_pr_point'), 
                    POPULATION_ESTABLISH_PR_CSV_COLNAME,
                    series = (
                        !is.null(params$establish_pr_type) 
                        && params$establish_pr_type == "spatiotemporal")
                    )
                if (is.list(establish_pr_data)) {
                    as.matrix(as.data.frame(establish_pr_data))
                } else if (!is.null(establish_pr_data)) {
                    as.matrix(establish_pr_data)
                }
            },
        none = 
            {
                NULL
            }
    )

    if (!is.null(params$growth_matrix)){
        growth_matrix = matrix(
            unlist(params$growth_matrix), 
            nrow = length(params$growth_matrix), 
            byrow = TRUE
        )
    } else {
        growth_matrix = NULL
    }

    if (!is.null(params$growth_matrix_mask)){
        mask = matrix(
            unlist(params$growth_matrix_mask), 
            nrow = length(params$growth_matrix_mask), 
            byrow = TRUE
        )
        attr(growth_matrix, "survivals") <- mask
    }

    if (!is.null(params$growth_type) 
        && params$growth_type %in% c("spatial", "spatiotemporal")
    ) {
        growth <- switch(params$region_type,
            raster = 
                {
                    rast = get_rast(
                        params, 'growth_rast', 
                        region$get_rast(1),
                        na_strategy = region$na_strategy
                    )
                    if (!is.null(rast)){
                        as.matrix(rast[region$get_indices()])
                    }
                },
            network = 
                {
                    growth_data <- get_df_col(
                        get_csv(params, 'growth_point'), 
                        POPULATION_GROWTH_CSV_COLNAME,
                        series = (params$growth_type == "spatiotemporal")
                    )
                    if (is.list(growth_data)) {
                        as.matrix(as.data.frame(growth_data))
                    } else {
                        as.matrix(growth_data)
                    }
                }
        )
    } else if (
        !is.null(params$growth_type) 
        && params$growth_type == "temporal" 
        && !is.null(params$growth_temp)
    ) {
        growth <- get_df_col(get_csv(params, 'growth_temp'), POPULATION_GROWTH_CSV_COLNAME)
        if (!is.null(growth)) {
            growth <- t(as.matrix(growth))
        }
    } else if (is.null(params$growth_type) || (!is.null(params$growth_type) && params$growth_type == "single")) {
        growth <- params[["growth"]]
    } else {
        growth <- NULL
    }

    if (!is.null(params$growth_mult)){
      if (!is.null(params$growth_mult$dim_type) 
        && params$growth_mult$dim_type %in% c("spatial", "spatiotemporal")
        ) {
          growth_mult <- switch(params$region_type,
              raster = 
                  {
                      rast = get_rast(
                          params, 'growth_mult', 
                          region$get_rast(1),
                          na_strategy = region$na_strategy
                      )
                      if (!is.null(rast)){
                          as.matrix(rast[region$get_indices()])
                      }
                  },
              network = 
                  {
                      growth_mult_data <- get_df_col(get_csv(params, 'growth_mult'), POPULATION_GROWTH_MULT_CSV_COLNAME,
                                                     series = (params$growth_mult$dim_type == "spatiotemporal"))
                      if (is.list(growth_mult_data)) {
                          as.matrix(as.data.frame(growth_mult_data))
                      } else {
                          as.matrix(growth_mult_data)
                      }
                  }
              )
        } else if ((!is.null(params$growth_mult$dim_type) 
            && params$growth_mult$dim_type == "temporal") 
            || params$region_type == "none"
        ) {
            growth_mult <- get_df_col(
                get_csv(params, 'growth_mult'), 
                POPULATION_GROWTH_MULT_CSV_COLNAME
            )
            if (!is.null(growth_mult)) {
                growth_mult <- t(as.matrix(growth_mult))
            }
        } else {
          growth_mult <- NULL
        }
        if (!is.null(growth_mult) 
            && !is.null(params$growth_mult$apply_to) 
            && params$growth_mult$apply_to %in% c("reproduction", "survival")
        ) {
            attr(growth_mult, "apply_to") <- params$growth_mult$apply_to
        }
        if (!is.null(growth_mult) && !is.null(params$growth_mult$stages)) {
            attr(growth_mult, "stages") <- unlist(params$growth_mult$stages)
        }
    } else {
        growth_mult <- NULL
    }
    
    if (!is.null(params$capacity_area)) {
      capacity_area <- switch(
          params$capacity_area,
          m   = 1,
          km  = 1000000,
          max = region$get_max_implicit_area()
      )
    } else {
        capacity_area <- NULL
    }

    population = switch(params$population_type,
        presence_only   = bsspread::PresencePopulation(region,
                            establish_pr = establish_pr,
                            spread_delay = params$spread_delay),
        unstructured    = bsspread::UnstructPopulation(region,
                            growth = growth,
                            capacity = capacity,
                            capacity_area = capacity_area,
                            establish_pr = establish_pr,
                            incursion_mean = params$incursion_mean),
        stage_structured = bsspread::StagedPopulation(region,
                            growth = growth_matrix,
                            growth_mult = growth_mult,
                            capacity = capacity,
                            capacity_area = capacity_area,
                            capacity_stages = unlist(params$capacity_stages),
                            establish_pr = establish_pr,
                            incursion_mean = params$incursion_mean, # moved to initialiser/incursions
                            incursion_stages = unlist(params$incursion_stages)) # TODO move to initialiser/incursions
        )

    if (is.null(population)){ 
        stop("bsspread_populaton_model_factory: model not initialized")
    }

    return(population)
}

#' bsspread::Attractor object factory
#'
#' @param params A JSON object described by 
#'               (Legacy) https://biosecuritycommons.org.au/bsspread/bsspread_dispersal_model_attractors.json
#'               (v1.19) https://biosecuritycommons.org.au/bsspread/bsspread_dispersal_model_items$definitions.attractor_rast.json
#'               (v1.19) https://biosecuritycommons.org.au/bsspread/bsspread_dispersal_model_items$definitions.attractor_point.json
#' @param region A \code{bsspread::Region} or inherited class object
bsspread_attractor_factory <- function(params, region) {
    if (is.null(params$type)) {
        params$type <- "attractor" # else retro source/destination
    }
    print(paste("bsspread_attractor_factory:", params$type))

    if (!is.null(region)){
        attractor_data = switch(region$get_type(),
            grid = {
                if (!"attractor_rast" %in% names(params)) {
                    params$attractor_rast <- params # v1.19
                }
                get_rast(
                    params, 'attractor_rast',
                    region$get_rast(1), 
                    na_strategy = region$na_strategy
                )
            },
            patch = {
                if (!"attractor_point" %in% names(params)) {
                  params$attractor_point <- params # v1.19
                }
                get_df_col(get_csv(params, 'attractor_point'), params$type)
            }
        )
    } else {
        attractor_data = NULL
    }

    attractor = bsspread::Attractor(attractor_data, region, is_dynamic = params$is_dynamic)
}

#' bsspread::Kernels object factory
#'
#' @param params A JSON object
bsspread_kernel_factory <- function(params) {
    print(paste("bsspread_kernel_factory:", params$function_name))

    if (!is.null(params$table)) {
        table = get_csv(params, 'table')
        
        # Check column names are present
        if (params$type == "direction" && !all(c("direction", "probability") %in% names(table))) {
            stop("Dispersal direction kernel lookup table should have columns 'direction' and 'probability'.",
                 call. = FALSE)
        }
        if (params$type == "distance" && !all(c("distance", "probability") %in% names(table))) {
          stop("Dispersal distance kernel lookup table should have columns 'distance' and 'probability'.",
               call. = FALSE)
        }
        
        # Squash duplicate directions via probability averaging and order as per bsspread package
        if (params$type == "direction") {
            if (any(duplicated(table[,"direction"]))) {
                table <- aggregate(probability ~ direction, data = table, FUN = mean)
            }
            table <- table[,c("direction", "probability")]
        }

        # Squash duplicate distances via probability averaging and order as per bsspread package
        if (params$type == "distance") {
            if (any(duplicated(table[,"distance"]))) {
                table <- aggregate(probability ~ distance, data = table, FUN = mean)
            }
            table <- table[,c("distance", "probability")]
        }
        
    } else {
        table = NULL
    }

    if (is.null(params$shift)){
        params$shift = 0
    }
    if (is.null(params$lower)){
        params$lower = 0
    }

    if (params$type == "direction") {
        if (is.null(params$direction_type)) {
            params$direction_type <- "navigational"
        }
        if (is.null(params$direction_orient)) {
            params$direction_orient <- "to"
        }
        kernels <- DirectionKernel(direction_type = params$direction_type,
                                   orientation = params$direction_orient,
                                   multiplier = 1)
        return(switch(params$function_name,
            Beta = kernels$get_beta_function(
                alpha=params$alpha,
                beta=params$beta,
                shift=params$shift),
            Lookup = kernels$get_lookup_function(table=table)
        ))
    } else { # distance
        kernels <- Kernels(multiplier = 1)
        return(switch(params$function_name,
            Beta = kernels$get_beta_function(
                alpha=params$alpha,
                beta=params$beta,
                upper=params$upper,
                shift=params$shift),
            Cauchy = kernels$get_cauchy_function(scale=params$scale),
            Exp = kernels$get_exp_function(mean=params$mean),
            Gaussian = kernels$get_gaussian_function(sd=params$sd),
            Lognormal = kernels$get_lognormal_function(mean=params$mean, sd=params$sd),
            Uniform = kernels$get_uniform_function(lower=params$lower, upper=params$upper),
            Weibull = kernels$get_weibull_function(shape=params$shape, scale=params$scale),
            Lookup = kernels$get_lookup_function(table=table)
        ))
    }
}

#' bsspread::DispersalModel object factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsspread/bsspread_dispersal_model.json
#' @param region A \code{bsspread::Region} or inherited class object
#' @param population_model A \code{bsspread::PopulationModel} or inherited class object
bsspread_dispersal_model_factory <- function(params, region, population_model) {
    print(paste("bsspread_dispersal_model_factory:", params$dispersal_type))

    # Resolve proportion
    if (!is.null(params$proportion_type) && params$proportion_type %in% c("spatial", "spatiotemporal")) {
        proportion <- switch(params$region_type,
            raster = 
                {
                    rast = get_rast(
                        params, 'proportion_rast', 
                        region$get_rast(1),
                        na_strategy = region$na_strategy
                    )
                    if (!is.null(rast)){
                        as.matrix(rast[region$get_indices()])
                    }
                },
            network = 
                {
                    proportion_data <- get_df_col(get_csv(params, 'proportion_point'), 'proportion', 
                                                  series = (params$proportion_type == "spatiotemporal"))
                    if (is.list(proportion_data)) {
                        as.matrix(as.data.frame(proportion_data))
                    } else {
                        as.matrix(proportion_data)
                    }
                }
        )
        if (params$region_type == "network" && !is.null(proportion) && !region$is_compatible(proportion[,1])) {
            stop("bsspread_dispersal_model_factory: Region not compatible with proportion data")
        }
    } else if (!is.null(params$proportion_type) && params$proportion_type == "temporal" && !is.null(params$proportion_temp)) {
        proportion <- get_df_col(get_csv(params, 'proportion_temp'), 'proportion')
        if (!is.null(proportion)) {
            proportion <- t(as.matrix(proportion))
        }
    } else if (is.null(params$proportion_type) || (!is.null(params$proportion_type) && params$proportion_type == "single")) {
        proportion <- params[["proportion"]]
    } else {
        proportion <- NULL
    }

    # Set proportion to 1 for presence-only models
    if (is.null(proportion) && params$population_type == "presence_only") {
        proportion <- 1
    }
  
    # Resolve events
    if (!is.null(params$events_type) && params$events_type %in% c("spatial", "spatiotemporal")) {
        is_temporal <- (params$events_type == "spatiotemporal")
        events <- switch(params$region_type,
            raster = 
                {
                    rast = get_rast(
                        params, 'events_rast', 
                        region$get_rast(1),
                        na_strategy = region$na_strategy
                    )
                    if (!is.null(rast)){
                        as.matrix(rast[region$get_indices()])
                    }
                },
            network = 
                {
                    events_data <- get_df_col(get_csv(params, 'events_point'), 'events', 
                                              series = (params$events_type == "spatiotemporal"))
                    if (is.list(events_data)) {
                        as.matrix(as.data.frame(events_data))
                    } else {
                        as.matrix(events_data)
                    }
                }
        )
        if (params$region_type == "network" && !is.null(events) && !region$is_compatible(events[,1])) {
            stop("bsspread_dispersal_model_factory: Region not compatible with events data")
        }
    } else if (!is.null(params$events_type) && params$events_type == "temporal" && !is.null(params$events_temp)) {
        events <- get_df_col(get_csv(params, 'events_temp'), 'events')
        if (!is.null(events)) {
            events <- t(as.matrix(events))
        }
    } else if (is.null(params$events_type) || (!is.null(params$events_type) && params$events_type == "single")) {
        events <- params[["events"]]
    } else {
        events <- NULL
    }

    # Resolve permeability
    if (!is.null(region)){
        permeability_data = switch(region$get_type(),
            grid = get_rast(
                params, 'permeability_rast', region$get_rast(1), 
                na_strategy = region$na_strategy
            ),
            patch = 
                { 
                    perm_point = get_csv(params, 'permeability_point')
                    if (!is.null(perm_point) && !region$is_compatible(perm_point)){
                        stop("bsspread_dispersal_model_factory: Region not compatible with 'permeability_point'")
                    }
                    perm_point
                }
        )
    } else {
        permeability_data = NULL
    }

    if (!is.null(permeability_data)) {
        print("Creating bsspread::Permeability")
        permeability = bsspread::Permeability(permeability_data, region)
    } else {
        permeability = NULL
    }

    attractors <- list()
    if (length(params$attractors) > 0) {
        params$attractors <- lapply(params$attractors, function(a) {
            if (is.logical(params$dynamic_attractors) && params$dynamic_attractors) {
                a$is_dynamic <- TRUE
            } else {
                a$is_dynamic <- FALSE
            }
            a
        })
        for (att in 1:length(params$attractors)) {
        
            # Load attractor object
            attractors[[att]] <-
              bsspread_attractor_factory(params$attractors[[att]], region)

            # Retrofit source/both types via proportion and/or events (< v1.19)
            if (!is.null(params$attractors[[att]]$type) &&
                params$attractors[[att]]$type %in% c("source", "both")){
              if (is.numeric(proportion)) {
                proportion <- proportion*attractors[[att]]$get_values()
              }
              if (is.numeric(events)) {
                events <- events*attractors[[att]]$get_values()
              }
              if (params$attractors[[att]]$type == "source") {
                attractors[[att]] <- NULL
              }
            }

        }
        attractors <- attractors[!sapply(attractors, is.null)] # remove NULLs
    }
    
    direction_function = NULL
    if ("direction_function" %in% names(params) && 
        "function_name" %in% attributes(params$direction_function)$names) {
        # Direction function always uses an upper of 360
        if (params$direction_function$function_name == 'Beta'){
            params$direction_function$upper = 360
        }
        params$direction_function$type <- "direction"
        direction_function = bsspread_kernel_factory(params$direction_function)
    }

    max_distance = NULL
    distance_adjust = NULL
    distance_function = NULL
    if ("distance_function" %in% names(params)) {
        if ("function_name" %in% names(params$distance_function)) {
            params$distance_function$type <- "distance"
            distance_function = bsspread_kernel_factory(params$distance_function)
        }

        if ("function_name" %in% names(params$distance_function) &&
            params$distance_function$function_name %in% c("Beta", "Uniform")) {
            max_distance = params$distance_function$upper
        } else if ("max_distance" %in% names(params$distance_function)) {
            max_distance = params$distance_function$max_distance
        }

        if ("distance_adjust" %in% names(params$distance_function) && region$get_type() == "grid"){
            distance_adjust = params$distance_function$distance_adjust
        }
    }

    # Set default distance scale to 1000 (i.e. kilometres) when absent (temporary)
    if (is.null(params$distance_scale)) {
      params$distance_scale <- 1000
    }

    switch(params$dispersal_type,
           diffusion = bsspread::Diffusion(
             region,
             population_model,
             diffusion_rate = params$diffusion_rate,
             proportion = proportion,
             direction_function = direction_function,
             attractors = attractors,
             permeability = permeability),
           
           radial_diffusion = bsspread::Diffusion(
             region,
             population_model,
             diffusion_rate = params$diffusion_rate),
           
           area_spread = bsspread::AreaSpread(
             region,
             population_model),
           
           reaction_diffusion = bsspread::Diffusion(
             region,
             population_model,
             diffusion_rate = params$diffusion_rate),
           
           gravity = bsspread::Gravity(
             region,
             population_model,
             attractors = attractors,
             attractor_function = params$attractor_function,
             beta = params$beta,
             distance_scale = params$distance_scale,
             dispersal_stages = unlist(params$dispersal_stages),
             proportion = proportion,
             events = events,
             direction_function = direction_function,
             permeability = permeability),
           
           dispersal = bsspread::Dispersal(
             region,
             population_model,
             dispersal_stages = unlist(params$dispersal_stages),
             proportion = proportion,
             events = events,
             distance_function = distance_function,
             direction_function = direction_function,
             distance_adjust = distance_adjust,
             attractors = attractors,
             permeability = permeability,
             max_distance = max_distance)
    )
}

#' bsspread::DispersalModel collection factory
#'
#' @param dispersal_models_params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsspread/bsspread_dispersal_models.json
#' @param region A \code{bsspread::Region} or inherited class object
#' @param population_model A \code{bsspread::PopulationModel} or inherited class object
bsspread_dispersal_models_factory <- function(dispersal_models_params, region, population_model, impact_params){
    dispersal_models = list()
    for (dm in 1:length(dispersal_models_params)) {
        
        # Dynamic attractors only when linked to impacts
        if (!any(sapply(impact_params, function(p) {
          if (is.logical(p$link_attractors) && p$link_attractors) TRUE else FALSE
        }))) {
          dispersal_models_params[[dm]]$dynamic_attractors <- NULL
        }
      
        # Add dispersal model
        dispersal_models[[dm]] <-
            bsspread_dispersal_model_factory(dispersal_models_params[[dm]], region, population_model)
    }
    return (dispersal_models)
}


## Platform wrapper for bsspread package functions
##
## https://github.com/cebra-analytics/bsspread/
## 
JOB_WEBLOG_URL=Sys.getenv("JOB_WEBLOG_URL")
JOB_WEBLOG_AUTH=Sys.getenv('NF_WEBLOG_BASIC_TOKEN')

# Number of cpu cores to optimise for
DEFAULT_CPUS=1

if (Sys.getenv("TASK_CPUS") > 0){
    cpus = as.numeric(Sys.getenv("TASK_CPUS"))
} else {
    cpus = DEFAULT_CPUS
}

paste("available cpus:", cpus)

library(httr2)
library(bsmanage)
paste("bsimpact:", packageVersion("bsimpact"))
paste("bsspread:", packageVersion("bsspread"))
paste("bsmanage:", packageVersion("bsmanage"))
paste("terra:", packageVersion("terra"))
paste("data.table:", packageVersion("data.table"))
paste("profvis:", packageVersion("profvis"))
paste("htmlwidgets:", packageVersion("htmlwidgets"))


# Load parameters
params = rjson::fromJSON(file="params.json")
input.params <- params
input.env <- params[['platform:env']]
rm(params)

if(!dir.exists(input.env$outputdir)) { dir.create(input.env$outputdir, recursive = TRUE) }

# Set terra tmpdir to job scratch space
terra::terraOptions(tempdir = input.env$workdir)

# Set working directory
setwd(input.env$workdir)

# Optional behavioural equivalence collector (manual workflow helper)
eq_enabled <- FALSE
eq_state <- NULL
eq_helper_paths <- c(
    "tools/equivalence_metrics.R",
    "/workdir/script/tools/equivalence_metrics.R",
    "/workdir/tools/equivalence_metrics.R"
)
eq_helper <- eq_helper_paths[file.exists(eq_helper_paths)][1]
if (!is.na(eq_helper)) {
    source(eq_helper)
    eq_enabled <- exists("eq_init") && exists("eq_capture_step") &&
        exists("eq_write_outputs")
}

# Common get parameter method
get_param <- function(params, name, region, use_aggr_fun = "mean", layers = "multiple", na_strategy = "zero") {
  if (is.null(params[[name]]$dim_type) || params[[name]]$dim_type == "spatial") {
    switch(params$region_type,
           raster = {
             param_rast <- get_rast(params, name)
             if (!is.null(param_rast)) {
               param_rast <- conform_layer(param_rast, region$get_template(), na_strategy = na_strategy, use_aggr_fun = use_aggr_fun)
               if (terra::nlyr(param_rast) > 1 && layers == "multiple") {
                 if (length(unique(names(param_rast))) < terra::nlyr(param_rast)) {
                   names(param_rast) <- paste0(names(param_rast), 
                                               1:terra::nlyr(param_rast))
                 }
                 param_series <- lapply(names(param_rast), function(n) {
                   param_rast[n][region$get_indices()][,1]
                 })
                 names(param_series) <- names(param_rast)
                 param_series
               } else {
                 param_rast[region$get_indices()][,1]
               }
             } else {
               NULL
             }
           },
           network = {
             if (!is.null(params[[name]]$csv) && params[[name]]$csv) {
               param_df <- region$get_data() # will need update to bsspread Region
             } else {
               param_df <- get_csv(params, name)
               if (is.data.frame(param_df) && all(c("lon", "lat") %in% names(param_df))) {
                 region_idx <- order(region$get_coords()[,"lon"], region$get_coords()[,"lat"])
                 param_idx <- order(param_df[,"lon"], param_df[,"lat"])
                 param_df <- param_df[param_idx[order(region_idx)],]
               }
             }
             get_df_col(param_df, name)
           },
           none = {
             if ("value" %in% names(params[[name]])) {
               for (n in names(params[[name]])[which(names(params[[name]]) != "value")]) {
                 attr(params[[name]]$value, n) <- params[[name]][[n]]
               }
               params[[name]]$value
             } else {
               params[[name]]
             }
           }
    )
  } else if (params[[name]]$dim_type == "single") {
    params[[name]]$value
  }
}

#' bsmanage::ManageImpacts collection factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsmanage/sim.json#input.properties.impacts
#' @param params A \code{bsspread::Region} or inherited class object
#' @param params A \code{bsspread::PopulationModel} or inherited class object
#' @param params A \code{bsmanage::ManageContext} or inherited class object
bsmanage_impacts_factory <- function(params, region, population_model, context, sim_params) {
    print(paste("bsmanage_impacts_factory:"))
    
    if (length(params)) {

        # Find groupings for impacts based on type, measure, population stages,
        # and recovery delay
        groups = c()
        for (i in 1:length(params)) {
            valuation_type <- params[[i]]$valuation_type
            if (valuation_type == "dynamic") {
                valuation_type <- paste(valuation_type, i)
            }
            if (population_model$get_type() == "stage_structured") {
                params[[i]]$impact_stages <- unlist(params[[i]]$impact_stages)
            } else {
                params[[i]]$impact_stages <- 1
            }
            impact_measure <- params[[i]]$impact_measure
            if (is.null(impact_measure)) {
                impact_measure <- paste("none", i)
                params[[i]]$impact_measure <- ""
            }
            if (is.null(params[[i]]$impact_type)) {
                if (region$spatially_implicit()) {
                    params[[i]]$impact_type <- "area"
                } else {
                    params[[i]]$impact_type <- "presence"
                }
            }
            recovery_delay <- params[[i]]$recovery_delay
            if (is.list(recovery_delay)) {
              recovery_delay <- paste(recovery_delay$value,
                                      recovery_delay$unit)
            } else {
              recovery_delay <- "0"
            }
            groups[i] <- paste(c(valuation_type, impact_measure, params[[i]]$impact_type,
                                 params[[i]]$impact_stages, recovery_delay),
                               collapse = " ")
        }
        
        # Only group when measure, stages, & recovery delay match for each
        # valuation type
        valuation_types = sapply(params, function(p) p$valuation_type)
        for (type in unique(valuation_types)) {
            idx <- which(valuation_types == type)
            if (length(unique(groups[idx])) == 1) {
                groups[idx] <- type
            } else {
                groups[idx] <- paste0(type, "_", idx)
            }
        }
      
        # Collect impacts for unique groupings
        unique_groups <- unique(groups)
        impacts <- list()
        for (i in 1:length(unique_groups)) {
            idx <- which(groups == unique_groups[i])
            impacts[[unique_groups[i]]] <- list(
                region_type = params[[idx[1]]]$region_type,
                valuation_type = params[[idx[1]]]$valuation_type,
                impact_measure = params[[idx[1]]]$impact_measure,
                asset_names = sapply(params[idx], function(p) p$asset_name),
                asset_values = lapply(params[idx], function(p) p$asset_value),
                loss_rates = sapply(params[idx], function(p) p$loss_rate),
                discount_rates = lapply(params[idx], function(p) p$discount_rate),
                impact_stages = params[[idx[1]]]$impact_stages,
                recovery_delay = params[[idx[1]]]$recovery_delay,
                impact_type = params[[idx[1]]]$impact_type,
                dynamic_links = lapply(params[idx], function(p) {
                  list(suitability = if (is.logical(p$dynamic_establish_pr) &&
                                         p$dynamic_establish_pr &&
                                         is.logical(p$link_establish_pr) &&
                                         p$link_establish_pr) TRUE else FALSE,
                       capacity = if (is.logical(p$dynamic_capacity) &&
                                      p$dynamic_capacity &&
                                      is.logical(p$link_capacity) &&
                                      p$link_capacity) TRUE else FALSE,
                       attractors = if (is.logical(p$link_attractors) &&
                                        p$link_attractors) TRUE else FALSE)
                }),
                impact_id = i,
                simulator = sim_params
            )
            impacts[[i]] <- bsmanage_manage_impacts_factory(impacts[[i]], region, population_model, context)
        }
        return(impacts)
    } else {
        return(list())
    }
}

#' bsmanage::ManageImpacts collection factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsmanage/sim.json#input.properties.impacts
#' @param params A \code{bsspread::Region} or inherited class object
#' @param params A \code{bsspread::PopulationModel} or inherited class object
#' @param params A \code{bsmanage::ManageContext} or inherited class object
bsmanage_manage_impacts_factory <- function(params, region, population_model, context) {
    print(paste("bsmanage_manage_impacts_factory:"))
    
    # Set impact context
    impact_types <- context$get_impact_types()
    if (length(impact_types) == 1 && impact_types == "none") { # ignore "none", there must be some, so select all
        impact_types <- c("social", "economic", "ecological")
    }
    if (params$valuation_type == "dynamic") { # assumed non-monetary
        valuation_type <- "non-monetary"
        is_dynamic <- TRUE
    } else {
        valuation_type <- params$valuation_type
        is_dynamic <- FALSE
    }
    impact_context <- bsimpact::Context(species_names = context$get_species_names(),
                                        species_types = context$get_species_types(),
                                        impact_types = impact_types,
                                        impact_scope = params$asset_names,
                                        valuation_type = valuation_type,
                                        impact_measures = params$impact_measure)

    # Set impact region and incursion objects
    if (region$spatially_implicit()) {
        impact_region <- bsimpact::Region()
        incursion <- bsimpact::Incursion(0, region, type = params$impact_type)
    } else {
        impact_region <- region # use bsspread::Region object
        incursion <- bsimpact::Incursion(rep(0, region$get_locations()), region, type = params$impact_type)
    }
    
    # Process asset_values into 'impact_layers'
    impact_layers <- lapply(params$asset_values, function(x) {
        time_mult <- 1
        if (params$valuation_type == "monetary" && !is.null(x$unit) && x$unit == "year") {
            time_mult <- switch(params$simulator$step_units,
                                years = 1,
                                months = 1/12,
                                weeks = 7/365.25,
                                days = 1/365.25)*params$simulator$step_duration
        }
        time_mult*get_param(list(value = x, region_type = params$region_type),
                            "value", region, layers = "single", na_strategy = region$na_strategy)
    })
    names(impact_layers) <- params$asset_names
    
    # Spatially implicit area-based impact values per metres squared
    if (region$spatially_implicit()) {
        for (i in 1:length(impact_layers)) {
            if ("area" %in% names(attributes(impact_layers[[i]]))) {
                if (attr(impact_layers[[i]], "area") == "max") {
                    if (is.numeric(region$get_max_implicit_area())) {
                        impact_layers[[i]] <- as.numeric(impact_layers[[i]])/region$get_max_implicit_area()
                    } else {
                        stop(paste("No maximum area is provided, though impacted", params$asset_names[i], "asset values are configured for maximum area."), call. = FALSE)
                    }
                } else if (attr(impact_layers[[i]], "area") == "km") {
                    impact_layers[[i]] <- as.numeric(impact_layers[[i]])/1000000
                } else {
                    impact_layers[[i]] <- as.numeric(impact_layers[[i]])
                }
            } else {
              impact_layers[[i]] <- as.numeric(impact_layers[[i]])
              print(paste("No area units are given for impacted", params$asset_names[i], "asset values. Assuming metres squared."))
            }
        }
    }

    # Extract loss rates
    loss_rates <- params$loss_rates
    names(loss_rates) <- params$asset_names
    
    # Extract and transform discount rates
    discount_rates <- sapply(params$discount_rates, function(rate) {
        if (!is.null(rate) && !is.null(rate$value)) {
            if (rate$unit == "year") {
                (1 + rate$value)^(
                    switch(params$simulator$step_units,
                           years = 1,
                           months = 1/12,
                           weeks = 7/365.25,
                           days = 1/365.25)*params$simulator$step_duration) - 1
            } else {
                rate$value
            }
        } else { # zero at first
            0
        }
    })
    if (any(discount_rates > 0)) {
        names(discount_rates) <- params$asset_names
    } else { # all zeros
        discount_rates <- NULL
    }

    # Only combine if multiple aspects and impact measure is defined
    if (length(params$asset_names) > 1 && !is.null(params$impact_measure) && params$impact_measure != "") {
        combine_function <- "sum"
    } else {
        combine_function <- "none"
    }
    
    # Value impacts object
    impacts <- bsimpact::ValueImpacts(
        impact_context,
        region = impact_region,
        incursion,
        impact_layers = impact_layers,
        loss_rates = loss_rates,
        discount_rates = discount_rates,
        is_dynamic = is_dynamic,
        combine_function = combine_function)
    impacts$set_id(params$impact_id)
    
    # Recovery delay (in time steps)
    recovery_delay <- NULL
    if (is.list(params$recovery_delay) && is.numeric(params$recovery_delay$value)) {
        recovery_delay <- params$recovery_delay$value
        if (!is.null(params$recovery_delay$unit) && params$recovery_delay$unit == "year") {
            recovery_delay <- recovery_delay*switch(params$simulator$step_units,
                                                    years = 1,
                                                    months = 12,
                                                    weeks = 365.25/7,
                                                    days = 365.25)/params$simulator$step_duration
        }
    }
    
    # Resolve dynamic links
    dynamic_links <- unlist(params$dynamic_links)
    if (is.logical(dynamic_links) && length(which(dynamic_links))) {
      dynamic_links <- unique(names(dynamic_links[which(dynamic_links)]))
    } else {
      dynamic_links <- NULL
    }

    return(ManageImpacts(impacts,
                         population_model,
                         impact_stages = params$impact_stages,
                         calc_total = TRUE,
                         dynamic_links = dynamic_links,
                         recovery_delay = recovery_delay))
}

#' actions collection factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsmanage/sim.json#input.properties.actions
#' @param params A \code{bsspread::Region} or inherited class object
#' @param params A \code{bsspread::PopulationModel} or inherited class object
#' @param params A \code{bsmanage::ManageContext} or inherited class object
bsmanage_actions_factory <- function(params, region, population_model, context) {
    print(paste("bsmanage_actions_factory:"))
  
    # Create list of actions
    actions = list()
    if (length(params)){
        for (i in 1:length(params)) {
            
            # Transform stages and schedules into vectors
            if (population_model$get_type() == "stage_structured") {
                params[[i]]$stages <- unlist(params[[i]]$stages)
            }
            if (!is.null(params[[i]]$schedule)) {
                params[[i]]$schedule <- unlist(params[[i]]$schedule)
            }
            
            # Generate manage detection, control, or removal action objects
            actions[[i]] <-
                switch(params[[i]]$action_type,
                    detection = bsmanage_detection_factory(params[[i]], region, population_model, context),
                    control = bsmanage_control_factory(params[[i]], region, population_model, context),
                    removal = bsmanage_removal_factory(params[[i]], region, population_model, context)
                )
        }
      
      # Set action ids (allows duplicate action types)
      for (i in 1:length(actions)) {
        actions[[i]]$set_id(i)
      }
    }
    
    # # Ensure detection actions precede other actions that apply to detected
    # # (allows users to place actions in unexpected order)
    # action_labels <- sapply(actions, function(a) a$get_label(include_id = FALSE))
    # if ("detected" %in% action_labels && any(!action_labels %in% c("detected", "control_search_destroy"))) {
    #   if (which("detected" == action_labels)[1] > which(!action_labels %in% c("detected", "control_search_destroy"))[1]) {
    #     # Place detections first
    #     actions <- actions[c(which(action_labels == "detected"), which(action_labels != "detected"))]
    #   }
    # }

    return (actions)
}

#' bsmanage::ManageDetection factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsmanage/sim.json#input.properties.actions
#' @param params A \code{bsspread::Region} or inherited class object
#' @param params A \code{bsspread::PopulationModel} or inherited class object
#' @param params A \code{bsmanage::ManageContext} or inherited class object
bsmanage_detection_factory <- function(params, region, population_model, context) {
    print("bsmanage_detection_factory")
    
    # Place detection sensitivity into a surveillance design object
    sensitivity <- get_param(params, "sensitivity", region, layers = "single", na_strategy = region$na_strategy)
    surv_context <- bsdesign::Context(species_name = paste(context$get_species_names(), collapse = " "),
                                      species_type = context$get_species_types()[1])
    divisions <- bsdesign::Divisions(matrix(0, c(region$get_locations(), 1)))
    surveillance <- bsdesign::SurveillanceDesign(surv_context,
                                                 divisions,
                                                 optimal = "none",
                                                 exist_sens = sensitivity)
    
    # Surveillance cost & unit
    surv_cost <- get_param(params, "surv_cost", region, layers = "single", na_strategy = region$na_strategy)
    if (!is.null(surv_cost) && !is.null(params$surv_cost$unit)) {
      attr(surv_cost, "unit") <- params$surv_cost$unit
    }
    
    # Sensitivity type and threshold
    if (!is.null(params$sensitivity_type) && params$sensitivity_type %in% c("individual", "population", "presence")) {
      sensitivity_type <- params$sensitivity_type
    } else { # default
      sensitivity_type <- "individual"
    }
    if (sensitivity_type == "population") {
      sensitivity_threshold <- params$sensitivity_threshold
    } else {
      sensitivity_threshold <- NULL
    }
    
    return(ManageDetection(region,
                           population_model,
                           surveillance,
                           sensitivity_type = sensitivity_type,
                           sensitivity_threshold = sensitivity_threshold,
                           surv_cost = surv_cost,
                           stages = params$stages,
                           schedule = params$schedule))
}

#' bsmanage::ManageControls factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsmanage/sim.json#input.properties.actions
#' @param params A \code{bsspread::Region} or inherited class object
#' @param params A \code{bsspread::PopulationModel} or inherited class object
#' @param params A \code{bsmanage::ManageContext} or inherited class object
bsmanage_control_factory <- function(params, region, population_model, context) {
    print(paste("bsmanage_control_factory:", params$control_type))
    
    
    # Divisions for control design object
    divisions <- bsdesign::Divisions(matrix(0, c(region$get_locations(), 1)))
    
    # Collect parameters for control
    if (params$control_type == "search_destroy") {
        
        manage_pr <- get_param(params, "manage_pr", region, layers = "single", na_strategy = region$na_strategy)
        control_design <- ManageDesign(context,
                                       divisions,
                                       dim_type = "spatial",
                                       optimal = "none",
                                       exist_manage_pr = manage_pr)
        
        if (!is.null(params$manage_pr_type) && params$manage_pr_type %in% c("individual", "population")) {
          manage_pr_type <- params$manage_pr_type
        } else { # default
          manage_pr_type <- "individual"
        }
        control_mult <- NULL
        apply_to <- NULL
        
    } else if (params$control_type %in% c("growth", "spread", "establishment")) {
        
        # Existing control passed via a manage design object
        exist_control <- get_param(params, "exist_control", region, layers = "single", na_strategy = region$na_strategy)
        if (!is.null(exist_control)) {
            control_design <- ManageDesign(context,
                                           divisions,
                                           dim_type = "spatial",
                                           optimal = "none",
                                           exist_alloc = +(exist_control > 0))
        } else {
            control_design <- NULL
        }
        manage_pr_type <- "individual" # not used
        
        # Control multiplier
        control_mult <- get_param(params, "control_mult", region, layers = "single", na_strategy = region$na_strategy)
        
        # Apply stage-based growth multiplier applied to reproduction or survival
        if (population_model$get_type() == "stage_structured" && params$control_type == "growth" &&
            !is.null(params$apply_to) && params$apply_to %in% c("reproduction", "survival")) {
            apply_to <- params$apply_to
        } else {
            apply_to <- NULL
        }
    }
    
    # Control cost & unit
    control_cost <- get_param(params, "control_cost", region, layers = "single", na_strategy = region$na_strategy)
    if (!is.null(control_cost) && !is.null(params$control_cost$unit)) {
      attr(control_cost, "unit") <- params$control_cost$unit
    }
    
    # Spatially implicit area-based cost values per metres squared for treatment controls
    if (region$spatially_implicit() && !is.null(control_cost) &&
        params$control_type %in% c("growth", "spread", "establishment")) {
        if ("area" %in% names(attributes(control_cost))) {
          if (attr(control_cost, "area") == "max") {
            if (is.numeric(region$get_max_implicit_area())) {
              control_cost <- control_cost/region$get_max_implicit_area()
            } else {
              stop(paste("No maximum area is provided, though", params$control_type, "control is configured for maximum area."), call. = FALSE)
            }
          } else if (attr(control_cost, "area") == "km") {
            control_cost <- control_cost/1000000
          }
        } else {
          print(paste("No area units are given for", params$control_type, "control. Assuming metres squared."))
        }
        attr(control_cost, "area") <- NULL
    }

    return(ManageControls(region,
                          population_model,
                          control_type = params$control_type,
                          control_design = control_design,
                          manage_pr_type = manage_pr_type,
                          control_mult = control_mult,
                          control_cost = control_cost,
                          radius = params$radius,
                          stages = params$stages,
                          apply_to = apply_to,
                          schedule = params$schedule))
}

#' bsmanage::ManageRemovals factory
#'
#' @param params A JSON object described by 
#'               https://biosecuritycommons.org.au/bsmanage/sim.json#input.properties.actions
#' @param params A \code{bsspread::Region} or inherited class object
#' @param params A \code{bsspread::PopulationModel} or inherited class object
#' @param params A \code{bsmanage::ManageContext} or inherited class object
bsmanage_removal_factory <- function(params, region, population_model, context) {
    print("bsmanage_removal_factory")
    removal_pr <- get_param(params, "removal_pr", region, layers = "single", na_strategy = region$na_strategy)
    if (!is.null(params$removal_pr_type) && params$removal_pr_type %in% c("individual", "population")) {
      removal_pr_type <- params$removal_pr_type
    } else { # default
      removal_pr_type <- "individual"
    }
    remove_always <- FALSE
    if (is.logical(params$remove_always)) {
        remove_always <- params$remove_always
    }
    detected_only <- FALSE
    if (is.logical(params$detected_only)) {
        detected_only <- params$detected_only
    }
    
    # Removal cost & unit
    removal_cost <- get_param(params, "removal_cost", region, layers = "single", na_strategy = region$na_strategy)
    if (!is.null(removal_cost) && !is.null(params$removal_cost$unit)) {
      attr(removal_cost, "unit") <- params$removal_cost$unit
    }
    
    # Spatially implicit area-based cost values per metres squared
    if (region$spatially_implicit() && !is.null(removal_cost)) {
      if ("area" %in% names(attributes(removal_cost))) {
        if (attr(removal_cost, "area") == "max") {
          if (is.numeric(region$get_max_implicit_area())) {
            removal_cost <- removal_cost/region$get_max_implicit_area()
          } else {
            stop("No maximum area is provided, though removal is configured for maximum area.", call. = FALSE)
          }
        } else if (attr(removal_cost, "area") == "km") {
          removal_cost <- removal_cost/1000000
        }
      } else {
        print("No area units are given for removal. Assuming metres squared.")
      }
      attr(removal_cost, "area") <- NULL
    }

    return(ManageRemovals(region,
                          population_model,
                          removal_pr,
                          removal_pr_type = removal_pr_type,
                          remove_always = remove_always,
                          detected_only = detected_only,
                          removal_cost = removal_cost,
                          radius = params$radius,
                          stages = params$stages,
                          schedule = params$schedule))
}


#' Generate animations from result layers
#'
#' @param region A \code{bsspread::Region} or inherited class object
#' @param population_model A \code{bsspread::PopulationModel} or inherited class object
#' @param result_layers A list of raster layers to animate
generate_result_animations <- function(region, population_model, result_layers) {
  tryCatch({
    # Process habitat suitability raster
    hs <- region$get_rast(population_model$get_establish_pr()) |>
      terra::crop(result_layers[[1]]) |>
      terra::mask(result_layers[[1]][[1]])
    
    hs_gt0 <- hs > 0
    
    # Set tmap options
    tmap::tmap_options(raster.max_cells = 20e6)
  }, error = function(e) {
    warning(
      sprintf(
        "Failed to initialise animation setup: %s. Animations skipped.", 
        e$message
      ),
      call. = FALSE
    )
    return(invisible(NULL))
  })
  
  for (layer in names(result_layers)) {
    tryCatch({
      r <- result_layers[[layer]]
      metadata <- attr(r, "metadata")
      r[r == 0] <- NA
      
      message(sprintf("Animating %s", metadata$label))
      
      suppressWarnings(ranges <- unlist(terra::global(r, range, na.rm = TRUE)))
      suppressWarnings(zlim <- range(ranges[!is.infinite(ranges)])) # this handles NaN occupancy sd at t = 0
      no_values <- all(is.infinite(zlim))
      is_constant <- no_values || length(unique(zlim)) == 1
      
      if(no_values) {
        message(sprintf("No layers have values for layer %s.", layer))
      }
      
      names(r) <- sprintf(
        sprintf(
          "%%s\nt = %% %ss",
          max(nchar(names(r))) - nchar(names(r)) + 1
        ),
        metadata$label,
        as.integer(names(r))
      )
      
      scl <- switch(
        metadata$scale_type,
        
        continuous_log10 = {
          zlim_log10 <- c(floor(log10(zlim[1])), ceiling(log10(zlim[2])))
          tmap::tm_scale_continuous_log10(
            limits = 10^zlim_log10,
            ticks = 10^seq.int(zlim_log10[1], zlim_log10[2]),
            label.format = scales::label_number(
              scale_cut = scales::cut_short_scale()
            )
          )
        },
        
        continuous = {
          tmap::tm_scale_continuous(limits = zlim)
        },
        
        percent = {
          tmap::tm_scale_continuous(
            limits = zlim,
            label.format = scales::label_percent()
          ) 
        },
        
        # default scale
        tmap::tm_scale_continuous(limits = zlim) 
      )
      
      if(
        (metadata$category == "occupancy" && metadata$summary == "mean") && 
        !is_constant
      ) {
        scl$ticks <- seq(0, 1, 0.2)
        scl$limits <- c(0, 1)
      }
      
      scl$values <- "matplotlib.rainbow"
      
      # Override scl with discrete scale if values are constant across time 
      # and space
      if(is_constant) {
        scl <- tmap::tm_scale_discrete(
          ticks = unique(zlim),
          label.format = scl$label.format,
          values = "firebrick"
        )
      }
      
      p <- tmap::tm_shape(hs_gt0) +
        tmap::tm_raster(
          col.legend = tmap::tm_legend_hide(),
          col.scale = tmap::tm_scale_categorical(
            levels.drop = FALSE, 
            levels = c(FALSE, TRUE),
            values = c("gray55", "gray85"),
            labels = c("Unsuitable", "Unoccupied")
          )
        ) +
        tmap::tm_shape(r) +
        tmap::tm_raster(
          col.legend = tmap::tm_legend(
            title = "",  reverse = TRUE, 
            frame = FALSE, text.size = 1
          ),
          col.scale = scl,
          col.free = FALSE
        ) +
        tmap::tm_layout(
          panel.show = TRUE,
          panel.label.fontface = "bold",
          panel.label.size = 1.5,
          panel.label.height = 3.5,
          panel.label.bg.color = "white",
          panel.label.frame = FALSE,
          legend.position = tmap::tm_pos_out("right", "center", pos.v = "center")
        )
      
      if(metadata$summary != "sd") {
        p <- p + tmap::tm_add_legend( # so that this legend appears second
          labels = c("Unoccupied", "Unsuitable"),
          fill = c("gray85", "gray55"), type = "polygons",
          text.size = 1
        ) 
      }
      
      p <- p + tmap::tm_animate(fps = 1, play = "once")
      
      asp <- tmaptools::get_asp_ratio(p) 
      # ^ can potentially optimise this step by checking a single-layer map
      
      tmap::tmap_animation(
        p,
        filename = paste0(layer, ".mp4"),
        height = if(asp < 0) 1440 else NA,
        width = if(asp >= 0) 1440 else NA,
        vfilter = "null",
        scale = 1.5,
        audio = NULL,
        verbose = FALSE
      )
    }, error = function(e) {
      warning(
        sprintf(
          "Failed to generate animation for layer '%s': %s. Continuing with other layers.", 
          layer, e$message
        ),
        call. = FALSE
      )
    })
  }
}

# Checkpoint helpers (paths live in Region's closure, not on the region list)
get_region_paths <- function(region) {
    environment(region$calculate_paths)$paths
}

set_region_paths <- function(region, paths) {
    environment(region$calculate_paths)$paths <- paths
}

# Terra/permeability pointers in paths do not survive saveRDS; drop before save.
prepare_paths_for_checkpoint <- function(paths) {
    paths$perms <- NULL
    if (is.list(paths$graphs) && !is.null(paths$graphs$poly)) {
        paths$graphs$poly <- NULL
    }
    paths
}

restore_region_paths <- function(region, saved_paths) {
    env <- environment(region$calculate_paths)
    live_perms <- env$paths$perms
    set_region_paths(region, saved_paths)
    if (!is.null(live_perms)) {
        environment(region$calculate_paths)$paths$perms <- live_perms
    }
}

default_checkpoint_file <- function(tm, r, outputdir) {
    file.path(outputdir, sprintf("checkpoint_tm%d_r%d.rds", tm, r))
}

checkpoint_file_in <- function() {
    explicit <- Sys.getenv("CHECKPOINT_FILE_IN", unset = "")
    if (nzchar(explicit)) {
        return(explicit)
    }
    Sys.getenv("CHECKPOINT_FILE", unset = "")
}

checkpoint_file_out <- function(tm, r, outputdir) {
    explicit <- Sys.getenv("CHECKPOINT_FILE_OUT", unset = "")
    if (nzchar(explicit)) {
        return(explicit)
    }
    default_checkpoint_file(tm, r, outputdir)
}

save_simulation_checkpoint <- function(file, tm, r, n, region,
                                       calc_impacts = NULL) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    payload <- list(
        rng = .Random.seed,
        n = n,
        paths = prepare_paths_for_checkpoint(get_region_paths(region)),
        tm_saved = tm,
        tm_next = tm + 1L,
        r = r,
        calc_impacts = calc_impacts
    )
    saveRDS(payload, file)
    message(sprintf(
        "Wrote checkpoint (tm_next=%d, r=%d, %d path origins) -> %s",
        payload$tm_next,
        r,
        length(payload$paths$idx),
        file
    ))
}

load_simulation_checkpoint <- function(file) {
    if (!file.exists(file)) {
        stop("Checkpoint file not found: ", file, call. = FALSE)
    }
    bench <- readRDS(file)
    required <- c("rng", "n", "paths", "tm_next", "r")
    missing <- setdiff(required, names(bench))
    if (length(missing)) {
        stop("Checkpoint missing fields: ", paste(missing, collapse = ", "),
            call. = FALSE)
    }
    tm_saved_label <- if (is.null(bench$tm_saved)) "?" else bench$tm_saved
    message(sprintf(
        "Loaded checkpoint (saved tm=%s, resume tm_next=%d, r=%d, %d path origins)",
        tm_saved_label,
        bench$tm_next,
        bench$r,
        length(bench$paths$idx)
    ))
    bench
}

# Deep copy for repeated benchmark reps without re-readRDS (serialize round-trip).
duplicate_simulation_checkpoint <- function(bench) {
    unserialize(serialize(bench, NULL))
}

duplicate_checkpoint_paths <- function(paths) {
    unserialize(serialize(paths, NULL))
}

# Reset state from an immutable scaling master (paths deep-copied; n/rng from master).
apply_scaling_checkpoint_restore <- function(region, population_model, master) {
    .Random.seed <- master$rng
    restore_region_paths(region, duplicate_checkpoint_paths(master$paths))
    population_model$set_capacity_mult(master$n)
    n <<- master$n
    master$calc_impacts
}

# Create objects
if (!is.null(input.params$context) && !is.null(input.params$context$impacts) && input.params$context$impacts$include) {
    impact_types <- unlist(input.params$context$impacts[c("social", "economic", "ecological")])
    if (any(impact_types)) {
        impact_types <- names(impact_types[which(impact_types)])
    } else {
        impact_types <- "none"
    }
} else {
    impact_types <- "none"
}
if (!is.null(input.params$context) && !is.null(input.params$context$actions) && input.params$context$actions$include) {
    action_types <- unlist(input.params$context$actions[c("surveillance", "control", "removal")])
    if (any(action_types)) {
        action_types <- names(action_types[which(action_types)])
    } else {
        action_types <- "none"
    }
} else {
    action_types <- "none"
}

context <- bsmanage::ManageContext(
    species_names = "not available",  
    species_types = c("pest", "weed", "disease"),
    impact_types = impact_types,
    action_types = action_types,
    manage_scope = "scenarios"
)

region <- bsspread_region_factory(input.params$region)

population_model <- bsspread_populaton_model_factory(
    input.params$population_model,
    region
)

initializer <- bsspread_initializer_factory(
    input.params$initializer,
    region,
    population_model,
    input.env$outputdir
)

dispersal_models <- bsspread_dispersal_models_factory(
    input.params$dispersal_models,
    region,
    population_model,
    impact_params = input.params$impacts
)

impacts <- bsmanage_impacts_factory(
    input.params$impacts,
    region,
    population_model,
    context, 
    sim_params = input.params$simulator
)

actions <- bsmanage_actions_factory(
    input.params$actions,
    region,
    population_model,
    context
)

# Progress weblog
if (JOB_WEBLOG_URL != '' && JOB_WEBLOG_AUTH != ''){
    progress_weblog_enabled <- TRUE
    message('Using progress webhook ', JOB_WEBLOG_URL)
} else {
    progress_weblog_enabled <- FALSE
}

weblog <- function(log) {
    message(log)
    if (progress_weblog_enabled){
        tryCatch(
            expr = {
                auth = strsplit(JOB_WEBLOG_AUTH, split = ':')[[1]]
                weblog <- list(
                    event="progress",
                    metadata=list(message=log)
                )
                req <- httr2::request(JOB_WEBLOG_URL)
                req |>
                    httr2::req_body_json(weblog) |>
                    httr2::req_auth_basic(auth[1], auth[2]) |>
                    httr2::req_throttle(capacity = 5, fill_time_s = 60) |>
                    httr2::req_retry(max_tries = 10, retry_on_failure = TRUE) |>
                    httr2::req_perform()
            },
            error = function(e){          
                message(e)
                progress_weblog_enabled <<- FALSE
            }
        )
    }
}

apply_rate_limit <- function(func, n_seconds) {
  # Create a closure that stores the time of the last hit internally
  # This decouples logging from execution frequency.
  last_hit <- Sys.time() - n_seconds
  
  function(...) {
    time_since_last_hit <- as.numeric(
      difftime(Sys.time(), last_hit, units = "secs")
    )
    if (time_since_last_hit >= n_seconds) {
      last_hit <<- Sys.time()
      func(...)
    } 
  }
}

rate_limited_weblog <- apply_rate_limit(weblog, n_seconds = 30)

progress_function <- function(n, r, tm) {
    time_steps <- input.params$simulator$time_steps
    replicates <- input.params$simulator$replicates
    progress_pct <-
        (((r - 1) * time_steps + tm) / (replicates * time_steps)) * 100
    rate_limited_weblog(sprintf(
        "Completed timestep %d of replicate %d (%.2f%%)",
        tm, r, progress_pct))
    return(n)
}

# Use serial (not parallel) when there is a random seed (unless benchmark override)
if (is.numeric(input.params$random_seed)) {
    set.seed(input.params$random_seed)
    if (identical(Sys.getenv("BENCHMARK_ALLOW_PARALLEL"), "1")) {
        message(
            "random_seed set but BENCHMARK_ALLOW_PARALLEL=1: using TASK_CPUS/DEFAULT_CPUS"
        )
    } else {
        cpus <- 1
        print("Using seed, this will run without parallel threads!")
    }
}

simulator <- bsmanage::ManageSimulator(
    region = region,
    time_steps = input.params$simulator$time_steps,
    step_units = input.params$simulator$step_units,
    step_duration = input.params$simulator$step_duration,
    collation_steps = input.params$simulator$collation_steps,
    replicates = input.params$simulator$replicates,
    result_stages = unlist(input.params$simulator$result_stages),
    parallel_cores = cpus,
    initializer = initializer,
    population_model = population_model,
    dispersal_models = dispersal_models,
    impacts = impacts,
    actions = actions,
    user_function = progress_function
)

#### DEBUG

list2env(as.list(environment(simulator$run)), envir = .GlobalEnv)

object_absence <- c(is.null(initializer), is.null(population_model))
if (any(object_absence)) {
    stop(sprintf("The simulator requires the %s to be set.", 
        paste(c("initializer", "population model")[object_absence], 
            collapse = " and ")), call. = FALSE)
}
continued_incursions <- initializer$continued_incursions()
results <<- ManageResults(region, population_model, impacts = impacts, 
    actions = actions, time_steps = time_steps, step_duration = step_duration, 
    step_units = step_units, collation_steps = collation_steps, 
    replicates = replicates, combine_stages = result_stages)

if (eq_enabled) {
    eq_run_label <- Sys.getenv("EQ_RUN_LABEL")
    if (!nzchar(eq_run_label)) {
        eq_run_label <- paste0(
            "run_",
            if (!is.null(input.params$job_uuid)) input.params$job_uuid else
                format(Sys.time(), "%Y%m%d_%H%M%S")
        )
    }
    eq_git_sha <- Sys.getenv("EQ_GIT_SHA")
    if (!nzchar(eq_git_sha)) eq_git_sha <- NA_character_
    eq_branch <- Sys.getenv("EQ_BRANCH")
    if (!nzchar(eq_branch)) eq_branch <- NA_character_
    eq_scenario <- Sys.getenv("EQ_SCENARIO")
    if (!nzchar(eq_scenario)) eq_scenario <- "manual_docker_run"
    eq_out_dir <- Sys.getenv("EQ_OUT_DIR")
    if (!nzchar(eq_out_dir)) {
        eq_out_dir <- file.path(input.env$outputdir, "equivalence")
    }
    eq_state <- eq_init(
        run_label = eq_run_label,
        git_sha = eq_git_sha,
        branch = eq_branch,
        scenario = eq_scenario,
        out_dir = eq_out_dir,
        save_occupancy = TRUE
    )
}

checkpoint_data <- NULL
if (checkpoint_resume) {
    ckpt_path <- checkpoint_file_in()
    if (!nzchar(ckpt_path)) {
        stop(
            "CHECKPOINT_RESUME=1 requires CHECKPOINT_FILE_IN or CHECKPOINT_FILE",
            call. = FALSE
        )
    }
    checkpoint_data <- load_simulation_checkpoint(ckpt_path)
    if (as.integer(Sys.getenv("BENCHMARK_SCALING_REPS", unset = "0")) > 0L) {
        checkpoint_data <- duplicate_simulation_checkpoint(checkpoint_data)
    }
    .Random.seed <- checkpoint_data$rng
    restore_region_paths(region, checkpoint_data$paths)
    population_model$set_capacity_mult(checkpoint_data$n)
    message("Checkpoint: region path cache restored (fresh permeability rasters)")
}

gc_time_prev <- sum(gc.time())

benchmark_replicate <- Sys.getenv("BENCHMARK_REPLICATE", unset = "")
checkpoint_warmup_reps <- Sys.getenv("CHECKPOINT_WARMUP_REPS", unset = "")
if (nzchar(benchmark_replicate)) {
    replicate_seq <- as.integer(benchmark_replicate)
} else if (nzchar(checkpoint_warmup_reps)) {
    warmup_reps <- as.integer(strsplit(
        checkpoint_warmup_reps, ",", fixed = TRUE)[[1]])
    replicate_seq <- c(warmup_reps, checkpoint_replicate)
    message(sprintf(
        paste0(
            "Checkpoint create: warmup replicate(s) %s (shared region path cache), ",
            "then replicate %d"
        ),
        paste(warmup_reps, collapse = ", "),
        checkpoint_replicate
    ))
} else if (checkpoint_resume && !is.null(checkpoint_data)) {
    # Resume only the replicate stored in the checkpoint (e.g. tm=29 r=2 -> run r=2)
    replicate_seq <- checkpoint_data$r
    message(sprintf(
        "Checkpoint resume: running replicate %d only (set BENCHMARK_REPLICATE to override)",
        checkpoint_data$r
    ))
} else if (benchmark_parallel_replicates) {
    replicate_seq <- if (nzchar(benchmark_replicate)) {
        as.integer(benchmark_replicate)
    } else {
        seq_len(input.params$simulator$replicates)
    }
} else {
    replicate_seq <- seq_len(2)
}
benchmark_time_end <- Sys.getenv("BENCHMARK_TIME_END", unset = "")
tm_run_end <- if (nzchar(benchmark_time_end)) {
    as.integer(benchmark_time_end)
} else {
    time_steps
}
checkpoint_warmup_time_end <- Sys.getenv("CHECKPOINT_WARMUP_TIME_END", unset = "")
warmup_tm_end <- if (nzchar(checkpoint_warmup_time_end)) {
    as.integer(checkpoint_warmup_time_end)
} else {
    time_steps
}
warmup_reps_vec <- if (nzchar(checkpoint_warmup_reps)) {
    as.integer(strsplit(checkpoint_warmup_reps, ",", fixed = TRUE)[[1]])
} else {
    integer(0)
}
if (length(warmup_reps_vec) && nzchar(checkpoint_warmup_reps)) {
    message(sprintf(
        "Checkpoint create: warmup replicate(s) run through tm=%d; checkpoint replicate through tm=%d",
        warmup_tm_end,
        tm_run_end
    ))
}

force_serial_inner_parallel <- function() {
    region$set_cores(1)
    if (length(dispersal_models)) {
        for (i in seq_along(dispersal_models)) {
            dispersal_models[[i]]$set_cores(1)
        }
    }
}

run_one_replicate_parallel <- function(r) {
    collations <- list()
    if (is.numeric(input.params$random_seed)) {
        if (benchmark_reseed_per_replicate) {
            set.seed(input.params$random_seed)
            message(sprintf(
                "BENCHMARK_RESEED_PER_REPLICATE: set.seed for replicate %d (pid %d)",
                r, Sys.getpid()))
        } else {
            # Forked workers inherit identical RNG; per-rep seed (not serial chaining).
            rep_seed <- input.params$random_seed + as.integer(r) - 1L
            set.seed(rep_seed)
            message(sprintf(
                "Parallel replicate %d: set.seed(%d) (pid %d)",
                r, rep_seed, Sys.getpid()))
        }
    }
    force_serial_inner_parallel()

    n <- initializer$initialize()
    if (region$spatially_implicit()) {
        if (any(sapply(dispersal_models, function(dm) inherits(dm,
            "Diffusion")))) {
            idx <- which(sapply(dispersal_models, function(dm) inherits(dm,
                "Diffusion")))[1]
            attr(n, "initial_n") <- n
            attr(n, "diffusion_rate") <-
                dispersal_models[[idx]]$get_diffusion_rate()
            attr(n, "diffusion_radius") <- 0
        }
        if (any(sapply(dispersal_models, function(dm) inherits(dm,
            "AreaSpread")))) {
            capacity <- population_model$get_capacity()
            capacity_area <- attr(capacity, "area")
            if (population_model$get_type() == "stage_structured") {
                stages <- population_model$get_capacity_stages()
                attr(n, "spread_area") <- sum(n[, stages]) /
                    as.numeric(capacity) * capacity_area
            } else {
                attr(n, "spread_area") <- n / as.numeric(capacity) *
                    capacity_area
            }
        }
    }
    if (length(impacts)) {
        calc_impacts <- vector("list", length(impacts))
        for (i in seq_along(impacts)) {
            n <- impacts[[i]]$calculate(n, 0)
            calc_impacts[[i]] <- attr(n, "impacts")
        }
        attr(n, "impacts") <- NULL
        population_model$set_capacity_mult(n)
    } else {
        calc_impacts <- NULL
    }
    if (length(actions)) {
        for (i in seq_along(actions)) {
            n <- actions[[i]]$apply(n, 0)
        }
    }
    collations[[1L]] <- list(r = r, tm = 0L, n = n, calc_impacts = calc_impacts)

    gc_time_prev <- sum(gc.time())
    for (tm in seq_len(tm_run_end)) {
        t0 <- Sys.time()
        n <- population_model$grow(n, tm)
        t1 <- Sys.time()

        if (length(dispersal_models)) {
            n <- dispersal_models[[1]]$pack(n)
            for (i in seq_along(dispersal_models)) {
                t_dm <- Sys.time()
                n <- dispersal_models[[i]]$disperse(n, tm)
                message(sprintf(
                    "    dispersal model %d: %.2fs (r=%d pid %d)",
                    i,
                    as.numeric(Sys.time() - t_dm, units = "secs"),
                    r, Sys.getpid()))
            }
            n <- dispersal_models[[1]]$unpack(n)
        }
        t2 <- Sys.time()

        if (length(impacts)) {
            calc_impacts <- vector("list", length(impacts))
            for (i in seq_along(impacts)) {
                n <- impacts[[i]]$calculate(n, tm)
                calc_impacts[[i]] <- attr(n, "impacts")
            }
            attr(n, "impacts") <- NULL
            population_model$set_capacity_mult(n)
        }
        t3 <- Sys.time()

        if (length(actions)) {
            for (i in seq_along(actions)) {
                n <- actions[[i]]$clear_attributes(n)
            }
            for (i in seq_along(actions)) {
                n <- actions[[i]]$apply(n, tm)
            }
        }
        t4 <- Sys.time()

        if (is.function(user_function)) {
            n_attr <- attributes(n)
            if (length(formals(user_function)) == 3) {
                n <- user_function(n, r, tm)
            } else {
                n <- user_function(n)
            }
            if (length(n_attr)) {
                for (i in seq_along(n_attr)) {
                    if (!names(n_attr[i]) %in% names(attributes(n))) {
                        attr(n, names(n_attr[i])) <- n_attr[[i]]
                    }
                }
            }
        }
        if (is.function(continued_incursions)) {
            n <- continued_incursions(tm, n)
        }
        t5 <- Sys.time()

        collations[[length(collations) + 1L]] <- list(
            r = r, tm = tm, n = n, calc_impacts = calc_impacts)

        elapsed <- function(a, b) as.numeric(b - a, units = "secs")
        mem_info <- ps::ps_system_memory()
        proc_mem <- ps::ps_memory_info(ps::ps_handle(Sys.getpid()))
        gc_time_now <- sum(gc.time())
        gc_step_s <- gc_time_now - gc_time_prev
        gc_time_prev <- gc_time_now
        rng_check <- runif(1)
        message(sprintf(
            paste0(
                "tm=%s r=%s | grow=%.2fs  dispersal=%.2fs  ",
                "impacts=%.2fs  actions=%.2fs  other=%.2fs | total=%.2fs | ",
                "%s sys avail | %s rss | gc=%.3fs | rng check=%.4f | pid %d"
            ),
            tm, r,
            elapsed(t0, t1),
            elapsed(t1, t2),
            elapsed(t2, t3),
            elapsed(t3, t4),
            elapsed(t4, t5),
            elapsed(t0, t5),
            prettyunits::pretty_bytes(mem_info$avail),
            prettyunits::pretty_bytes(proc_mem["rss"]),
            gc_step_s,
            rng_check,
            Sys.getpid()
        ))
    }
    list(r = r, collations = collations)
}

if (benchmark_parallel_replicates) {
    if (.Platform$OS.type != "unix") {
        stop("BENCHMARK_PARALLEL_REPLICATES requires Unix (mclapply)",
            call. = FALSE)
    }
    if (checkpoint_resume || checkpoint_tm > 0L || nzchar(checkpoint_warmup_reps)) {
        stop(
            "BENCHMARK_PARALLEL_REPLICATES is incompatible with checkpoint ",
            "resume/create (use serial replicate loop)",
            call. = FALSE)
    }
    if (as.integer(Sys.getenv("BENCHMARK_SCALING_REPS", unset = "0")) > 0L) {
        stop(
            "BENCHMARK_PARALLEL_REPLICATES is incompatible with ",
            "BENCHMARK_SCALING_REPS",
            call. = FALSE)
    }
    if (isTRUE(profile_timestep)) {
        stop(
            "BENCHMARK_PARALLEL_REPLICATES is incompatible with profile_timestep",
            call. = FALSE)
    }
    if (eq_enabled) {
        message(
            "BENCHMARK_PARALLEL_REPLICATES: equivalence capture disabled ",
            "(per-worker merge not implemented)")
        eq_enabled <- FALSE
    }
    force_serial_inner_parallel()
    parallel_rep_cores <- min(as.integer(cpus), length(replicate_seq))
    message(sprintf(
        paste0(
            "BENCHMARK_PARALLEL_REPLICATES: %d replicate(s), ",
            "mc.cores=%d, inner dispersal/paths serial"
        ),
        length(replicate_seq), parallel_rep_cores))
    rep_wall_start <- Sys.time()
    rep_outputs <- parallel::mclapply(
        replicate_seq,
        run_one_replicate_parallel,
        mc.cores = parallel_rep_cores)
    rep_wall_s <- as.numeric(Sys.time() - rep_wall_start, units = "secs")
    for (out in rep_outputs) {
        for (col in out$collations) {
            results$collate(col$r, col$tm, col$n, col$calc_impacts)
        }
    }
    message(sprintf(
        "BENCHMARK_PARALLEL_REPLICATES complete: wall=%.2fs (%d reps, mc.cores=%d)",
        rep_wall_s, length(replicate_seq), parallel_rep_cores))
}

if (!benchmark_parallel_replicates) {
for (r in replicate_seq) {
    if (identical(Sys.getenv("PATHS_EQ_DUMP", unset = ""), "1")) {
        Sys.setenv(PATHS_EQ_REPLICATE = as.character(r))
    }

    resume_this_r <- !is.null(checkpoint_data) && r == checkpoint_data$r
    tm_start <- if (resume_this_r) checkpoint_data$tm_next else 1L

    if (resume_this_r) {
        n <- checkpoint_data$n
        calc_impacts <- checkpoint_data$calc_impacts
        message(sprintf(
            "Checkpoint resume: replicate %d from tm=%d (skipping initialise and tm 0:%d)",
            r, tm_start, tm_start - 1L))
    } else {
    if (benchmark_reseed_per_replicate &&
        is.numeric(input.params$random_seed)) {
        set.seed(input.params$random_seed)
        message(sprintf(
            "BENCHMARK_RESEED_PER_REPLICATE: set.seed for replicate %d",
            r))
    }
    n <- initializer$initialize()
    if (region$spatially_implicit()) {
        if (any(sapply(dispersal_models, function(dm) inherits(dm, 
            "Diffusion")))) {
            idx <- which(sapply(dispersal_models, function(dm) inherits(dm, 
                "Diffusion")))[1]
            attr(n, "initial_n") <- n
            attr(n, "diffusion_rate") <- dispersal_models[[idx]]$get_diffusion_rate()
            attr(n, "diffusion_radius") <- 0
        }
        if (any(sapply(dispersal_models, function(dm) inherits(dm, 
            "AreaSpread")))) {
            capacity <- population_model$get_capacity()
            capacity_area <- attr(capacity, "area")
            if (population_model$get_type() == "stage_structured") {
                stages <- population_model$get_capacity_stages()
                attr(n, "spread_area") <- sum(n[, stages])/as.numeric(capacity) * 
                capacity_area
            } else {
                attr(n, "spread_area") <- n/as.numeric(capacity) * 
                capacity_area
            }
        }
    }
    if (length(impacts)) {
        calc_impacts <- lapply(impacts, function(impacts_i) {
            n <<- impacts_i$calculate(n, 0)
            attr(n, "impacts")
        })
        attr(n, "impacts") <- NULL
        population_model$set_capacity_mult(n)
    } else {
        calc_impacts <- NULL
    }
    if (length(actions)) {
        for (i in 1:length(actions)) {
            n <- actions[[i]]$apply(n, 0)
        }
    }
    results$collate(r, 0, n, calc_impacts)
    if (eq_enabled) {
        eq_state <- eq_capture_step(
            state = eq_state,
            n = n,
            timestep = 0,
            replicate = r
        )
    }
    } # end !resume_this_r initialise block

    tm_end_r <- if (length(warmup_reps_vec) && r %in% warmup_reps_vec) {
        warmup_tm_end
    } else {
        tm_run_end
    }

    benchmark_scaling_reps <- as.integer(Sys.getenv(
        "BENCHMARK_SCALING_REPS", unset = "0"))
    scaling_ckpt_master <- NULL
    if (benchmark_scaling_reps > 0L) {
        if (is.null(checkpoint_data) || !resume_this_r) {
            stop(
                "BENCHMARK_SCALING_REPS requires CHECKPOINT_RESUME for the checkpoint replicate",
                call. = FALSE)
        }
        scaling_ckpt_master <- checkpoint_data
        message(sprintf(
            paste0(
                "BENCHMARK_SCALING_REPS=%d: one process, scenario built once; ",
                "frozen checkpoint + paths copy per rep (no readRDS)"
            ),
            benchmark_scaling_reps))
    }
    scaling_rep_count <- if (benchmark_scaling_reps > 0L) {
        benchmark_scaling_reps
    } else {
        1L
    }

    for (scaling_rep_i in seq_len(scaling_rep_count)) {
        if (!is.null(scaling_ckpt_master)) {
            calc_impacts <- apply_scaling_checkpoint_restore(
                region, population_model, scaling_ckpt_master)
            # n assigned in parent via <<- inside apply_scaling_checkpoint_restore
            message(sprintf(
                "BENCHMARK scaling rep %d/%d (paths copy from frozen checkpoint)",
                scaling_rep_i, scaling_rep_count))
            if (identical(Sys.getenv("BENCHMARK_SCALING_GC"), "1") &&
                scaling_rep_i < scaling_rep_count) {
                gc(FALSE)
            }
        }

        for (tm in tm_start:tm_end_r) {
            t0 <- Sys.time()
            profile_this <- isTRUE(profile_timestep) &&
                r == profile_replicate &&
                tm == profile_timestep_n

            if (length(dispersal_models)) {
                run_dispersal <- function() {
                    n <<- dispersal_models[[1]]$pack(n)
                    for (i in seq_along(dispersal_models)) {
                        t_dm <- Sys.time()
                        n <<- dispersal_models[[i]]$disperse(n, tm)
                        message(sprintf(
                            "    dispersal model %d: %.2fs",
                            i,
                            as.numeric(Sys.time() - t_dm, units = "secs")))
                    }
                    n <<- dispersal_models[[1]]$unpack(n)
                }
            }

            if (profile_this &&
                requireNamespace("profvis", quietly = TRUE) &&
                requireNamespace("htmlwidgets", quietly = TRUE)) {
                profvis_out <- profile_output
                if (!nzchar(profvis_out)) {
                    profvis_out <- file.path(
                        input.env$outputdir,
                        sprintf("profvis_timestep_r%d_tm%d.html", r, tm))
                }
                prof_dir <- dirname(profvis_out)
                if (nzchar(prof_dir)) {
                    dir.create(prof_dir, recursive = TRUE, showWarnings = FALSE)
                }
                message(sprintf(
                    "Profiling timestep at r=%d tm=%d -> %s",
                    r, tm, profvis_out))
                prof_obj <- profvis::profvis({
                    n <<- population_model$grow(n, tm)
                    t1 <<- Sys.time()
                    if (length(dispersal_models)) {
                        run_dispersal()
                    }
                    t2 <<- Sys.time()
                    if (length(impacts)) {
                        calc_impacts <<- lapply(impacts, function(impacts_i) {
                            n <<- impacts_i$calculate(n, tm)
                            attr(n, "impacts")
                        })
                        attr(n, "impacts") <- NULL
                        population_model$set_capacity_mult(n)
                    }
                    t3 <<- Sys.time()
                    if (length(actions)) {
                        for (i in 1:length(actions)) {
                            n <<- actions[[i]]$clear_attributes(n)
                        }
                        for (i in 1:length(actions)) {
                            n <<- actions[[i]]$apply(n, tm)
                        }
                    }
                    t4 <<- Sys.time()
                    if (is.function(user_function)) {
                        n_attr <- attributes(n)
                        if (length(formals(user_function)) == 3) {
                            n <<- user_function(n, r, tm)
                        } else {
                            n <<- user_function(n)
                        }
                        if (length(n_attr)) {
                            for (i in 1:length(n_attr)) {
                                if (!names(n_attr[i]) %in% names(attributes(n))) {
                                    attr(n, names(n_attr[i])) <- n_attr[[i]]
                                }
                            }
                        }
                    }
                    results$collate(r, tm, n, calc_impacts)
                    if (is.function(continued_incursions)) {
                        n <<- continued_incursions(tm, n)
                    }
                })
                prof_ok <- tryCatch({
                    htmlwidgets::saveWidget(prof_obj, profvis_out, selfcontained = TRUE)
                    TRUE
                }, error = function(e) {
                    message("profvis saveWidget failed: ", conditionMessage(e))
                    FALSE
                })
                if (prof_ok && file.exists(profvis_out)) {
                    message("Wrote profvis output to ", profvis_out)
                } else {
                    message("Profvis output not found at ", profvis_out)
                }
            } else {
                n <- population_model$grow(n, tm)
                t1 <- Sys.time()

                if (length(dispersal_models)) {
                    run_dispersal()
                }
                t2 <- Sys.time()

                if (length(impacts)) {
                    calc_impacts <- lapply(impacts, function(impacts_i) {
                        n <<- impacts_i$calculate(n, tm)
                        attr(n, "impacts")
                    })
                    attr(n, "impacts") <- NULL
                    population_model$set_capacity_mult(n)
                }
                t3 <- Sys.time()

                if (length(actions)) {
                    for (i in 1:length(actions)) {
                        n <- actions[[i]]$clear_attributes(n)
                    }
                    for (i in 1:length(actions)) {
                        n <- actions[[i]]$apply(n, tm)
                    }
                }
                t4 <- Sys.time()

                if (is.function(user_function)) {
                    n_attr <- attributes(n)
                    if (length(formals(user_function)) == 3) {
                        n <- user_function(n, r, tm)
                    } else {
                        n <- user_function(n)
                    }
                    if (length(n_attr)) {
                        for (i in 1:length(n_attr)) {
                            if (!names(n_attr[i]) %in% names(attributes(n))) {
                                attr(n, names(n_attr[i])) <- n_attr[[i]]
                            }
                        }
                    }
                }
                results$collate(r, tm, n, calc_impacts)
                if (is.function(continued_incursions)) {
                    n <- continued_incursions(tm, n)
                }
            }
            t5 <- Sys.time()

            elapsed <- function(a, b) as.numeric(b - a, units = "secs")
            mem_info <- ps::ps_system_memory()
            proc_mem <- ps::ps_memory_info(ps::ps_handle(Sys.getpid()))
            gc_time_now <- sum(gc.time())
            gc_step_s <- gc_time_now - gc_time_prev
            gc_time_prev <- gc_time_now
            rng_check <- runif(1)
            if (eq_enabled) {
                eq_state <- eq_capture_step(
                    state = eq_state,
                    n = n,
                    timestep = tm,
                    replicate = r,
                    grow_s = elapsed(t0, t1),
                    dispersal_s = elapsed(t1, t2),
                    impacts_s = elapsed(t2, t3),
                    actions_s = elapsed(t3, t4),
                    other_s = elapsed(t4, t5),
                    total_s = elapsed(t0, t5),
                    rng_check = rng_check
                )
            }
            bench_rep_prefix <- if (benchmark_scaling_reps > 0L) {
                sprintf("BENCHMARK_REP=%d ", scaling_rep_i)
            } else {
                ""
            }
            message(sprintf(
                paste0(
                    "%stm=%s r=%s | grow=%.2fs  dispersal=%.2fs  ",
                    "impacts=%.2fs  actions=%.2fs  other=%.2fs | total=%.2fs | ",
                    "%s sys avail | %s rss | gc=%.3fs | rng check=%.4f"
                ),
                bench_rep_prefix, tm, r,
                elapsed(t0, t1),
                elapsed(t1, t2),
                elapsed(t2, t3),
                elapsed(t3, t4),
                elapsed(t4, t5),
                elapsed(t0, t5),
                prettyunits::pretty_bytes(mem_info$avail),
                prettyunits::pretty_bytes(proc_mem["rss"]),
                gc_step_s,
                rng_check
            ))

            if (checkpoint_tm > 0L && tm == checkpoint_tm && r == checkpoint_replicate) {
                ckpt_out <- checkpoint_file_out(
                    checkpoint_tm, r, input.env$outputdir)
                save_simulation_checkpoint(
                    file = ckpt_out,
                    tm = tm,
                    r = r,
                    n = n,
                    region = region,
                    calc_impacts = calc_impacts
                )
                if (checkpoint_stop) {
                    message("CHECKPOINT_STOP=1: exiting after checkpoint save")
                    quit(save = "no", status = 0)
                }
            }
        }
    } # end scaling_rep_i
}
} # !benchmark_parallel_replicates

benchmark_scaling_reps_done <- as.integer(Sys.getenv(
    "BENCHMARK_SCALING_REPS", unset = "0"))
if (benchmark_scaling_reps_done > 0L) {
    message("BENCHMARK_SCALING_REPS complete")
    quit(save = "no", status = 0)
}

if (eq_enabled) {
    eq_files <- eq_write_outputs(eq_state)
    message(sprintf(
        "Equivalence outputs saved: metrics=%s occupancy=%s",
        eq_files$metrics_path,
        eq_files$occupancy_path
    ))
}

if (benchmark_skip_save) {
    message("BENCHMARK_SKIP_SAVE=1: skipping results finalize and output writes")
} else {
    results$finalize()
    rate_limited_weblog(sprintf(
        "Completed timestep %s of replicate %s (100%%)",
        input.params$simulator$time_steps,
        input.params$simulator$replicates
    ))
    message("Simulation complete. Saving outputs..")

    setwd(input.env$outputdir)

    result_layers <- NULL
    if (region$get_type() == "grid") {
        result_layers <- results$save_rasters(
            gdal = c("COMPRESS=LZW", "TILED=YES"),
            overwrite = TRUE
        )
    }

    results$save_csv()
    tryCatch(
        results$save_plots(width = 1000, height = 500),
        error = function(e) {
            print(e)
        }
    )

    if (region$get_type() == "grid" && !is.null(result_layers)) {
        tryCatch(
            generate_result_animations(region, population_model, result_layers),
            error = function(e) {
                warning(
                    sprintf("Animation generation failed: %s.", e$message),
                    call. = FALSE
                )
            }
        )
    }
}
