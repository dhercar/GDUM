#' Format environmental data
#'
#' @description
#' Formats species x site table into a long data frame of community dissimilarities
#'
#' @param env community matrix (sites x environmental variables)
#' @param id Optional: Unique name for each location
#' @param collapse Whether environmental distances should be calculated for each variable (`collapse = FALSE`) or as a whole
#' @param scale Whether environmental variables should be scaled before calculating environmental distances. This is particularly relevant when collapse = TRUE
#' @param method Distance method passed to `stats::dist`
#' @param var_type Single value or vector indicating variable type. For categorical variables, the distance between samples is 0 when samples belong to the same category and 1 otherwise.
#' @param var_trans Transformation to be applied to raw data before calculating distance. Single function or vector with ether functions or NA for variables that should not be transformed
#' @param check_trans Check for non-finite values after transformation(s)
#' @param na.rm Whether rows with NAs should be removed?
#' @param na.replace If not NULL, NAs are substituted with value provided
#'
#' @return A vector or data frame with environmental distances
#' @export
#'
#' @examples
#'
#' n_sites = 50
#' data <- data.frame( x1 = rnorm(n_sites),
#'                    x2 = rpois(n_sites, 10),
#'                   x3 = sample(c('1','2','3'),
#'                               size = n_sites,
#'                               replace = TRUE))
#'
# default: Individual distances with no transformation
#'make_x_df(data)
#'
# Transform data
#'make_x_df(data,
#'          var_type = 'numeric',
#'          var_trans = c(NA, log1p, NA))
#'
#'# Collapse into one environmental distance
#'make_x_df(data,
#'         collapse  = TRUE,
#'          var_type = 'numeric')


make_x_df <-
  function(env = NULL,
           id = NULL,
           collapse = FALSE,
           scale = TRUE,
           method = 'euclidean',
           var_type = NULL,
           var_trans = NULL,
           check_trans = TRUE,
           na.rm = FALSE,
           na.replace = NULL
  ){

    # ---- Check env class ----
    if(is.null(env)){
      cli::cli_abort(c("{.var env} must be supplied",
                       "x" = "You have not provided an environmental matrix"))
    }
    # ---- Check data class ----
    if(!(is.matrix(env) | is.data.frame(env))) {
      cli::cli_abort(c("{.var env} must be a site x variable matrix or a data frame",
                       "x" = "You have supplied a {.cls {class(env)}}"))
    }

    try(env <- data.frame(env))

    # ---- Transform matrix ----
    if(is.null(var_type)){
      var_type = unlist(sapply(env, class))
    }

    if(any(!(var_type %in% c('numeric', 'factor', 'character', 'integer')))){
      cli::cli_abort(c('{.var var_type} contains unknown variable types',
                       'i' = '{.var var_type} must be either "character", "factor", "interger" or "numeric"',
                       'i' = '{.var var_type} provided: {unique(var_type)}'))
    }

    if(length(var_type) == 1){
      var_type = rep(var_type, ncol(env))

    }else if(length(var_type) > 1){
      if(length(var_type) != ncol(env)){
        cli::cli_abort(c('length of {.var var_type} does not match the dimensions of the environmental matrix',
                         'i' = 'environmental matrix {.var env} has {ncol(env)} element{?s}, but {.var var_type} has {length(var_type)}',
                         'i' = '{.var var_type} must be either "numeric", "factor", or a vector with as many elements as environmental variables in {.var env}'))
      }
    }

    # Ensure variables are in the right class
    tryCatch({
      env[,var_type == 'factor'] <- sapply(env[,var_type == 'factor'], as.factor)
      env[,var_type == 'character'] <- sapply(env[,var_type == 'character'], as.factor)
      env[,var_type == 'numeric'] <- sapply(env[,var_type == 'numeric'], as.numeric)
      env[,var_type == 'integer'] <- sapply(env[,var_type == 'integer'], as.numeric)

    }, error = function(cond) {
      cli::cli_abort(c('Variable type in {.var} not compatible with data provided in{.var env}',
                     'i' = '{.var var_type} provided: {unique(var_type)}'))
    })

    # ---- NA handling ----

    if(!is.null(na.replace)){
      tryCatch({
        n_na <- sum(is.na(env))
        if(n_na> 0){
          env[is.na(env)] <- na.replace
          cli::cli_warn(c('NAs have been produced!',
                          'i' = '{n_na} {?has/have} been replaced with {na.replace}'))
        }
      }, error = function(cond){
        cli::cli_abort(c('Could not substiture NAs with values supplied in {.var replace}',
                       'i' = '{.var na.replace}: {na.replace}'))
      })

    }else if(na.rm == TRUE){

      # rows containing missing values
      NA_row <- apply(env, 1, function(x) any(is.na(x)))

      # remove samples from id and env
      id <- id[!NA_row]
      env <- env[!NA_row,]

      if (sum(NA_row) > 0) {
      cli::cli_warn(c('Rows with NAs removed!',
                      'i'  = '{(sum(NA_row))} row(s) containing NAs have been removed'))
      }
    }

    if(any(is.na(env))){
      cli::cli_abort(c('NAs in environmental matrix!',
                       'x' = 'NAs are not allowed in {.var env}',
                       'i' = '{sum(is.na(env))} NA{?s} found in {.var env}',
                       'i' = 'This error can also occour when a variable class is misspecified'))
    }


    # Transform numerical variables
    if(!is.null(var_trans)){

      if(length(var_trans) == 1){

        var_trans = rep(list(var_trans), ncol(env))

      }else if(length(var_trans) != ncol(env)){

        cli::cli_abort(c('Uncompatible dimensions!',
                         'x' = '{.var var_trans} must be single function or a list of functions of the same length as the numnber of environmental variables',
                         'i' = '{.var var_trans} contains {length(var_trans)} element{?s} and {.var env} has {ncol(env)} variable{?s}'))

      }

      trans.na <- is.na(var_trans)

      if(any(lapply(var_trans[!trans.na],class) != 'function')){

        cli::cli_abort(c('Uncompatible transormation!',
                       'x' = "{.var var_trans}' must be the name of a function or an ordered vector of functions",
                       'i' = 'For instance: `var_trans = sqrt`, not `var_trans = "sqrt"` or `var_trans = sqrt()`',
                       'i' = 'You can also provide custom functions like this: `var_trans = function(x) log(x + 1)`',
                       'i' = 'If a vector is provided, variables that should not be transformed must be indicated with an NA in {.var var_trans}'))
      }

      for(i in 1:ncol(env)){
        if(is.numeric(env[,i])){

          if(trans.na[i] == TRUE) trans_i = identity else trans_i = var_trans[[i]]

            tryCatch({
              env[, i] <- trans_i(env[, i])

            }, error = function(cond) {
              cli::cli_abort(c('Unable to transform variables',
                               'x' = 'Variable {names(env)[i]} could not be transformed with function provided'))
            })

          if(any(!is.finite(env[,i])) & check_trans){
            cli::cli_abort(c('Transformation produced non-finite values',
                             'x' = 'Transformation of variable {names(env)[i]} produced non-infinite values',
                             'i' = '{sum(!is.finite(env[,i]))} non-infinite value{?s} {?was/were} produced'))
          }
        }
      }
    }

    # ---- Calculate distance matrix ----

    row_col <- expand.grid(s1 =1:nrow(env), s2 = 1:nrow(env))
    row_col <- row_col[row_col$s1 > row_col$s2,]

    if(collapse){ #One single euclidean distance

      # Hard encode categorical variables
      env <- stats::model.matrix(~., data = env)[,-1]
      #Scale all variables
      if(scale){
        env <- scale(env)
      }
      #Euclidian distance
      dist <- as.matrix(stats::dist(env, method = method))
      #Data
      dist.data <- data.frame(row_col,
                              env_dist = dist[cbind(row_col$s1, row_col$s2)])

    }else{ #Distance for each variable

      env.dist <- list()

      for(i in 1:ncol(env)){

        if(is.numeric(env[,i])){
          if(scale){
            env[,i] <- scale(env[,i])
          }

          dist <- as.matrix(stats::dist(env[,i], method = method))
          dist <- dist[cbind(row_col$s1, row_col$s2)]

        }else{
          dist <- as.numeric(env[row_col$s1,i] != env[row_col$s2,i])
        }

        env.dist[[length(env.dist)+1]] <- dist

      }

      env.dist <- data.frame(do.call(cbind,env.dist))
      names(env.dist) <- paste0("dist_", names(env))

      dist.data <- data.frame(row_col, env.dist)

    }

    if(!is.null(id)){
      dist.data$s1 <- id[dist.data$s1]
      dist.data$s2 <- id[dist.data$s2]
    }

    return(dist.data)
  }



