#' Format community data
#' @description
#' Formats species x site table into a long data frame of community dissimilarities
#'
#' @param com community matrix (sites x species)
#' @param id Optional: Unique name for each location
#' @param method Method used to calculate community dissimilarity between pairs of sites. Either a distance metric supported by vegan, 'abcd' to obtain each component of the binary contingecy matrix, 'decomp1' or 'decomp2' for Pierre Legendre's or Andres Baselga's turnover and nestedness decomposition
#' @param num_den For some dissimilarity indices, whether function should return the numerator and denominator in two different columns
#' @param trans Transformation to be applied to raw data before calculating distance matrix. Either a function or 'binary' to transform to presence-absence)
#' @param drop_empty_rows Should empty rows (i.e., sites containing 0s only) be removed?
#' @param drop_empty_cols Should empty columns (i.e., species without observations) be removed?
#' @param binary  Should abundance data be transform into presence absence?
#' @param check_trans Check for non-finite values in transformation
#' @param na.rm Remove rows with NA?
#' @param na.replace Replace NAs with specific value
#'
#' @return A data frame the following elements
#' \itemize{
#'   \item s1, s2  -  A combination of sites (rows).
#'   \item Additional columns indicating the dissimilarity between sites. If `method = 'abcd'` or `num_den = TRUE`, the data frame will contain these individual elements as different columns.
#' }
#'
#' @examples
#'
#' n_sites <- 10
#' n_sp <- 10
#' data <- matrix(rpois(n_sites*n_sp, 10), ncol = n_sp)
#'
#' ID <- paste0('sample_', 1:nrow(data))
#'
#' com <- make_y_df(com = data,
#'             id = ID,
#'             trans = log1p)
#'
#' head(com)
#'
#'@export

make_y_df <-
  function(com = NULL,
           id = NULL,
           method = 'bray',
           num_den = FALSE,
           trans = NULL,
           drop_empty_rows = FALSE,
           drop_empty_cols = TRUE,
           binary = FALSE,
           check_trans = TRUE,
           na.rm = FALSE,
           na.replace = NULL){

    if(is.null(com)){
      cli::cli_abort(c("{.var com} must be supplied",
                     "x" = "You have not provided a community matrix"))
    }
    # ---- Check data class ----
    if(!(is.matrix(com) | is.data.frame(com))) {
      cli::cli_abort(c("{.var com} must be a site x species matrix or a data frame",
                     "x" = "You have supplied a {.cls {class(com)}}"))
    }

    #Convert to matrix
    try(com <- as.matrix(com, rownames =TRUE, colnames = TRUE))

    tryCatch({
      class(com) <- 'numeric'
    },
    # Error message
    error = function(cond) {
      cli::cli_abort('x' = '{.var com} Community matrix "com" could not be transformed into numeric format')
    })

    if (method == 'sorensen') {
      binary = TRUE
    }


    # ---- NA handling ----

    if(!is.null(na.replace)){
      tryCatch({
        n_na <- sum(is.na(com))
        if(n_na> 0){
          com[is.na(com)] <- na.replace
          cli::cli_warn(c("Values replaced!",
                          "i" = "{n_na} NA{?s} {?has/have} been replaced with {na.replace}"))
        }
      }, error = function(cond){
        cli::cli_abort('x' = 'Could not substiture NAs with {na.replace}')
      })

    }else if(na.rm == TRUE){

      # rows containing missing values
      NA_row <- apply(com, 1, function(x) any(is.na(x)))

      # remove samples from id and env
      id <- id[!NA_row]
      com <- com[!NA_row,]

      cli::cli_warn(c("Rows removed!",
                      "i" = "{sum(NA_row)} row{?s} containing NAs {?has/have} been removed"))
    }

    # ---- Transform matrix ----
    # Remove empty row (sites)
    if(drop_empty_rows){
      row_drop <- rowSums(com) == 0 #Check for empty rows

      if(any(row_drop)){
        com <- com[!row_drop,]

        if(is.null(id)){
          site_char <- paste(which(row_drop), sep = ",", collapse = ",")
        }else{
          site_char <- paste(id[which(row_drop)], sep = ",", collapse = ",")
          id <- id[!row_drop]
        }
        cli::cli_warn(c("Sites removed!",
                        "i" = "{sum(row_drop)} empty row{?s} (all 0s) {?has/have} been removed",
                        "i" = 'Row{?s}: {site_char}'))
      }
    }

    # Remove empty col (sp)
    if(drop_empty_cols){
      col_drop <- colSums(com) == 0 #Check for empty cols

      if(any(col_drop)){
        com <- com[,!col_drop]

        sp_char <- paste(names(col_drop)[(which(col_drop))], sep =",", collapse = ",")
        cli::cli_warn(c("Species removed!",
                        "i" = "{sum(col_drop)} empty column{?s} (all 0s) {?has/have} been removed",
                        "i" = 'Column{?s}: {sp_char}'))
      }
    }

    if (!is.null(trans)) {
      #Apply provided function
      if (is.function(trans)) {
        # Try function
        tryCatch({
          com <- trans(com)
        },
        # Error message
        error = function(cond) {
          cli::cli_end('Unable to transform data with function provided')
        })

        if (all(!is.finite(com)) & check_trans) {
          cli::cli_abort(c('Non-infinite values!',
                          'x' = 'transformation produced non-infinite values only'))

        }else if (any(!is.finite(com)) & check_trans) {
          cli::cli_warn(c('Non-infinite values created!',
                          'i' = 'transformation produced {sum(!is.finite(com))} non-infinite value{?s}'))
        }
      } else {
        cli::cli_abort(
            c('{.var trans} must be a function',
              'i' = '{.var trans} must be the name of a function such as log, log1p or sqrt',
              'i' = 'For instnace: trans = sqrt, not trans = "sqrt" or trans = sqrt())'
              ))

      }
    }

    if(binary){
      com = (com>0)
    }

    # ---- Calculate distance matrix ----

    row_col <- expand.grid(s1 =1:nrow(com), s2 = 1:nrow(com))
    row_col <- row_col[row_col$s1 > row_col$s2,]

    # This could be parallelised in the future´
    if(num_den){
      if(method == 'sorensen'){

        betas <- lapply(c('a', 'b', 'c'),
                        function(x) {
                          as.matrix(vegan::designdist(com,method = x,abcd = TRUE))
                        })

        a = betas[[1]][cbind(row_col$s1, row_col$s2)]
        b = betas[[2]][cbind(row_col$s1, row_col$s2)]
        c = betas[[3]][cbind(row_col$s1, row_col$s2)]

        dist.data <- data.frame(row_col,
                                num_sor = (b+c),
                                den_sor = (2*a+b+c))

      }else if(method == 'jaccard'){

        betas <- lapply(c('a', 'b', 'c'),
                        function(x) {
                          as.matrix(vegan::designdist(com,method = x,abcd = TRUE))
                        })

        a = betas[[1]][cbind(row_col$s1, row_col$s2)]
        b = betas[[2]][cbind(row_col$s1, row_col$s2)]
        c = betas[[3]][cbind(row_col$s1, row_col$s2)]

        dist.data <- data.frame(row_col,
                                num_jac = (b+c),
                                den_jac = (a+b+c))

      }else if(method == 'bray'){

        bray1 <- data.frame(t(apply(row_col,
                                    MARGIN = 1,
                                    function(x){
                                      bray1(com[x[1],], com[x[2],])
                                    })))

        names(bray1) <- c('num_bra','den_bra')
        dist.data <- data.frame(row_col,bray1)

      }else{
        cli::cli_abort(c('{.var num_den} not supported for {.var method} provided!',
                         'i' = '{.var num_den} only supported for methods "bray", "sorensen" and "jaccard"',
                         'i' = '{.var method} provided: {method}'))
      }
    }else{
      if (method == 'abcd') {
        # List of a, b, c, d, components
        betas <- lapply(c('a', 'b', 'c', 'd'), function(x) {
          as.matrix(vegan::designdist(com, method = x, abcd = TRUE))
        })

        # Create data frame
        dist.data <- data.frame(
          row_col,
          a = betas[[1]][cbind(row_col$s1, row_col$s2)],
          b = betas[[2]][cbind(row_col$s1, row_col$s2)],
          c = betas[[3]][cbind(row_col$s1, row_col$s2)],
          d = betas[[4]][cbind(row_col$s1, row_col$s2)]
        )
      } else if (method %in% c('decomp1', 'decomp2')) {
        betas <- lapply(c('a', 'b', 'c'), function(x) {
          as.matrix(vegan::designdist(com, method = x, abcd = TRUE))
        })

        a = betas[[1]][cbind(row_col$s1, row_col$s2)]
        b = betas[[2]][cbind(row_col$s1, row_col$s2)]
        c = betas[[3]][cbind(row_col$s1, row_col$s2)]

        if (method == 'decomp1') {
          #Legendre's apprach
          dist.data <- data.frame(
            row_col,
            sim = 2 * pmin(b, c) / (2 * a + b + c),
            sne = abs(b - c) / (2 * a + b + c),
            sor = (b + c) / (2 * a + b + c)
          )

        } else if (method == 'decomp2') {
          #Baselga's approach
          dist.data <- data.frame(
            row_col,
            sim = pmin(b, c) / (a + pmin(b + c)),
            sor = (b + c) / (2 * a + b + c)
          )
            dist.data$sne <- dist.data$sor - dist.data$sim
        }
      } else {

        if(method == 'sorensen'){
          method = 'bray'
        }

        tryCatch({
          betas <- as.matrix(vegan::vegdist(com, method = method))
        }, error = function(cond) {
          cli::cli_abort(c('Unable to calcualte community dissimilarity matrix!',
                           'i' = '{.var method} should be one of:',
                           '- A method copatible compatible with vegan `vegdist`',
                           '- "decomp1" or "decomp2" for Sorensen-based decomposition of beta-diversity into species replacement and richness differences',
                           '- "abcd" for individual elements of the binary contingency table'))
        })

        dist.data <- data.frame(row_col, diss = betas[cbind(row_col$s1, row_col$s2)])
      }
    }

    if (!is.null(id)) {
      dist.data$s1 = id[dist.data$s1]
      dist.data$s2 = id[dist.data$s2]
      row.names(com) <- id
    }

    return(dist.data)
  }
