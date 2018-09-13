#' Generate Cache Key
#'
#' Generates a unique hash key from a list of objects.
#'
#' @param ... Any number of objects used to make a key.
#' @return A unique key.
#'
#' @seealso \code{\link{cache.save}}, \code{\link{cache.load}}
#' @export
cache.key <- function(...){

  if(!requireNamespace("digest", quietly = TRUE)){
    stop("Uh oh! This method depends on digest. ",
         "Try running: install.packages('digest')")
  }

  args <- as.list(eval(list(...)))
  key <- paste0("cache", digest::digest(args, serialize = TRUE))
  return(key)
}

#' Save Value to Cache
#'
#' Saves a value in raw to an SQLite database under a specified key.
#'
#' @param key A unique key for the cache.
#' @param value The value to associate with this key.
#' @param where The cache file.
#'
#' @seealso \code{\link{cache.key}}, \code{\link{cache.load}}
#' @export
cache.save <- function(key, value, where = paste0(getwd(), "/.cache")){

  if(!requireNamespace("RMySQL", quietly = TRUE)){
    stop("Uh oh! This method depends on RMySQL. ",
         "Try running: install.packages('RMySQL')")
  }

  if(!requireNamespace("RSQLite", quietly = TRUE)){
    stop("Uh oh! This method depends on RSQLite. ",
         "Try running: install.packages('RSQLite')")
  }

  con <- RMySQL::dbConnect(RSQLite::SQLite(), where)
  value.serial <- data.frame("value" = serialize(value, connection = NULL))
  RMySQL::dbWriteTable(con, name = key, value = value.serial, overwrite = TRUE)
  RMySQL::dbDisconnect(con)
  return(TRUE)
}

#' Load Value from Cache
#'
#' Loads a value in raw from an SQLite database under a specified key.
#'
#' @param key A unique key for the cache.
#' @param where The cache file.
#' @return The imported value.
#'
#' @seealso \code{\link{cache.key}}, \code{\link{cache.save}}
#' @export
cache.load <- function(key, where = paste0(get(), "/.cache")){

  if(!requireNamespace("RMySQL", quietly = TRUE)){
    stop("Uh oh! This method depends on RMySQL. ",
         "Try running: install.packages('RMySQL')")
  }

  if(!requireNamespace("RSQLite", quietly = TRUE)){
    stop("Uh oh! This method depends on RSQLite. ",
         "Try running: install.packages('RSQLite')")
  }

  con <- RMySQL::dbConnect(RSQLite::SQLite(), where)
  value.serial <- RMySQL::dbReadTable(con, key)
  RMySQL::dbDisconnect(con)

  value.raw <- as.raw(as.hexmode(value.serial$value))
  value <- unserialize(value.raw)
  return(value)
}
