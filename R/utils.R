#' @keywords internal
#' @noRd
#' @importFrom methods is
#' @importFrom cli cli_abort
.check_datasources <- function(datasources) {
  if (!(is.list(datasources) &&
        all(unlist(lapply(datasources, function(d) methods::is(d, "DSConnection")))))) {
    cli::cli_abort("'datasources' must be a list of DSConnection objects.")
  }
}


#' @keywords internal
#' @noRd
#' @importFrom cli cli_abort
.check_string <- function(x, arg_name) {
  if (is.null(x)) {
    cli::cli_abort("`{arg_name}` must not be NULL")
  }
  if (!is.character(x) || length(x) != 1L) {
    cli::cli_abort(
      c(
        "x" = "`{arg_name}` must be a single character string",
        "i" = "You provided an object of class {class(x)}"
      )
    )
  }
}
