#' Create Export Tables
#'
#' @description
#' Create Export Tables for object of class "saegamlss"
#'
#'
#'
#' @param x An R object of class "saegamlss"
#' @param caption Character vector of length 1 or 2 containing the table's caption or title. If length is 2, the second item is the "short caption" used when LaTeX generates a "List of Tables". Set to NULL to suppress the caption. Default value is NULL
#' @param label Character vector of length 1 containing the LaTeX label or HTML anchor. Set to NULL to suppress the label. Default value is NULL
#' @param align Character vector of length equal to the number of columns of the resulting table, indicating the alignment of the corresponding columns
#' @param digits Numeric vector of length equal to one (in which case it will be replicated as necessary) or to the number of columns of the resulting table or matrix of the same size as the resulting table, indicating the number of digits to display in the corresponding columns
#' @param display Character vector of length equal to the number of columns of the resulting table, indicating the format for the corresponding columns
#' @param auto Logical, indicating whether to apply automatic format when no value is passed to align, digits, or display. This ‘autoformat’ (based on xalign, xdigits, and xdisplay) can be useful to quickly format a typical matrix or data.frame. Default value is FALSE
#' @param ... Other parameters
#'
#' @return an xtable
#'
#' @export
#'
#' @examples
#'
#' est <-  sa_p_index(sample = s_data, y = "y",
#'                    sa = "sa", fdis = "LOGNO",
#'                    sigma.f = TRUE, index = "all")
#'
#' xtable(x = est)
#'

xtable.saegamlss <- function(x, caption = NULL, label = NULL, align = NULL, digits = NULL,
                             display = NULL, auto = FALSE, ...){

  y <- as.data.frame.saegamlss(x)

  xtab <- xtable::xtable(x = y, caption = caption, label = label, align = align,
                         digits = digits, display = display, auto = auto)

  return(xtab)

}
