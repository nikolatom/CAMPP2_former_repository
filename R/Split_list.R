#' @title Split list
#' @description Spliting list of files using ","
#' @param my.list list
#' @export
#' @return split list
#' @examples \dontrun{
#' ...
#' }

SplitList <- function(my.list) {
#    my.list <- as.character(unlist(strsplit(my.list, split=",")))
    my.list <- as.character(unlist(my.list))

    return(my.list)
}
