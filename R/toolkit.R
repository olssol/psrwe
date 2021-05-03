#' @title Export results into a template file
#'
#' @param numbers Vector of results.
#' @param template.f Template file name.
#' @param out.f Output file name.
#' @param sub.str Pattern of string to be replaced.
#'
#' @noRd
tkExpRst <- function(numbers, template.f,  
                     out.f = "rst.txt", sub.str = "AA") {
    if (!file.exists(template.f)) {
        return(NULL)
    }
    
    # Read template
    tpla <- readChar(template.f, file.info(template.f)$size)

    # Substitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    # Write out
    write(tpla, file = out.f)
}


#' @title Import objects in a list into a designated environment
#'
#' @param alist List of objects.
#' @param dest.env Designated environment.
#'
#' @noRd
tkMakeLocal <- function(alist, dest.env) {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env)
    }
}


#' @title Call function by its name organized as a vector
#'
#' @param vec Function names as a vector.
#' @param ... Parameters needed for the actual function.
#'
#' @noRd
tkCallFun <- function(vec, ...) {
    rst <- NULL
    eval(parse(text = paste("rst <- ",
                            paste0(vec, collapse = ""),
                            "(...)",
                            sep = "")
               )
         )
    rst
}
