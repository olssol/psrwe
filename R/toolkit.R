#' @title Export results into a template file
#'
#' @param numbers Vector of results.
#' @param template.f Template file name.
#' @param out.f Output file name.
#' @param sub.str Pattern of string to be replaced.
#'
#' @noRd
tkExpRst <- function(numbers, template_f,
                     out_f = "rst.txt", sub_str = "AA") {
    if (!file.exists(template_f)) {
        return(NULL)
    }
    
    # Read template
    tpla <- readChar(template_f, file.info(template_f)$size)

    # Substitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub_str, numbers[i], tpla);
    }

    # Write out
    write(tpla, file = out_f)
}


#' @title Import objects in a list into a designated environment
#'
#' @param alist List of objects.
#' @param dest_env Designated environment.
#'
#' @noRd
tkMakeLocal <- function(alist, dest_env) {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest_env)
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
