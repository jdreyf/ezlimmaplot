#' Lengthen output table from \code{ezlimma::roast_*} to be plotted with \code{ggplot2}
#'
#' Lengthen output table from \code{ezlimma::roast_*} to be plotted with \code{ggplot2} in \code{ezlimmaplot}.
#' Name is based on \code{tidyr::pivot_longer} function.
#'
#' @param tab Table of output from \code{ezlimma::roast_*}.
#' @param prefix.v Character vector of prefixes.
#' @param direction.v Character vector of directions to include in output.
#' @return Table of pathway statistics to be plotted using \code{ggplot2}.
#' @export

# don't need to import !!: https://stackoverflow.com/questions/55383205/how-to-use-rlang-operators-in-a-package
pivot_roast_longer <- function(tab, prefix.v, direction.v = c("Up", "Down", "Mixed")){
  stopifnot(direction.v %in% c("Up", "Down", "Mixed"), "NGenes" %in% colnames(tab))
  ds <- data.frame()
  for (prefix in prefix.v){
    for (dir.tmp in direction.v){
      col.lab.nm <- paste(prefix, dir.tmp, sep=".")
      if (dir.tmp != "Mixed"){
        prop.suffix <- paste0("Prop", dir.tmp, "P05")
        dir.col.tmp <- paste0(prefix, ".Direction")
        rows.tmp <- tab |> dplyr::pull(!!rlang::sym(dir.col.tmp)) == dir.tmp
        if (!any(rows.tmp)) next
        tab.tmp <- tab[rows.tmp,]
        prop.v <- tab.tmp |> dplyr::pull(!!rlang::sym(paste(prefix, prop.suffix, sep=".")))
      } else {
        tab.tmp <- tab
        prop.suffix.v <- paste0("Prop", c("Up", "Down"), "P05")
        # stopifnot(paste(prefix, prop.suffix.v, sep=".") %in% colnames(tab.tmp))
        # prop.v <- rowSums(tab.tmp[, paste(prefix, prop.suffix.v, sep=".")])
        prop.v <- tab.tmp |> dplyr::select(!!!rlang::syms(paste(prefix, prop.suffix.v, sep="."))) |> rowSums()
        prefix <- paste(prefix, "Mixed", sep=".")
      }
      ds <- dplyr::bind_rows(ds, dplyr::bind_cols(Comparison = col.lab.nm, Pwy = rownames(tab.tmp), Prop_p05 = prop.v,
                                                  Count = round(tab.tmp$NGenes * prop.v),
                                                  p=tab.tmp |> dplyr::pull(!!rlang::sym(paste0(prefix, ".p"))),
                                                  FDR=tab.tmp |> dplyr::pull(!!rlang::sym(paste0(prefix, ".FDR"))) ))
    }
  }
  ds
}
