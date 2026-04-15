#' Parse UniProt Names
#' @importFrom purrr map_chr
#' @importFrom stringr str_split str_extract
#' @noRd
parse_uniprot_names <- function(seq_names) {
  purrr::map_chr(seq_names, function(nm) {
    if (nm %in% c("Consensus", " Sequence Numbering")) return(nm)

    parts <- stringr::str_split(nm, '\\|| ')[[1]]

    if (length(parts) >= 3) {
      org <- stringr::str_extract(nm, 'OS.*OX')
      org <- if (!is.na(org)) gsub('OS\\=| OX', '', org) else "Unknown"

      gene <- gsub('_.*', '', parts[3])
      id <- parts[2]
      return(paste(org, gene, id, sep = ' | '))
    } else {
      return(nm)
    }
  })
}

#' Format Positions
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @noRd
format_positions <- function(df) {
  df %>% dplyr::mutate(
    group = substr(formatC(position, width = 5, format = "d", flag = "0"), 1, 3),
    Position = substr(formatC(position, width = 5, format = "d", flag = "0"), 4, 5)
  )
}

#' Clean List Columns
#' @noRd
clean_list_col <- function(x, extract_col = NULL) {
  if (is.list(x)) {
    sapply(x, function(v) {
      if (is.null(v) || length(v) == 0) return(NA_character_)
      if (is.data.frame(v)) {
        if (!is.null(extract_col) && extract_col %in% names(v)) {
          return(paste(v[[extract_col]], collapse = ", "))
        } else {
          return(paste(v[[1]], collapse = ", "))
        }
      }
      paste(v, collapse = ", ")
    })
  } else {
    as.character(x)
  }
}
