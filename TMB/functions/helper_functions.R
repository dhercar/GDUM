print_title <- function(text, width = getOption("width"), symb = '-') {
  # Full-width line of dashes
  separator_line <- strrep(symb, width)
  
  # Centered text
  padding <- floor((width - nchar(text)) / 2)
  if (padding < 0) padding <- 0
  centered_text <- paste0(strrep(" ", padding), text)
  
  # Print all
  cat(separator_line, "\n")
  cat(centered_text, "\n")
  cat(separator_line, "\n\n")
}
